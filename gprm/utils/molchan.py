import numpy as np
import pandas as pd
import pygmt
from scipy.interpolate import RegularGridInterpolator
from .proximity import contour_proximity, polyline_proximity, polygons_buffer, points_proximity, boundary_proximity, reconstruct_and_rasterize_polygons
from .create_gpml import gpml2gdf
import xarray as xr
import shapely
import pygplates
from .geometry import apply_nearest_feature, apply_reconstruction
from .spatial import topology_lookup


from collections import OrderedDict


DEFAULT_DISTANCE_MAX = 1e7
DEFAULT_DISTANCE_STEP = 2e4
DEFAULT_GEOGRAPHIC_EXTENT = [-180.,180.,-90.,90.]
DEFAULT_GEOGRAPHIC_SAMPLING = 0.25


def scipy_interpolater(da, points):
    f = RegularGridInterpolator((da.x.data,da.y.data), da.data.T, method='linear')
    return f(points)


def molchan_test(grid, 
                 points, 
                 distance_max = DEFAULT_DISTANCE_MAX, 
                 distance_step = DEFAULT_DISTANCE_STEP,
                 buffer_radius=1,
                 interpolater='pygmt'):
    """
    Molchan test for a set of points against one grid
    """
    
    earth_area = (6371.**2)*(4*np.pi)
    
    # distance of each point to the target
    if interpolater=='scipy':
        points['distance'] = scipy_interpolater(grid, points[['Longitude','Latitude']])
        
    else:
        points['distance'] = pygmt.grdtrack(points=pd.DataFrame(data=points[['Longitude','Latitude']]), 
                                            grid=grid, 
                                            no_skip=False, 
                                            interpolation='l',
                                            radius=buffer_radius,
                                            newcolname='dist')['dist']
    
    # percentage of target distance grid within each contour level
    grid_histogram = pygmt.grdvolume(grid, contour=[0, distance_max, distance_step], f='g', unit='k')
    
    (point_histogram,
     bin_edges) = np.histogram(points['distance'], 
                               bins=np.arange(-distance_step, distance_max+distance_step, distance_step))
    
    permissible_area = grid_histogram.iloc[0,1]
    print('Total permissible area is {:0.1f}% of total Earth surface'.format(100*permissible_area/earth_area))
    
    grid_fraction = 1-grid_histogram.iloc[:,1]/permissible_area 
    # Note that this computation will penalize models where there are lots of 
    # invalid points (since they are 'missed' at any grid fraction)
    points_fraction = 1-np.cumsum(point_histogram)/len(points)
    
    Skill = 0.5+np.trapz(grid_fraction, points_fraction)
    
    return grid_fraction[::-1], points_fraction[::-1], Skill



def molchan_point(grid, 
                  points, 
                  distance_max = DEFAULT_DISTANCE_MAX, 
                  distance_step = DEFAULT_DISTANCE_STEP, 
                  buffer_radius=1, 
                  interpolater='pygmt',
                  verbose=False,
                  return_fraction=True):
    """
    Molchan test for a single point
    """
    
    earth_area = (6371.**2)*(4*np.pi)
    
    # distance of each point to the target
    if interpolater=='scipy':
        points['distance'] = scipy_interpolater(grid, points[['Longitude','Latitude']])
    else:
        points['distance'] = pygmt.grdtrack(grid=grid,
                                            points=points[['Longitude','Latitude']], 
                                            no_skip=False, 
                                            interpolation='l',
                                            radius=buffer_radius,
                                            newcolname='dist')['dist']
        
    if not return_fraction:
        return float(points['distance'].values)
    
    # percentage of target distance grid within each contour level
    grid_histogram = pygmt.grdvolume(grid, contour=[0, distance_max, distance_step], 
                                     f='g', unit='k', verbose='e')
    
    permissible_area = grid_histogram.iloc[0,1]
    if verbose:
        print('Total permissible area is {:0.1f}% of total Earth surface'.format(100*permissible_area/earth_area))
    
    grid_fraction = 1-grid_histogram.iloc[:,1]/permissible_area 
    
    area_better_than_points = np.interp(points['distance'], grid_histogram.loc[:,0], grid_histogram.loc[:,1])
    
    return float(points['distance'].values), float(1-area_better_than_points/permissible_area)
    
    

def space_time_molchan_test(raster_dict, 
                            point_distances,
                            healpix_resolution=128,
                            distance_max=DEFAULT_DISTANCE_MAX, 
                            distance_step=DEFAULT_DISTANCE_STEP,
                            interpolater='pygmt'):
    """
    Given a raster sequence and a set of point distances already extracted from them, 
    compute a molchan test result where the grid fraction is summed over all rasters
    in the sequence
    """
    from gprm import PointDistributionOnSphere
    hp = PointDistributionOnSphere(distribution_type='healpix', N=healpix_resolution)
    hp_dataframe = pd.DataFrame(data={'x':hp.longitude, 'y':hp.latitude})

    space_time_distances = []

    for reconstruction_time in raster_dict.keys():
        if interpolater=='scipy':
            smpl = scipy_interpolater(raster_dict[reconstruction_time], 
                                      hp_dataframe)
            space_time_distances.extend(smpl[np.isfinite(smpl)].tolist())
        else:
            smpl = pygmt.grdtrack(
                grid=raster_dict[reconstruction_time], 
                points=hp_dataframe, 
                no_skip=False,
                interpolation='l',
                newcolname='distance'
            )
            space_time_distances.extend(smpl['distance'].dropna().tolist())

    # Determine for both the grids and the points, the fraction of overall
    # points within each distance contour
    (hp_histogram,
     bin_edges) = np.histogram(space_time_distances, 
                               bins=np.arange(-distance_step, distance_max+distance_step, distance_step))

    (pm_histogram,
     bin_edges) = np.histogram(point_distances,
                               bins=np.arange(-distance_step, distance_max+distance_step, distance_step))

    grid_fraction = np.cumsum(hp_histogram)/len(space_time_distances)
    point_fraction = 1-np.cumsum(pm_histogram)/len(point_distances)
    
    Skill = np.trapz(1-grid_fraction, 1-point_fraction) - 0.5

    return grid_fraction, point_fraction, Skill



def combine_raster_sequences(raster_dict1, raster_dict2):
    """
    Given two raster sequences (dictionaries with coincident keys), 
    generate a new raster sequence that multiplies the coincident rasters
    from each sequence
    """
    
    raster_dict3 = OrderedDict()

    for key in raster_dict1.keys():
        raster_dict3[key] = raster_dict1[key] * raster_dict2[key]
        
    return raster_dict3

    

def space_time_distances(raster_dict, gdf, age_field_name='age', 
                         distance_max=DEFAULT_DISTANCE_MAX, 
                         distance_step=DEFAULT_DISTANCE_STEP, 
                         buffer_radius=1,
                         interpolater='pygmt'):
    """
    Computes the distances to targets rconstructed to their time of appearance 
    from a raster sequence of raster grids
    
    The input gdf is assumed to have reconstructed coordinates in its geometry
    """
    
    results = []

    for i,row in gdf.iterrows():
        reconstruction_time = row[age_field_name]
        result = molchan_point(raster_dict[reconstruction_time],
                               pd.DataFrame(data={'Longitude': [row.geometry.x], 
                                                  'Latitude': [row.geometry.y]}),
                               distance_max=distance_max, 
                               distance_step=distance_step, 
                               buffer_radius=buffer_radius, 
                               interpolater=interpolater,
                               )
        results.append(result)

    return pd.DataFrame(data=results, 
                        columns=['distance', 'area_fraction'])

'''

def generate_raster_sequence_from_polygons(features,
                                           rotation_model,
                                           reconstruction_times,
                                           sampling=DEFAULT_GEOGRAPHIC_SAMPLING,
                                           buffer_distance=None):
    """
    Given some reconstrutable polygon features, generates a series of 
    """
    
    raster_dict = OrderedDict()

    for reconstruction_time in reconstruction_times:  
        
        tmp = reconstruct_and_rasterize_polygons(features,
                                                 rotation_model,
                                                 reconstruction_time,
                                                 sampling=sampling)

        tmp = tmp.where(tmp!=0, np.nan)
    
        if buffer_distance is not None:
            bn = boundary_proximity(tmp)
            tmp.data[bn.data<=buffer_distance] = 1
    
        raster_dict[reconstruction_time] = tmp
        
    return raster_dict



def generate_distance_raster_sequence(target_features, 
                                      reconstruction_model,
                                      reconstruction_times,
                                      sampling=DEFAULT_GEOGRAPHIC_SAMPLING,
                                      region=DEFAULT_GEOGRAPHIC_EXTENT):
    
    
    
    prox_grid_sequence = OrderedDict()
    
    for reconstruction_time in reconstruction_times:
        if isinstance(target_features, dict):
        
            r_target_features = gpml2gdf(target_features[reconstruction_time])
    
        else:
            #generate distance raster, masked against the permissive area
            r_target_features = reconstruction_model.reconstruct(target_features, 
                                                                 reconstruction_time, 
                                                                 use_tempfile=False)

        if r_target_features is not None:
            prox_grid = polyline_proximity(r_target_features,
                                           spacing=sampling, 
                                           region=region)
        else:
            prox_grid = np.ones_like(tmp) * np.nan

        prox_grid_sequence[reconstruction_time] = prox_grid
        
    return prox_grid_sequence

'''

from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from tqdm import tqdm
from collections import OrderedDict
import numpy as np


def generate_raster_sequence_from_polygons(features,
                                           rotation_model,
                                           reconstruction_times,
                                           sampling=DEFAULT_GEOGRAPHIC_SAMPLING,
                                           buffer_distance=None,
                                           max_workers=None):
    """
    Given some reconstrutable polygon features, generates a series of rasterized outputs
    using multithreading with progress bar.
    
    Parameters:
    -----------
    features : object
        Polygon features for reconstruction
    rotation_model : object
        Rotation model for reconstruction
    reconstruction_times : list
        List of reconstruction times
    sampling : float
        Geographic sampling parameter
    buffer_distance : float, optional
        Buffer distance for boundary proximity (default: None)
    max_workers : int, optional
        Maximum number of worker threads (default: None uses ThreadPoolExecutor default)
    
    Returns:
    --------
    OrderedDict : Ordered dictionary of rasterized data keyed by reconstruction time
    """
    
    def process_single_time(reconstruction_time):
        """Process a single reconstruction time."""
        try:
            # Reconstruct and rasterize polygons
            tmp = reconstruct_and_rasterize_polygons(features,
                                                   rotation_model,
                                                   reconstruction_time,
                                                   sampling=sampling)

            # Replace 0 values with NaN
            tmp = tmp.where(tmp != 0, np.nan)
        
            # Apply buffer distance if specified
            if buffer_distance is not None:
                bn = boundary_proximity(tmp)
                tmp.data[bn.data <= buffer_distance] = 1
            
            return reconstruction_time, tmp
            
        except Exception as e:
            print(f"Error processing time {reconstruction_time}: {str(e)}")
            return reconstruction_time, None
    
    # Initialize results storage
    results = {}
    
    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_time = {
            executor.submit(process_single_time, time): time 
            for time in reconstruction_times
        }
        
        # Process completed tasks with progress bar
        with tqdm(total=len(reconstruction_times), desc="Processing polygon rasterization") as pbar:
            for future in as_completed(future_to_time):
                reconstruction_time, raster_data = future.result()
                results[reconstruction_time] = raster_data
                pbar.update(1)
    
    # Create OrderedDict maintaining the original order of reconstruction_times
    raster_dict = OrderedDict()
    for reconstruction_time in reconstruction_times:
        raster_dict[reconstruction_time] = results[reconstruction_time]
    
    return raster_dict


def generate_distance_raster_sequence(target_features, 
                                      reconstruction_model,
                                      reconstruction_times,
                                      sampling=DEFAULT_GEOGRAPHIC_SAMPLING,
                                      region=DEFAULT_GEOGRAPHIC_EXTENT,
                                      max_workers=None):
    """
    Generate distance raster sequence using multithreading with progress bar.
    
    Parameters:
    -----------
    target_features : dict or other
        Target features for reconstruction
    reconstruction_model : object
        Model for reconstruction
    reconstruction_times : list
        List of reconstruction times
    sampling : float
        Geographic sampling parameter
    region : object
        Geographic region
    max_workers : int, optional
        Maximum number of worker threads (default: None uses ThreadPoolExecutor default)
    
    Returns:
    --------
    OrderedDict : Ordered dictionary of proximity grids keyed by reconstruction time
    """
    
    def process_single_time(reconstruction_time):
        """Process a single reconstruction time."""
        try:
            if isinstance(target_features, dict):
                if not target_features[reconstruction_time]:
                    prox_grid = zeros_grid_like(sampling=sampling, region=region) * np.nan
                    return reconstruction_time, prox_grid
                elif isinstance(target_features[reconstruction_time][0], pygplates.Feature):
                    r_target_features = gpml2gdf(target_features[reconstruction_time])
                elif isinstance(target_features[reconstruction_time], pd.GeoDataFrame):
                    r_target_features = target_features[reconstruction_time]

            else:
                r_target_features = reconstruction_model.reconstruct(target_features, 
                                                                     reconstruction_time, 
                                                                     use_tempfile=False)

            # Generate distance raster, masked against the permissive area
            if r_target_features is not None:

                if isinstance(r_target_features.geometry.iloc[0], shapely.geometry.point.Point):
                    prox_grid = points_proximity(r_target_features.geometry.x,
                                                 r_target_features.geometry.y,
                                                 spacing=sampling, 
                                                 region=region)

                elif isinstance(r_target_features.geometry.iloc[0], shapely.geometry.linestring.LineString):
                    prox_grid = polyline_proximity(r_target_features,
                                                   spacing=sampling, 
                                                   region=region)

                elif isinstance(r_target_features.geometry.iloc[0], shapely.geometry.polygon.Polygon):    
                    prox_grid = polygons_buffer(r_target_features,
                                                spacing=sampling, 
                                                region=region)
                    
                else:
                    raise ValueError("Unsupported geometry type in target features.")
                    
            else:
                # Note: You'll need to define what 'tmp' should be in this context
                # For now, assuming you want a grid of the same shape as expected output
                #prox_grid = np.ones_like(tmp) * np.nan
                #prox_grid = xr.DataArray(
                #    np.full((len(np.arange(-90, 90+sampling, sampling)), len(np.arange(-180, 180+sampling, sampling))), np.nan),
                #    coords={'y': np.arange(-90, 90+sampling, sampling), 'x': np.arange(-180, 180+sampling, sampling)},
                #   dims=['y', 'x']
                #    )
                prox_grid = zeros_grid_like(sampling=sampling, region=region) * np.nan

            return reconstruction_time, prox_grid
            
        except Exception as e:
            print(f"Error processing time {reconstruction_time}: {str(e)}")
            return reconstruction_time, None
    
    # Initialize results storage
    results = {}
    
    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_time = {
            executor.submit(process_single_time, time): time 
            for time in reconstruction_times
        }
        
        # Process completed tasks with progress bar
        with tqdm(total=len(reconstruction_times), desc="Processing reconstruction times") as pbar:
            for future in as_completed(future_to_time):
                reconstruction_time, prox_grid = future.result()
                results[reconstruction_time] = prox_grid
                pbar.update(1)
    
    # Create OrderedDict maintaining the original order of reconstruction_times
    prox_grid_sequence = OrderedDict()
    for reconstruction_time in reconstruction_times:
        prox_grid_sequence[reconstruction_time] = results[reconstruction_time]
    
    return prox_grid_sequence


def generate_masked_distance_raster_sequence(
        reconstruction_model, boundary_lookup,
        reconstruction_times, sampling=DEFAULT_GEOGRAPHIC_SAMPLING, 
        polygon_buffer_distance=None):
    # Make two raster sequences, where:
    # 1. Mask rasters where the pixels lying within continents (from those lying outside, therefore unreconstructable)
    # 2. Distance rasters from subduction zone geometries

    reconstruction_raster_dict = generate_raster_sequence_from_polygons(
        reconstruction_model.continent_polygons[0],
        reconstruction_model.rotation_model,
        reconstruction_times,
        sampling=sampling,
        buffer_distance=polygon_buffer_distance
    )

    target_distance_dict = generate_distance_raster_sequence(
        boundary_lookup,
        reconstruction_model,
        reconstruction_times,
        sampling=sampling
    )
    
    #Combine the distance rasters with the mask rasters
    target_distance_dict_mask = combine_raster_sequences(target_distance_dict, 
                                                         reconstruction_raster_dict)
    
    return target_distance_dict_mask


def generate_random_distance_sequence(target_distance_dict_mask, healpix_N=64):
    
    from gprm import PointDistributionOnSphere
    hp = PointDistributionOnSphere(distribution_type='healpix', N=healpix_N)
    hp_dataframe = pd.DataFrame(data={'x':hp.longitude, 'y':hp.latitude})

    space_time_distances = []
    for reconstruction_time in target_distance_dict_mask.keys():
        smpl = pygmt.grdtrack(
            grid=target_distance_dict_mask[reconstruction_time], 
            points=hp_dataframe, 
            no_skip=False,
            interpolation='l',
            newcolname='distance'
        )
        smpl['reconstruction_time'] = reconstruction_time
        space_time_distances.append(smpl.dropna(subset=['distance']))

    #space_time_distances
    return pd.concat(space_time_distances)


def zeros_grid_like(sampling=DEFAULT_GEOGRAPHIC_SAMPLING, 
                    region=DEFAULT_GEOGRAPHIC_EXTENT):
    """
    Helper function to get grid of zeros with expected grid shape and sampling
    """
    # This should return a grid with the expected shape based on sampling and region
    # You'll need to replace this with your actual implementation
    # For example:
    # return create_grid_from_region(region, sampling)
    return xr.DataArray(
        np.full((len(np.arange(region[2], region[3]+sampling, sampling)), len(np.arange(region[0], region[1]+sampling, sampling))), 0.),
        coords={'y': np.arange(region[2], region[3]+sampling, sampling), 'x': np.arange(region[0], region[1]+sampling, sampling)},
        dims=['y', 'x']
        )
    #pass



def sample_distance_analysis(data_df, reconstruction_model, 
                             age_field='age', time_min=0, time_max=1000., 
                             reconstruction_time_step=1, targets='subduction'):
    """
    Perform nearest distance analysis between a set of samples and a set of target features,
    using measurements from the geometries themselves. 
    This is more accurate, but possibly slower for larger data sets, compared to using a distance
    raster.
    'targets' is a dictionary of 
    """
    
    # If not provided, create a lookup table for the target features
    if isinstance(targets, dict):
        target_lookup = targets
    elif targets in ['subduction', 'midoceanridge', 'other']:
        target_lookup = topology_lookup(reconstruction_model, 
                                        np.arange(time_min, time_max+reconstruction_time_step, reconstruction_time_step),
                                        boundary_types=[targets])
    else:
        raise ValueError("Unsupported input for targets...")
        
        
    # From the input data, extract the data within the determined age range and assign plateids
    data_select = data_df[(data_df[age_field]<=time_max) & (data_df[age_field]>=time_min)]
    data_select = reconstruction_model.assign_plate_ids(data_select, 
                                                        keep_unpartitioned_features=False,
                                                        copy_valid_times=True)
    # TODO this will hit problems if FROMAGE is already assigned in the input
    data_select = data_select[data_select[age_field] <= data_select['FROMAGE']].reset_index(drop=True)

    # assign a reconstruction time which is the nearest time step to the age associated with the data point
    data_select['reconstruction_time'] = np.round(data_select[age_field]/reconstruction_time_step)*reconstruction_time_step
    
    # reconstruct the points
    data_select['rgeometry'] = data_select.apply(lambda x: apply_reconstruction(x, reconstruction_model.rotation_model), axis=1)
    
    # Determine the shortest distance to the target features at the associated time
    data_select['distance_to_target'] = data_select.apply(
        lambda x: apply_nearest_feature(x, 
                                        target_lookup, 
                                        geometry_field='rgeometry',
                                        age_field='reconstruction_time'), axis=1)

    return data_select

