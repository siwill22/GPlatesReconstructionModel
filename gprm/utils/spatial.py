"""Spatial analysis utilities: polygon rasterisation, plate partitioning, and topology queries."""
import pygplates
from ptt.utils import points_in_polygons
from ptt.utils import points_spatial_tree
from ptt.utils.proximity_query import *
from .create_gpml import create_gpml_regular_long_lat_mesh, create_gpml_healpix_mesh
from .create_gpml import gpml2gdf, gdf2gpml
import numpy as np

import geopandas as gpd
from shapely.geometry import Polygon
from shapely.validation import make_valid


def merge_polygons(polygons, rotation_model,
                   reconstruction_time=0, sampling=1., area_threshold=None, filename=None,
                   return_raster=False, anchor_plate_id=0):
    """Merge reconstructed polygons into contour polygons via grid point-in-polygon test and image contouring.

    :param polygons: pygplates FeatureCollection of polygon features.
    :param rotation_model: pygplates RotationModel.
    :param reconstruction_time: Age in Ma (default 0).
    :param sampling: Grid spacing in degrees used for the point-in-polygon test (default 1).
    :param area_threshold: If set, discard output polygons smaller than this area in steradians.
    :param filename: If set, write the resulting feature collection to this file and return None.
    :param return_raster: If True, return the binary numpy grid instead of polygon features (default False).
    :param anchor_plate_id: Plate ID used as the fixed reference frame (default 0).
    :returns: List of pygplates features, or a 2-D numpy binary array if return_raster=True; None if filename is given.
    """
    from skimage import measure

    multipoints = create_gpml_regular_long_lat_mesh(sampling)
    grid_dims = (int(180/sampling)+1,int(360/sampling)+1)

    for multipoint in multipoints:
        for mp in multipoint.get_all_geometries():
            points = mp.to_lat_lon_point_list()

    bi = run_grid_pip(reconstruction_time,points,polygons,rotation_model,grid_dims,anchor_plate_id=anchor_plate_id)
    
    if return_raster:
        return bi
    
    else:
        # To handle edge effects, pad grid before making contour polygons  
        ## --- start
        pad_hor = np.zeros((1,bi.shape[1]))
        pad_ver = np.zeros((bi.shape[0]+2,1))
        pad1 = np.vstack((pad_hor,bi,pad_hor))      # add row of zeros to top and bottom
        pad2 = np.hstack((pad_ver,pad1,pad_ver))    # add row of zeros to left and right
        #pad3 = np.hstack((pad2,pad_ver))
        contours = measure.find_contours(pad2, 0.5, fully_connected='low')
        ## --- end
    
        contour_polygons = []
        contour_features = []
    
        for n,cp in enumerate(contours):
        
            # To handle edge effects again - strip off parts of polygon
            # due to padding, and adjust from image coordinates to long/lat
            # --- start
            cp[:,1] = (cp[:,1]*sampling)-sampling
            cp[:,0] = (cp[:,0]*sampling)-sampling
            cp[np.where(cp[:,0]<0.),0] = 0
            cp[np.where(cp[:,0]>180.),0] = 180
            cp[np.where(cp[:,1]<0.),1] = 0
            cp[np.where(cp[:,1]>360.),1] = 360
            ## --- end
        
            cpf = pygplates.PolygonOnSphere(zip(cp[:,0]-90,cp[:,1]-180))
            contour_polygons.append(cpf)
        
            feature = pygplates.Feature()
            feature.set_geometry(cpf)
            contour_features.append(feature)

        if filename is not None:
            pygplates.FeatureCollection(contour_features).write(filename)

        else:
            return contour_features


def rasterise_polygons(polygon_features, rotation_model, reconstruction_time, raster_domain_points=None,
					   sampling=0.5, meshtype='LongLatGrid', masking=None):
    """Partition a regular grid or healpix mesh by polygon features, mapping shapefile attributes to each point.

    :param polygon_features: pygplates FeatureCollection of polygon features.
    :param rotation_model: pygplates RotationModel.
    :param reconstruction_time: Age in Ma.
    :param raster_domain_points: Pre-built mesh of domain points; generated from sampling/meshtype if None.
    :param sampling: Grid spacing in degrees, or healpix nSide integer if meshtype='healpix' (default 0.5).
    :param meshtype: 'LongLatGrid' (default) or 'healpix'.
    :param masking: None = return all points; 'outside' = points inside polygons only; 'inside' = points outside polygons only.
    :returns: List of pygplates point features with polygon attributes copied.
    """

    if not raster_domain_points:
        if meshtype=='healpix':
            raster_domain_points = create_gpml_healpix_mesh(sampling,filename=None,feature_type='MeshNode')
        else:
            raster_domain_points = create_gpml_regular_long_lat_mesh(sampling,filename=None,feature_type='MeshNode')

    plate_partitioner = pygplates.PlatePartitioner(polygon_features, rotation_model, reconstruction_time=reconstruction_time)

    if masking is not None:
        pg_points = plate_partitioner.partition_features(raster_domain_points,
														 partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned,
                                                         properties_to_copy=[pygplates.PropertyName.gpml_shapefile_attributes])
        if masking == 'outside':
            pg_points = pg_points[0]
        elif masking == 'inside':
            pg_points = pg_points[1]

    else:
        pg_points = plate_partitioner.partition_features(raster_domain_points,
                                                         properties_to_copy=[pygplates.PropertyName.gpml_shapefile_attributes])

    return pg_points


def force_polygon_geometries(input_features):
    """Convert all feature geometries to closed polygons, preserving plate IDs and valid times."""

    polygons = []
    for feature in input_features: 
        for geom in feature.get_all_geometries():
            polygon = pygplates.Feature(feature.get_feature_type())
            polygon.set_geometry(pygplates.PolygonOnSphere(geom))
            polygon.set_reconstruction_plate_id(feature.get_reconstruction_plate_id())
            # some features in COBTerranes had invalid time ranges - these with throw an error if 
            # we try to create a new feature with same times
            if feature.get_valid_time()[0]>=feature.get_valid_time()[1]:
                polygon.set_valid_time(feature.get_valid_time()[0],feature.get_valid_time()[1])
                polygons.append(polygon)
    polygon_features = pygplates.FeatureCollection(polygons)

    return polygon_features


def polygon_area_threshold(polygons, area_threshold):
    """Return only the polygons whose area exceeds area_threshold (in steradians)."""
    polygons_larger_than_threshold = []
    for polygon in polygons:
        if polygon.get_geometry() is not None:
            if polygon.get_geometry().get_area()>area_threshold:
                polygons_larger_than_threshold.append(polygon)

    return polygons_larger_than_threshold


def run_grid_pip(time, points, polygons, rotation_model, grid_dims, anchor_plate_id=0):
    """Point-in-polygon test on a regular grid of points using reconstructed polygons; returns binary array.

    :param time: Reconstruction age in Ma.
    :param points: List of pygplates PointOnSphere objects covering the grid.
    :param polygons: pygplates FeatureCollection of polygon features.
    :param rotation_model: pygplates RotationModel.
    :param grid_dims: Tuple (n_rows, n_cols) giving the shape of the output array.
    :param anchor_plate_id: Plate ID used as the fixed reference frame (default 0).
    :returns: 2-D numpy binary array (1 = inside a polygon, 0 = outside) with shape grid_dims.
    """
    reconstructed_polygons = []
    pygplates.reconstruct(polygons,rotation_model,reconstructed_polygons,time,anchor_plate_id=anchor_plate_id)

    rpolygons = []
    for polygon in reconstructed_polygons:
        if polygon.get_reconstructed_geometry():
            rpolygons.append(polygon.get_reconstructed_geometry())

    polygons_containing_points = points_in_polygons.find_polygons(points, rpolygons)

    lat = []
    lon = []
    zval = []
    for pcp,point in zip(polygons_containing_points,points):
        lat.append(point.get_latitude())
        lon.append(point.get_longitude())
        if pcp is not None:
            zval.append(1)
        else:
            zval.append(0)

    bi = np.array(zval).reshape(grid_dims[0],grid_dims[1])

    return bi


# Function to run efficient point in/near polygons
# returns two numbers - one is distance to polygon edge,
def run_grid_pnp(recon_time,
                 points,
                 spatial_tree_of_uniform_recon_points,
                 polygons,
                 rotation_model,
                 distance_threshold_radians=2,
                 anchor_plate_id=0):
    """Distance to polygon boundaries on a grid using a spatial tree; returns inside and boundary distance arrays.

    :param recon_time: Reconstruction age in Ma.
    :param points: List of pygplates PointOnSphere objects (uniform grid).
    :param spatial_tree_of_uniform_recon_points: Pre-built spatial tree over the same point list (from points_spatial_tree).
    :param polygons: pygplates FeatureCollection of polygon features.
    :param rotation_model: pygplates RotationModel.
    :param distance_threshold_radians: Search radius for the spatial tree query in radians (default 2).
    :param anchor_plate_id: Plate ID used as the fixed reference frame (default 0).
    :returns: Tuple (distance_to_polygon, distance_to_polygon_boundary) as 1-D numpy arrays in radians;
        distance_to_polygon is 0 for points inside a polygon.
    """
    reconstructed_polygons = []
    pygplates.reconstruct(polygons, rotation_model, reconstructed_polygons, recon_time, anchor_plate_id=anchor_plate_id)
    rpolygons = []
    for polygon in reconstructed_polygons:
        if polygon.get_reconstructed_geometry():
            rpolygons.append(polygon.get_reconstructed_geometry())
                
    res1 = find_closest_geometries_to_points_using_points_spatial_tree(points,
                                                                    spatial_tree_of_uniform_recon_points,
                                                                    rpolygons,
                                                                    distance_threshold_radians = distance_threshold_radians,
                                                                    geometries_are_solid = False)

    distance_to_polygon_boundary = np.array(list(zip(*res1))[0])

    # Make a copy of list of distances.
    distance_to_polygon = list(distance_to_polygon_boundary)

    # Set distance to zero for any points inside a polygon (leave other points unchanged).
    res2 = points_in_polygons.find_polygons_using_points_spatial_tree(points,
                                                                     spatial_tree_of_uniform_recon_points,
                                                                     rpolygons)
    for point_index, rpolygon in enumerate(res2):
        # If not inside any polygons then result will be None.
        if rpolygon:
            distance_to_polygon[point_index] = 0.0

    distance_to_polygon = np.array(distance_to_polygon)

    return distance_to_polygon,distance_to_polygon_boundary


def get_merged_cob_terrane_polygons(COBterrane_file, rotation_model, reconstruction_time,
                                    sampling, area_threshold=None, return_raster=False):
    """Reconstruct COB terrane polygons and merge them into contour polygons, optionally filtered by area."""

    polygon_features = pygplates.FeatureCollection(COBterrane_file)

    cobter = force_polygon_geometries(polygon_features)

    cf = merge_polygons(cobter, rotation_model, reconstruction_time=reconstruction_time, sampling=sampling)
    
    if area_threshold is not None:
        sieve_polygons = polygon_area_threshold(cf, area_threshold)
        return sieve_polygons

    else:
        return cf

def get_merged_cob_terrane_raster(COBterrane_file, rotation_model, reconstruction_time,
                                  sampling, method='pygplates', anchor_plate_id=0):
    """Rasterize reconstructed COB terrane polygons to a binary mask array.

    :param COBterrane_file: Path to a GPlates-compatible COB terrane polygon file.
    :param rotation_model: pygplates RotationModel.
    :param reconstruction_time: Age in Ma.
    :param sampling: Grid spacing in degrees.
    :param method: 'pygplates' (default) uses contour-based merging; 'rasterio' burns polygons directly.
    :param anchor_plate_id: Plate ID used as the fixed reference frame (default 0).
    :returns: 2-D numpy array (1 = inside terrane, 0 = outside).
    """

    if method == 'pygplates':
        polygon_features = pygplates.FeatureCollection(COBterrane_file)

        cobter = force_polygon_geometries(polygon_features)

        mask = merge_polygons(cobter, rotation_model, reconstruction_time=reconstruction_time,
                                sampling=sampling, return_raster=True, anchor_plate_id=anchor_plate_id)

    elif method=='rasterio':
        import tempfile
        import geopandas as gpd
        from rasterio.features import rasterize, Affine

        polygon_features = pygplates.FeatureCollection(COBterrane_file)

        polygon_features = force_polygon_geometries(polygon_features)

        with tempfile.TemporaryDirectory() as temporary_directory:
            pygplates.reconstruct(polygon_features, rotation_model, '{:s}/masking_temp.shp'.format(temporary_directory), reconstruction_time, anchor_plate_id=anchor_plate_id)

            gdf = gpd.read_file('{:s}/masking_temp.shp'.format(temporary_directory))

        dims = (int(180./sampling)+1, int(360./sampling)+1)
        transform = Affine(sampling, 0.0, -180.-sampling/2., 0.0, sampling, -90.-sampling/2.)

        geometry_zval_tuples = [(x.geometry, 1) for i, x in gdf.iterrows()]

        mask = rasterize(
            geometry_zval_tuples,
            transform=transform,
            out_shape=dims)

        # the first and last columns should match, but may not due to the imposed dateline
        mask[:,0] = mask[:,-1]

    return mask


def get_dynamic_polygon_raster(dynamic_polygons, rotation_model, reconstruction_time,
                               sampling):
    """Resolve topological plate polygons at a given time and rasterize them by plate ID.

    :param dynamic_polygons: List of pygplates FeatureCollection objects with topological polygon features.
    :param rotation_model: pygplates RotationModel.
    :param reconstruction_time: Age in Ma.
    :param sampling: Grid spacing in degrees.
    :returns: 2-D numpy array where each cell contains the plate ID of the enclosing polygon (0 = unassigned).
    """

    import tempfile
    import geopandas as gpd
    from rasterio.features import rasterize, Affine

    with tempfile.TemporaryDirectory() as temporary_directory:
        pygplates.resolve_topologies(dynamic_polygons, 
                                     rotation_model, 
                                     '{:s}/masking_temp.shp'.format(temporary_directory), 
                                     reconstruction_time)

        gdf = gpd.read_file('{:s}/masking_temp.shp'.format(temporary_directory))

    dims = (int(180./sampling)+1, int(360./sampling)+1)
    transform = Affine(sampling, 0.0, -180.-sampling/2., 0.0, sampling, -90.-sampling/2.)

    geometry_zval_tuples = [(x.geometry, x.PLATEID1) for i, x in gdf.iterrows()]
    
    #with rasterio.open(raster_file) as src:
        # iterate over features to get (geometry, id value) pairs
    mask = rasterize(
        geometry_zval_tuples,
        transform=transform,
        out_shape=dims)

    # the first and last columns should match, but may not due to the imposed dateline
    mask[:,0] = mask[:,-1]

    return mask



def plate_boundary_intersections(cross_section_geometry, shared_boundary_sections, ProfileX_kms):
    """Find subduction, ridge, and other plate boundary crossings along a great-circle cross-section.

    :param cross_section_geometry: pygplates PolylineOnSphere defining the cross-section path.
    :param shared_boundary_sections: List of resolved topological boundary sections from pygplates.resolve_topologies.
    :param ProfileX_kms: Cumulative distance array (km) along the cross-section, one value per segment start.
    :returns: Tuple (subduction_intersections, ridge_intersections, other_intersections); each is a list of
        [PointOnSphere, distance_along_profile_km, polarity] for subduction, or [PointOnSphere, distance_km] for others.
    """
    subduction_intersections = []
    ridge_intersections = []
    other_intersections = []

    for shared_boundary_section in shared_boundary_sections:

        for shared_subsegment in shared_boundary_section.get_shared_sub_segments():

            (min_distance_to_feature,
            closest_point_on_section,
            closest_point_on_topology,
            section_index,
            topology_index) = pygplates.GeometryOnSphere.distance(
                cross_section_geometry,
                shared_subsegment.get_resolved_geometry(),
                return_closest_positions=True,
                return_closest_indices=True)

            if min_distance_to_feature == 0:

                # find the distance along the section profile
                #print section_index
                cross_section_segment = cross_section_geometry.get_segments()[section_index]
                distance_along_segment = pygplates.GeometryOnSphere.distance(closest_point_on_section,cross_section_segment.get_start_point())
                distance_long_profile = distance_along_segment*pygplates.Earth.mean_radius_in_kms + ProfileX_kms[section_index]

                if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('SubductionZone'):
                    polarity = get_subduction_polarity(shared_subsegment,topology_index,cross_section_segment,distance_along_segment)
                    subduction_intersections.append([closest_point_on_section,distance_long_profile,polarity])

                elif shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('MidOceanRidge'):
                    ridge_intersections.append([closest_point_on_section,distance_long_profile])
                else:
                    other_intersections.append([closest_point_on_section,distance_long_profile])

    return subduction_intersections,ridge_intersections,other_intersections


def get_subduction_polarity(shared_subsegment, topology_index, cross_section_segment, distance_along_segment):
    """Return the cross-section dip direction (True = left-to-right) for a subduction segment intersection."""
    topology_section_segment = shared_subsegment.get_resolved_geometry().get_segments()[topology_index]

    # Get cross-section segment direction at point of intersection.
    if not cross_section_segment.is_zero_length():
        cross_section_segment_direction = cross_section_segment.get_arc_direction(
                distance_along_segment / cross_section_segment.get_arc_length())
    else:
        cross_section_segment_direction = cross_section_segment.get_arc_direction(0)

    # Imagine topology going South to North. If cross section crosses it going from West to East then
    # cross section segment and normal to topology segment will point 'away' from each other (ie, < 0).
    # Note that topology segment normal points to the left side (West).
    # Note: Same sort of thing applies when imagining topology going North to South (hence reference to left/right instead).
    cross_section_left_to_right = pygplates.Vector3D.dot(
            cross_section_segment_direction,
            topology_section_segment.get_great_circle_normal()) < 0

    # If overriding plate on left (West) and cross-section goes left to right (West to East) then
    # subduction dips from left to right along cross-section (same for overriding on right and cross-section going right to left).
    # The other two cases (out of four cases total) result in subduction dipping right to left along cross-section.
    subduction_polarity = shared_subsegment.get_feature().get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)
    if ((subduction_polarity == 'Left' and cross_section_left_to_right) or
        (subduction_polarity == 'Right' and not cross_section_left_to_right)):
        cross_section_dips_left_to_right = True
    else:
        # NOTE: We'll also get here if (subduction_polarity == 'Unknown').
        cross_section_dips_left_to_right = False

    return cross_section_dips_left_to_right
    

def polyline_zonal_lengths(gdf, binsize=10):
    '''
    Given a geodataframe containing polylines, divides the lines into
    latitude bands of a user-defined width and returns the length of the lines
    within each band 
    '''

    bin_lengths = []
    
    for bin_edge in np.arange(-90.,90.,binsize):
        
        polygon = Polygon([(-180, bin_edge), (-180, bin_edge+binsize), (180, bin_edge+binsize), (180, bin_edge), (-180, bin_edge)])
        poly_gdf = gpd.GeoDataFrame([1], geometry=[polygon], crs=4326)
        
        # Fix geometries that are invalid (e.g. near the poles)
        # https://gis.stackexchange.com/questions/430384/using-shapely-methods-explain-validity-and-make-valid-on-shapefile
        gdf.geometry = gdf.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)

        poly_clip = gpd.clip(gdf, poly_gdf)
        
        if poly_clip.empty:
            bin_length = 0
        else:
            bin_length = 0
            for i,feature in poly_clip.explode().iterrows():
                if len(feature.geometry.xy[0])>1:
                    bin_length += pygplates.PolylineOnSphere(
                        [(lat,lon) for lat,lon in zip(feature.geometry.xy[1], 
                                                      feature.geometry.xy[0])]).get_arc_length()
                    
        bin_lengths.append(bin_length * pygplates.Earth.mean_radius_in_kms)
        
    return bin_lengths


def polygon_zonal_areas(gdf, binsize=10, method='polygon', raster_sampling=1):
    '''
    Given a geodataframe containing polygons, divides the polygons into
    latitude bands of a user-defined width and returns the area of the polygon
    within each band 
    '''
    if method=='polygon':
        bin_areas = []
        
        for bin_edge in np.arange(-90.,90.,binsize):
            
            polygon = Polygon([(-180, bin_edge), (-180, bin_edge+binsize), (180, bin_edge+binsize), (180, bin_edge), (-180, bin_edge)])
            poly_gdf = gpd.GeoDataFrame([1], geometry=[polygon], crs=4326)
            
            # Fix geometries that are invalid (e.g. near the poles)
            # https://gis.stackexchange.com/questions/430384/using-shapely-methods-explain-validity-and-make-valid-on-shapefile
            gdf.geometry = gdf.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)

            poly_clip = gpd.clip(gdf, poly_gdf)
            
            if poly_clip.empty:
                bin_area = 0
            else:
                bin_area = 0
                tmp = gdf2gpml(poly_clip)
                for f in tmp:
                    if f.get_geometry():
                        bin_area += f.get_geometry().get_area()
                        
            bin_areas.append(bin_area * pygplates.Earth.mean_radius_in_kms**2)
        
        return bin_areas

    elif method=='rasterize':
        import xarray as xr
        from rasterio.features import rasterize, Affine

        #TODO change the raster definition to be pixel-registered
        dims = (int(180./raster_sampling)+1, int(360./raster_sampling)+1)
        transform = Affine(raster_sampling, 0.0, -180.-raster_sampling/2., 0.0, raster_sampling, -90.-raster_sampling/2.)

        geometry_zval_tuples = [(x.geometry, 1) for i, x in gdf.iterrows()]

        #with rasterio.open(raster_file) as src:
            # iterate over features to get (geometry, id value) pairs
        mask = rasterize(
            geometry_zval_tuples,
            transform=transform,
            out_shape=dims)

        # the first and last columns should match, but may not due to the imposed dateline
        mask[:,0] = mask[:,-1]

        lats = np.arange(-90,90+raster_sampling,raster_sampling)

        bin_areas = raster_zonal_areas(mask, lats, binsize)

        return bin_areas
    
    else:
        raise ValueError('Unknown value {} for method parameter'.format(method))


def raster_zonal_areas(grd, lats, binsize):
    """Sum grid-cell areas within latitude bins from a binary mask, accounting for spherical area weighting.

    :param grd: 2-D numpy binary array (1 = region of interest, 0 = elsewhere).
    :param lats: 1-D array of latitude values corresponding to grd rows (uniform spacing assumed).
    :param binsize: Width of each latitude bin in degrees.
    :returns: List of areas in km² for each latitude bin from -90 to 90.
    """
    raster_sampling = np.abs(lats[1]-lats[0])

    area_weights = pygplates.Earth.mean_radius_in_kms**2 * np.sin(np.radians(90-lats)) * np.radians(raster_sampling)**2

    weighted_areas = grd.sum(axis=1) * area_weights

    bin_areas = []
    for bin_edge in np.arange(-90.,90.,binsize):
        ind = np.logical_and(lats>=bin_edge, lats<bin_edge+binsize)
        bin_areas.append(np.sum(weighted_areas[ind]))

    return bin_areas


def topology_lookup(reconstruction_model,
                    reconstruction_times=np.arange(0,1001,10),
                    boundary_types=['subduction']):
    """Build a time-keyed dict of reconstructed plate boundary features for fast per-time queries.

    :param reconstruction_model: ReconstructionModel instance with dynamic polygons loaded.
    :param reconstruction_times: Iterable of ages in Ma to pre-compute (default 0–1000 Ma in 10 Ma steps).
    :param boundary_types: List of boundary types to include; valid values are 'subduction', 'midoceanridge', 'other'.
    :returns: Dict mapping each reconstruction time (float, Ma) to a list of resolved boundary features.
    """
    lookup_dict = {}
    for reconstruction_time in reconstruction_times:
        snapshot = reconstruction_model.plate_snapshot(reconstruction_time=reconstruction_time)
        lookup_dict[reconstruction_time] = snapshot.get_boundary_features(boundary_types=boundary_types)

    return lookup_dict





