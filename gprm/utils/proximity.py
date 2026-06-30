"""Distance and proximity rasters from polygon, polyline, and point features."""
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

import warnings
try:
    from xrspatial import viewshed, proximity
except Exception:
    warnings.warn('gprm.utils.proximity functions based on xrspatial not available')
from rasterio.features import rasterize, Affine
from .spatial import get_merged_cob_terrane_raster


def mask_to_da(mask, sampling=1):
    """Convert a 2D numpy mask to a global lat/lon xarray DataArray.

    :param mask: 2-D numpy array with shape (n_lats, n_lons).
    :param sampling: Grid spacing in degrees used to build the coordinate axes (default 1).
    :returns: xarray DataArray with 'x' (longitude) and 'y' (latitude) coordinates.
    """

    # the first and last columns should match, but may not due to the imposed dateline
    mask[:,0] = mask[:,-1]
    
    coords = [('y',np.arange(-90,90+sampling,sampling)), 
              ('x',np.arange(-180,180+sampling,sampling))]

    return xr.DataArray(mask,
                        coords=coords,
                        name='z')


def rasterize_polygons(gdf, sampling=1, region=[-180, 180, -90, 90], zval_field=None):
    """Rasterize geodataframe polygons to a DataArray, optionally using a per-polygon attribute as the cell value.

    :param gdf: GeoDataFrame containing polygon geometries.
    :param sampling: Grid spacing in degrees (default 1).
    :param region: Bounding box [xmin, xmax, ymin, ymax] (default global).
    :param zval_field: Column name to use as the raster value; if None, all polygons are burned as 1.
    :returns: xarray DataArray.
    """
    
    dims = (int((region[3]-region[2])/sampling)+1, 
            int((region[1]-region[0])/sampling)+1)
    transform = Affine(sampling, 0.0, region[0]-sampling/2., 0.0, sampling, region[2]-sampling/2.)

    if zval_field is not None:
        geometry_zval_tuples = [(x.geometry, x[zval_field]) for i, x in gdf.iterrows()]
    else:
        geometry_zval_tuples = [(x.geometry, 1) for i, x in gdf.iterrows()]

    #with rasterio.open(raster_file) as src:
        # iterate over features to get (geometry, id value) pairs
    mask = rasterize(
        geometry_zval_tuples,
        transform=transform,
        out_shape=dims)

    return mask_to_da(mask, sampling=sampling)


def reconstruct_and_rasterize_polygons(features, rotation_model, reconstruction_time, sampling=1, anchor_plate_id=0):
    """Reconstruct polygon features to a specified time and return a rasterized DataArray.

    :param features: pygplates FeatureCollection or path to a GPlates-compatible polygon file.
    :param rotation_model: pygplates RotationModel.
    :param reconstruction_time: Age in Ma.
    :param sampling: Grid spacing in degrees (default 1).
    :param anchor_plate_id: Plate ID used as the fixed reference frame (default 0).
    :returns: xarray DataArray with 1 inside reconstructed polygons and 0 outside.
    """

    mask = get_merged_cob_terrane_raster(features, rotation_model, reconstruction_time,
                                         sampling=sampling, method='rasterio',
                                         anchor_plate_id=anchor_plate_id)

    return mask_to_da(mask, sampling=sampling)



def polygons_buffer(gdf, sampling=1, region=[-180, 180, -90, 90], inside=False):
    """Return a great-circle distance-to-polygon-boundary raster from geodataframe polygons.

    :param gdf: GeoDataFrame containing polygon geometries.
    :param sampling: Grid spacing in degrees (default 1).
    :param region: Bounding box [xmin, xmax, ymin, ymax] (default global).
    :param inside: If False (default), measure distance from outside; if True, from inside. See boundary_proximity for 'both'/'boundary' options.
    :returns: xarray DataArray of great-circle distances in metres.
    """

    ds = rasterize_polygons(gdf, sampling=sampling, region=region)
    
    return boundary_proximity(ds, inside=inside)


def raster_buffer(ds, clipval=0, inside=False):
    """Return a distance-to-boundary raster from a continuous DataArray thresholded at clipval.

    :param ds: xarray DataArray (continuous values).
    :param clipval: Threshold: cells >= clipval are treated as 'inside' (default 0).
    :param inside: Distance direction; see boundary_proximity for accepted values.
    :returns: xarray DataArray of great-circle distances.
    """
    
    ds_binary = ds.where(ds>=clipval, other=0)
    ds_binary = ds_binary.where(ds_binary<=0, other=1)
    
    return boundary_proximity(ds_binary, inside=inside)

    
def handle_da_coordinates(da):
    """Normalise DataArray coordinate names to 'x' and 'y' for xrspatial compatibility."""
    
    coord_keys = [key for key in da.coords.keys()]  # updated for python3 compatibility

    if 'x' in coord_keys:
        return da
    else:
        if 'lon' in coord_keys[0].lower():
            latitude_key=1; longitude_key=0   
        else:
            latitude_key=0; longitude_key=1

        da = da.rename({coord_keys[longitude_key]:'x', 
                        coord_keys[latitude_key]:'y'})
        
        return da

    
def boundary_proximity(da, inside=False):
    """Compute great-circle distance to polygon boundaries from a binary DataArray.

    :param da: Binary xarray DataArray (1 = inside polygon, 0 = outside); coords must be 'x' and 'y'.
    :param inside: Controls what distance is returned: False = distance from outside to polygon edge;
        True = distance from inside to polygon edge; 'both' = returns (outside_dist, inside_dist) tuple;
        'boundary' = sum of outside and inside distances (distance to the boundary from either side).
    :returns: xarray DataArray of distances, or a tuple of two DataArrays if inside='both'.
    """
    
    da = handle_da_coordinates(da)

    if inside:
        return proximity(da, target_values=[0], distance_metric='GREAT_CIRCLE')
    elif not inside:
        return proximity(da, target_values=[1], distance_metric='GREAT_CIRCLE')
    elif inside in ['both', 'boundary']:
        prox_outside = proximity(da, target_values=[1], distance_metric='GREAT_CIRCLE')
        prox_inside = proximity(da, target_values=[0], distance_metric='GREAT_CIRCLE')
        if inside=='both':
            return prox_outside, prox_inside
        elif inside=='boundary':
            return prox_outside + prox_inside


def contour_proximity(da, target_value=0, inside='boundary'):
    """Compute distance to a contour level in a continuous raster by thresholding at target_value.

    :param da: Continuous xarray DataArray.
    :param target_value: Contour level; cells below this value become the 'inside' region (default 0).
    :param inside: Distance mode passed to boundary_proximity (default 'boundary').
    :returns: xarray DataArray of distances.
    """
    da2 = da.copy(deep=True)
    da2.data[da.data>=target_value] = 0
    da2.data[da.data<target_value] = 1
    return boundary_proximity(da2, inside=inside)



def points_proximity(x, y, spacing=1, region=[-180, 180, -90, 90]):
    """Compute a great-circle distance-to-nearest-point raster from arrays of point coordinates.

    :param x: Longitude values of the input points.
    :param y: Latitude values of the input points.
    :param spacing: Output grid spacing in degrees (default 1).
    :param region: Bounding box [xmin, xmax, ymin, ymax] (default global).
    :returns: xarray DataArray of great-circle distances to the nearest input point.
    """
    from datashader import Canvas

    df = pd.DataFrame({"x": x, "y": y,})

    dims = (int((region[3]-region[2])/spacing)+1, int((region[1]-region[0])/spacing)+1)

    # Note the creation of a canvas slightly larger than what may appear needed to force
    # the grid to conform to the desired gridline registered coordinates - but there
    # may be some issue here with pixel versus gridline concepts??
    cvs = Canvas(plot_width=dims[1], plot_height=dims[0], 
                 x_range=(region[0]-spacing/2., region[1]+spacing/2.), 
                 y_range=(region[2]-spacing/2., region[3]+spacing/2.))
    #cvs = Canvas(plot_width=dims[1], plot_height=dims[0], 
    #             x_range=(-180, 180.), 
    #             y_range=(-90, 90))

    points_agg = cvs.points(df, x="x", y="y")
    points_agg.data[~np.isfinite(points_agg.data)] = 0

    target_proximity_agg = proximity(
        points_agg, distance_metric="GREAT_CIRCLE"
    )
    return target_proximity_agg


def polyline_proximity(features, spacing=1, region=[-180, 180, -90, 90]):
    """Compute a great-circle distance-to-polyline raster from GPlates features or a GeoDataFrame.

    :param features: pygplates FeatureCollection (with polyline geometries) or a GeoDataFrame of LineStrings.
    :param spacing: Output grid spacing in degrees; polylines are tessellated to spacing/5 (default 1).
    :param region: Bounding box [xmin, xmax, ymin, ymax] (default global).
    :returns: xarray DataArray of great-circle distances to the nearest polyline.
    """

    import pygplates
    tesselation_spacing=spacing/5

    features_as_points = []
    if isinstance(features, gpd.GeoDataFrame):
        for i,feature in features.explode(index_parts=False).iterrows():
            geometry = pygplates.PolylineOnSphere(zip(feature.geometry.coords.xy[1], feature.geometry.coords.xy[0]))
            features_as_points.extend(geometry.to_tessellated(np.radians(tesselation_spacing)).to_lat_lon_list())
    else:
        for f in features:
            if f.get_geometry():
                features_as_points.extend(f.get_geometry().to_tessellated(np.radians(tesselation_spacing)).to_lat_lon_list())

    return points_proximity(x=[lon for lat,lon in features_as_points],
                            y=[lat for lat,lon in features_as_points],
                            spacing=spacing, 
                            region=region)

    
def generate_shadows(da, x, y, observer_elev):
    """Compute a viewshed raster from an observer point (x, y) and elevation grid."""
    da = handle_da_coordinates(da)
    
    #gridc = pygmt.grdclip(grid, below=[-100,-100])
    #topo0 = topo_low.rename({'lon':'x', 'lat':'y'})
    #topo0.data[topo0.data<0] = 0
    res = viewshed(da, x=x, y=y, observer_elev=observer_elev)
    
    