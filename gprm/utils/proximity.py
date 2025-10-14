import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
#import xrspatial as xrs

import warnings
try:
    from xrspatial import viewshed, proximity
except:
    warnings.warn('gprm.utils.proximity functions based on xrspatial not available')
#import pygmt
from rasterio.features import rasterize, Affine
from .spatial import get_merged_cob_terrane_raster


def mask_to_da(mask, sampling=1):
    # given a numpy array, create a xr.dataarray assuming that the extent
    # is global in lat/long 

    # the first and last columns should match, but may not due to the imposed dateline
    mask[:,0] = mask[:,-1]
    
    coords = [('y',np.arange(-90,90+sampling,sampling)), 
              ('x',np.arange(-180,180+sampling,sampling))]

    return xr.DataArray(mask,
                        coords=coords,
                        name='z')


def rasterize_polygons(gdf, sampling=1, region=[-180, 180, -90, 90], zval_field=None):
    # given a geodataframe with some polygons, returns a raterized version
    # TODO add region option
    
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


def reconstruct_and_rasterize_polygons(features, rotation_model, reconstruction_time, sampling=1):
    # given a set of reconstructable polygon features, together with a rotation model and 
    # reconstruction time, returns a raster that is a rasterized version of the reconstructed 
    # polygon geometries 

    mask = get_merged_cob_terrane_raster(features, rotation_model, reconstruction_time,
                                         sampling=sampling, method='rasterio')

    return mask_to_da(mask, sampling=sampling)



def polygons_buffer(gdf, sampling=1, region=[-180, 180, -90, 90], inside=False):
    # given a geodataframe containing a set of polygons,
    # computes a raster where each mode is the distance to the polygon boundaries
    # Options are to compute distance to inside, outside, or edge

    # TODO rename to polygon_proximity??

    ds = rasterize_polygons(gdf, sampling=sampling, region=region)
    
    return boundary_proximity(ds, inside=inside)


def raster_buffer(ds, clipval=0, inside=False):
    # DUPLICATING contour_proximity??
    
    ds_binary = ds.where(ds>=clipval, other=0)
    ds_binary = ds_binary.where(ds_binary<=0, other=1)
    
    return boundary_proximity(ds_binary, inside=inside)

    
def handle_da_coordinates(da):
    # utility function to ensure the geographic coordinates from a xarray dataarray
    # work correctly when passed to xrspatial functions
    
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
    # given a dataarray assumed to contain ones and zeros, returns a raster 
    # of the same dimensions where each grid node contains the distance to some target values 
    # options are to compute 'inside', 'outside', both inside and outside (returning both distance 
    # arrays separately) or distance to the edge from both inside and outside
    
    da = handle_da_coordinates(da)

    if inside==True:
        return proximity(da, target_values=[0], distance_metric='GREAT_CIRCLE')
    elif inside==False:
        return proximity(da, target_values=[1], distance_metric='GREAT_CIRCLE')
    elif inside in ['both', 'boundary']:
        prox_outside = proximity(da, target_values=[1], distance_metric='GREAT_CIRCLE')
        prox_inside = proximity(da, target_values=[0], distance_metric='GREAT_CIRCLE')
        if inside=='both':
            return prox_outside, prox_inside
        elif inside=='boundary':
            return prox_outside + prox_inside


def contour_proximity(da, target_value=0, inside='boundary'):

    da2 = da.copy(deep=True)
    da2.data[da.data>=target_value] = 0
    da2.data[da.data<target_value] = 1
    return boundary_proximity(da2, inside=inside)



def points_proximity(x, y, spacing=1, region=[-180, 180, -90, 90]):
    '''
    Given a set of points, return a dataArray of distances to the nearest
    point. 
    Should work reasonable for lines given the points along the lines are tessellated
    to a close spacing, and avoids dateline issues
    Simplified from the example here:
    https://xarray-spatial.org/user_guide/proximity.html
    '''
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
    # compute raster of distances to a set of polylines

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
    
    da = handle_da_coordinates(da)
    
    #gridc = pygmt.grdclip(grid, below=[-100,-100])
    #topo0 = topo_low.rename({'lon':'x', 'lat':'y'})
    #topo0.data[topo0.data<0] = 0
    res = viewshed(da, x=x, y=y, observer_elev=observer_elev)
    
    