import pygplates
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import pygmt
from .create_gpml import geometries_to_geodataframe, geodataframe_to_geometries
from shapely.geometry import LineString, Polygon
from .raster import xyz2grd
# import litho1pt0 as litho


DEFAULT_COLLISION_PARAMETERS = (0.7, 10)
DEFAULT_RECONSTRUCTION_TIME_INCREMENT = 1.
DEFAULT_DEACTIVATE_POINTS = True
DEFAULT_RETURN_INACTIVE_POINTS = True


def find_deforming_mesh_geometry(polygons, names):
    """
    Given a list of polygons and a name or list of names,
    returns a list of polygons with names matching the list
    """

    if isinstance(names, str):
        names = [names]

    selected_polygons = []
    for name in names:
        for polygon in polygons:
            if polygon.get_feature().get_name() == name:
                selected_polygons.append(polygon.get_resolved_geometry())

    return selected_polygons


def points_within_mesh_geometry(pts, mesh_geoms):
    """
    Given a list of point geometries and a list of polygon geometries,
    returns a list of point geometries within the polygons
    """

    pts_in_poly = []
    for mesh_geom in mesh_geoms:
        for pt in pts.multipoint.get_points():
            if mesh_geom.is_point_in_polygon(pt):
                pts_in_poly.append(pt) 

    return pts_in_poly


def create_circles(centre_points, circle_radius_degrees=0.5, polygon=None):
    # create small circles centred on the provide points, with user-specified radii
    # optionally, limit the points to those which lie within a provided polygon geometry
    # requires pyshtools
    import pyshtools as pysh

    if polygon:
        pts_in_polygon = []
        for pt in centre_points.multipoint.get_points():
            if polygon.is_point_in_polygon(pt):
                pts_in_polygon.append(pt.to_lat_lon())
        #print(pts_in_polygon[0])
    
        circle_list = []
        for lat,lon in pts_in_polygon:
            circle_list.append(pysh.utils.MakeCircleCoord(lat,lon,circle_radius_degrees))
            
    else:
        circle_list = []
        for lon,lat in zip(centre_points.longitude, centre_points.latitude):
            circle_list.append(pysh.utils.MakeCircleCoord(lat,lon,circle_radius_degrees))

    return circle_list

    
def create_graticule(graticule_spacing=1, xlims=[-180,180], ylims=[-90,90], 
                     along_line_point_density=0.1, polygon=None):
    # create a latitude-longitude graticule, in the form of line geometries along each meridian and parallel 
    # optionally specify the graticule spacing, extent on the globe, sampling of the vertices along each line
    # optionally clip the graticule to a user-provided polygon geometry
    graticule_lines = []
    for xlevel in np.arange(xlims[0], xlims[1]+graticule_spacing, graticule_spacing):
        graticule_lines.append(np.array([(ylevel,xlevel) for ylevel in np.arange(ylims[0], ylims[1]+along_line_point_density, along_line_point_density)]))
        
    for ylevel in np.arange(ylims[0], ylims[1]+graticule_spacing, graticule_spacing):
        graticule_lines.append(np.array([(ylevel,xlevel) for xlevel in np.arange(xlims[0], xlims[1]+along_line_point_density, along_line_point_density)]))

    if polygon:
        graticule_lines_gdf = geometries_to_geodataframe(graticule_lines, geometry_type='PolyLine')
        polygon_gdf = geometries_to_geodataframe([polygon.to_lat_lon_array()])
        graticule_lines_in_polygon_gdf = gpd.overlay(graticule_lines_gdf, polygon_gdf, how='intersection', keep_geom_type=False)
        graticule_lines = geodataframe_to_geometries(graticule_lines_in_polygon_gdf)

    return graticule_lines


def get_crustal_thickness_points(points, grid=None, 
                                 top_name='CRUST1-TOP', 
                                 bottom_name='CRUST3-BOTTOM'):
    
    # if no grid is provided, we take the layer thickness from litho1.0
    if not grid:
        import litho1pt0 as litho
    
        ptlats = np.array([pt.to_lat_lon()[0] for pt in points])
        ptlons = np.array([pt.to_lat_lon()[1] for pt in points])
        top_depth = litho.layer_depth(ptlats, ptlons, top_name)
        bottom_depth = litho.layer_depth(ptlats, ptlons, bottom_name)
        layer_thickness = bottom_depth-top_depth

        return layer_thickness

    # if a grid is provided, sample using grdtrack
    else:

        ptlats = np.array([pt.to_lat_lon()[0] for pt in points])
        ptlons = np.array([pt.to_lat_lon()[1] for pt in points])
        layer_thickness = pygmt.grdtrack(points=np.vstack((ptlons, ptlats)), grid=grid, newcolname='z')
        
        return np.array(layer_thickness['z'])


def topological_reconstruction(topological_model, points, 
                               reconstruction_time, 
                               initial_time=0., 
                               oldest_time=None, 
                               youngest_time=0., 
                               time_increment=DEFAULT_RECONSTRUCTION_TIME_INCREMENT,
                               initial_scalars=None, 
                               return_inactive_points=DEFAULT_RETURN_INACTIVE_POINTS, 
                               deactivate_points=DEFAULT_DEACTIVATE_POINTS, 
                               collision_parameters=DEFAULT_COLLISION_PARAMETERS):

    if not oldest_time:
        oldest_time = reconstruction_time

    # If deactivate points is a boolean, we use it with the default thresholds
    if deactivate_points is not None:
        deactivate_points = pygplates.ReconstructedGeometryTimeSpan.DefaultDeactivatePoints(
            threshold_velocity_delta=collision_parameters[0], 
            threshold_distance_to_boundary=collision_parameters[1],
            deactivate_points_that_fall_outside_a_network = deactivate_points)

    # TODO determine default behaviour for optional arguments
    time_spans = topological_model.reconstruct_geometry(
        points,
        initial_time=initial_time,
        oldest_time=oldest_time,
        youngest_time=youngest_time,
        initial_scalars=initial_scalars,
        time_increment=time_increment,
        deactivate_points = deactivate_points)

    reconstructed_points = time_spans.get_geometry_points(reconstruction_time, return_inactive_points=return_inactive_points)

    #TODO iterate over scalar values and get reconstructed value

    #print(points,reconstructed_points)
    #TODO for cases where this could lead to an array of inconsistent length - maybe should allow points to be 'None'??
    valid_index = [reconstructed_point is not None for reconstructed_point in reconstructed_points]
    pts = list(zip(*[reconstructed_point.to_lat_lon() for reconstructed_point in reconstructed_points if reconstructed_point is not None]))

    return pts, valid_index


def geodataframe_topological_reconstruction(gdf, topological_model, 
                                            reconstruction_time, 
                                            initial_time=0, 
                                            oldest_time=None,
                                            youngest_time=0., 
                                            time_increment=DEFAULT_RECONSTRUCTION_TIME_INCREMENT,
                                            return_inactive_points=DEFAULT_RETURN_INACTIVE_POINTS, 
                                            deactivate_points=DEFAULT_DEACTIVATE_POINTS):

    # Given a geodataframe, will reconstruct using a topological model to a given reconstruction time   
    # TODO check if this is the default behaviour anyway??? 
    if not oldest_time:
        oldest_time = reconstruction_time
    
    # Preprocessing:
    gdf = gdf[gdf.geometry.is_valid]   # remove invalid geometries
    if all([item in gdf.columns for item in ['FROMAGE','TOAGE']]):
        gdf = gdf[(gdf.FROMAGE>=reconstruction_time) & (gdf.TOAGE<=reconstruction_time)]   # select points valid at reconstruction time
    gdf = gdf.explode()   # multipart to singlepart
    gdf.reset_index(inplace=True)   # reset index
    
     
    # It is quicker to reconstruct features with just one call to the topological model, but we need 
    # to check that this is possible for the input dataframe
    if len(gdf.geom_type.unique())==1 and gdf.geom_type.unique()[0]=='Point':  

        # Should work for points, but for polygon/polyline need to iterate over each feature 
        geometry_points = [(lat,lon) for lat,lon in zip(gdf.geometry.y,gdf.geometry.x)]

        (pts, 
         valid_index) = topological_reconstruction(topological_model, geometry_points, reconstruction_time, 
                                                   initial_time, oldest_time, youngest_time, time_increment,
                                                   return_inactive_points=return_inactive_points,
                                                   deactivate_points=deactivate_points)
        
        reconstructed_gdf = gdf.iloc[valid_index]
        reconstructed_gdf = reconstructed_gdf.set_geometry(gpd.points_from_xy(pts[1], pts[0]))

    # For all polylines and polygons, and points treated as individual features, we reconstruct by iterating over features
    else:
        
        reconstructed_gdf = gdf.copy()
        
        for i,feature in gdf.iterrows():
            
            # Point features not handled yet
            if feature.geometry.geom_type in ['LineString']:
                geometry_points = [(lat,lon) for lat,lon in zip(feature.geometry.xy[1], feature.geometry.xy[0])]
            elif feature.geometry.geom_type in ['Polygon']:
                # clearly this isn't handling interior rings
                geometry_points = [(lat,lon) for lat,lon in zip(feature.geometry.exterior.coords.xy[1], 
                                                                feature.geometry.exterior.coords.xy[0])]
            
            # TODO add option to tesselate features??

            (pts, 
             valid_index) = topological_reconstruction(topological_model, geometry_points, reconstruction_time, 
                                                       initial_time, oldest_time, youngest_time, time_increment,
                                                       return_inactive_points=return_inactive_points,
                                                       deactivate_points=deactivate_points)
            
            # TODO put something in here to deal with cases where the whole geometry has become invalid
            # 
            if feature.geometry.geom_type in ['LineString']:
                geom = LineString([tuple(coord) for coord in zip(pts[1], pts[0])])
            elif feature.geometry.geom_type in ['Polygon']:
                geom = Polygon([tuple(coord) for coord in zip(pts[1], pts[0])])
            #else:
            #    geom = Point(tuple(pts[1], pts[0]))
            reconstructed_gdf.loc[i, 'geometry'] = geom
            
    return reconstructed_gdf


def raster_topological_reconstruction(grid, topological_model, reconstruction_time, reverse=True,
                                      region=None, spacing=1.):
    """
    Generate a reconstructed raster based on a topological reconstruction
    
    By default, the operation will use the input raster to determine the extent and sampling
    of the output raster (Which would work for global grids).
    Optionally, a different region and sampling for the output grid can be specified.
    """
    
    if region:
        coords = [('lat',np.arange(region[2],region[3]+spacing, spacing)), ('lon',np.arange(region[0],region[1]+spacing, spacing))]
        XX,YY = np.meshgrid(np.arange(region[0],region[1]+spacing, spacing),
                            np.arange(region[2],region[3]+spacing, spacing))
    else:
        coords = [('lat',grid.lat.data), ('lon',grid.lon.data)]
        XX,YY = np.meshgrid(grid.lon.data, grid.lat.data)

    geometry_points = [(lat,lon) for lat,lon in zip(YY.flatten(),XX.flatten())]

    if reverse:
        pts, valid_index = topological_reconstruction(topological_model, geometry_points, 0., initial_time=reconstruction_time, oldest_time=reconstruction_time, 
                                                      return_inactive_points=True, deactivate_points=False)

        res = pygmt.grdtrack(grid=grid, 
                             points=pd.DataFrame(data={'x':pts[1],'y':pts[0]}), newcolname='z')
        x=XX.flatten()[valid_index] 
        y=YY.flatten()[valid_index]

        resg = xyz2grd(x,y,np.array(res['z']),XX,YY)
        
        reconstructed_raster = xr.DataArray(resg, coords=coords, name='z')
        
        return reconstructed_raster

    else:
        print('Forward reconstruction not yet implemented')
    # TODO implement the forward reconstruction case


#def feature_collection_topological_reconstruction():
#
#    
#    return


