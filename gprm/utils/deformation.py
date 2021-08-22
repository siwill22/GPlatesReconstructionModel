import pygplates
import numpy as np
import geopandas as gpd
import pygmt
from .create_gpml import geometries_to_geodataframe, geodataframe_to_geometries
from shapely.geometry import LineString, Polygon
# import litho1pt0 as litho


def find_deforming_mesh_geometry(polygons, name):
    for polygon in polygons:
        if polygon.get_feature().get_name() == name:
            return polygon.get_resolved_geometry()


def points_within_mesh_geometry(pts, geom):
    pts_in_poly = []
    for pt in pts.multipoint.get_points():
        if geom.is_point_in_polygon(pt):
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

    
def create_graticule(graticule_spacing=1, xlims=[-180,180], ylims=[-90,90], along_line_point_density=0.1, polygon=None):
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


def get_crustal_thickness_points(points, grid=None, top_name='CRUST1-TOP', bottom_name='CRUST3-BOTTOM'):
    
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


def topological_reconstruction(topological_model, points, reconstruction_time, 
                               initial_time=0, final_time=None, initial_scalars=None,
                               return_inactive_points=True, deactivate_points=True):

    if not final_time:
        final_time = reconstruction_time

    # If deactivate points is a boolean, we use it with the default thresholds
    if deactivate_points is not None:
        deactivate_points = pygplates.ReconstructedGeometryTimeSpan.DefaultDeactivatePoints(
            threshold_velocity_delta=0.7, 
            threshold_distance_to_boundary=10,
            deactivate_points_that_fall_outside_a_network = deactivate_points)

        
    # TODO determine default behaviour for optional arguments
    time_spans = topological_model.reconstruct_geometry(
        points,
        initial_time=initial_time,
        oldest_time=final_time,
        youngest_time=initial_time,
        initial_scalars=initial_scalars,
        deactivate_points = deactivate_points)

    reconstructed_points = time_spans.get_geometry_points(reconstruction_time, return_inactive_points=return_inactive_points)

    #TODO iterate over scalar values and get reconstructed value

    #TODO for cases where this could lead to an array of inconsistent length - maybe should allow points to be 'None'??
    valid_index = [reconstructed_point is not None for reconstructed_point in reconstructed_points]
    pts = list(zip(*[reconstructed_point.to_lat_lon() for reconstructed_point in reconstructed_points if reconstructed_point is not None]))

    return pts, valid_index


def geodataframe_topological_reconstruction(gdf, topological_model, 
                                            reconstruction_time, initial_time=0, final_time=None,
                                            return_inactive_points=True,
                                            deactivate_points_that_fall_outside_a_network=True):

    # Given a geodataframe, will reconstruct using a topological model to a given reconstruction time    
    if not final_time:
        final_time = reconstruction_time
        
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
                                                   initial_time, final_time, 
                                                   return_inactive_points=return_inactive_points,
                                                   deactivate_points_that_fall_outside_a_network=deactivate_points_that_fall_outside_a_network)
        
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
                                                       initial_time, final_time)
            
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



#def feature_collection_topological_reconstruction():
#
#    
#    return


