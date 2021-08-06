import pygplates
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString, Polygon
#import litho1pt0 as litho


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


def geometries_to_geodataframe(geometries, geometry_type='polygon'):
    gdf = gpd.GeoDataFrame()
    gdf['geometry'] = None
    for i,geometry in enumerate(geometries):
        if geometry_type in ['PolyLine','Polyline']:
            poly = LineString([tuple(coord) for coord in np.fliplr(geometry)])
        else:
            poly = Polygon([tuple(coord) for coord in np.fliplr(geometry)])
        gdf.loc[i, 'geometry'] = poly

    return gdf


def geodataframe_to_geometries(gdf):
    geometry_list = []
    gdf = gdf.explode()
    for i,row in gdf.iterrows():
        geometry_list.append([(lat,lon) for lat,lon in zip(row.geometry.xy[1], row.geometry.xy[0])])
    return geometry_list


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
        graticule_lines = dataframe_to_geometries(graticule_lines_in_polygon_gdf)

    return graticule_lines