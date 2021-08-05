import pygplates
import litho1pt0 as litho


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





