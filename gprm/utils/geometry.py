import pygplates
import numpy as np
from shapely.geometry import Point, LineString, Polygon
import geopandas as _gpd
import sys

def apply_reconstruction(feature, rotation_model, 
                         reconstruction_time_field='reconstruction_time',
                         reconstruction_plate_id_field='PLATEID1',
                         anchor_plate_id=0, reverse=False):
    '''
    function that can be used within the 'apply' method of
    a geodataframe to return a reconstructed geometry 
    '''

    rotation_pole = rotation_model.get_rotation(
                        feature[reconstruction_time_field],
                        feature[reconstruction_plate_id_field],
                        anchor_plate_id=anchor_plate_id)

    if reverse:
        rotation_pole = rotation_pole.get_inverse()
    
    if feature.geometry.geom_type=='Point':
        rp = rotation_pole * pygplates.PointOnSphere(feature.geometry.y, feature.geometry.x)
        return Point(rp.to_lat_lon()[::-1])

    elif feature.geometry.geom_type in ['LineString']:
        rp = rotation_pole * pygplates.PolylineOnSphere([(lat,lon) for lat,lon in zip(feature.geometry.xy[1], 
                                                                                      feature.geometry.xy[0])])
        return LineString([tuple(point.to_lat_lon()[::-1]) for point in rp.get_points()])
    elif feature.geometry.geom_type in ['Polygon']:
        rp = rotation_pole * pygplates.PolygonOnSphere([(lat,lon) for lat,lon in zip(feature.geometry.exterior.coords.xy[1], 
                                                                                     feature.geometry.exterior.coords.xy[0])])
        return Polygon([tuple(point.to_lat_lon()[::-1]) for point in rp.get_points()])


def apply_nearest_feature(point, lookup_dict, geometry_field='geometry', age_field='age'):
    # function to apply the nearest feature function assuming we have points in geodataframe
    # that also contains an 'age' field, and an existing lookup
    # dictionary of reconstructed features
    
    d = nearest_feature(pygplates.PointOnSphere(point[geometry_field].y, point[geometry_field].x), 
                           lookup_dict[point[age_field]])
    if d is None:
        return np.nan
    else:
        return d*pygplates.Earth.mean_radius_in_kms


def nearest_feature(point, features, return_nearest_feature=False):
    # The minimum distance to all features and the nearest feature.
    min_distance_to_all_features = None
    nearest_feature_ = None

    for feature in features:
        for geometry in feature.get_geometries():

            # Get the minimum distance from point to the current reconstructed geometry.
            min_distance_to_feature = pygplates.GeometryOnSphere.distance(
                    point,
                    geometry,
                    min_distance_to_all_features)

            # If the current geometry is nearer than all previous geometries then
            # its associated feature is the nearest feature so far.
            if min_distance_to_feature is not None:
                min_distance_to_all_features = min_distance_to_feature
                nearest_feature_ = feature

    if return_nearest_feature:
        return min_distance_to_all_features, nearest_feature_
    else:
        return min_distance_to_all_features


def distance_between_reconstructed_points_and_features(reconstructed_point_features,features):
    reconstructed_lat = []
    reconstructed_lon = []
    distances = []
    for point in reconstructed_point_features:
        reconstructed_lat.append(point.get_reconstructed_geometry().to_lat_lon()[0])
        reconstructed_lon.append(point.get_reconstructed_geometry().to_lat_lon()[1])

        dist = nearest_feature(point.get_reconstructed_geometry(),
                               features)
        distances.append(dist*pygplates.Earth.mean_radius_in_kms)

    return reconstructed_lon, reconstructed_lat, distances


def wrap_polyline_feature(polyline_feature, date_line_wrapper=None):
    
    if not date_line_wrapper:
        date_line_wrapper = pygplates.DateLineWrapper(0.0)

    polyline = pygplates.PolylineOnSphere(
        [(lat,lon) for lat,lon in zip(polyline_feature.geometry.xy[1], 
                                      polyline_feature.geometry.xy[0])])
    wrapped_polyline = date_line_wrapper.wrap(polyline)
    return LineString([tuple(point.to_lat_lon()[::-1]) for point in wrapped_polyline[0].get_points()])


def wrap_polygon_feature(polygon_feature, date_line_wrapper=None):
    
    if not date_line_wrapper:
        date_line_wrapper = pygplates.DateLineWrapper(0.0)

    polygon = pygplates.PolygonOnSphere(
        [(lat,lon) for lat,lon in zip(polygon_feature.geometry.exterior.coords.xy[1], 
                                      polygon_feature.geometry.exterior.coords.xy[0])])
    wrapped_polygon = date_line_wrapper.wrap(polygon)
    if len(wrapped_polygon)>1:
        print("Warning: polygon was split by dateline wrapping")
    return Polygon([tuple(point.to_lat_lon()[::-1]) for point in wrapped_polygon[0].get_points()])


def wrap_polygon_features(polygon_features, date_line_wrapper=None):
    
    if not date_line_wrapper:
        date_line_wrapper = pygplates.DateLineWrapper(0.0)

    results = []
    for idx, polygon_feature in polygon_features.iterrows():

        polygon = pygplates.PolygonOnSphere(
            [(lat,lon) for lat,lon in zip(polygon_feature.geometry.exterior.coords.xy[1], 
                                        polygon_feature.geometry.exterior.coords.xy[0])])
        wrapped_polygon_list = date_line_wrapper.wrap(polygon)
        wrapped_polygon_features = []
        for row in wrapped_polygon_list:
            wrapped_polygon_feature = polygon_feature.copy()
            wrapped_polygon_feature['geometry'] = Polygon([tuple(point.to_lat_lon()[::-1]) for point in row.get_points()])
            wrapped_polygon_features.append(wrapped_polygon_feature) 

        results.extend(wrapped_polygon_features)

    return _gpd.GeoDataFrame(results, crs=polygon_features.crs)


# Determine the overriding and subducting plates of the subduction shared sub-segment.
def find_overriding_and_subducting_plates(subduction_shared_sub_segment, time=-999):
    
    # Get the subduction polarity of the nearest subducting line.
    subduction_polarity = subduction_shared_sub_segment.get_feature().get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)
    if (not subduction_polarity) or (subduction_polarity == 'Unknown'):
        print('Unable to find the overriding plate of the subducting shared sub-segment "{0}"'.format(
            subduction_shared_sub_segment.get_feature().get_name()), file=sys.stderr)
        print('    subduction zone feature is missing subduction polarity property or it is set to "Unknown".', file=sys.stderr)
        return

    # There should be two sharing topologies - one is the overriding plate and the other the subducting plate.
    sharing_resolved_topologies = subduction_shared_sub_segment.get_sharing_resolved_topologies()
    if len(sharing_resolved_topologies) != 2:
        print('Unable to find the overriding and subducting plates of the subducting shared sub-segment "{0}" at {1}Ma'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time), file=sys.stderr)
        print('    there are not exactly 2 topologies sharing the sub-segment.', file=sys.stderr)
        print(str(sharing_resolved_topologies[0].get_resolved_feature().get_reconstruction_plate_id()), file=sys.stderr)
        return

    overriding_plate = None
    subducting_plate = None
    
    geometry_reversal_flags = subduction_shared_sub_segment.get_sharing_resolved_topology_geometry_reversal_flags()
    for index in range(2):

        sharing_resolved_topology = sharing_resolved_topologies[index]
        geometry_reversal_flag = geometry_reversal_flags[index]

        if sharing_resolved_topology.get_resolved_boundary().get_orientation() == pygplates.PolygonOnSphere.Orientation.clockwise:
            # The current topology sharing the subducting line has clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line not reversed in the topology.
            if ((subduction_polarity == 'Left' and geometry_reversal_flag) or
                (subduction_polarity == 'Right' and not geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
            else:
                subducting_plate = sharing_resolved_topology
        else:
            # The current topology sharing the subducting line has counter-clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is not reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line reversed in the topology.
            if ((subduction_polarity == 'Left' and not geometry_reversal_flag) or
                (subduction_polarity == 'Right' and geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
            else:
                subducting_plate = sharing_resolved_topology
    
    if overriding_plate is None:
        print('Unable to find the overriding plate of the subducting shared sub-segment "{0}" at {1}Ma'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time), file=sys.stderr)
        print('    both sharing topologies are on subducting side of subducting line.', file=sys.stderr)
        return
    
    if subducting_plate is None:
        print('Unable to find the subducting plate of the subducting shared sub-segment "{0}" at {1}Ma'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time), file=sys.stderr)
        print('    both sharing topologies are on overriding side of subducting line.', file=sys.stderr)
        return
    
    return (overriding_plate, subducting_plate, subduction_polarity)


   