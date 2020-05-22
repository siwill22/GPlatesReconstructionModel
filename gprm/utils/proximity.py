import pygplates

def nearest_feature(point,features):
    # The minimum distance to all features and the nearest feature.
    min_distance_to_all_features = None
    nearest_feature = None

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
                #nearest_feature = feature

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
