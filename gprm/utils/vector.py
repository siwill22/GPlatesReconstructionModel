import pygplates
from .create_gpml import create_gpml_healpix_mesh

def get_velocities(rotation_model,topology_features,time,velocity_domain_features=None,delta_time=1,velocity_type='MagAzim'):

    if velocity_domain_features is None:
        velocity_domain_features = create_gpml_healpix_mesh(32,feature_type='MeshNode')

    # All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
    all_domain_points = []
    all_velocities = []
    plate_ids = []

    # Partition our velocity domain features into our topological plate polygons at the current 'time'.
    plate_partitioner = pygplates.PlatePartitioner(topology_features, rotation_model, time)

    for velocity_domain_feature in velocity_domain_features:

        # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
        # Iterate over them all.
        for velocity_domain_geometry in velocity_domain_feature.get_geometries():

            for velocity_domain_point in velocity_domain_geometry.get_points():

                all_domain_points.append(velocity_domain_point)

                partitioning_plate = plate_partitioner.partition_point(velocity_domain_point)
                if partitioning_plate:

                    # We need the newly assigned plate ID to get the equivalent stage rotation of that tectonic plate.
                    partitioning_plate_id = partitioning_plate.get_feature().get_reconstruction_plate_id()

                    # Get the stage rotation of partitioning plate from 'time + delta_time' to 'time'.
                    equivalent_stage_rotation = rotation_model.get_rotation(time, partitioning_plate_id, time + delta_time)

                    # Calculate velocity at the velocity domain point.
                    # This is from 'time + delta_time' to 'time' on the partitioning plate.
                    velocity_vectors = pygplates.calculate_velocities(
                        [velocity_domain_point],
                        equivalent_stage_rotation,
                        delta_time)

                    if velocity_type=='east_north':
                        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
                        velocities = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                                [velocity_domain_point],
                                velocity_vectors)
                        all_velocities.append(velocities[0])

                    else:
                        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
                        velocities = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                                [velocity_domain_point],
                                velocity_vectors)
                        all_velocities.append(velocities[0])

                    plate_ids.append(partitioning_plate_id)

                else:
                    # If point is not within a polygon, set velocity and plate_id to zero
                    all_velocities.append((0,0,0))
                    plate_ids.append(0)

    pt_vel1=[]
    pt_vel2=[]
    if velocity_type=='east_north':
        for velocity_vector in all_velocities:
            if getattr(velocity_vector,'get_x',None) is not None:
                pt_vel1.append(velocity_vector.get_y())
                pt_vel2.append(velocity_vector.get_x())
            else:
                pt_vel1.append(0.)
                pt_vel2.append(0.)
    else:
        for velocity_vector in all_velocities:
            pt_vel1.append(velocity_vector[0])
            pt_vel2.append(velocity_vector[1])

    pt_lon = []
    pt_lat = []
    for pt in all_domain_points:
        pt_lon.append(pt.to_lat_lon()[1])
        pt_lat.append(pt.to_lat_lon()[0])

    return pt_lat,pt_lon,pt_vel1,pt_vel2,plate_ids



def reconstruct_vector_points(point_tuples, vector_tuples, point_rotation):
    # given lists of points and vectors, a rotation pole, returns 
    # a list of reconstructed points and vector directions
   
    if len(point_tuples) != len(vector_tuples):
        raise ValueError('expected same number of points and vectors')
        
    reconstructed_points = [point_rotation * pygplates.PointOnSphere(point_tuple) for point_tuple in point_tuples]
   
    reconstructed_vectors = [point_rotation * pygplates.LocalCartesian.convert_from_magnitude_azimuth_inclination_to_geocentric(point_tuple, vector_tuple) for (point_tuple,vector_tuple) in zip(point_tuples,vector_tuples)]
   
    reconstructed_vector_magnitude_azimuth_inclinations = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(reconstructed_points,
                                                                                                                                            reconstructed_vectors)

    return (
        [reconstructed_point.to_lat_lon() for reconstructed_point in reconstructed_points],
        reconstructed_vector_magnitude_azimuth_inclinations)


def reconstruct_vectors(gdf, rotation_model, reconstruction_time, direction_col_name='Direction', anchor_plate_id=0):

    point_tuples = list(zip(gdf.geometry.y,gdf.geometry.x))

    vector_tuples = [(1,direction,0) for direction in gdf[direction_col_name]]

    point_rotations = [rotation_model.get_rotation(reconstruction_time,
                                                   reconstruction_plate_id,
                                                   anchor_plate_id=anchor_plate_id) for reconstruction_plate_id in gdf.PlateID1]

    reconstruct_vector_points(point_tuples, vector_tuples, point_rotations)


