import pygplates
import pandas as pd

def generate_rotation_feature(rotation_filenames):

    rotation_features = pygplates.FeatureCollection()
    for rotation_filename in rotation_filenames:
        rotation_features.add(pygplates.FeatureCollection(rotation_filename))

    return rotation_features


def get_rotation_table(rotation_features, plate_id_list=None, asdataframe=False):

    list_of_rotations = []
    for feature in rotation_features:

        (fixed_plate_id, 
         moving_plate_id, 
         total_reconstruction_pole) = feature.get_total_reconstruction_pole()

        # Ignore moving plate IDs equal to 999 since these are commented lines in the PLATES4 rotation format.
        if moving_plate_id == 999:
            continue

        if plate_id_list is not None:
            if moving_plate_id not in plate_id_list:
                continue

        # Print the time period of the rotation feature.
        # This is the times of the first and last enabled rotation time samples.
        enabled_time_samples = total_reconstruction_pole.get_enabled_time_samples()

        # Print the rotation pole information from the enabled rotation time samples.
        #print '  time samples:'
        for time_sample in enabled_time_samples:

            # Get the finite rotation from the GpmlFiniteRotation property value instead the GpmlTimeSample.
            finite_rotation = time_sample.get_value().get_finite_rotation()

            # Extract the pole and angle (in degrees) from the finite rotation.
            pole_lat, pole_lon, pole_angle = finite_rotation.get_lat_lon_euler_pole_and_angle_degrees()

            # The time and optional description come from the GpmlTimeSample.
            time = time_sample.get_time()
            description = time_sample.get_description()

            #list_of_rotation_descriptions.append(description)
            list_of_rotations.append([moving_plate_id, time, pole_lat, pole_lon, pole_angle, fixed_plate_id, description])
            
    if asdataframe:
        return pd.DataFrame(list_of_rotations, 
                            columns=['MovingPlateID','Time','Lat','Long','Angle','FixedPlateID','Description'])
    else:
        return list_of_rotations

