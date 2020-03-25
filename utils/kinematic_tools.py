import pandas as pd
import numpy as np
import pygplates

def subduction_parameters(seed_point_features,rotation_model):

    df_AllTimes = pd.read_hdf('/Users/Simon/GIT/pygplates-alpha/subduction/SubductionTable_0_200Ma.h5')
    time_step = 1

    DataFrameTemplate = df_AllTimes.columns

    df_OreDepositBirthTimeStats = pd.DataFrame(columns = DataFrameTemplate)

    for seed_point in seed_point_features:
        
        nearest_time = np.round(seed_point.get_valid_time()[0]/time_step) * time_step
        
        ore_deposit_rotation = rotation_model.get_rotation(nearest_time, 
                                                           seed_point.get_reconstruction_plate_id(), 
                                                           anchor_plate_id=0)
        
        reconstructed_point = ore_deposit_rotation * seed_point.get_geometry()
        reconstructed_point_degrees = reconstructed_point.to_lat_lon_list()[0]
        
        subset = df_AllTimes[df_AllTimes.time==nearest_time]
        
        if subset.empty is False:
            ore_deposit_subduction_values = get_nearest_subduction_point(
                subset,reconstructed_point_degrees)  
    
            ore_deposit_subduction_values['lat0'] = seed_point.get_geometry().to_lat_lon()[0]
            ore_deposit_subduction_values['lon0'] = seed_point.get_geometry().to_lat_lon()[1]
            #ore_deposit_subduction_values['index'] = np.int(seed_point.get_name())
            ore_deposit_subduction_values['DepositAge'] = seed_point.get_valid_time()[0]
    
            df_OreDepositBirthTimeStats = df_OreDepositBirthTimeStats.append(ore_deposit_subduction_values)
    
    return df_OreDepositBirthTimeStats


################################################

def get_nearest_subduction_point(subduction_data,seed_point):

    subduction_stats_at_seed_points = []
    
    min_distance_to_all_features = None
    nearest_subduction_point = None

    for index,row in subduction_data.iterrows():
        min_distance_to_feature = pygplates.GeometryOnSphere.distance(
            pygplates.PointOnSphere(row['lat'], row['lon']),
            pygplates.PointOnSphere(seed_point),
            min_distance_to_all_features)

        # If the current geometry is nearer than all previous geometries then
        # its associated feature is the nearest feature so far.
        if min_distance_to_feature is not None:
            min_distance_to_all_features = min_distance_to_feature
            nearest_subduction_point = row

    return nearest_subduction_point

###############################################
