import pygplates
import sys
import glob, re
import numpy as np
import matplotlib.pyplot as plt
try:
    from mpl_toolkits.basemap import Basemap
except:
    print 'failed to load plotting dependencies'
import xarray as xr
from . import points_spatial_tree

import .polygon_processing as pp
import .paleogeography as pg




def get_paleogeography_time_list(basedir):
    # make a sorted list of the (midpoint) times for paleogeography polygons
    tmp = glob.glob('%s/*/' % basedir)

    time_list = []
    for tm in tmp:
        time_list.append(float(re.findall(r'\d+Ma+',tm)[1][:-2]))

    return time_list.sort()


def get_masked_multipoint(coords,masking_array,plate_partitioner,valid_time=None):
    # Inputs:
    # a list of coordinates (typically a regular grid of Lat/Long points),
    # an array of indices into that list of coordinates,
    # a set of polygons to use for cookie-cutting
    # a valid time to assign to the output features
    # Returns:
    # a multipoint feature that

    multipoint_feature = pygplates.Feature()
    multipoint_feature.set_geometry(pygplates.MultiPointOnSphere(zip(np.array(coords[0])[masking_array],
                                                                     np.array(coords[1])[masking_array])))
    (pg_points_masked,
     dummy) = plate_partitioner.partition_features(multipoint_feature,
                                                   partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned)

    if valid_time is not None:
        for feature in pg_points_masked:
            feature.set_valid_time(valid_time[0],valid_time[1])

    return pg_points_masked


def get_change_masks(t1,points,spatial_tree_of_uniform_recon_points,psl_t1,psl_t2,rotation_model):

    distance_to_land_t1,distance_to_psl_t1 = pp.run_grid_pnp(t1,
                                                             points,
                                                             spatial_tree_of_uniform_recon_points,
                                                             psl_t1,
                                                             rotation_model)

    distance_to_land_t2,distance_to_psl_t2 = pp.run_grid_pnp(t1,
                                                             points,
                                                             spatial_tree_of_uniform_recon_points,
                                                             psl_t2,
                                                             rotation_model)

    # this mask is true where one (and only one) of the indicators is not zero
    # --> delineates points that change environment between time steps
    #msk = np.logical_xor(distance_to_land_t1,distance_to_land_t2)

    regression_msk = np.logical_and(distance_to_land_t1==0,distance_to_land_t2>0)

    transgression_msk = np.logical_and(distance_to_land_t1>0,distance_to_land_t2==0)

    always_land_msk = np.logical_and(distance_to_land_t1==0,distance_to_land_t2==0)


    return distance_to_land_t1,distance_to_psl_t1,distance_to_land_t2,distance_to_psl_t2,\
           regression_msk,transgression_msk,always_land_msk


# function to run point in/near polygon test for two successive time slices
def get_change_mask_multipoints(pg_features,t1,t2,psl_t1,psl_t2,
                                points,spatial_tree_of_uniform_recon_points,
                                rotation_model,plot=False):

    print 'Working on interpolation from %0.2f Ma to %0.2f Ma .....' % (t1,t2)

    plate_partitioner = pygplates.PlatePartitioner(pg_features, rotation_model, reconstruction_time=t1)

    (distance_to_land_t1,
     distance_to_psl_t1,
     distance_to_land_t2,
     distance_to_psl_t2,
     regression_msk,
     transgression_msk,
     always_land_msk) = get_change_masks(t1,points,spatial_tree_of_uniform_recon_points,
                                         psl_t1,psl_t2,rotation_model)

    coords = zip(*[point.to_lat_lon() for point in points])

    pg_points_regression = get_masked_multipoint(coords,regression_msk,plate_partitioner,
                                                 valid_time=[t2,t1+0.01])
    pg_points_transgression = get_masked_multipoint(coords,transgression_msk,plate_partitioner,
                                                    valid_time=[t2,t1+0.01])
    pg_points_always_land = get_masked_multipoint(coords,always_land_msk,plate_partitioner,
                                                  valid_time=[t2,t1])

    pygplates.FeatureCollection(pg_points_regression).write('./tween_feature_collections/mountain_regression_%0.2fMa_%0.2fMa.gpmlz' % (t1,t2))
    pygplates.FeatureCollection(pg_points_transgression).write('./tween_feature_collections/mountain_transgression_%0.2fMa_%0.2fMa.gpmlz' % (t1,t2))
    pygplates.FeatureCollection(pg_points_always_land).write('./tween_feature_collections/mountain_stable_%0.2fMa_%0.2fMa.gpmlz' % (t1,t2))



def get_vertical_change_multipoints(pg_features,t1,t2,psl_t1,psl_t2,
                                    points,spatial_tree_of_uniform_recon_points,
                                    rotation_model,plot=False):
    # NOT WORKING DUE TO LACK OF SUPPORT FOR SCALAR COVERAGES

    print 'Working on interpolation from %0.2f Ma to %0.2f Ma .....' % (t1,t2)

    plate_partitioner = pygplates.PlatePartitioner(pg_features, rotation_model, reconstruction_time=t1)

    (distance_to_land_t1,
     distance_to_psl_t1,
     distance_to_land_t2,
     distance_to_psl_t2,
     regression_msk,
     transgression_msk,
     always_land_msk) = get_change_masks(t1,points,spatial_tree_of_uniform_recon_points,
                                         psl_t1,psl_t2,rotation_model)

    coords = zip(*[point.to_lat_lon() for point in points])

    pg_points_regression = get_masked_multipoint(coords,regression_msk,plate_partitioner,
                                                 valid_time=[t2,t1+0.01])
    pg_points_transgression = get_masked_multipoint(coords,transgression_msk,plate_partitioner,
                                                    valid_time=[t2,t1+0.01])
    pg_points_always_land = get_masked_multipoint(coords,always_land_msk,plate_partitioner,
                                                  valid_time=[t2,t1])


    # make a scalar coverage
    multi_point = pygplates.MultiPointOnSphere(points)
    scalar_coverages = {
        pygplates.ScalarType.create_gpml('distance_to_land_t1'): distance_to_land_t1,
        pygplates.ScalarType.create_gpml('distance_to_psl_t1'): distance_to_psl_t1,
        pygplates.ScalarType.create_gpml('distance_to_land_t2'): distance_to_land_t2,
        pygplates.ScalarType.create_gpml('distance_to_psl_t2'): distance_to_psl_t2}

    sc_feature = pygplates.Feature()
    sc_feature.set_geometry((multi_point,scalar_coverages))
    sc_feature.set_name('Paleotopography Test Points')

    (cc_sc_features,
     dummy) = plate_partitioner.partition_features(sc_feature,
                                                  partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned)

    pygplates.FeatureCollection(cc_sc_features).write('./tween_feature_collections/mountain_scalar_coverages_%0.2fMa_%0.2fMa.gpmlz' % (t1,t2))




# function to run point in/near polygon test for two successive time slices
def interpolate_paleoshoreline_for_stage(pg_features,t1,t2,psl_t1,psl_t2,time_step,
                                         points,spatial_tree_of_uniform_recon_points,
                                         rotation_model,plot=False):

    print 'Working on interpolation from %0.2f Ma to %0.2f Ma .....' % (t1,t2)

    plate_partitioner = pygplates.PlatePartitioner(pg_features, rotation_model, reconstruction_time=t1)

    (distance_to_land_t1,
     distance_to_psl_t1,
     distance_to_land_t2,
     distance_to_psl_t2,
     regression_msk,
     transgression_msk,
     always_land_msk) = get_change_masks(t1,points,spatial_tree_of_uniform_recon_points,
                                         psl_t1,psl_t2,rotation_model)


    # normalised distance derivation
    # for each point, divide the distance to shoreline at t0 by the total distance to both shorelines
    # --> if the point is halfway between the shorelines, value will be 0.5
    #     if the point is closer to the t1 shoreline, the value will be less than 0.5
    #     all values will be between 0 and 1
    psl_dist_norm = np.divide(distance_to_psl_t1,(distance_to_psl_t1+distance_to_psl_t2))


    t_diff = (t2-t1)



    # before looping over the time steps, create an empty list to put the features in
    pg_points_land_list = []
    pg_points_marine_list = []

    # don't need to do t2 itself, since this will be first step in next iteration
    for reconstruction_time in np.arange(t1,t2,time_step):

        if reconstruction_time==t1:
            land_points = np.where(distance_to_land_t1==0)[0]

        else:
            # normalised time, in range 0 to 1 between start and end of stage
            t_norm = (reconstruction_time-t1)/t_diff

            is_transgressing_land_msk = np.less_equal(psl_dist_norm,t_norm)

            is_regressing_land_msk = np.greater_equal(psl_dist_norm,t_norm)

            land_points = np.where(
                np.logical_or(
                    np.logical_or(
                        np.logical_and(is_regressing_land_msk,regression_msk),
                        np.logical_and(is_transgressing_land_msk,transgression_msk)
                    ),
                always_land_msk))[0]



        marine_mask = np.ones(distance_to_land_t1.shape,dtype=bool)
        marine_mask[land_points] = False
        marine_points = np.where(marine_mask)

        coords = zip(*[point.to_lat_lon() for point in points])


        pg_points_land = get_masked_multipoint(coords,land_points,plate_partitioner,
                                               valid_time=[reconstruction_time+(time_step/2.),
                                               reconstruction_time-(time_step/2.)+0.01])

        pg_points_marine = get_masked_multipoint(coords,marine_points,plate_partitioner,
                                                 valid_time=[reconstruction_time+(time_step/2.),
                                                 reconstruction_time-(time_step/2.)+0.01])


        # append the point features for this time to the overall list
        pg_points_land_list+=pg_points_land
        pg_points_marine_list+=pg_points_marine


    pygplates.FeatureCollection(pg_points_land_list).write('./tween_feature_collections/tweentest_land_%0.2fMa_%0.2fMa.gpmlz' % (t1,t2))
    pygplates.FeatureCollection(pg_points_marine_list).write('./tween_feature_collections/tweentest_ocean_%0.2fMa_%0.2fMa.gpmlz' % (t1,t2))
