import pygplates
import glob, re
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import xarray as xr

from . import paleogeography_tweening as pgt

from .proximity_query import *
from .create_gpml import create_gpml_regular_long_lat_mesh
from . import points_in_polygons
from .sphere_tools import sampleOnSphere
from . import points_spatial_tree
from .fileio import load_netcdf
from .spatial import polygon_area_threshold, merge_polygons

from ptt.utils.call_system_command import call_system_command
import tempfile



# define a function that loads paleogeography multipoints at a specified time
# NOTE this time can be anything, not a time where the multipoints fit nicely together,
# hence the gaps and overlaps will be present
def add_reconstructed_points_to_xyz(points_file,rotation_model,reconstruction_time,zval,mask_file=None):
    
    reconstructed_points = []
    pygplates.reconstruct(points_file,rotation_model,reconstructed_points,reconstruction_time)
    
    #mask_file = './OrogenMasks.gpml'
    mask_file = None
    if mask_file is not None:
        reconstructed_masks = []
        pygplates.reconstruct(mask_file,rotation_model,reconstructed_masks,reconstruction_time)

        for reconstructed_point in reconstructed_points:
            for reconstructed_mask in reconstructed_masks:
                if reconstructed_mask.get_reconstructed_geometry().is_point_in_polygon(reconstructed_point.get_reconstructed_geometry().get_centroid()):
                    reconstructed_point.get_feature().set_name(reconstructed_mask.get_feature().get_name())
        
    point_array = []
    zval_array = []
    for reconstructed_point in reconstructed_points:
        feature_coordinates = reconstructed_point.get_reconstructed_geometry().to_lat_lon_array()
        point_array.append(feature_coordinates)
        #if reconstructed_point.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('OrogenicBelt'):
        if reconstructed_point.get_feature().get_name() == 'OrogenicBelt':
            zval_array.append(np.ones((feature_coordinates.shape[0],1))*zval*3)
        elif reconstructed_point.get_feature().get_name() == 'Cordillera':
            zval_array.append(np.ones((feature_coordinates.shape[0],1))*zval*2)
        else:
            zval_array.append(np.ones((feature_coordinates.shape[0],1))*zval)
      
    xy_array = np.vstack(point_array)
    xyz_array = np.hstack((xy_array,np.vstack(zval_array)))
    
    return xyz_array


# function to facilitate the smoothing of topography
# at the edge of mountain range polygons
def get_distance_to_mountain_edge(point_array,reconstruction_basedir,rotation_model,time,area_threshold):
    
    distance_threshold_radians=None
    env_list = ['m']

    if time==0:
        pg_dir = './present_day_paleogeography.gmt'
        pg_features = pg.load_paleogeography(pg_dir,env_list,single_file=True,env_field='Layer')
    else:
        pg_dir = '%s/PresentDay_Paleogeog_Matthews2016_%dMa/' % (reconstruction_basedir,time)
        pg_features = pg.load_paleogeography(pg_dir,env_list)

    cf = merge_polygons(pg_features, rotation_model, time=time, sampling=0.25)
    sieve_polygons_t1 = polygon_area_threshold(cf, area_threshold)

    polygons_as_list = []
    for feature in sieve_polygons_t1:
        polygons_as_list.append(feature.get_geometry())
        
    res1 = find_closest_geometries_to_points([pygplates.PointOnSphere(point) for point in zip(point_array[:,0],point_array[:,1])],
                                             polygons_as_list,
                                             distance_threshold_radians = distance_threshold_radians)
    
    distance_to_polygon_boundary = np.degrees(np.array(zip(*res1)[0]))

    # Make a copy of list of distances.
    distance_to_polygon = list(distance_to_polygon_boundary)

    # Set distance to zero for any points inside a polygon (leave other points unchanged).
    res2 = points_in_polygons.find_polygons([pygplates.PointOnSphere(point) for point in zip(point_array[:,0],point_array[:,1])],
                                            polygons_as_list)

    for point_index, rpolygon in enumerate(res2):
        # If not inside any polygons then result will be None.
        if rpolygon is None:
            distance_to_polygon[point_index] = 0.0
            
    return distance_to_polygon


# use merged seive_polygons to get a regular lat-long multipoint that will contain points
# only within the COB Terranes (ie not within the 'deep ocean')
def get_land_sea_multipoints(sieve_polygons, sampling, depth_for_unknown_ocean, subdivision_depth=4):

    multipoints = create_gpml_regular_long_lat_mesh(sampling)
    grid_dims = (int(180/sampling)+1,int(360/sampling)+1)

    for multipoint in multipoints:
        for mp in multipoint.get_all_geometries():
            points = mp.to_lat_lon_point_list()

    rpolygons = []
    for polygon in sieve_polygons:
        if polygon.get_geometry():
            rpolygons.append(polygon.get_geometry())

    polygons_containing_points = points_in_polygons.find_polygons(points, rpolygons, subdivision_depth=subdivision_depth)

    lat = []
    lon = []
    zval = []

    lat_deep = []
    lon_deep = []
    zval_deep = []

    for pcp,point in zip(polygons_containing_points, points):
        if pcp is not None:
            lat.append(point.get_latitude())
            lon.append(point.get_longitude())
        else:
            lat_deep.append(point.get_latitude())
            lon_deep.append(point.get_longitude())
            zval_deep.append(depth_for_unknown_ocean)
            
    return lat,lon,zval,lat_deep,lon_deep,zval_deep



def paleotopography_job(reconstruction_time, paleogeography_timeslice_list, 
                        tween_basedir, reconstruction_basedir, output_dir, 
                        file_format, rotation_file, COBterrane_file, agegrid_file_template,
                        lowland_elevation, shallow_marine_elevation, max_mountain_elevation, depth_for_unknown_ocean, 
                        sampling, mountain_buffer_distance_degrees, area_threshold,
                        grid_smoothing_wavelength_kms, merge_with_bathymetry, 
                        land_or_ocean_precedence='land', netcdf3_output=False, subdivision_depth=4):


    print('Working on Time %0.2fMa\n' % reconstruction_time)
        
    rotation_model = pygplates.RotationModel(rotation_file)
                        
    # find times that bracket the selected exact time in the paleogeography source files
    time_stage_max = paleogeography_timeslice_list[np.where(paleogeography_timeslice_list>reconstruction_time)[0][0]]
    time_stage_min = paleogeography_timeslice_list[np.where(paleogeography_timeslice_list<=reconstruction_time)[0][-1]]

    # Note the logic for selecting the times:
    # The main issue is that each set of paleogeography polygons has a defined 'midpoint' time
    # --> if the reconstruction time is between these, the choice of t1 and t2 is obvious
    # --> if the reconstruction time matches one of these times, then we can work directly on
    #     the geometries that match this time - hence the two routes through the if statement below
    
    print('Selected Time is in the stage %0.2fMa to %0.2fMa' % (time_stage_min,time_stage_max))

    land_points_file = '%s/tweentest_land_%0.2fMa_%0.2fMa.%s' % (tween_basedir,time_stage_min,time_stage_max,file_format)
    marine_points_file = '%s/tweentest_ocean_%0.2fMa_%0.2fMa.%s' % (tween_basedir,time_stage_min,time_stage_max,file_format)
    mountains_going_up_file = '%s/mountain_transgression_%0.2fMa_%0.2fMa.%s' % (tween_basedir,time_stage_min,time_stage_max,file_format)
    mountains_going_down_file = '%s/mountain_regression_%0.2fMa_%0.2fMa.%s' % (tween_basedir,time_stage_min,time_stage_max,file_format)
    mountains_stable_file = '%s/mountain_stable_%0.2fMa_%0.2fMa.%s' % (tween_basedir,time_stage_min,time_stage_max,file_format)

    
    # get a nx3 array defining the points above sea-level, reconstructed to time of interest
    # columns are [lat, long, elevation assigned for lowland]
    land_point_array = add_reconstructed_points_to_xyz(land_points_file,
                                                       rotation_model,
                                                       reconstruction_time,
                                                       lowland_elevation)
    
    # get a nx3 array defining shallow marine areas, reconstructed to time of interest
    # columns are [lat, long, elevation assigned for shallow marine]
    marine_point_array = add_reconstructed_points_to_xyz(marine_points_file,
                                                         rotation_model,
                                                         reconstruction_time,
                                                         shallow_marine_elevation)
    
    # Note that the two arrays just created are based on 'regular' lat/long grids, but 
    # are not aligned with the regular lat/long grid that we want to output
    # since they are (usually) reconstructed to a different time from the one at which they
    # were created (and anyway may be at a different resolution to the grid sampling specified
    # here)

    # combine the previous two arrays
    pg_point_array = np.vstack((land_point_array, marine_point_array))

    # get a merged version of COB terranes, optionally excluding polygons that are small in area
    # TODO deal with donut polygons better
    sieve_polygons = get_merged_cob_terrane_polygons(COBterrane_file, rotation_model,
                                                     reconstruction_time, sampling)

    # get arrays defining the land and sea based on which points fall within the COB terranes
    # NOTE this step is where we create the points that ARE on the regular lat/long grid we 
    # will ultimately output
    (lat,lon,zval,
     lat_deep,lon_deep,zval_deep) = get_land_sea_multipoints(sieve_polygons, sampling, depth_for_unknown_ocean,
                                                             subdivision_depth=subdivision_depth)


    # sample the land/marine points onto the points within the COB Terranes
    # This will fill the gaps that exist within continents, and average out overlaps
    d,l = sampleOnSphere(pg_point_array[:,1],pg_point_array[:,0],pg_point_array[:,2],
                         np.array(lon),np.array(lat),n=1)

    land_marine_interp_points = pg_point_array[:,2].ravel()[l]

    # At this point, the land points are all considered to be 'lowland'......
    
    ####################################
    # Deal with the mountains
    if np.equal(reconstruction_time, time_stage_min):
        print('Temporary fix for valid time')
        #dat3 = add_reconstructed_points_to_xyz(mountains_going_up_file,rotation_model,reconstruction_time,3)
        dat4 = add_reconstructed_points_to_xyz(mountains_going_down_file, rotation_model, reconstruction_time+0.01,1)
        dat5 = add_reconstructed_points_to_xyz(mountains_stable_file, rotation_model, reconstruction_time+0.01,1)
        mountains_tr_point_array = np.vstack((dat4, dat5))
        
        dist_tr = get_distance_to_mountain_edge(mountains_tr_point_array, reconstruction_basedir,
                                                rotation_model, reconstruction_time, area_threshold)
        dist_tr_cap = np.array(dist_tr)
        dist_tr_cap[np.array(dist_tr)>mountain_buffer_distance_degrees] = mountain_buffer_distance_degrees

        dist_tr_cap = dist_tr_cap*mountains_tr_point_array[:,2]

        normalized_mountain_elevation = dist_tr_cap
        
    else:
        # load in the mountain points but at three different times: t1 and t2, and the reconstruction time
        # note that these three arrays should all be identical in size, since they are the same multipoints
        # just reconstructed to three slightly different times
        dat3 = add_reconstructed_points_to_xyz(mountains_going_up_file, rotation_model, time_stage_max, 1)
        dat4 = add_reconstructed_points_to_xyz(mountains_going_down_file, rotation_model, time_stage_max, 1)
        dat5 = add_reconstructed_points_to_xyz(mountains_stable_file, rotation_model, time_stage_max, 1)
        mountains_t2_point_array = np.vstack((dat3, dat4, dat5))

        dat3 = add_reconstructed_points_to_xyz(mountains_going_up_file, rotation_model, time_stage_min+0.01, 1)
        dat4 = add_reconstructed_points_to_xyz(mountains_going_down_file, rotation_model, time_stage_min+0.01, 1)
        dat5 = add_reconstructed_points_to_xyz(mountains_stable_file, rotation_model, time_stage_min+0.01, 1)
        mountains_t1_point_array = np.vstack((dat3, dat4, dat5))

        dat3 = add_reconstructed_points_to_xyz(mountains_going_up_file, rotation_model, reconstruction_time, 1)
        dat4 = add_reconstructed_points_to_xyz(mountains_going_down_file, rotation_model, reconstruction_time, 1)
        dat5 = add_reconstructed_points_to_xyz(mountains_stable_file, rotation_model, reconstruction_time, 1)
        mountains_tr_point_array = np.vstack((dat3, dat4, dat5))

        # calculate distances of the mountain points to the edge of the mountain region at t1 and t2,
        # using the pg polygons that they should exactly correspond to 
        dist_t1 = get_distance_to_mountain_edge(mountains_t1_point_array, reconstruction_basedir,
                                                rotation_model, time_stage_min, area_threshold)
        dist_t2 = get_distance_to_mountain_edge(mountains_t2_point_array, reconstruction_basedir,
                                                rotation_model, time_stage_max, area_threshold)
        
        #is_in_orogeny_index = find_mountain_type(mountains_tr_point_array,
        #                                         orogeny_feature_filename,
        #                                         reconstruction_time)


        # cap the distances at some arbitrary value defined earlier
        dist_t1_cap = np.array(dist_t1)
        dist_t1_cap[np.array(dist_t1)>mountain_buffer_distance_degrees] = mountain_buffer_distance_degrees

        dist_t1_cap = dist_t1_cap*mountains_t1_point_array[:,2]
        
        dist_t2_cap = np.array(dist_t2)
        dist_t2_cap[np.array(dist_t2)>mountain_buffer_distance_degrees] = mountain_buffer_distance_degrees
        
        dist_t2_cap = dist_t2_cap*mountains_t2_point_array[:,2]

        # get the normalised time within this time stage
        # for example we are at 0.25 between the t1 and t2
        t_diff = (time_stage_max-time_stage_min)
        t_norm = (reconstruction_time-time_stage_min)/t_diff

        # use 1d interpolation to get the 'normalized' height of the mountains at the preceding
        # and subsequent times to the specific reconstruction time
        # [note this is not spatial interpolation - rather it is interpolation at each individual point
        # between the heights at earlier and later times]
        tmp = np.vstack((dist_t1_cap, dist_t2_cap))
        f = interpolate.interp1d([0,1],tmp.T)
        normalized_mountain_elevation = f(t_norm)
    
    
    # interpolate the elevations at tr onto the regular long lat points that we will ultimately use 
    # for the grid output
    # note the k value here controls number of neighbouring points used in inverse distance average
    d,l = sampleOnSphere(mountains_tr_point_array[:,1], mountains_tr_point_array[:,0], normalized_mountain_elevation,
                         np.array(lon), np.array(lat), k=4)
    w = 1./d**2
    normalized_mountain_elevation_interp_points = np.sum(w * normalized_mountain_elevation.ravel()[l],axis=1) / np.sum(w,axis=1)

    # this index isolates only those points that are within a certain distance of the mountain range
    # (since the interpolation will give values everywhere in the 'land', so we want to re-isolate only 
    # those points that fall within the 'mountain' regions)
    # TODO this should be set to the sampling??
    mountain_proximity_index = np.degrees(np.min(d, axis=1))< sampling*2 #mountain_buffer_distance_degrees

    
    #####################################
    # Put the grid together
    #####################################
    
    # write the land/marine points to a file
    land_marine_xyz_file = tempfile.NamedTemporaryFile()
    mountain_xyz_file = tempfile.NamedTemporaryFile()
    land_marine_nc_file = tempfile.NamedTemporaryFile()
    mountain_nc_file = tempfile.NamedTemporaryFile()
    
    write_xyz_file(land_marine_xyz_file.name, zip(lon+lon_deep,
                                                  lat+lat_deep,
                                                  np.hstack((land_marine_interp_points, zval_deep))))

    # convert the normalized mountain elevations to metres, then write to file
    mountain_elevation_factor = max_mountain_elevation / mountain_buffer_distance_degrees
    mountain_elevation_array = normalized_mountain_elevation_interp_points[mountain_proximity_index]*mountain_elevation_factor
    write_xyz_file(mountain_xyz_file.name, zip(np.array(lon)[mountain_proximity_index],
                                               np.array(lat)[mountain_proximity_index],
                                               mountain_elevation_array))

    # all the points are already on the same regular lat/long grid (but with gaps) - just 
    # need to piece them all together and combine.
    # Note we assume the the mountain elevation is the height IN ADDITION to the lowland elevation
    # so that we can simply add them
    call_system_command(['gmt', 'xyz2grd', land_marine_xyz_file.name, '-Rd', '-I%0.8f' % sampling, 
                         '-G%s' % land_marine_nc_file.name])
    call_system_command(['gmt', 'xyz2grd', mountain_xyz_file.name, '-Rd', '-I%0.8f' % sampling, '-di0', 
                         '-G%s' % mountain_nc_file.name])
    call_system_command(['gmt', 'grdmath', mountain_nc_file.name, land_marine_nc_file.name, 'ADD', '=', 
                         '%s/paleotopo_%0.2fd_%0.2fMa.nc' % (output_dir, sampling, reconstruction_time)])

    # clean-up temp files
    land_marine_xyz_file.delete
    mountain_xyz_file.delete
    land_marine_nc_file.delete
    mountain_nc_file.delete
    
    # load result back into python
    topoX,topoY,topoZ = load_netcdf('%s/paleotopo_%0.2fd_%0.2fMa.nc' % (output_dir, sampling, reconstruction_time))

    
    if merge_with_bathymetry:
    
    # TODO create seperate function for this step
    
        # PALEOBATHYMETRY based on age grids
        # load age grid for this time and calculate paleobathymetry
        agegrid_file = agegrid_file_template % reconstruction_time

        ageX,ageY,ageZ = load_netcdf(agegrid_file)

        paleodepth = pg.age2depth(ageZ,model='GDH1')

        # get index for grid nodes where age grid is nan, replace values with topography/shallow bathymetry
        land_or_ocean_precedence = 'land'
        if land_or_ocean_precedence is 'ocean':
            not_bathy_index = np.isnan(paleodepth)
            paleodepth[not_bathy_index] = topoZ[not_bathy_index]
        else:
            not_bathy_index = np.greater(topoZ,depth_for_unknown_ocean)
            paleodepth[not_bathy_index] = topoZ[not_bathy_index]
            leftover_nans = np.isnan(paleodepth)
            paleodepth[leftover_nans] = depth_for_unknown_ocean
            
        
        plt.figure()
        plt.imshow(not_bathy_index)
        plt.show()

        paleotopobathy_nc_file = tempfile.NamedTemporaryFile()
        paleotopobathy_smooth_nc_file = tempfile.NamedTemporaryFile()
        
        # save the merged grid (forcing compatibility with GPlates-readable netCDF in case it helps)
        ds = xr.DataArray(paleodepth,
                          coords=[('lat',topoY), ('lon',topoX)],
                          name='elevation')
        ds.to_netcdf(paleotopobathy_nc_file.name, format='NETCDF3_CLASSIC')

        # smooth the grid using GMT [wavelength is optional
        #pg.smooth_topography_grid('paleotopobathy.nc','paleotopobathy_smooth_%0.2fMa.nc' % reconstruction_time,400.)
        # TODO skip this step if grid_smoothing_wavelength_kms set to zero
        call_system_command(['gmt', 'grdfilter', paleotopobathy_nc_file.name, '-G%s' % paleotopobathy_smooth_nc_file.name, 
                             '-Fg%0.2f' % grid_smoothing_wavelength_kms, '-fg', '-D4', '-Vl'])

        if netcdf3_output:
            # finally, once again force GPlates-readable netCDF (ie netCDF v3) and put the 
            # grid in the output folder with a filename containing the age
            call_system_command(['gmt', 'grdconvert', paleotopobathy_smooth_nc_file.name, 
                                 '-G%s=cf' % '%s/paleotopobathy_buffer%0.2dd_filter_%0.2fkm_%0.2fMa.nc' % (output_dir,
                                                                                              mountain_buffer_distance_degrees,
                                                                                              grid_smoothing_wavelength_kms,
                                                                                              sampling, 
                                                                                              reconstruction_time)])
        else:
            call_system_command(['cp', paleotopobathy_smooth_nc_file.name, 
                                '%s/paleotopobathy_buffer%0.2dd_filter%0.2fkm_%0.2fd_%0.2fMa.nc' % (output_dir,
                                                                                              mountain_buffer_distance_degrees,
                                                                                              grid_smoothing_wavelength_kms,
                                                                                              sampling, 
                                                                                              reconstruction_time)])
        
        # load and plot the result
        topo_smoothX,topo_smoothY,topo_smoothZ = load_netcdf(paleotopobathy_smooth_nc_file.name)
        #
        plt.figure(figsize=(25,11))
        plt.imshow(topo_smoothZ,origin='lower',
                   extent=[-180,180,-90,90],
                   cmap=plt.cm.terrain,vmin=-5000,vmax=5000)
        plt.title('%0.2fMa' % reconstruction_time)
        plt.colorbar()
        plt.savefig('%s/paleotopobathy_buffer%0.2dd_filter%0.2fkm_%0.2fd_%0.2fMa.png' % (output_dir,
                                                                     mountain_buffer_distance_degrees,
                                                                     grid_smoothing_wavelength_kms,
                                                                     sampling,
                                                                     reconstruction_time))
        plt.close()

        paleotopobathy_nc_file.delete
        paleotopobathy_smooth_nc_file.delete
