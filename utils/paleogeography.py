import pygplates
import glob
import numpy as np
import os
import sys
import xarray as xr
import scipy.interpolate as spi
from create_gpml import create_gpml_regular_long_lat_mesh, create_gpml_healpix_mesh

try:
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
except:
    print('Failed to load plotting dependencies')


def load_paleogeography(pg_dir,env_list=None,
                        single_file=False,env_field='ENV'):

    # default environment_list is the format used for Cao++ 2017
    if env_list is None:
        env_list = ['lm','m','sm','i']

    if single_file:
        print(pg_dir)
        features = pygplates.FeatureCollection(pg_dir)
        pg_features = []
        for feature in features:
            if feature.get_shapefile_attribute(env_field) in env_list:
                feature.set_shapefile_attribute('Layer',feature.get_shapefile_attribute(env_field))
                pg_features.append(feature)

    else:
        pg_features = []
        for env in env_list:
            try:
                filename = glob.glob('%s/%s_*.shp' % (pg_dir,env))
                print(filename)
                features = pygplates.FeatureCollection(filename[0])
                for feature in features:
                    feature.set_shapefile_attribute('Layer',env)
                    pg_features.append(feature)

            except:
                print('no features of type %s' % env)

    return pg_features


def rasterise_paleogeography(pg_features,rotation_model,time,
							 sampling=0.5,env_list=None,meshtype='LongLatGrid',
							 masking=None):
    # takes paleogeography polygons like those from Cao++ 2017 and converts them
    # into a raster
    # if meshtype is set to 'healpix', sampling should be set to an integer defining nSide

    #pg_features = load_paleogeography(pg_dir,env_list)
    if meshtype=='healpix':
        raster_domain = create_gpml_healpix_mesh(sampling,filename=None,feature_type='MeshNode')
    else:
        raster_domain = create_gpml_regular_long_lat_mesh(sampling,filename=None,feature_type='MeshNode')

    plate_partitioner = pygplates.PlatePartitioner(pg_features, rotation_model, reconstruction_time=time)

    if masking is not None:
        pg_points = plate_partitioner.partition_features(raster_domain,
														 partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned,
                                                         properties_to_copy=[pygplates.PropertyName.gpml_shapefile_attributes])
        if masking == 'Outside':
            pg_points = pg_points[0]
        elif masking == 'Inside':
            pg_points = pg_points[1]

    else:
        pg_points = plate_partitioner.partition_features(raster_domain,
                                                         properties_to_copy=[pygplates.PropertyName.gpml_shapefile_attributes])

    return pg_points


def paleogeography2topography_xyz(pg_points,topo_dict,sampling,
                                  fill_value=-3000,bathymetry_points=None,grdfile=None):
    # given some point features, and a dictionary relating the 'Layer'
    # to height, returns a grid of points with height values mapped from
    # the disctionary relating environment to height

    Xr = []
    Yr = []
    Zr = []

    for feature in pg_points:
        env = feature.get_shapefile_attribute('Layer')
        if env is not None:
            height = topo_dict[env]
        elif bathymetry_points is None:
            height = fill_value
        else:
            continue
        for geometry in feature.get_geometries():
            for point in geometry.get_points():
                Xr.append(point.to_lat_lon()[1])
                Yr.append(point.to_lat_lon()[0])
                Zr.append(height)

    # if bathymetry points are provided, append them point-by-point
    if bathymetry_points is not None:
        for bathymetry_point in bathymetry_points:
            Xr.append(bathymetry_point.get_geometry().to_lat_lon()[1])
            Yr.append(bathymetry_point.get_geometry().to_lat_lon()[0])
            Zr.append(bathymetry_point.get_shapefile_attribute('depth'))

    if grdfile is None:
        return Xr,Yr,Zr
    else:
        tmp = np.vstack((Xr,Yr,Zr)).T
        np.savetxt('test.asc',tmp,fmt='%0.4f,%0.4f,%0.4f')

        os.system('gmt xyz2grd test.asc -Rd -I%0.6f -G%s' % (sampling,grdfile))



def smooth_topography_grid(grdfile,filt_grdfile,wavelength):
    # smooths a GMT grid using Gaussian filter of specified size (in kms)

    os.system('gmt grdfilter %s -G%s -Fg%0.2f -D4 -Vl' % (grdfile,filt_grdfile,wavelength))


def load_netcdf(grdfile,z_field_name='z'):

    ds_disk = xr.open_dataset(grdfile)

    data_array = ds_disk[z_field_name]
    coord_keys = data_array.coords.keys()

    if 'lat' in coord_keys[0].lower():
        latitude_key=0; longitude_key=1
    else:
        latitude_key=1; longitude_key=0

    try:
        gridX = data_array.coords[coord_keys[longitude_key]].data
        gridY = data_array.coords[coord_keys[latitude_key]].data
        gridZ = data_array.data
    except:
        # attempt to handle old-school GMT netcdfs (e.g. produced by grdconvert)
        gridX = np.linspace(ds_disk.data_vars['x_range'].data[0],
                            ds_disk.data_vars['x_range'].data[1],
                            ds_disk.data_vars['dimension'].data[0])
        gridY = np.linspace(ds_disk.data_vars['y_range'].data[0],
                            ds_disk.data_vars['y_range'].data[1],
                            ds_disk.data_vars['dimension'].data[1])
        gridZ = np.flipud(ds_disk.data_vars[z_field_name].data.reshape(ds_disk.data_vars['dimension'].data[1],
                                                                       ds_disk.data_vars['dimension'].data[0]))

    ds_disk.close()

    return gridX,gridY,gridZ


def create_slice(gridX,gridY,gridZ,GCPts,ProfilePoints):
    # make a cross-section across a grid, given (two or more) points
    # defined in lat/long. Profiles are defined as great-circle paths between
    # defined points

    f = spi.RectBivariateSpline(gridX,gridY,gridZ.T)
    XVals = f.ev(GCPts[:,1], GCPts[:,0])
    Zval = XVals.flatten()

    return Zval


def create_profile_points(PtLons,PtLats,PointSpacing = 0.5):

    polyline_features = []
    polyline = pygplates.PolylineOnSphere(zip(PtLats,PtLons))
    polyline_feature = pygplates.Feature()
    polyline_feature.set_geometry(polyline)
    polyline_features.append(polyline_feature)

    # Define point spacing in arc-degrees
    PointSpacing = 0.5

    for feature in polyline_features:
        geometry = feature.get_geometry()
        arc_distance = np.degrees(geometry.get_arc_length())
        tesselated_polyline = geometry.to_tessellated(np.radians(PointSpacing))
        GCPts = tesselated_polyline.to_lat_lon_array()

        # Actually this is wrong - since it will only give 'near to' 15 degrees, not exact
        label_points = geometry.to_tessellated(np.radians(15)).to_lat_lon_array()

        arc_distance = np.degrees(geometry.get_arc_length())
        ProfilePoints = np.linspace(-arc_distance/2,arc_distance/2,GCPts.shape[0])

    return GCPts,ProfilePoints,arc_distance


def profile_plate_ids(resolved_topologies,rotation_model,GreatCirclePoints):

    partitioner = pygplates.PlatePartitioner(resolved_topologies,rotation_model)

    plate_ids = []
    for point in GreatCirclePoints:
        partitioned_point = partitioner.partition_point(pygplates.PointOnSphere(point))
        plate_ids.append(partitioned_point.get_feature().get_reconstruction_plate_id())

    return plate_ids


def plate_boundary_intersections(cross_section_geometry,shared_boundary_sections,ProfileX_kms):

    # Given a polyline, and the subduction boundary sections, finds places where the cross-section
    # intersects a plate boundary
    # returns the Lat/Long coordinates and the distance along profile

    subduction_intersections = []
    ridge_intersections = []
    other_intersections = []

    for shared_boundary_section in shared_boundary_sections:

        for shared_subsegment in shared_boundary_section.get_shared_sub_segments():

            (min_distance_to_feature,
            closest_point_on_section,
            closest_point_on_topology,
            section_index,
            topology_index) = pygplates.GeometryOnSphere.distance(
                cross_section_geometry,
                shared_subsegment.get_resolved_geometry(),
                return_closest_positions=True,
                return_closest_indices=True)

            if min_distance_to_feature == 0:

                # find the distance along the section profile
                #print section_index
                cross_section_segment = cross_section_geometry.get_segments()[section_index]
                distance_along_segment = pygplates.GeometryOnSphere.distance(closest_point_on_section,cross_section_segment.get_start_point())
                distance_long_profile = distance_along_segment*pygplates.Earth.mean_radius_in_kms + ProfileX_kms[section_index]

                if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('SubductionZone'):
                    polarity = get_subduction_polarity(shared_subsegment,topology_index,cross_section_segment,distance_along_segment)
                    subduction_intersections.append([closest_point_on_section,distance_long_profile,polarity])

                elif shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('MidOceanRidge'):
                    ridge_intersections.append([closest_point_on_section,distance_long_profile])
                else:
                    other_intersections.append([closest_point_on_section,distance_long_profile])

    return subduction_intersections,ridge_intersections,other_intersections


def get_subduction_polarity(shared_subsegment,topology_index,cross_section_segment,distance_along_segment):
# gets the subduction polarity at locations where a subduction segment intersects
# another line segment

    topology_section_segment = shared_subsegment.get_resolved_geometry().get_segments()[topology_index]

    # Get cross-section segment direction at point of intersection.
    if not cross_section_segment.is_zero_length():
        cross_section_segment_direction = cross_section_segment.get_arc_direction(
                distance_along_segment / cross_section_segment.get_arc_length())
    else:
        cross_section_segment_direction = cross_section_segment.get_arc_direction(0)

    # Imagine topology going South to North. If cross section crosses it going from West to East then
    # cross section segment and normal to topology segment will point 'away' from each other (ie, < 0).
    # Note that topology segment normal points to the left side (West).
    # Note: Same sort of thing applies when imagining topology going North to South (hence reference to left/right instead).
    cross_section_left_to_right = pygplates.Vector3D.dot(
            cross_section_segment_direction,
            topology_section_segment.get_great_circle_normal()) < 0

    # If overriding plate on left (West) and cross-section goes left to right (West to East) then
    # subduction dips from left to right along cross-section (same for overriding on right and cross-section going right to left).
    # The other two cases (out of four cases total) result in subduction dipping right to left along cross-section.
    subduction_polarity = shared_subsegment.get_feature().get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)
    if ((subduction_polarity == 'Left' and cross_section_left_to_right) or
        (subduction_polarity == 'Right' and not cross_section_left_to_right)):
        cross_section_dips_left_to_right = True
    else:
        # NOTE: We'll also get here if (subduction_polarity == 'Unknown').
        cross_section_dips_left_to_right = False

    return cross_section_dips_left_to_right


def topo2moho(topo_profile,ref_depth=20000,rhoM = 3300.,rhoC = 2700.):
    # TODO handle both air-loaded and water-loaded

    base_surface = (topo_profile*rhoC)/(rhoM-rhoC)
    moho_depth = -base_surface-ref_depth

    return moho_depth



########################
# PALEOBATHYMETRY
def find_distance_to_nearest_ridge(resolved_topologies,shared_boundary_sections,
                                   point_features,fill_value=5000.):

    all_point_distance_to_ridge = []
    all_point_lats = []
    all_point_lons = []

    for topology in resolved_topologies:
        plate_id = topology.get_resolved_feature().get_reconstruction_plate_id()
        print('Generating distances for Plate %d ...' % plate_id)

        # Section to isolate the mid-ocean ridge segments that bound the current plate
        mid_ocean_ridges_on_plate = []
        for shared_boundary_section in shared_boundary_sections:

            if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('MidOceanRidge'):
                for shared_subsegment in shared_boundary_section.get_shared_sub_segments():
                    sharing_resolved_topologies = shared_subsegment.get_sharing_resolved_topologies()
                    for resolved_polygon in sharing_resolved_topologies:
                        if resolved_polygon.get_feature().get_reconstruction_plate_id() == plate_id:
                            mid_ocean_ridges_on_plate.append(shared_subsegment.get_resolved_geometry())

        point_distance_to_ridge = []
        point_lats = []
        point_lons = []

        for point_feature in point_features:

            for points in point_feature.get_geometries():
                for point in points:

                    if topology.get_resolved_geometry().is_point_in_polygon(point):

                        if len(mid_ocean_ridges_on_plate)>0:

                            min_distance_to_ridge = None

                            for ridge in mid_ocean_ridges_on_plate:
                                distance_to_ridge = pygplates.GeometryOnSphere.distance(point,ridge,min_distance_to_ridge)

                                if distance_to_ridge is not None:
                                    min_distance_to_ridge = distance_to_ridge

                            point_distance_to_ridge.append(min_distance_to_ridge*pygplates.Earth.mean_radius_in_kms)
                            point_lats.append(point.to_lat_lon()[0])
                            point_lons.append(point.to_lat_lon()[1])

                        else:

                            point_distance_to_ridge.append(fill_value)
                            point_lats.append(point.to_lat_lon()[0])
                            point_lons.append(point.to_lat_lon()[1])

        all_point_distance_to_ridge.extend(point_distance_to_ridge)
        all_point_lats.extend(point_lats)
        all_point_lons.extend(point_lons)


    return all_point_lons,all_point_lats,all_point_distance_to_ridge


#
def age2depth(age_array,model='GDH1'):

    if model is 'GDH1':
        paleodepth = 2600. + 365. * np.sqrt(age_array)
        paleodepth[age_array>=20.] = 5651 - 2473*np.exp(-0.0278*age_array[age_array>=20.])
        paleodepth = -paleodepth

    elif model is 'Crosby':
        paleodepth = 2652. + (324. * np.sqrt(age_array))
        paleodepth[age_array>75.] = 5028. + 5.26*age_array[age_array>75.] - 250.*np.sin((age_array[age_array>75.]-75.)/30.)
        paleodepth[age_array>160.] = 5750.
        paleodepth = -paleodepth

    else:
        print('unknown depth model')

    return paleodepth


def paleobathymetry_from_topologies(resolved_topologies,shared_boundary_sections,
                                    deep_ocean_features,
                                    model='GDH1',half_spreading_rate=50.):

    # Approximation of paleobathymetry based on distance to MORs
    # given some resolved topologies, and some point features (typically in the deep ocean),
    # calculates the distance of each point to the nearest mid-ocean ridge segment that
    # forms part of the boundary that the point is located within - then, determines
    # the implied age assuming a constant spreading rate and given age-depth model

    pX,pY,pZ = find_distance_to_nearest_ridge(resolved_topologies,shared_boundary_sections,
                                              deep_ocean_features)

    age = np.array(pZ) / half_spreading_rate

    pdepth = age2depth(age,model=model)

    pdepth_points = []
    for (lon,lat,depth) in zip(pX,pY,pdepth):
        point_feature = pygplates.Feature()
        point_feature.set_geometry(pygplates.PointOnSphere(lat,lon))
        point_feature.set_shapefile_attribute('depth',np.float(depth))
        pdepth_points.append(point_feature)

    return pdepth_points


########################
# PLOTTING
def paleogeography_points_basemap(pg_points,env_color_dict,fill_color='darkblue',markersize=2,alpha=1):
    m = Basemap(projection='robin', lon_0=0, resolution='c')
    m.drawmapboundary(fill_color='white')
    for feature in pg_points:
        env = feature.get_shapefile_attribute('Layer')
        if env is not None:
            color=env_color_dict[env]
        else:
            color=fill_color

        for geometry in feature.get_geometries():
            x,y = m(geometry.to_lat_lon_array()[:,1],geometry.to_lat_lon_array()[:,0])
            plt.plot(x,y,'.',color=color,markersize=markersize)

    return m


def paleogeography_cross_section(ProfileX_kms,topo_profile,moho_profile,
                                 subduction_intersections,ridge_intersections,
                                 vertical_exaggeration=20.):

    plt.plot(ProfileX_kms,topo_profile,'k')
    plt.plot(ProfileX_kms,moho_profile,'r')
    #plt.plot([0,ProfileX_kms[-1]],[0,0],'lightblue',linewidth=3,zorder=1)

    plt.fill_between(ProfileX_kms,topo_profile,moho_profile,color='pink',zorder=2)
    plt.fill_between(ProfileX_kms,0,-7000,color='lightblue')
    plt.fill_between(ProfileX_kms,-7000,-1e7,color='magenta')

    for point in subduction_intersections:
        plt.arrow(point[1],5000, 0.0, -4000, fc="b", ec="b",head_width=40, head_length=1000, linewidth=5,zorder=2)
        if point[2]:
            plt.plot([point[1]+25,point[1]-250],[-8000,-50000],linewidth=12,color='pink',zorder=1)
        else:
            plt.plot([point[1]-25,point[1]+250],[-8000,-50000],linewidth=12,color='pink',zorder=1)
    for point in ridge_intersections:
        plt.arrow(point[1],5000, 0.0, -4000, fc="r", ec="r",head_width=40, head_length=1000, linewidth=5,zorder=2)
    plt.gca().axis('tight')
    plt.gca().set_aspect(vertical_exaggeration/1000.)  # 1000 because intended units are km for distance, but meters for depth
    plt.ylim(-65000,5000)


def paleo_age_grid_cross_section(ProfileX_kms, profile_plate_ids, seafloor_age_profile,
                                 subduction_intersections, daspect = 50, smoothing_iterations=20,
                                 age_min = -50, age_max = 250, cmap=plt.cm.plasma_r):

    seafloor_depth_profile = age2depth(seafloor_age_profile)/1000

    subduction_indices = []
    for point in subduction_intersections:

        # the index will always be the nearest profile point on the updip side
        temp_array = np.array(point[1] - ProfileX_kms)
        if point[2]:
            temp_array = temp_array*-1.
        temp_array[temp_array<0] = 9e20

        index_of_trench = temp_array.argmin()
        subduction_indices.append(index_of_trench)

        subducting_plate = profile_plate_ids[index_of_trench]

        # get an array with the depths just within one plate polygon
        depths_in_plate = seafloor_depth_profile[np.equal(profile_plate_ids,subducting_plate)]

        #print depths_in_plate, depths_in_plate.shape
        if ~point[2]:
            depths_in_plate = np.flip(depths_in_plate) #[-1:0:-1]
        if np.any(np.isnan(depths_in_plate)):
            index = np.where(~np.isnan(depths_in_plate))[0].min()

            if ~point[2]:
                seafloor_depth_profile[index_of_trench-index:index_of_trench+1] = depths_in_plate[index]
            else:
                seafloor_depth_profile[index_of_trench:index_of_trench+index] = depths_in_plate[index]

    # index
    land_ocean_index = ~np.isnan(seafloor_depth_profile)

    smooth_topo = profile_smoothing(seafloor_depth_profile, n_iter=smoothing_iterations)
    moho = topo2moho(smooth_topo*1000,ref_depth=22000, rhoC=2200)/1000

    moho[land_ocean_index] = seafloor_depth_profile[land_ocean_index]-6


    # set up coloured ocean crust cells along profiles
    xgrid = np.vstack((ProfileX_kms,ProfileX_kms))
    zgrid = np.vstack((smooth_topo,moho))
    tmp=np.vstack((seafloor_age_profile,seafloor_age_profile))

    # PLOTTING
    norm = matplotlib.colors.Normalize(vmin=age_min, vmax=age_max)

    plt.figure(figsize=(20,10))

    # first, fill in all crust with a grey background
    plt.fill_between(ProfileX_kms,smooth_topo,moho,color='grey')

    for point,subduction_index in zip(subduction_intersections,subduction_indices):
        trench_depth = seafloor_depth_profile[subduction_index]
        if point[2]:
            slab_top_X = np.array([point[1],(point[1]-(10*daspect))])
            slab_top_Y = np.array([trench_depth,-30])
            slab_base_X = slab_top_X+150
            slab_base_Y = slab_top_Y-6
        else:
            slab_top_X = np.array([point[1],(point[1]+(10*daspect))])
            slab_top_Y = np.array([trench_depth,-30])
            slab_base_X = slab_top_X-150
            slab_base_Y = slab_top_Y-6

        plt.fill_betweenx(slab_top_Y,slab_top_X,slab_base_X,color=cmap(norm(seafloor_age_profile[subduction_index])),zorder=5)

        plt.pcolormesh(xgrid,zgrid,tmp,cmap=plt.cm.inferno_r,vmin=-10,vmax=150)


    plt.pcolormesh(xgrid,zgrid,tmp,cmap=cmap,vmin=age_min,vmax=age_max)

    plt.ylim(-50,5)
    plt.xlim((ProfileX_kms.min(),ProfileX_kms.max()))
    plt.gca().set_aspect(daspect)
    #plt.show()


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def profile_smoothing(profileZ,n_iter=5):

    is_ocean_depth = np.copy(profileZ)
    for i in range(n_iter):
        #print i
        profileZ[np.isnan(profileZ)] = 2.0
        profileZ = smooth(profileZ,3)
        profileZ[~np.isnan(is_ocean_depth)] = is_ocean_depth[~np.isnan(is_ocean_depth)]

    return profileZ
