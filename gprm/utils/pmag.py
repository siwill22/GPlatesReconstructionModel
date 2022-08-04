from pandas.core.indexing import is_nested_tuple
import pygplates
import numpy as np
import pandas as _pd
import geopandas as _gpd
#from shapely.geometry.point import Point
import os


def vgp_to_dataframe(vgp_feature_collection, as_geodataframe=True, return_feature_id=False):
    '''
    Read a gpml file containing virtual geomagnetic poles to a (geo)pandas dataframe
    options: 
    input can be either a pygplates FeatureCollection or a 
    as_geodataframe [Bool]: If True (default), returns geopandas geodataframe. If False, 
                            returns a pandas dataframe
    return_feature_id [Bool]: If True, include a column with GPlates unique feature ids 
                              (default False)                       
    TODO allow .vgp files to also be read
    TODO handle cases where mean age or age ranges could be inferred from one another
    '''
    
    if os.path.isfile(vgp_feature_collection):
        feature_collection = pygplates.FeatureCollection(vgp_feature_collection)
    elif isinstance(vgp_feature_collection, pygplates.FeatureCollection):
        feature_collection = vgp_feature_collection
    else:
        raise ValueError('Unable to load {:s} as vgp input'.format(vgp_feature_collection))


    DataFrameTemplate = ['AverageSampleSiteLongitude','AverageSampleSiteLatitude',
                         'Name','Description','PoleLongitude','PoleLatitude','PoleA95',
                         'AverageAge','MaximumAge','MinimumAge','PLATEID1']
    if return_feature_id:
        DataFrameTemplate.append('Feature_ID')
    
    # Get attribute (other than coordinate) names from first feature
    for feature in feature_collection: 
        if feature.get_shapefile_attributes() is not None:
            for attribute in feature.get_shapefile_attributes():
                DataFrameTemplate.append(attribute) 
            break

    vgps = []
    for feature in feature_collection:
        vgp = []
        
        # default columns
        average_sample_site_position = feature.get(pygplates.PropertyName.create_gpml('averageSampleSitePosition')).get_value().get_geometry().to_lat_lon()
        #print(type(average_sample_site_position.to_lat_lon()))
        vgp.append(float(average_sample_site_position[1]))
        vgp.append(float(average_sample_site_position[0]))
        vgp.append(str(feature.get_name()))
        vgp.append(str(feature.get_description()))
        pole_position = feature.get_geometry().to_lat_lon()
        vgp.append(float(pole_position[1]))
        vgp.append(float(pole_position[0]))
        if feature.get(pygplates.PropertyName.create_gpml('poleA95')):
            vgp.append(float(feature.get(pygplates.PropertyName.create_gpml('poleA95')).get_value().get_double()))
        if feature.get(pygplates.PropertyName.create_gpml('averageAge')):
            vgp.append(float(feature.get(pygplates.PropertyName.create_gpml('averageAge')).get_value().get_double()))
        feature_valid_time = feature.get_valid_time()
        vgp.append(float(feature_valid_time[0]))
        vgp.append(float(feature_valid_time[1]))
        vgp.append(int(feature.get_reconstruction_plate_id()))
        
        # optional
        if return_feature_id:
            vgp.append(str(feature.get_feature_id()))
            
        # depending on input file
        if feature.get_shapefile_attributes() is not None:
            for attribute in feature.get_shapefile_attributes():
                vgp.append(feature.get_shapefile_attribute(attribute))
            
        vgps.append(vgp)
        
    if as_geodataframe:
        df = _pd.DataFrame(vgps,columns=DataFrameTemplate)
        return _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.AverageSampleSiteLongitude, 
                                                                  df.AverageSampleSiteLatitude), crs=4326)
    
    else:
        return _pd.DataFrame(vgps,columns=DataFrameTemplate)


def assign_plate_ids(vgps, reconstruction_model):
    '''
    assign plate ids to Virtual Geomagnetic Poles (vgps), which is a special case
    of plate partitioning where we must use the 'AverageSampleSitePosition' rather than
    the feature geometry
    The input type can be a geodataframe or a pygplates FeatureCollection. The output type
    will match the input 
    '''

    plate_partitioner = pygplates.PlatePartitioner(reconstruction_model.static_polygons, 
                                                   reconstruction_model.rotation_model)

    if isinstance(vgps, _gpd.GeoDataFrame):
        partition_plate_ids = []
        for i,row in vgps.iterrows():
            partition_polygon = plate_partitioner.partition_point(pygplates.PointOnSphere(row.geometry.y,
                                                                                          row.geometry.x))
            partition_plate_ids.append(partition_polygon.get_feature().get_reconstruction_plate_id())

        vgps['PLATEID1'] = partition_plate_ids

        return vgps

    elif isinstance(vgps, (pygplates.FeatureCollection, list)):
        if isinstance(vgps, list):
            vgps = pygplates.FeatureCollection(vgps)
        partitioned_vgps = []
        for vgp in vgps:
            partition_polygon = plate_partitioner.partition_point(vgp.get(pygplates.PropertyName.gpml_average_sample_site_position).get_value().get_geometry())
            vgp.set_reconstruction_plate_id(partition_polygon.get_feature().get_reconstruction_plate_id())
            partitioned_vgps.append(vgp)

        return pygplates.FeatureCollection(partitioned_vgps)

    else:
        raise TypeError('Unexpected type {:} for vgp input'.format(type(vgps)))


def rotate_to_common_reference(vgps, reconstruction_model, reference_plate_id=701):
    '''
    Rotate a collection of vgps to a common reference plate
    '''

    if isinstance(vgps, _gpd.GeoDataFrame):
        rotated_vgps = []
        for i,row in vgps.iterrows():
            vgp_geometry = pygplates.PointOnSphere(row.PoleLatitude,row.PoleLongitude)
            feature_rotation = reconstruction_model.rotation_model.get_rotation(row.AverageAge, 
                                                                                row.PLATEID1, 
                                                                                anchor_plate_id=reference_plate_id)

            reconstructed_geometry = feature_rotation * vgp_geometry
            rotated_vgps.append(reconstructed_geometry.to_lat_lon())

        vgps.PoleLatitude = list(zip(*rotated_vgps))[0]
        vgps.PoleLongitude = list(zip(*rotated_vgps))[1]

        return vgps
        

    elif isinstance(vgps, pygplates.FeatureCollection):
        rotated_vgps = []
        for vgp in vgps:
            feature_rotation = reconstruction_model.rotation_model.get_rotation(vgp.get(pygplates.PropertyName.gpml_average_age).get_value().get_double(), 
                                                                                vgp.get_reconstruction_plate_id(), 
                                                                                anchor_plate_id=reference_plate_id)
            reconstructed_geometry = feature_rotation * vgp.get_geometry()
            vgp.set_geometry(reconstructed_geometry)
            vgp.set_reconstructed_plate_id(reference_plate_id)
            rotated_vgps.append(vgp)

        return pygplates.FeatureCollection(rotated_vgps)


def generate_running_mean_path(vgps,time_list,time_window=20,right=True):

    import pmagpy.ipmag as ipmag

    running_mean_path = []

    if isinstance(vgps, _gpd.GeoDataFrame):
        vgps_df = vgps[['PoleLatitude', 'PoleLongitude', 'AverageAge', 'PoleA95']]
    elif isinstance(vgps, pygplates.FeatureCollection):
        vgp_list = []
        for vgp in vgps:
            vgp_list.append((vgp.get_geometry().to_lat_lon()[0],
                             vgp.get_geometry().to_lat_lon()[1],
                             float(vgp.get(pygplates.PropertyName.create_gpml('averageAge')).get_value().get_double()),
                             float(vgp.get(pygplates.PropertyName.create_gpml('poleA95')).get_value().get_double())))
        vgps_df = _pd.DataFrame(vgp_list, columns=['PoleLatitude', 'PoleLongitude', 'AverageAge','PoleA95'])
    else:
        raise TypeError('Unexpected type {:s} for vgp input'.format(type(vgps)))
        
    for mean_pole_age in time_list:
        if right:
            vgps_window = vgps_df[(vgps_df['AverageAge']>=mean_pole_age-time_window/2.) 
                                & (vgps_df['AverageAge']<=mean_pole_age+time_window/2.)]
        else:
            vgps_window = vgps_df[(vgps_df['AverageAge']>=mean_pole_age-time_window/2.) 
                                & (vgps_df['AverageAge']<mean_pole_age+time_window/2.)]
        
        if vgps_window.empty:
            running_mean_path.append((mean_pole_age,np.nan,np.nan,np.nan,0))

        elif len(vgps_window)==1:
            running_mean_path.append((mean_pole_age, np.array(vgps_window['PoleLongitude'])[0], 
                                      np.array(vgps_window['PoleLatitude'])[0], np.array(vgps_window['PoleA95'])[0], 1))

        else:
            #print(vgps_window)
            mean_pole = ipmag.fisher_mean(np.array(vgps_window.PoleLongitude), 
                                          np.array(vgps_window.PoleLatitude))
            #print(mean_pole)

            running_mean_path.append((mean_pole_age, mean_pole['dec'],
                                      mean_pole['inc'],mean_pole['alpha95'],len(vgps_window)))
        

    return _pd.DataFrame(running_mean_path, columns=['Age','PoleLongitude','PoleLatitude','PoleA95','N'])


def write_vgp_feature(vgp, mapping, half_time_range = 10.):
    '''
    Create a vgp feature from one row of a dataframe
    TODO handle cases where some fields (e.g. description) are not present
    '''
    other_properties = [(pygplates.PropertyName.create_gpml('poleA95'), pygplates.XsDouble(vgp[mapping['PoleA95']])),
                        (pygplates.PropertyName.create_gpml('averageAge'), pygplates.XsDouble(vgp[mapping['AverageAge']]))]
    if 'geometry' in vgp:
        other_properties.append(
            (pygplates.PropertyName.create_gpml('averageSampleSitePosition'),
            pygplates.GmlPoint(pygplates.PointOnSphere([float(float(vgp.geometry.y)), 
                                                        float(float(vgp.geometry.x))])))
        )

    vgpFeature = pygplates.Feature.create_reconstructable_feature(
                 pygplates.FeatureType.create_gpml('VirtualGeomagneticPole'),
                 pygplates.PointOnSphere([vgp[mapping['PoleLatitude']], vgp[mapping['PoleLongitude']]]),
                 name = str(vgp[mapping['Name']]),
                 description = str(vgp[mapping['Description']]),
                 valid_time=(float(vgp[mapping['AverageAge']])+half_time_range, float(vgp[mapping['AverageAge']])-half_time_range),
                 other_properties = other_properties)

    if 'ReconstructionPlateID' in mapping:
        vgpFeature.set_reconstruction_plate_id(int(vgp[mapping['ReconstructionPlateID']]))

    return vgpFeature


def dataframe_to_vgps(gdf, mapping={'Name':'Name',
                                    'Description':'Description',
                                    'PoleLongitude':'PoleLongitude',
                                    'PoleLatitude':'PoleLatitude',
                                    'PoleA95':'PoleA95',
                                    'AverageAge':'AverageAge'}):
    
    vpgFeatureCollection = []

    for i,row in gdf.iterrows():

        vgpFeature = write_vgp_feature(row, mapping)
        
        # Add newly created feature to existing Feature Collection
        vpgFeatureCollection.append(vgpFeature)
    
    return pygplates.FeatureCollection(vpgFeatureCollection)

