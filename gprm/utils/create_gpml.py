import pygplates
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString, Polygon
from .sphere import healpix_mesh
import os, tempfile


def create_gpml_crustal_thickness(longitude_array,latitude_array,thickness,filename=None):

    multi_point = pygplates.MultiPointOnSphere(zip(latitude_array,longitude_array))

    scalar_coverages = {
        pygplates.ScalarType.create_gpml('CrustalThickness'): thickness}

    ct_feature = pygplates.Feature()
    ct_feature.set_geometry((multi_point,scalar_coverages))
    ct_feature.set_name('Crustal Thickness')

    output_feature_collection = pygplates.FeatureCollection(ct_feature)

    if filename is not None:
        output_feature_collection.write(filename)
    else:
        return output_feature_collection


def create_gpml_velocity_feature(longitude_array,latitude_array,filename=None,feature_type=None):
# function to make a velocity mesh nodes at an arbitrary set of points defined in Lat
# Long and Lat are assumed to be 1d arrays.

    multi_point = pygplates.MultiPointOnSphere(zip(latitude_array,longitude_array))

    # Create a feature containing the multipoint feature.
    # optionally, define as 'MeshNode' type, so that GPlates will recognise it as a velocity layer
    if feature_type=='MeshNode':
        meshnode_feature = pygplates.Feature(pygplates.FeatureType.create_from_qualified_string('gpml:MeshNode'))
        meshnode_feature.set_name('Velocity Mesh Nodes')
    else:
        meshnode_feature = pygplates.Feature()
        meshnode_feature.set_name('Multipoint Feature')

    meshnode_feature.set_geometry(multi_point)

    output_feature_collection = pygplates.FeatureCollection(meshnode_feature)

    if filename is not None:
        output_feature_collection.write(filename)
    else:
        return output_feature_collection


def create_gpml_healpix_mesh(nSide,filename=None,feature_type=None):

    # call the function to create a healpix array
    longitude_array,latitude_array = healpix_mesh(nSide)

    # call the function to create a multipoint feature, with user-defined type
    output_feature_collection = create_gpml_velocity_feature(longitude_array,latitude_array,filename,feature_type)

    if filename is not None:  # This is superfluous, since file has already been written in previous line???
        output_feature_collection.write(filename)
    else:
        return output_feature_collection


def create_gpml_regular_long_lat_mesh(Sampling=1,filename=None,feature_type=None):

    # call the function to create a healpix array
    longitude_array,latitude_array = np.meshgrid(np.arange(-180.,180.001,Sampling),np.arange(-90.,90.001,Sampling))
    longitude_array = longitude_array.flatten()
    latitude_array = latitude_array.flatten()

    # call the function to create a multipoint feature, with user-defined type
    output_feature_collection = create_gpml_velocity_feature(longitude_array,latitude_array,filename,feature_type)

    if filename is not None:
        output_feature_collection.write(filename)
    else:
        return output_feature_collection


def geometries_to_geodataframe(geometries, geometry_type='polygon'):
    gdf = gpd.GeoDataFrame()
    gdf['geometry'] = None
    for i,geometry in enumerate(geometries):
        if geometry_type in ['PolyLine','Polyline']:
            poly = LineString([tuple(coord) for coord in np.fliplr(geometry)])
        else:
            poly = Polygon([tuple(coord) for coord in np.fliplr(geometry)])
        gdf.loc[i, 'geometry'] = poly

    return gdf


def geodataframe_to_geometries(gdf):
# from a geopandas geodataframe, return a list of the geometries (as type pygplates.GeometryOnSphere)
    
    geometry_list = []
    gdf = gdf.explode()
    for i,row in gdf.iterrows():
        geometry_list.append([(lat,lon) for lat,lon in zip(row.geometry.xy[1], row.geometry.xy[0])])
    return geometry_list 


def gdf2gpml(gdf):
    """ 
    Given a geopandas geodataframe, returns a gplates feature collection
    """

    temporary_file = tempfile.NamedTemporaryFile(delete=True, suffix='.geojson')
    temporary_file.close()

    gdf.to_file(temporary_file.name, driver='GeoJSON')
    feature_collection = pygplates.FeatureCollection(temporary_file.name)

    os.unlink(temporary_file.name)
    
    return feature_collection


def gpml2gdf(features):
    """ 
    Given a gplates feature collection, or a list of features, or a list
    of reconstructed features, returns a geopandas geodataframe containing 
    the same features
    """

    if isinstance(features, pygplates.FeatureCollection):
        pass
    elif isinstance(features, list):
        if isinstance(features[0], pygplates.Feature):
            features = pygplates.FeatureCollection(features)
        elif isinstance(features[0], pygplates.ReconstructedFeatureGeometry):
            features = _reconstructed_features_to_features(features)
            features = pygplates.FeatureCollection(features)
        else:
            raise TypeError('Unexpected list item of type {:s} for gpml2gdf input'.format(type(features[0])))
    else:
        raise TypeError('Unexpected type {:s} for gpml2gdf input'.format(type(features)))


    temporary_file = tempfile.NamedTemporaryFile(delete=True, suffix='.geojson')
    temporary_file.close()

    features.write(temporary_file.name)
    gdf = gpd.read_file(temporary_file.name)

    gdf['NAME'] = gdf['NAME'].astype(str)

    os.unlink(temporary_file.name)
    
    return gdf

def _reconstructed_features_to_features(reconstructed_features):

    features = []
    for feature in reconstructed_features:
        f = feature.get_feature().clone()
        f.set_geometry(feature.get_reconstructed_geometry())
        features.append(f)

    return features
