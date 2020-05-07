import pygplates
import numpy as np
from .sphere_tools import healpix_mesh

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
