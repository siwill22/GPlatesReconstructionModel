import pygplates
import litho1pt0 as litho
from pprint import pprint 
#import requests
#from io import StringIO

from reconstruction_classes import *

class litho1_scalar_coverage(object):

    def __init__(self, distribution_type='healpix', N=32):
        self.points = PointDistributionOnSphere(distribution_type,N)
        self.litho = litho
        self.layer_keys = litho.l1_layer_decode.items()
        self.value_keys = litho.l1_data_decode.items()

    def to_scalar_coverage(self, filename, layer_names='All'):
        
        if layer_names is 'All':
            layer_names = [name[0] for name in self.layer_keys]

        scalar_coverages = {}
        for layer_name in layer_names:
            
            layerZ = litho.layer_depth(self.points.latitude, self.points.longitude, layer_name)
            scalar_coverages[pygplates.ScalarType.create_gpml(layer_name)] = layerZ

        ct_feature = pygplates.Feature()
        ct_feature.set_geometry((self.points.multipoint,scalar_coverages))
        ct_feature.set_name('litho1.0 layers')

        pygplates.FeatureCollection(ct_feature).write(filename)
