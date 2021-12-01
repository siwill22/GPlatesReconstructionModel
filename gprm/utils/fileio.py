'''
MIT License

Copyright (c) 2017-2021 Simon Williams

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import numpy as np
import xarray as xr
import pygplates
import pandas as _pd
import geopandas as _gpd
import os


def write_xyz_file(output_filename, output_data):
    """
    write data arrays to an xyz file
    :param filename: (string) name of output ascii file
    :param output_data: name of array containing data to be written
    """
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def load_netcdf(grdfile,z_field_name='z'):

    ds_disk = xr.open_dataset(grdfile)

    data_array = ds_disk[z_field_name]
    coord_keys = [key for key in data_array.coords.keys()]  # updated for python3 compatibility

    if 'lon' in coord_keys[0].lower():
        latitude_key=1; longitude_key=0
    elif 'x' in coord_keys[0].lower():
        latitude_key=1; longitude_key=0
    else:
        latitude_key=0; longitude_key=1

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

    if gridZ.shape[0]==gridX.shape[0]:
        gridZ = gridZ.T

    return gridX,gridY,gridZ


def write_netcdf_grid(filename, x, y, z, xname='x', yname='y', zname='z', format='NETCDF4_CLASSIC'):
    '''
    write grid to a netcdf file
    '''
    ds = xr.DataArray(z, coords=[(yname,y), (xname,x)], name=zname)
    ds.to_netcdf(filename, format=format)
    ds.close()



def gpml_to_dataframe(feature_collection, as_geodataframe=True):
    '''
    function to read in any gplates-compatible feature collection and 
    place it into a pandas dataframe
    # TODO handle polylines and polygons
    '''

    if os.path.isfile(feature_collection):
        feature_collection = pygplates.FeatureCollection(feature_collection)
    elif isinstance(feature_collection, pygplates.FeatureCollection):
        pass
    else:
        raise ValueError('Unable to load {:s} as vgp input'.format(feature_collection))


    DataFrameTemplate = ['Longitude','Latitude','Name','Description',
                         'MaximumAge','MinimumAge','PlateID']

    # Get attribute (other than coordinate) names from first feature
    for feature in feature_collection: 
        for attribute in feature.get_shapefile_attributes():
            DataFrameTemplate.append(attribute) 
        break

    fs = []
    for feature in feature_collection:
        f = []
        f.append(feature.get_geometry().to_lat_lon()[1])
        f.append(feature.get_geometry().to_lat_lon()[0])
        f.append(str(feature.get_name()))
        f.append(str(feature.get_description()))
        geom = feature.get_geometry().to_lat_lon()
        f.append(float(geom[1]))
        f.append(float(geom[0]))
        feature_valid_time = feature.get_valid_time()
        f.append(float(feature_valid_time[0]))
        f.append(float(feature_valid_time[1]))
        f.append(int(feature.get_reconstruction_plate_id()))
        f.append(str(feature.get_feature_type()))
        f.append(str(feature.get_feature_id()))

        for attribute in feature.get_shapefile_attributes():
            f.append(feature.get_shapefile_attribute(attribute))
        fs.append(f)

    if as_geodataframe:
        df = _pd.DataFrame(fs,columns=DataFrameTemplate)
        return _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.AverageSampleSiteLongitude, 
                                                                  df.AverageSampleSiteLatitude), crs=4326)
    
    else:
        return _pd.DataFrame(fs,columns=DataFrameTemplate)


