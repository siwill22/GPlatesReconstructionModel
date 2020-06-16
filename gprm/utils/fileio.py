import numpy as np
import xarray as xr


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


def write_netcdf_grid(filename, x, y, z, xname='x', yname='y', zname='z', format='NETCDF4'):
    ds = xr.DataArray(z, coords=[(yname,y), (xname,x)], name=zname)
    ds.to_netcdf(filename, format=format)
    ds.close()

