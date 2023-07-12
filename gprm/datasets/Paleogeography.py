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

from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
#import pandas as _pd
#import geopandas as _gpd
import os as _os
import collections
import xarray as _xr

def fetch_Paleomap(resolution='01d'):
    """
    PaleoDEM rasters from Scotese and Wright (2018)
    
    resolution can be '01d' (default) or '06m'
    """

    if resolution=='01d':
        fnames = _retrieve(
            url="https://zenodo.org/record/5460860/files/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc.zip?download=1",
            known_hash="md5:77147998623ab039d86ff3e0b5e40344",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='Paleomap_01d'),
        )

        dirname = '{:s}/Paleomap_01d/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2'.format(fnames[0].split('Paleomap_01d')[0])
        #dirname = '{:s}/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2'.format(_os.path.split(fnames[0])[0])

        # if downloading for first time, remove the unwanted cache files
        for file in _os.listdir(dirname):
            if file.endswith(".cache"):
                _os.remove('{:s}/{:s}'.format(dirname,file))

        raster_dict = {}
        for file in _os.listdir(dirname):
            if file.endswith(".nc"):
                # Replace whitespace with underscore to help pygmt plotting
                if ' ' in file:
                    _os.rename('{:s}/{:s}'.format(dirname,file), '{:s}/{:s}'.format(dirname,file.replace(' ','_')))
                raster_dict[float(file.split('_')[-1][:-5])] = '{:s}/{:s}'.format(dirname,file.replace(' ','_'))

        ordered_raster_dict = collections.OrderedDict(sorted(raster_dict.items()))

        return ordered_raster_dict


    elif resolution=='06m':
        fnames = _retrieve(
            url="https://zenodo.org/record/5460860/files/Scotese_Wright_2018_Maps_1-88_6minX6min_PaleoDEMS_nc.zip?download=1",
            known_hash="md5:89eb50d8645707ab221b023078535bda",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='Paleomap_06m'),
        )

        dirname = '{:s}/Paleomap_06m/Scotese_Wright_2018_Maps_1-88_6minX6min_PaleoDEMS_nc'.format(fnames[0].split('Paleomap_06m')[0])
        #dirname = '{:s}/Scotese_Wright_2018_Maps_1-88_6minX6min_PaleoDEMS_nc'.format(_os.path.split(fnames[0])[0])

        raster_dict = {}
        for file in _os.listdir(dirname):
            if file.endswith(".nc"):
                # Replace whitespace with underscore to help pygmt plotting
                if ' ' in file:
                    _os.rename('{:s}/{:s}'.format(dirname,file), '{:s}/{:s}'.format(dirname,file.replace(' ','_')))
                raster_dict[float(file.split('_')[-1][:-5])] = '{:s}/{:s}'.format(dirname,file.replace(' ','_'))

        ordered_raster_dict = collections.OrderedDict(sorted(raster_dict.items()))

        return ordered_raster_dict

    else:
        ValueError('Spacing for source grids must be either 01d (for 1 degree version) or 06m (for 6 minute version)')


def fetch_Pohl2022(return_xarray=False, value=None):
    """
    valid value names are:
    'area' [grid point area]
    'evp' [evaporation]
    'koppen' [Koppen-Geiger climatic zones]
    'PmE' [precipitation minus evaporation balance]
    'precip' [precipitation]
    'rnf' [runoff]
    'topo' [topography]
    'tssub1' [top soil layer temp]
    """

    fnames = _retrieve(
            url="https://zenodo.org/record/6620748/files/All_NC_files.zip?download=1",
            known_hash="md5:b0b8bf04647f3f084d282d106fa52a20",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='Pohl2022'),
        )

    dirname = '{:s}/Pohl2022/All_NC_files'.format(fnames[0].split('Pohl2022')[0])

    raster_dict = {}
    for file in _os.listdir(dirname):
        if file.endswith(".nc"):
            raster_dict[float(file.split('Ma')[0])] = '{:s}/{:s}'.format(dirname,file)

    ordered_raster_dict = collections.OrderedDict(sorted(raster_dict.items()))

    return ordered_raster_dict

    #return xr.DataArray()

