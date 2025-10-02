'''
MIT License

Copyright (c) 2017-2023 Simon Williams

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
import geopandas as _gpd
#import os as _os
#import collections
#import xarray as _xr

def fetch_GlobalTectonicMap(load=True):
    """
    Polygons from Hasterok et al (2022) ESR
    
    """

    # NB This is a specific commit from the Github repo from around the time the files were commited to zenodo
    # BUT avoids downloading ~1GB of data not needed here 
    
    '''
    fnames = _retrieve(
        url="https://github.com/dhasterok/global_tectonics/archive/ac1b9af78122d6f673f53f2eb4eeb2783b887f05.zip",
        known_hash="sha256:50ffc7f999a630dc8b8923d45dd9f1b493280528f999a4deeb2546edf96a5ed8", # Not sure why SHA of filename doesn't work? 
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='GlobalTectonics'),
    )

    dirname = '{:s}/GlobalTectonics/global_tectonics-ac1b9af78122d6f673f53f2eb4eeb2783b887f05/plates&provinces/shp/'.format(str(_os_cache('gprm')))

    return _gpd.read_file('{:s}/global_gprv.shp'.format(dirname)).set_crs('EPSG:4326')
    '''

    # Slightly older version, but includes the shapefile containing ages
    fnames = _retrieve(
        url="https://github.com/dhasterok/global_tectonics/archive/766e485af4b63c34c88a555621541f64fd7e68d2.zip",
        known_hash="sha256:bcf565c9eba08c7e684de5e88f77a4dd496e3b960a0e45cbf2b7dafde10612d0", # Not sure why SHA of filename doesn't work? 
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='GlobalTectonics'),
    )

    dirname = '{:s}/GlobalTectonics/global_tectonics-766e485af4b63c34c88a555621541f64fd7e68d2/plates&provinces/'.format(str(_os_cache('gprm')))

    if load:
        return _gpd.read_file('{:s}/global_gprv_wage.shp'.format(dirname)).set_crs('EPSG:4326')
    else:
        return '{:s}/global_gprv_wage.shp'.format(dirname)


def fetch_SurfaceGeology(load=True):

    fnames = _retrieve(
        url="https://github.com/siwill22/global-geology/raw/d213c8696ad21bce851a2ce0a4d8351f62187503/global_geology_shapefile.zip",
        known_hash="sha256:a366392c9f4116e1f55b8d2ff05ab70cfd299d46b6a8904b0f21eb3e009053a1", 
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='SurfaceGeology'),
    )

    dirname = '{:s}/SurfaceGeology/'.format(str(_os_cache('gprm')))

    if load:
        return _gpd.read_file('{:s}/nrcan_geology_with_ages.shp'.format(dirname))
    else:
        return '{:s}/nrcan_geology_with_ages.shp'.format(dirname)