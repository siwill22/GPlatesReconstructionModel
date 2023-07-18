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

def fetch_GlobalTectonicMap():
    """
    Polygons from Hasterok et al (2022) ESR
    
    """

    # NB This is a specific commit from the Github repo from around the time the files were commited to zenodo
    # BUT avoids downloading ~1GB of data not needed here 
    fnames = _retrieve(
        url="https://github.com/dhasterok/global_tectonics/archive/ac1b9af78122d6f673f53f2eb4eeb2783b887f05.zip",
        known_hash="sha256:50ffc7f999a630dc8b8923d45dd9f1b493280528f999a4deeb2546edf96a5ed8", # Not sure why SHA of filename doesn't work? 
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='GlobalTectonics'),
    )

    dirname = '{:s}/GlobalTectonics/global_tectonics-ac1b9af78122d6f673f53f2eb4eeb2783b887f05/plates&provinces/shp/'.format(str(_os_cache('gprm')))

    return _gpd.read_file('{:s}/global_gprv.shp'.format(dirname))

