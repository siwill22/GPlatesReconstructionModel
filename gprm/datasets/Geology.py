"""
Loaders for global tectonic map and surface geology datasets.

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
"""

from pooch import os_cache as _os_cache
from ._fetch import retrieve as _retrieve
from ._ages import stamp as _stamp
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
        return _stamp(_gpd.read_file('{:s}/global_gprv_wage.shp'.format(dirname)).set_crs('EPSG:4326'),
                      'Geology.GlobalTectonicMap')
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
        return _stamp(_gpd.read_file('{:s}/nrcan_geology_with_ages.shp'.format(dirname)),
                      'Geology.SurfaceGeology')
    else:
        return '{:s}/nrcan_geology_with_ages.shp'.format(dirname)

# Artemieva (2006) TC1, as distributed in C. Brenhin Keller's StatisticalGeochemistry repository
# (MIT licence, Copyright (c) 2016 C. Brenhin Keller), pinned to one commit
_TC1_URL = ('https://raw.githubusercontent.com/brenhinkeller/StatisticalGeochemistry/'
            '040aaace092cd599fc1e1097c93bf604b1b71e88/resources/tc1/')

# legend of tc1age.png: RGB colour, then the age bin it stands for (midpoint, min, max; Ma).
# Taken from findtc1age.m in the same repository. White is no data (oceans).
_TC1_LEGEND = [
    ((171,  45,  22),   25,    0,   50),
    ((230, 129,  55),  150,   50,  250),
    ((255, 248,  54),  395,  250,  540),
    ((171, 245,  49),  695,  540,  850),
    (( 94, 183,  47),  975,  850, 1100),
    ((114, 250, 253), 1400, 1100, 1700),
    (( 96, 138, 249), 2100, 1700, 2500),
    (( 15,  45, 138), 2750, 2500, 3000),
    ((120,  43, 183), 3250, 3000, 3500),
]


def fetch_TC1(load=True):
    """
    Tectono-thermal age of the continental lithosphere from the TC1 model of Artemieva (2006).

    Artemieva, I. M. (2006). Global 1x1 thermal model TC1 for the continental lithosphere:
    implications for lithosphere secular evolution. Tectonophysics, 416, 245-277,
    doi:10.1016/j.tecto.2005.11.022

    The ages come from the map image tc1age.png in C. Brenhin Keller's StatisticalGeochemistry
    repository (https://github.com/brenhinkeller/StatisticalGeochemistry, resources/tc1,
    MIT licence, Copyright (c) 2016 C. Brenhin Keller), downloaded from a fixed commit.
    The image is a global lon/lat map at 1/3 degree, and each pixel is decoded to the age bin
    of the nearest legend colour, as findtc1age.m in that repository does. The image was drawn
    from the 1 degree model, so bin boundaries are only good to about a pixel or two, and a few
    pixels on boundaries between two colours may decode to a third bin.

    Returns an xarray Dataset on cell centres (lat and lon ascending, -180 to 180) with
    'age' (bin midpoint, Ma), 'age_min' and 'age_max'; NaN where there is no data (oceans).
    With load=False, returns the path to the cached image.
    """
    fname = _retrieve(
        url=_TC1_URL + 'tc1age.png',
        known_hash="sha256:190124c07983c3c90768afa67391dc8ca2c27ed30a509bd225377fe84c58a208",
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        fname='TC1_tc1age.png',
    )
    if not load:
        return fname

    import numpy as np
    import xarray as xr
    from PIL import Image

    rgb = np.asarray(Image.open(fname).convert('RGB')).astype(int)     # (540, 1080, 3), north at the top
    colours = np.array([c for c, *_ in _TC1_LEGEND] + [(255, 255, 255)])
    nearest = ((rgb[:, :, None, :] - colours[None, None])**2).sum(-1).argmin(-1)

    n_lat, n_lon = nearest.shape
    lat = 90. - (np.arange(n_lat) + 0.5)*180./n_lat
    lon = -180. + (np.arange(n_lon) + 0.5)*360./n_lon

    data_vars = {}
    for i, name in enumerate(['age', 'age_min', 'age_max']):
        values = np.array([bins[i] for _, *bins in _TC1_LEGEND] + [np.nan], dtype=float)
        data_vars[name] = (('lat', 'lon'), values[nearest][::-1])
    ds = xr.Dataset(data_vars, coords={'lat': lat[::-1], 'lon': lon})
    ds['age'].attrs = {'long_name': 'tectono-thermal age, bin midpoint', 'units': 'Ma'}
    ds['age_min'].attrs = {'long_name': 'tectono-thermal age, bin minimum', 'units': 'Ma'}
    ds['age_max'].attrs = {'long_name': 'tectono-thermal age, bin maximum', 'units': 'Ma'}
    ds.attrs = {'source': 'Artemieva (2006) TC1, doi:10.1016/j.tecto.2005.11.022; age map image '
                          'from github.com/brenhinkeller/StatisticalGeochemistry (MIT licence, '
                          'Copyright (c) 2016 C. Brenhin Keller), commit 040aaac'}
    return ds


def fetch_TC1Thermal(load=True):
    """
    The 1x1 degree thermal model TC1 for the continental lithosphere (Artemieva, 2006).

    Artemieva, I. M. (2006). Global 1x1 thermal model TC1 for the continental lithosphere:
    implications for lithosphere secular evolution. Tectonophysics, 416, 245-277,
    doi:10.1016/j.tecto.2005.11.022

    The file TC1-thermal-0606-1x1.dat from C. Brenhin Keller's StatisticalGeochemistry
    repository (https://github.com/brenhinkeller/StatisticalGeochemistry, resources/tc1,
    MIT licence, Copyright (c) 2016 C. Brenhin Keller), downloaded from a fixed commit.

    Returns a DataFrame with one row per continental 1 degree node: lon, lat, z1300 (thermal
    lithosphere thickness, depth to the 1300 C isotherm, km), T50 (temperature at 50 km, C),
    gradTlith (mean lithospheric temperature gradient, C/km), z550 (depth to 550 C, km) and
    region (the regional geotherm used). With load=False, returns the path to the cached file.
    """
    fname = _retrieve(
        url=_TC1_URL + 'TC1-thermal-0606-1x1.dat',
        known_hash="sha256:a624560635c64f1e196384d0f566e7b207541b8111910760eef53bfe36c16b1d",
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        fname='TC1-thermal-0606-1x1.dat',
    )
    if not load:
        return fname

    import pandas as pd
    return pd.read_csv(fname, dtype={'region': str})
