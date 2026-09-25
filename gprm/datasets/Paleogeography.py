"""
Loaders for paleogeographic raster reconstructions.

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
"""

from pooch import os_cache as _os_cache
from ._fetch import retrieve as _retrieve
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
        raise ValueError('Spacing for source grids must be either 01d (for 1 degree version) or 06m (for 6 minute version)')


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


def fetch_LiHu2022(return_xarray=False):
    """
    Load the Li et al. (2022) high-resolution climate simulations for the past 540 Myr.

    Li, X., Hu, Y. et al. (2022). A high-resolution climate simulation dataset for the past
    540 million years. Scientific Data 9, 371.
    Data: figshare, doi:10.6084/m9.figshare.19920662.v1 (CC BY 4.0).

    55 CESM1.2.2 snapshot simulations on Scotese paleogeography, all in one netCDF file
    (754 MB). Downloaded once and cached; the checksum is figshare's own published md5.

    :param return_xarray: if True, return the dataset opened lazily with xarray instead of
        the path to the cached file.
    :returns: path to the cached netCDF file, or an xarray.Dataset.
    """
    fname = _retrieve(
        url="https://ndownloader.figshare.com/files/35396138",
        known_hash="md5:0cc0b557e6337369131c82eb8ca568e8",
        fname="High_Resolution_Climate_Simulation_Dataset_540_Myr.nc",
        downloader=_HTTPDownloader(progressbar=True),
        path=_os.path.join(str(_os_cache('gprm')), 'LiHu2022'),
    )
    if return_xarray:
        return _xr.open_dataset(fname, decode_times=False)
    return fname


def fetch_Valdes2021():
    """
    Load the annual-mean atmosphere fields of the Valdes, Scotese & Lunt (2021) BRIDGE
    simulations: the HadCM3 `scotese_02` run sequence, 109 time slices from 541 Ma to present.

    Valdes, P.J., Scotese, C.R. & Lunt, D.J. (2021). Deep ocean temperatures through time.
    Climate of the Past 17, 1483-1506, doi:10.5194/cp-17-1483-2021.
    Data: https://www.paleo.bristol.ac.uk/ummodel/data/<run>/climate/ (public, no login).

    One file per run, <run>a.pdclann.nc (~1.3 MB each, 140 MB in total), downloaded on first
    use and checked against checksums pinned in gprm. Among its 39 variables are
    ``temp_mm_1_5m`` (air temperature at 1.5 m, K) and ``precip_mm_srf`` (precipitation,
    kg m-2 s-1 despite its units attribute). There is no single total-evaporation field:
    evaporation is split across ``evapsea_mm_srf``, ``soilEvap_mm_srf``, ``canopyEvap_mm_can``,
    ``transpiration_mm_srf`` and ``srfSublim_mm_srf``.

    Only the annual atmosphere mean is fetched. The same server also has monthly means and
    two ocean streams, which are not pinned here.

    :returns: OrderedDict mapping age (Ma) to the path of that run's cached file, youngest first.
    """
    from ._valdes2021_runs import RUNS

    downloader = _HTTPDownloader(progressbar=False, headers={'User-Agent': 'gprm (pooch)'})
    cache = _os.path.join(str(_os_cache('gprm')), 'Valdes2021')

    raster_dict = {}
    for index, (run, age, sha256) in enumerate(RUNS):
        name = '{:s}a.pdclann.nc'.format(run)
        raster_dict[float(age)] = _retrieve(
            url='https://www.paleo.bristol.ac.uk/ummodel/data/{:s}/climate/{:s}'.format(run, name),
            known_hash='sha256:{:s}'.format(sha256),
            # Prefixed with the index: run codes differ only in case, and would collide on
            # a case-insensitive filesystem.
            fname='{:03d}_{:s}'.format(index, name),
            downloader=downloader,
            path=cache,
        )

    return collections.OrderedDict(sorted(raster_dict.items()))

