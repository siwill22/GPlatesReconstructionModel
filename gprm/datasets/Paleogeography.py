from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
#import pandas as _pd
#import geopandas as _gpd
import os as _os
import collections

def fetch_Paleomap(resolution='01d'):
    """
    PaleoDEM rasters from Scotese and Wright (2018)
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
                raster_dict[float(file.split('_')[-1][:-5])] = '{:s}/{:s}'.format(dirname,file)

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
                raster_dict[float(file.split('_')[-1][:-5])] = '{:s}/{:s}'.format(dirname,file)

        ordered_raster_dict = collections.OrderedDict(sorted(raster_dict.items()))

        return ordered_raster_dict

    else:
        ValueError('Spacing for source grids must be either 01d (for 1 degree version) or 06m (for 6 minute version)')

