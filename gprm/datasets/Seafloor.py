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
from pooch import Untar as _Untar
from pooch import Unzip as _Unzip
import pandas as _pd
import geopandas as _gpd
import os as _os


def MagneticPicks(load=True):
    '''
    Magnetic Picks from the 'Global Seafloor Fabric (and) Magnetic Linations' 
    database, returned as a geopandas dataframe
    Alternatively, select 'load=False' to return filname of '.gmt' file in 
    cache folder

    '''
    fname = _retrieve(
        url="http://www.soest.hawaii.edu/PT/GSFML/ML/DATA/GSFML.global.picks.gmt",
        known_hash="sha256:0895b76597f600a6c6184a7bec0edc0df5ca9234255f3f7bac0fe944317caf65",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
    )
    
    if load:
        return _gpd.read_file(fname)
    else:
        return fname


def SeafloorFabric(feature_type='FZ', load=True):
    '''
    Seafloor fabric from the 'Global Seafloor Fabric (and) Magnetic Linations'
    database, returned as a geopandas dataframe. 
    Alternatively, select 'load=False' to return filname of '.gmt' file in 
    cache folder

    Parameters
    ----------
    feature_type (str), choose from one of:
    
    'FZ': Fracture Zones (Kara Matthews)
    'FZLC': Fracture Zones, Less Certainty (Kara Matthews)
    'DZ': Discordant Zones
    'VANOM': V-Shaped Structures
    'PR': Propagating Ridges
    'ER': Extinct Ridge
    'UNCV': Unclassified V-Anomalies

    Additonal Fracture Zone interpretations:
    'FZ_JW': from Jo Whittaker
    'FZ_RM': from Robert Myhill
    'FZ_MC': from Michael Chandler
    '''
    fnames = _retrieve(
        url="http://www.soest.hawaii.edu/PT/GSFML/SF/DATA/GSFML_SF.tbz",
        known_hash="sha256:e27a73dc544611685144b4587d17f03bde24438ee4646963f10761f8ec2e6036",
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Untar(extract_dir='SeafloorFabric'),
    )
    
    FABRIC_TYPE = {
        "FZ": "GSFML_SF_FZ_KM.gmt",
        "FZLC": "GSFML_SF_FZLC_KM.gmt",
        "UNCV": "GSFML_SF_UNCV_KM.gmt",
        "DZ": "GSFML_SF_DZ_KM.gmt",
        "PR": "GSFML_SF_PR_KM.gmt",
        "VANOM": "GSFML_SF_VANOM_KM.gmt",
        "ER": "GSFML_SF_ER_KM.gmt",
        "FZ_JW": "GSFML_SF_FZ_JW.gmt",
        "FZ_RM": "GSFML_SF_FZ_RM.gmt",
        "FZ_MC": "GSFML_SF_FZ_MC.gmt",
    }

    if feature_type not in FABRIC_TYPE.keys():
        raise ValueError('Unknown feature type {:s}'.format(feature_type))

    for fname in fnames:
        if _os.path.split(fname)[1]==FABRIC_TYPE[feature_type]:
            if load:
                return _gpd.read_file(fname)
            else:
                return fname


def PacificSeamountAges(catalogue='2021', load=True):
    '''
    Pacific Seamount Age compilation from GMT website
    The options for 'catalogue' are:
    2021 [default] - data from Chase and Wessel (2021)
    2013 - data from GMT website, data from Clouard and Bonneville (2005) with updates by Wessel up to 2013
    '''
    if catalogue=='2013':
        fname = _retrieve(
            url="https://www.earthbyte.org/webdav/gmt_mirror/gmt/data/cache/Pacific_Ages.txt",
            known_hash="sha256:8c5e57b478c2c2f5581527c7aea5ef282e976c36c5e00452210885a92e635021",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )
        
        if load:
            df = _pd.read_csv(fname, comment='#', delim_whitespace=True,
                            names=['Long', 'Lat', 'Average_Age_Ma', 'Average_Age_Error_Ma', 'Tag', 'SeamountName', 'SeamountChain'])
            return _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Long, df.Lat))
        else:
            return fname
        
    elif catalogue=='2021':
        fname = _retrieve(
            url="https://zenodo.org/record/6558676/files/Pacific_Hotspot_Trails_Datasets.zip?download=1",
            known_hash="md5:97bec4ebde94694698077eb527bc1ef4",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='PHT2021')
        )
        dirname = '{:s}/PHT2021/'.format(str(_os_cache('gprm')))

        if load:
            df = _pd.read_csv('{:}/PHT2021/PHT2021_raw_data/PHT2021_pacific_ages.txt'.format(dirname), 
                              comment='#', delim_whitespace=True,
                              names=['Long', 'Lat', 'Average_Age_Ma', 'Average_Age_Error_Ma', 'Type', 'Ref', 'SampleName', 'Tag', 'SeamountChain'])
            return _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Long, df.Lat))
        else:
            return fname


def Seamounts(catalogue='KimWessel', load=True):
    '''
    Seamount Census from Kim and Wessel

    '''
    if catalogue in ['KimWessel', 'KW']:
        fname = _retrieve(
            url="http://www.soest.hawaii.edu/PT/SMTS/kwsmts/KWSMTSv01.txt",
            known_hash="sha256:91c93302c44463a424835aa4051b7b2a1ea04d6675d928ca8405b231ae7cea9a",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )
        
        if load:
            df = _pd.read_csv(fname, delim_whitespace=True, skiprows=17, comment='>', 
                    names=['Long', 'Lat', 'Azimuth', 'Major', 'Minor', 'Height', 'FAA', 'VGG', 'Depth', 'CrustAge', 'ID'])
            return _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Long, df.Lat))
        else:
            return fname
        
    if catalogue in ['SIO_all', 'SIO_good', 'SIO_shallow', 'SIO_short', 'SIO_tall']:
        fname = _retrieve(
            url="https://zenodo.org/record/7718512/files/SIO_Seamounts.zip?download=1",
            known_hash="md5:efe6f739d34391f68b568a17eac7fee7",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='seamounts')
        )
        
        dirname = '{:s}/seamounts/'.format(str(_os_cache('gprm')))

        if load:
            df = _pd.read_csv('{:}/SIO_Seamounts/Seamounts_Modeled/{:s}.xyhrdnc'.format(dirname, catalogue[4:]), 
                              delim_whitespace=True, skiprows=17, comment='>', 
                    names=['Long', 'Lat', 'Height', 'Radius', 'Base_Depth', 'Name', 'Charted'])
            return _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Long, df.Lat))
        else:
            return fname

    elif catalogue in ['HillierWatts', 'HW']:
        fname = _retrieve(
            url="https://www.wattsgeophysics.co.uk/downloadfile/5616459",
            known_hash="sha256:d0b9aa7d15754ad9aabecfedf881005d22254e79183af8edf0806be840a549ac",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )

        if load:
            df = _pd.read_csv(fname, delim_whitespace=True, names=['Long', 'Lat', 'Height'])
            return _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Long, df.Lat))
        else:
            return fname

    else:
        raise ValueError('Unknown catalogue {:s}'.format(catalogue))


def LargeIgneousProvinces(catalogue='Whittaker', load=True):
    '''
    (Large) Igneous Province polygons included in GPlates sample data:
    - 'Whittaker' [default], from Whittaker et al (2015)
    - 'Johansson' from Johansson et al (2018)
    and also
    - 'UTIG' from the 2011 version of the UTIG LIP compilation

    '''
    if catalogue in ['Whittaker', 'Johansson']:
        fnames = _retrieve(
                url="https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/SampleData_GPlates2.2/Individual/FeatureCollections/LargeIgneousProvinces_VolcanicProvinces.zip",
                known_hash="sha256:8f86ab86a12761f5534beaaeaddbed5b4e3e6d3d9b52b0c87ee9b15af2a797cd",  
                downloader=_HTTPDownloader(progressbar=True),
                path=_os_cache('gprm'),
                processor=_Unzip(extract_dir='LIPs'),
            )

        for fname in fnames:
            if _os.path.split(fname)[1] == 'License.txt':
                dirname = _os.path.split(fname)[0]

        if catalogue=='Whittaker':
            fname='{:s}/LargeIgneousProvinces_VolcanicProvinces/Whittaker_etal_2015_LargeIgneousProvinces/SHP/Whittaker_etal_2015_LIPs.shp'.format(dirname)
        elif catalogue=='Johansson':
            fname='{:s}/LargeIgneousProvinces_VolcanicProvinces/Johansson_etal_2018_VolcanicProvinces/SHP/Johansson_etal_2018_VolcanicProvinces_v2.shp'.format(dirname)

    elif catalogue=='UTIG':
        fname = _retrieve(
                url="http://www-udc.ig.utexas.edu/external/plates/data/LIPS/Data/LIPS.2011.gmt",
                known_hash="sha256:11cd037382c518ec0b54b93728fef5e476ec3d8d57e5c433a1ccf14420ee99dd",  
                downloader=_HTTPDownloader(progressbar=True),
                path=_os_cache('gprm'),
            )

    else:
        raise ValueError('Unknown catalogue {:s}'.format(catalogue))

    if load:
        return _gpd.read_file(fname)
    else:
        return fname

