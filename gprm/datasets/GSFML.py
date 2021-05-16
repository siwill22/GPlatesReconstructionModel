from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Untar as _Untar
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


def SeafloorFabric(feature_type, load=True):
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
        processor=_Untar(),
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
        raise ValueError('Unknown feature type {}'.format(feature_type))

    for fname in fnames:
        if _os.path.split(fname)[1]==FABRIC_TYPE[feature_type]:
            if load:
                return _gpd.read_file(fname)
            else:
                return fname

