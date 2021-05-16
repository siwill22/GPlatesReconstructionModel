from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Untar
import geopandas as gpd



def MagneticPicks(load=True):
    '''
    Magnetic Picks from the 'Global Seafloor Fabric (and) Magnetic Linations' database,
    returned as a geopandas dataframe

    '''
    fname = _retrieve(
        url="http://www.soest.hawaii.edu/PT/GSFML/ML/DATA/GSFML.global.picks.gmt",
        known_hash="sha256:0895b76597f600a6c6184a7bec0edc0df5ca9234255f3f7bac0fe944317caf65",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    
    data = gpd.read_file(fname)
    
    return data


def SeafloorFabric(feature_type, load=True):
    '''
    MarsTopo2600 is a 2600 degree and order spherical harmonic model of the
    shape of the planet Mars. The coefficients are in units of meters.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    '''
    fnames = _retrieve(
        url="http://www.soest.hawaii.edu/PT/GSFML/SF/DATA/GSFML_SF.tbz",
        known_hash="sha256:e27a73dc544611685144b4587d17f03bde24438ee4646963f10761f8ec2e6036",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
        processor=Untar(),
    )
    
    print(fnames)
    data = [gpd.read_file(fname) for fname in fnames if '.gmt' in fname]
    return data
