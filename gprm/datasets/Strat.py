from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
import pandas as _pd
import geopandas as _gpd
import os as _os


def PaleoCurrents():

    fname = _retrieve(
        # URL to one of Pooch's test files
        url="https://datadryad.org/stash/downloads/file_stream/89328",
        known_hash="c438bf62eb5169a22409cb8080a2665a2806127733a735adf908a9e0572ca899",
        downloader=_HTTPDownloader(progressbar=True),
    )

    df = _pd.read_excel(fname, sheet_name='Main')
    #if remove_invalid_coordinates:
    #    df = df.dropna(subset=['Longitude','Latitude'])
    #    df.reset_index(inplace=True)

    Paleocurrent_Indicator_Dict = {1: 'crossbedding',
                                   2: 'ripple marks', 
                                   3: 'paleocurrent indicator',
                                   4: 'sole marks',
                                   5: 'fossil orientation',
                                   6: 'wind direction',
                                   7: 'current direction',
                                   8: 'turbidity currents',
                                   9: 'topography',
                                   10: 'miscellaneous',
                                   11: 'slumps and folds',
                                   12: 'flute/grooves',
                                   13: 'imbrication',
                                   14: 'channel axes',
                                   15: 'parting lineations',
                                   16: 'model',
                                   17: 'provenance',
                                   18: 'sed thickening',
                                   19: 'electric log/dip log',
                                   20: 'grain orientation'}

    Environment_Dict = {1: 'marine general',
                        2: 'marine shallow', 
                        3: 'marine deep',
                        4: 'lacustrine',
                        5: 'fluvial deltaic',
                        6: 'fluviatile',
                        7: 'alluvial',
                        8: 'subaerial (eolian)'}

    Lithology_Dict = {1: 'sandstone',
                      2: 'shale',
                      3: 'siltstone or turbidites',
                      4: 'conglomerate',
                      5: 'limestone',
                      6: 'carbonate sand',
                      7: 'volcanic or glacial'}

    df = df.replace({'Paleocurrent Indicator': Paleocurrent_Indicator_Dict,
                    'Environment': Environment_Dict,
                    'Lithology': Lithology_Dict})


    return _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

