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


def pbdb(path_to_pbdb_data=None):
    """
    Load data from pbdb downloaded file (not retrieved by pooch)

    The function assumes that a file has already been downloaded using 
    the pbdb navigator interface to the pbdb web service (version 1.2)
    """

    if not path_to_pbdb_data:
        path_to_pbdb_data = '{:s}/pbdb/pbdb_data.csv'.format(str(_os_cache('gprm')))

    df = _pd.read_csv(path_to_pbdb_data, delimiter=',', skiprows=14)

    df.rename(columns={'lng':'Longitude', 'lat':'Latitude'}, inplace=True)

    gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    return gdf


def pbdb_elevation_mapping(pbdb):

    # first we define the mapping dictionary. Ultimately this should be moved somewhere 
    # outside this function
    #
    # first the bathymetry mapping, based on the work of:
    # Fernandes and Roberts (2020), GSA Bulletin
    # doi:
    marine_env_dict = {
        'basin reef': (-1000, -500),
        'basinal (carbonate)': (-3000, -500),
        'basinal (siliciclastic)': (-4000, -500),
        'basinal (siliceous)': (-4000, -500),
        'carbonate indet.': (-3000, 0),
        'coastal indet.': (-50, 0),
        'deep subtidal indet.': (-50, -15),
        'deep subtidal ramp': (-50, -15),
        'deep subtidal shelf': (-50, -15),
        'deep-water indet.': (-4000, -250),
        'delta front': (-15, -1),
        'estuary/bay': (-15, 0),
        'foreshore': (-1, 0),
        'interdistributary bay': (-2, 0),
        'intrashelf/intraplatform reef': (-100, 0),
        'lagoonal/restricted shallow subtidal': (-2, 0),
        'lagoonal': (-5, 0),
        'marginal marine indet.': (-15, 0),
        'marine indet.': (-4000, 0),
        'offshore': (-250, -50),
        'offshore indet.': (-250, -50), 
        'offshore ramp': (-250, -50), 
        'offshore shelf': (-250, -50), 
        'open shallow subtidal': (-15, -1), 
        'paralic indet.': (-15, 0), 
        'perireef or subreef': (-250, -1), 
        'peritidal': (-2, 0),
        'platform/shelf-margin reef': (-250, 0),
        'prodelta': (-100, -15),
        'reef, buildup or bioherm': (-250, 0), 
        'sand shoal': (-15, 0), 
        'shallow subtidal indet.': (-15, -1), 
        'shoreface': (-50, -1),
        'slope': (-4000, -250), 
        'slope/ramp reef': (-250, -50),
        'submarine fan': (-4000, -1500),
        'transition zone/lower shoreface': (-70, -15), 
    }

    # Second the terrestrial environments
    # No meaningful values, just a way to map non-marine enviroments to 
    # have land elevations
    terrestrial_env_dict = {
        '"channel"': (0, 1000), 
        '"floodplain"': (0, 1000), 
        'alluvial fan': (0, 1000), 
        'cave': (0, 1000), 
        'channel lag': (0, 1000), 
        'coarse channel fill': (0, 1000), 
        'crater lake': (0, 1000), 
        'crevasse splay': (0, 1000), 
        'delta plain': (0, 100), 
        'deltaic indet.': (0, 100), 
        'dry floodplain': (0, 1000), 
        'dune': (0, 1000), 
        'eolian indet.': (0, 1000), 
        'fine channel fill': (0, 1000), 
        'fissure fill': (0, 1000), 
        'fluvial indet.': (0, 1000), 
        'fluvial-deltaic indet.': (0, 1000), 
        'fluvial-lacustrine indet.': (0, 1000), 
        'glacial': (0, 1000), 
        'interdune': (0, 1000), 
        'lacustrine - large': (0, 1000), 
        'lacustrine - small': (0, 1000),  
        'lacustrine delta front': (0, 1000), 
        'lacustrine delta plain': (0, 1000), 
        'lacustrine deltaic indet.': (0, 1000), 
        'lacustrine indet.': (0, 1000), 
        'lacustrine interdistributary bay': (0, 1000), 
        'levee': (0, 1000), 
        'loess': (0, 1000), 
        'mire/swamp': (0, 1000), 
        'pond': (0, 1000), 
        'sinkhole': (0, 1000), 
        'spring': (0, 1000), 
        'tar': (0, 1000), 
        'terrestrial indet.': (0, 1000), 
        'wet floodplain': (0, 1000), 
    }

    marine_env_dict.update(terrestrial_env_dict)

    elevation_ranges = _pd.DataFrame(pbdb['environment'].map(marine_env_dict).to_list(), 
                                     index=pbdb.index, 
                                     columns=['elevation_min', 'elevation_max'])

    return pbdb.join(elevation_ranges)

