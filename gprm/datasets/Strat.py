"""
Loaders for stratigraphic and fossil occurrence datasets (paleo-currents, PBDB).

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
from ._columns import add_aliases as _add_aliases
from ._ages import stamp as _stamp
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
import pandas as _pd
import geopandas as _gpd
import numpy as _np
import os as _os

_DATA_DIR = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), 'Data')


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


def pbdb(path_to_pbdb_data=None, usecols=None):
    """
    Load data from pbdb downloaded file (not retrieved by pooch)

    The function assumes that a file has already been downloaded using 
    the pbdb navigator interface to the pbdb web service (version 1.2)
    """

    if not path_to_pbdb_data:
        path_to_pbdb_data = '{:s}/pbdb/pbdb_data.csv'.format(str(_os_cache('gprm')))

    if usecols is None:
        df = _pd.read_csv(path_to_pbdb_data, delimiter=',', skiprows=14)
    else:
        df = _pd.read_csv(path_to_pbdb_data, delimiter=',', skiprows=14, usecols=usecols)

    df = _add_aliases(df, {'lng': 'Longitude', 'lat': 'Latitude'})

    gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    return _stamp(gdf, 'Strat.pbdb')


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

    nan_dict = {_np.nan: (_np.nan, _np.nan)}

    marine_env_dict.update(terrestrial_env_dict)
    marine_env_dict.update(nan_dict)

    # TODO change this so that the eleva
    #pbdb = pbdb.dropna(subset=['environment']).reset_index(drop=True)

    elevation_ranges = _pd.DataFrame(pbdb['environment'].map(marine_env_dict).to_list(),
                                     index=pbdb.index,
                                     columns=['elevation_min', 'elevation_max'])

    return pbdb.join(elevation_ranges)


# LithCode -> indicator name, from Boucot, Chen & Scotese (2013). 'M' is Boucot's own
# code for two distinct indicators (mostly mangroves, some lateritic manganese) and is
# left conflated here rather than split, since the source data doesn't disambiguate it.
BOUCOT_INDICATORS = {
    'C': 'Coal',
    'E': 'Evaporite',
    'B': 'Bauxite',
    'K': 'Kaolinite',
    'CA': 'Calcrete',
    'T': 'Tillite',
    'CR': 'Crocodilian',
    'L': 'Laterite',
    'D': 'Dropstone',
    'PA': 'Palm',
    'M': 'Mangrove or lateritic manganese',
    'G': 'Glendonite',
    'O': 'Oolitic ironstone',
    'I': 'Ice crystal',
    'LF': 'Lungfish burrow',
    'H': 'Humid soil',
}

# LithCode -> climate group, following the 3-way scheme of Cao et al. (2018, Geol. Mag.)
# doi:10.1017/S0016756818000110. Coals indicate terrestrial humidity, evaporites indicate
# aridity, and tillites/dropstones/glendonites indicate glacial/cold conditions. Cao et al.
# deliberately excluded the remaining indicators as unreliable latitude proxies (palms,
# mangroves and crocodilians are sampling-biased toward mid-high latitudes; laterites and
# oolitic ironstones have too few occurrences) -- codes not listed here are left unmapped.
BOUCOT_CLIMATE_GROUPS = {
    'C': 'Humid',
    'E': 'Arid',
    'T': 'Glacial',
    'D': 'Glacial',
    'G': 'Glacial',
}


def PaleoLithology(lithology=None, reconstruction_time=None):
    """
    Load the Boucot, Chen & Scotese (2013) compilation of climate-sensitive
    palaeolithologic and biotic indicators (coals, evaporites, bauxites, calcretes,
    tillites, dropstones, glendonites, kaolinites, laterites, oolitic ironstones,
    palms, mangroves, crocodilians, and a few singletons), spanning the Cambrian
    to the Miocene.

    Boucot, A.J., Chen, X., Scotese, C.R. & Morley, R.J. (2013). Phanerozoic
    Paleoclimate: An Atlas of Lithologic Indicators of Climate. SEPM Concepts in
    Sedimentology and Paleontology, 11.

    The returned GeoDataFrame carries no plate id: the source shapefile's plate ids
    were assigned by partitioning against an unrecorded static polygon set, so they
    are dropped rather than risk misleading anyone reconstructing against a
    different plate model. Before calling ``ReconstructionModel.reconstruct()``, run
    the result through ``ReconstructionModel.assign_plate_ids()`` against your
    chosen static polygons.

    :param lithology: a single indicator code or name (or a list of them) to select,
        e.g. 'C', 'Coal', ['T', 'D', 'G']. Case-insensitive. Default: all indicators.
    :param reconstruction_time: if given, keep only points valid at this age (Ma),
        i.e. where TOAGE <= reconstruction_time <= FROMAGE.
    :returns: GeoDataFrame with columns including LithCode, Indicator, GeogComm,
        Continent, Country, Stage, FROMAGE, TOAGE, ReconstructionAge, Lithology,
        Formation, LithComm, PrimRef, SeeAlso, geometry.
    """
    fname = '{:s}/boucot_paleolithology.gpkg'.format(_DATA_DIR)
    gdf = _gpd.read_file(fname)
    gdf['Indicator'] = gdf['LithCode'].map(BOUCOT_INDICATORS)

    if lithology is not None:
        if isinstance(lithology, str):
            lithology = [lithology]
        indicator_names_lower = {name.lower(): code for code, name in BOUCOT_INDICATORS.items()}
        codes = set()
        for item in lithology:
            item = item.strip()
            if item.upper() in BOUCOT_INDICATORS:
                codes.add(item.upper())
            elif item.lower() in indicator_names_lower:
                codes.add(indicator_names_lower[item.lower()])
            else:
                raise ValueError(f"'{item}' is not a recognised LithCode or indicator name")
        gdf = gdf[gdf['LithCode'].isin(codes)]

    if reconstruction_time is not None:
        gdf = gdf[(gdf['TOAGE'] <= reconstruction_time) & (reconstruction_time <= gdf['FROMAGE'])]

    return _stamp(gdf.reset_index(drop=True), 'Strat.PaleoLithology')


def paleolithology_climate_mapping(gdf):
    """
    Add a 'ClimateGroup' column (Humid/Arid/Glacial) to a GeoDataFrame returned by
    PaleoLithology(), following the 3-way scheme of Cao et al. (2018, Geol. Mag.)
    doi:10.1017/S0016756818000110. Indicators outside that scheme are left as NaN --
    this is not a judgement that they carry no climate signal, only that Cao et al.
    considered them unreliable latitude proxies (see BOUCOT_CLIMATE_GROUPS).

    :param gdf: a GeoDataFrame with a LithCode column, as returned by PaleoLithology().
    :returns: the same GeoDataFrame with a ClimateGroup column added.
    """
    return gdf.assign(ClimateGroup=gdf['LithCode'].map(BOUCOT_CLIMATE_GROUPS))

