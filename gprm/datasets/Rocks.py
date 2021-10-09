from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
import pandas as _pd
import geopandas as _gpd
import os as _os


def Geochem(usecols=None, return_column_names=False, remove_invalid_coordinates=True):
    '''
    Load the geochemistry database of Gard et al (2019)
    doi: https://doi.org/10.5194/essd-11-1553-2019

    Options:
    usecols: optionally define a list of columns to load (rather than the full table) [default=None]
    remove_invalid_coordinates: specify whether to remove rows from the table for which the latitude
                                and/or longitide are invalid [default=True] 
    return_column_names: instead of loading table into memory, return a list of column names

    '''
    fname = _retrieve(
        url="https://zenodo.org/record/3359791/files/complete.zip",
        known_hash="md5:9b97b54887ee7184c6650c845b4e92d4",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='GeochemDB'),
    )[0]
    
    if usecols:
        # TODO if 'Long and Lat fields are not included, the attempt to create a gdf will throw an error - fix this
        usecols = ['longitude' if x=='Longitude' else x for x in usecols]
        usecols = ['latitude' if x=='Latitude' else x for x in usecols]

    # TODO include some shorthands for the usecols, e.g. 'Major', 'REE', etc

    if return_column_names:
        return _pd.read_csv(fname, index_col=0, nrows=0, engine='python').columns.tolist()

    else:
        df = _pd.read_csv(fname, usecols=usecols, engine='python', encoding="ISO-8859-1")
        if remove_invalid_coordinates:
            df = df.dropna(subset=['longitude','latitude'])
            df.reset_index(inplace=True)
        df.rename(columns={'longitude':'Longitude', 'latitude':'Latitude'}, inplace=True)
        return _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)



def BaseMetalDeposits(deposit_type, keep_unknown_age_samples=False):
    '''
    Load the base metal deposit compilation from Hoggard et al (2020)
    doi: https://doi.org/10.1038/s41561-020-0593-2

    deposit_type must be one of: 'PbZn-CD', 'PbZn-MVT', 'Cu-sed', 'Magmatic Ni', 'VMS', 'Cu-por', 'IOCG'
    '''
    fname = _retrieve(
        url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41561-020-0593-2/MediaObjects/41561_2020_593_MOESM3_ESM.xls",
        known_hash="sha256:c1ddf941c490dcc55cce4ec5da40eac5cf2b88715941f6520235d3fbefb6de80",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
    )

    if deposit_type not in ['PbZn-CD', 'PbZn-MVT', 'Cu-sed', 'Magmatic Ni', 'VMS', 'Cu-por', 'IOCG']:
        raise ValueError('Unknown deposit type {}'.format(deposit_type))

    df = _pd.read_excel(fname, sheet_name=deposit_type)
    df.rename(columns={'Lon.':'Longitude', 'Lat.':'Latitude',
                       'Lon':'Longitude', 'Lat':'Latitude'}, inplace=True)
    gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    if not keep_unknown_age_samples:
        return gdf[gdf['Age (Ga)'] != 'ND']
    else:
        return gdf


def Kimberlites():
    '''
    Load the Kimberlite compilation from Faure (2010)
    '''
    fname = _retrieve(
        url="https://consorem2.uqac.ca/production_scientifique/fiches_projets/world_kimberlites_and_lamproites_consorem_database_v2010.xls",
        known_hash="sha256:8d9d8d89afa9304b6494ad32b8f66f6838de5631287c94155d498b7f6d413ac4",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
    )

    df = _pd.read_excel(fname, skiprows=1)
    gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    return gdf


