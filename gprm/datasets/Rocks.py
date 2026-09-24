"""
Loaders for geochemical, mineral deposit, and rock-type datasets.

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
from ._columns import add_aliases as _add_aliases
from ._ages import stamp as _stamp
from ._gpml import read_gpml_points as _read_gpml_points
from ._remote_zip import retrieve_zip_member as _retrieve_zip_member
import numpy as _np
import pandas as _pd
import geopandas as _gpd
import os as _os


def Geochem(usecols=None, return_column_names=False, remove_invalid_coordinates=True):
    '''
    Load the geochemistry database of Gard et al (2019)
    doi: https://doi.org/10.5194/essd-11-1553-2019

    Options:
    usecols: optionally define a list of columns to load (rather than the full table) [default=None]
    remove_invalid_coordinates: specify whether to drop rows with a missing longitude/latitude
                                or a latitude outside [-90, 90] (which cannot be a real point).
                                A longitude outside [-180, 180] is not dropped -- a small number
                                of rows in the source data record it in 0-360 convention (e.g.
                                194.4 instead of -165.6), and the sphere does not care which
                                convention is used, so these are wrapped into [-180, 180) instead
                                of being discarded. The wrap is applied to the 'Longitude' column
                                gprm adds; the source's own 'longitude' is left as in the file.
                                [default=True]
    return_column_names: instead of loading table into memory, return a list of column names

    Ages: the source's 'age' column is returned unchanged, and gprm adds 'Age' (Ma), which is the
    same except that (a) the 2 values older than the Earth (4567 Ma) and the 3 negative values
    larger than 10 kyr (-7.1, -5.7, -0.2) are NaN, and (b) the 119 values between 0 and -0.01 Ma,
    all young volcanics and apparently historical eruption dates entered as negative years, are
    0. No rows are dropped for their age. gprm.datasets.age_description('Rocks.Geochem') points
    at 'Age'.

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
        # Source columns keep their names; gprm's conventional ones are added alongside
        df = _add_aliases(df, {'longitude': 'Longitude', 'latitude': 'Latitude'})
        if remove_invalid_coordinates:
            df = df.dropna(subset=['Longitude','Latitude'])
            # A latitude outside [-90, 90] cannot be wrapped into a real point and is dropped;
            # a longitude outside [-180, 180] can be (see docstring) and is wrapped, not dropped.
            df = df[df.Latitude.between(-90, 90)]
            df['Longitude'] = ((df.Longitude + 180.) % 360.) - 180.
            df.reset_index(inplace=True, drop=True)
        if 'age' in df.columns:
            age = _pd.to_numeric(df['age'], errors='coerce')
            historical = age.between(-0.01, 0, inclusive='left')    # eruptions of recent decades
            impossible = (age > 4567.) | (age < -0.01)
            df['Age'] = age.mask(historical, 0.).mask(impossible)
        gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)
        return _stamp(gdf, 'Rocks.Geochem')



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
    # Source columns keep their names ('Lon.' in some sheets, 'Lon' in others); gprm's are added
    df = _add_aliases(df, {'Lon.': 'Longitude', 'Lat.': 'Latitude',
                           'Lon': 'Longitude', 'Lat': 'Latitude'})
    gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    # 'Age (Ga)' stays as in the file, including its 'ND' (no data) entries; 'Age' is in Ma
    gdf['Age'] = _pd.to_numeric(gdf['Age (Ga)'], errors='coerce')*1000.

    if not keep_unknown_age_samples:
        gdf = gdf.dropna(subset=['Age'])
    return _stamp(gdf, 'Rocks.BaseMetalDeposits')


# sha256 of the Tappe et al (2018) kimberlite centroids as distributed inside the Flament et al
# (2022) supplement, one file per tectonic reconstruction. M21 and M21NNR share a file, as do
# Y19 and M16 -- removing net rotation changes the rotations, not which plate a point sits on.
_TAPPE2018 = {
    'M21':    'sha256:6cde325c5c75f06d5c7da7ff06a1617c86b84b3ae49fc9238045897c6c3672c1',
    'M21NNR': 'sha256:6cde325c5c75f06d5c7da7ff06a1617c86b84b3ae49fc9238045897c6c3672c1',
    'Y19':    'sha256:866a9eacee933782170a94e7d93528f00f8da71ad31764d06244197a24b7a2e0',
    'M16':    'sha256:866a9eacee933782170a94e7d93528f00f8da71ad31764d06244197a24b7a2e0',
}


def Kimberlites(catalogue='Faure2010', reconstruction='M21', load=True):
    '''
    Load a global kimberlite occurrence compilation.

    catalogue must be one of:
    - 'Faure2010' [default], the CONSOREM world kimberlites and lamproites database (Faure 2010).
      Ages are in the column 'Age1_Ma'; this table carries no plate ids.
    - 'Tappe2018', the compilation of Tappe et al (2018), doi:10.1016/j.epsl.2018.02.013, as
      distributed with the supplement to Flament et al (2022), doi:10.1038/s41586-022-04538-y
      (Zenodo record 6031641). Point locations with a recommended eruption age (0-2848 Ma), the
      geochronometer and method behind it, and a plate id.

    For 'Tappe2018', reconstruction selects which model's plate ids are attached: 'M21' [default,
    Merdith et al 2021, matching fetch_Merdith2021], 'M21NNR' (the same with net rotation
    removed), 'Y19' (Young et al 2019) or 'M16' (Matthews et al 2016). This also changes how many
    kimberlites are returned -- 1129 for M21/M21NNR against 1174 for Y19/M16 -- because the
    reconstructions' static polygons do not assign a plate id to the same set of points. The
    source columns are kept under the (shapefile-truncated) names Tappe et al used; 'Age' is
    added as the recommended age in Ma, following the convention of the other loaders here.

    Note that each feature's validity runs from its eruption age to the present ('FROMAGE' = the
    age, 'TOAGE' = -999, GPlates' value for the distant future), as in the published file. So reconstructing the whole table to a given
    time returns only the kimberlites that had already erupted by then; to follow a time slice,
    select on 'Age' and widen the validity yourself.

    load=False returns the path to the cached file instead of reading it.
    '''
    if catalogue == 'Faure2010':
        fname = _retrieve(
            url="https://consorem2.uqac.ca/production_scientifique/fiches_projets/world_kimberlites_and_lamproites_consorem_database_v2010.xls",
            known_hash="sha256:8d9d8d89afa9304b6494ad32b8f66f6838de5631287c94155d498b7f6d413ac4",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )

        if not load:
            return fname

        df = _pd.read_excel(fname, skiprows=1)
        gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)
        return _stamp(gdf, 'Rocks.Kimberlites:Faure2010')

    elif catalogue == 'Tappe2018':
        if reconstruction not in _TAPPE2018:
            raise ValueError('Unknown reconstruction {} (expected one of {})'.format(
                reconstruction, ', '.join(sorted(_TAPPE2018))))

        # 772 MB archive, 6.7 MB file: pulled out with range requests rather than downloaded whole
        fname = _retrieve_zip_member(
            url="https://zenodo.org/records/6031641/files/Assembly_African_basal_mantle_structure_supplement.zip",
            member=('Assembly_African_basal_mantle_structure_supplement/Volcanic_eruption_locations/'
                    '{0}/T18/T18_centroids_{0}_plateIDs.gpml'.format(reconstruction)),
            known_hash=_TAPPE2018[reconstruction],
            path=_os_cache('gprm'),
        )

        if not load:
            return fname

        gdf = _read_gpml_points(fname)
        gdf['Age'] = _pd.to_numeric(gdf['Age (recom'], errors='coerce')
        return _stamp(gdf, 'Rocks.Kimberlites:Tappe2018')

    else:
        raise ValueError('Unknown catalogue {}'.format(catalogue))


def Carbonatites(keep_unknown_age_samples=False):
    '''
    Load Carbonatite data from Humphreys-Williams and Zahirovic (2021)

    217 of the 593 carbonatites have Age_ma = 0 with Error_ma = 0, and all but one of those have
    no reference: 0 here means "no age determined", not 0 Ma. They are left out by default.
    keep_unknown_age_samples=True keeps them. Either way the source columns are unchanged, and
    the 'Age' column gprm adds (in Ma) is NaN for them, so a table kept whole cannot be
    reconstructed to 0 Ma by accident.

    Note that FROMAGE is Age_ma + Error_ma, not the age (see gprm.datasets.age_description).
    '''
    fnames = _retrieve(
        url="https://zenodo.org/record/5968095/files/1_CarbonatitesShapefile_WithAgeConstraints.zip?download=1",
        known_hash="md5:7f219044c7a1ea9d81fc3410b64b2876",
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Carbonatites')
    )

    for fname in fnames:
        if fname.endswith('carbonatites_gplates.shp'):
            gdf = _gpd.read_file(fname)
            unknown_age = (gdf['Age_ma'] == 0) & (gdf['Error_ma'] == 0)
            gdf['Age'] = gdf['Age_ma'].where(~unknown_age)
            if not keep_unknown_age_samples:
                gdf = gdf[~unknown_age].reset_index(drop=True)
            return _stamp(gdf, 'Rocks.Carbonatites')

    raise FileNotFoundError(
        'carbonatites_gplates.shp was not found in the downloaded Carbonatites archive. The '
        'download may be incomplete or the archive may have been repackaged upstream; clearing '
        'the gprm cache (see gprm.datasets.cache_path()) and retrying is the usual fix.')


def Metamorphism():
    '''
    Load the Metamorphims compilation from Brown and Johnson, as reported in the SM of Liu et al (2022)
    '''
    fname = _retrieve(
        url="https://gsapubs.figshare.com/ndownloader/files/33947312",
        known_hash="fa815bc1f07c347834dc4e9724285bfff12176bd42d5330cbc0cabeab1757fc4",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
    )

    df = _pd.read_excel(fname, skiprows=1)
    df = _add_aliases(df, {'LONGITUDE (˚E)': 'Longitude', 'LATITUDE (˚N)': 'Latitude'})
    df['Age'] = df['AGE(Ga)']*1000.
    gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    return _stamp(gdf, 'Rocks.Metamorphism')
