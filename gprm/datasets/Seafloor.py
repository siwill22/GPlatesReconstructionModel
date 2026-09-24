"""
Loaders for seafloor datasets: magnetic picks, fabric, seamounts, and LIPs.

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
from pooch import Untar as _Untar
from pooch import Unzip as _Unzip
from ._columns import add_aliases as _add_aliases
from ._ages import stamp as _stamp
from ._gpml import read_gpml_points as _read_gpml_points
from ._remote_zip import retrieve_zip_member as _retrieve_zip_member
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
        return _stamp(_gpd.read_file(fname), 'Seafloor.MagneticPicks')
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

    raise FileNotFoundError(
        '{:s} was not found in the downloaded seafloor fabric archive. The download may be '
        'incomplete or the archive may have been repackaged upstream; clearing the gprm cache '
        '(see gprm.datasets.cache_path()) and retrying is the usual fix.'.format(
            FABRIC_TYPE[feature_type]))


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
            # No header row, so the columns take the names in the file's commented header line:
            #   #Lon Lat Average_age(Ma) Average_error(Ma) Tag Name(island,seamount,plateau_or_sample) Island_or_seamount_chain
            # except that in every data row the name comes before the two-letter chain code
            # (e.g. 'Macdonald  AC  Austral'), so those two are named by what they hold.
            df = _pd.read_csv(fname, comment='#', sep=r'\s+',
                              names=['Lon', 'Lat', 'Average_age(Ma)', 'Average_error(Ma)',
                                     'Name(island,seamount,plateau_or_sample)', 'Tag',
                                     'Island_or_seamount_chain'])
            df = _add_aliases(df, {'Lon': 'Long',
                                   'Average_age(Ma)': 'Average_Age_Ma',
                                   'Average_error(Ma)': 'Average_Age_Error_Ma',
                                   'Name(island,seamount,plateau_or_sample)': 'SeamountName',
                                   'Island_or_seamount_chain': 'SeamountChain'})
            gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Lon, df.Lat))
            return _stamp(gdf, 'Seafloor.PacificSeamountAges:2013')
        else:
            return fname
        
    elif catalogue=='2021':
        fnames = _retrieve(
            url="https://zenodo.org/record/6558676/files/Pacific_Hotspot_Trails_Datasets.zip?download=1",
            known_hash="md5:97bec4ebde94694698077eb527bc1ef4",
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='PHT2021')
        )

        fname = None
        for candidate in fnames:
            if candidate.endswith('PHT2021_pacific_ages.txt'):
                fname = candidate
        if fname is None:
            raise FileNotFoundError(
                'PHT2021_pacific_ages.txt was not found in the downloaded Pacific Hotspot '
                'Trails archive. The download may be incomplete or the archive may have been '
                'repackaged upstream; clearing the gprm cache (see gprm.datasets.cache_path()) '
                'and retrying is the usual fix.')

        if load:
            # No header row, so the columns take the names in the file's commented header line:
            #   # Lon Lat Age Error Type Ref Name(Sample) Tag Chain
            df = _pd.read_csv(fname, comment='#', sep=r'\s+',
                              names=['Lon', 'Lat', 'Age', 'Error', 'Type', 'Ref', 'Name(Sample)',
                                     'Tag', 'Chain'])
            df = _add_aliases(df, {'Lon': 'Long',
                                   'Age': 'Average_Age_Ma',
                                   'Error': 'Average_Age_Error_Ma',
                                   'Name(Sample)': 'SampleName',
                                   'Chain': 'SeamountChain'})
            gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Lon, df.Lat))
            return _stamp(gdf, 'Seafloor.PacificSeamountAges:2021')
        else:
            return fname

    else:
        raise ValueError(
            "Unknown catalogue '{:s}'. Valid options are '2021' and '2013'.".format(catalogue))


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
            # 17 '#' header lines, the last of which names the columns; '>' lines separate the
            # ocean basins. Note CrustAge is the age of the seafloor beneath the seamount (from
            # the AGE 3.2 grid, per the header), not the age of the seamount.
            df = _pd.read_csv(fname, sep=r'\s+', skiprows=17, comment='>',
                    names=['Longitude', 'Latitude', 'Azimuth', 'Major', 'Minor', 'Height', 'FAA',
                           'VGG', 'Depth', 'CrustAge', 'ID'])
            df = _add_aliases(df, {'Longitude': 'Long', 'Latitude': 'Lat'})
            gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude))
            return _stamp(gdf, 'Seafloor.Seamounts:KimWessel')
        else:
            return fname
        
    if catalogue in ['SIO_all', 'SIO_good', 'SIO_shallow', 'SIO_short', 'SIO_tall']:
        fnames = _retrieve(
            url="https://zenodo.org/record/7718512/files/SIO_Seamounts.zip?download=1",
            known_hash="md5:efe6f739d34391f68b568a17eac7fee7",
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='seamounts')
        )

        target = '{:s}.xyhrdnc'.format(catalogue[4:])
        fname = None
        for candidate in fnames:
            if _os.path.split(candidate)[1] == target:
                fname = candidate
        if fname is None:
            raise FileNotFoundError(
                '{:s} was not found in the downloaded SIO Seamounts archive. The download may '
                'be incomplete or the archive may have been repackaged upstream; clearing the '
                'gprm cache (see gprm.datasets.cache_path()) and retrying is the usual fix.'.format(target))

        if load:
            # No header row: every line is a seamount (the README gives 39399 for good.xyhrdnc).
            # Column names are the README's; its seventh column ('1 or 0 charted or uncharted')
            # is unnamed there, so it keeps gprm's name 'Charted'.
            df = _pd.read_csv(fname, sep=r'\s+', comment='>',
                    names=['longitude', 'latitude', 'height', 'radius', 'base_depth', 'name', 'Charted'])
            df = _add_aliases(df, {'longitude': 'Long', 'latitude': 'Lat', 'height': 'Height',
                                   'radius': 'Radius', 'base_depth': 'Base_Depth', 'name': 'Name'})
            gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.longitude, df.latitude))
            return _stamp(gdf, 'Seafloor.Seamounts:SIO')
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
            df = _pd.read_csv(fname, sep=r'\s+', names=['Long', 'Lat', 'Height'])
            gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Long, df.Lat))
            return _stamp(gdf, 'Seafloor.Seamounts:HillierWatts')
        else:
            return fname

    else:
        raise ValueError('Unknown catalogue {:s}'.format(catalogue))


# sha256 of the Johansson et al (2018) volcanic province centroids as distributed inside the
# Flament et al (2022) supplement, one file per tectonic reconstruction. All four hold the same
# 185 centroids; they differ in the plate ids attached (and in GPlates' internal feature ids).
_JOHANSSON_CENTROIDS = {
    'M21':    'sha256:17489d665ba216125c146560a24b4d60030f59430e776560936822f7dc56282a',
    'M21NNR': 'sha256:c7771492aa6e8826c8e78bda2ef43654150381df5cf1bf9b784d422ab4a1cde7',
    'Y19':    'sha256:34c5c5cdef0ed1a6dc79f99746cbd09ca07f61af61183e2bd3aad6483e92eca7',
    'M16':    'sha256:2a3202b428084bc33e1655aa90b67effc95fd4bc4b25f5d3dae0af886350891c',
}


def LargeIgneousProvinces(catalogue='Whittaker', reconstruction='M21', load=True,
                          keep_unknown_age_samples=False):
    '''
    (Large) Igneous Province polygons included in GPlates sample data:
    - 'Whittaker' [default], from Whittaker et al (2015)
    - 'Johansson' from Johansson et al (2018)
    and also
    - 'UTIG' from the 2011 version of the UTIG LIP compilation

    - 'Johansson_centroids' returns **points, not polygons**: the 185 centroids of the same
      Johansson et al (2018) catalogue, each carrying its emplacement age ('Age', in Ma, from
      the begin time of the feature's validity) and a plate id. Taken from the supplement to
      Flament et al (2022), doi:10.1038/s41586-022-04538-y (Zenodo record 6031641). Use this
      when you want ages and plate ids to reconstruct with; use 'Johansson' for the outlines.
      Each centroid is valid from its emplacement age to the present, so reconstructing the whole
      set to a given time returns only the provinces already emplaced by then.

    reconstruction applies to 'Johansson_centroids' only, and selects whose plate ids are
    attached: 'M21' [default, Merdith et al 2021, matching fetch_Merdith2021], 'M21NNR' (the same
    with net rotation removed), 'Y19' (Young et al 2019) or 'M16' (Matthews et al 2016).

    For both Johansson catalogues, an emplacement age (FROMAGE) of 0 is treated as an error
    rather than a real age: 144 of the 2526 polygons and 2 of the 185 centroids have one
    (e.g. the Tuamotu seamounts). These are left out by default; keep_unknown_age_samples=True
    keeps them. Either way FROMAGE is unchanged, and the 'Age' column gprm adds is NaN for them.

    '''
    if catalogue == 'Johansson_centroids':
        if reconstruction not in _JOHANSSON_CENTROIDS:
            raise ValueError('Unknown reconstruction {} (expected one of {})'.format(
                reconstruction, ', '.join(sorted(_JOHANSSON_CENTROIDS))))

        # 772 MB archive, 0.33 MB file: pulled out with range requests rather than downloaded whole
        fname = _retrieve_zip_member(
            url="https://zenodo.org/records/6031641/files/Assembly_African_basal_mantle_structure_supplement.zip",
            member=('Assembly_African_basal_mantle_structure_supplement/Volcanic_eruption_locations/'
                    '{0}/J18/J18_centroids_{0}_plateIDs.gpml'.format(reconstruction)),
            known_hash=_JOHANSSON_CENTROIDS[reconstruction],
            path=_os_cache('gprm'),
        )

        if not load:
            return fname

        gdf = _johansson_ages(_read_gpml_points(fname), keep_unknown_age_samples)
        return _stamp(gdf, 'Seafloor.LargeIgneousProvinces:Johansson_centroids')

    elif catalogue in ['Whittaker', 'Johansson']:
        fnames = _retrieve(
                url="https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/SampleData_GPlates2.2/Individual/FeatureCollections/LargeIgneousProvinces_VolcanicProvinces.zip",
                known_hash="sha256:8f86ab86a12761f5534beaaeaddbed5b4e3e6d3d9b52b0c87ee9b15af2a797cd",  
                downloader=_HTTPDownloader(progressbar=True),
                path=_os_cache('gprm'),
                processor=_Unzip(extract_dir='LIPs'),
            )

        dirname = None
        for fname in fnames:
            if _os.path.split(fname)[1] == 'License.txt':
                dirname = _os.path.split(fname)[0]
        if dirname is None:
            raise FileNotFoundError(
                'The Large Igneous Provinces archive did not contain the expected License.txt, '
                'so the location of the data files could not be determined. The download may '
                'be incomplete or the archive may have been repackaged upstream; clearing the '
                'gprm cache (see gprm.datasets.cache_path()) and retrying is the usual fix.')

        if catalogue=='Whittaker':
            fname='{:s}/LargeIgneousProvinces_VolcanicProvinces/Whittaker_etal_2015_LargeIgneousProvinces/SHP/Whittaker_etal_2015_LIPs.shp'.format(dirname)
        elif catalogue=='Johansson':
            fname='{:s}/LargeIgneousProvinces_VolcanicProvinces/Johansson_etal_2018_VolcanicProvinces/SHP/Johansson_etal_2018_VolcanicProvinces_v2.shp'.format(dirname)

    elif catalogue=='UTIG':
        import pygplates
        fname = _retrieve(
                url="http://www-udc.ig.utexas.edu/external/plates/data/LIPS/Data/LIPS.2011.gmt",
                known_hash="sha256:11cd037382c518ec0b54b93728fef5e476ec3d8d57e5c433a1ccf14420ee99dd",  
                downloader=_HTTPDownloader(progressbar=True),
                path=_os_cache('gprm'),
            )
        pygplates.FeatureCollection(fname).write('{:s}/LIPS_2011.gmt'.format(str(_os_cache('gprm'))))
        fname = '{:s}/LIPS_2011.gmt'.format(str(_os_cache('gprm')))

    else:
        raise ValueError('Unknown catalogue {:s}'.format(catalogue))

    if load:
        gdf = _gpd.read_file(fname)
        if catalogue == 'Johansson':
            gdf = _johansson_ages(gdf, keep_unknown_age_samples)
        return _stamp(gdf, 'Seafloor.LargeIgneousProvinces:' + catalogue)
    else:
        return fname


def _johansson_ages(gdf, keep_unknown_age_samples):
    """Add 'Age' (FROMAGE, NaN where FROMAGE is 0) and drop the zero-age rows unless asked not to.

    An emplacement age of 0 in the Johansson et al (2018) catalogue is taken to be an error, not a
    real age (it includes, e.g., the Cenozoic Tuamotu seamounts).
    """
    unknown_age = gdf['FROMAGE'] == 0
    gdf['Age'] = gdf['FROMAGE'].where(~unknown_age)
    if not keep_unknown_age_samples:
        gdf = gdf[~unknown_age].reset_index(drop=True)
    return gdf

