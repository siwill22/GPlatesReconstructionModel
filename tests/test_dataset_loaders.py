"""Correctness checks for individual dataset loaders, as opposed to test_network.py's URL
liveness checks -- these actually download, parse, and inspect the result. Marked ``network``
and excluded from the default run, same as test_network.py, because a failure here means an
upstream dataset changed shape, not that gprm's own logic is broken.
"""
import pytest

pytestmark = pytest.mark.network


def test_geochem_longitudes_are_wrapped_not_dropped():
    """remove_invalid_coordinates=True (the default) used to only drop rows with a missing
    longitude/latitude; 37 of the 1.006M rows record longitude in 0-360 convention (e.g. 194.4
    instead of -165.6), which survived that check unchanged. A sphere has no invalid longitude
    short of NaN, so these should be wrapped into [-180, 180), not left out of range."""
    from gprm.datasets.Rocks import Geochem

    gdf = Geochem(usecols=['Longitude', 'Latitude', 'sample_name'])

    assert len(gdf) > 1_000_000
    assert gdf.geometry.x.between(-180, 180, inclusive='left').all()
    assert gdf.geometry.y.between(-90, 90).all()
    assert 'index' not in gdf.columns  # reset_index(drop=True), not a stray column


def test_carbonatites_loads():
    from gprm.datasets.Rocks import Carbonatites

    gdf = Carbonatites()
    assert len(gdf) > 0
    assert gdf.geometry.notna().all()


def test_carbonatites_with_no_age_are_left_out_unless_asked_for():
    """217 rows have Age_ma = 0 and Error_ma = 0, which means 'no age determined'. Dropped by
    default; kept on request, with their source columns untouched but 'Age' NaN."""
    from gprm.datasets.Rocks import Carbonatites

    dated = Carbonatites()
    everything = Carbonatites(keep_unknown_age_samples=True)

    assert len(everything) == 593
    assert len(dated) == 593 - 217
    assert dated.Age.notna().all() and (dated.Age > 0).all()

    unknown = everything[everything.Age.isna()]
    assert len(unknown) == 217
    assert (unknown.Age_ma == 0).all() and (unknown.Error_ma == 0).all()
    assert dated.attrs['gprm_age']['field'] == 'Age'


def test_pacific_seamount_ages_2021_loads_and_load_false_returns_one_path():
    from gprm.datasets.Seafloor import PacificSeamountAges

    gdf = PacificSeamountAges(catalogue='2021')
    assert len(gdf) > 0

    path = PacificSeamountAges(catalogue='2021', load=False)
    assert isinstance(path, str)
    assert path.endswith('PHT2021_pacific_ages.txt')


def test_sio_seamounts_loads_and_load_false_returns_one_path():
    from gprm.datasets.Seafloor import Seamounts

    gdf = Seamounts(catalogue='SIO_good')
    assert len(gdf) > 0

    path = Seamounts(catalogue='SIO_good', load=False)
    assert isinstance(path, str)
    assert path.endswith('good.xyhrdnc')


def test_kimberlites_default_is_still_the_faure_table():
    """Kimberlites() grew a catalogue= argument; the no-argument call must not have changed."""
    from gprm.datasets.Rocks import Kimberlites

    gdf = Kimberlites()
    assert len(gdf) > 0
    assert 'Age1_Ma' in gdf.columns


def test_tappe2018_kimberlites_extracted_from_the_flament_supplement():
    """1129 points for M21, against 1174 for Y19 -- the reconstructions' static polygons do not
    assign a plate id to the same set of kimberlites. Both counts read off the published files."""
    from gprm.datasets.Rocks import Kimberlites

    gdf = Kimberlites(catalogue='Tappe2018')

    assert len(gdf) == 1129
    assert gdf.geometry.x.between(-180, 180).all()
    assert gdf.geometry.y.between(-90, 90).all()
    assert gdf.Age.notna().all()
    assert gdf.Age.between(0, 2848).all()
    assert (gdf.PLATEID1 > 0).all()
    assert 'Kimberlite' in gdf.columns and 'Country' in gdf.columns

    assert len(Kimberlites(catalogue='Tappe2018', reconstruction='Y19')) == 1174


def test_tappe2018_only_caches_the_member_not_the_772mb_archive():
    from gprm.datasets.Rocks import Kimberlites

    path = Kimberlites(catalogue='Tappe2018', load=False)
    assert isinstance(path, str)
    assert path.endswith('T18_centroids_M21_plateIDs.gpml')

    import os
    assert os.path.getsize(path) < 20e6


def test_johansson_lip_centroids_carry_ages_and_plate_ids():
    from gprm.datasets.Seafloor import LargeIgneousProvinces

    gdf = LargeIgneousProvinces(catalogue='Johansson_centroids')

    assert len(gdf) == 183                      # 185 less the two with an emplacement age of 0
    assert gdf.geom_type.eq('Point').all()      # points, unlike the 'Johansson' polygons
    assert gdf.Age.gt(0).all() and gdf.Age.le(650).all()
    assert (gdf.PLATEID1 > 0).all()


@pytest.mark.parametrize('catalogue, total, zero_age', [('Johansson', 2526, 144),
                                                        ('Johansson_centroids', 185, 2)])
def test_johansson_zero_ages_are_left_out_unless_asked_for(catalogue, total, zero_age):
    """An emplacement age of 0 in Johansson et al (2018) is an error, not a real age."""
    from gprm.datasets.Seafloor import LargeIgneousProvinces

    everything = LargeIgneousProvinces(catalogue=catalogue, keep_unknown_age_samples=True)
    dated = LargeIgneousProvinces(catalogue=catalogue)

    assert len(everything) == total
    assert len(dated) == total - zero_age
    assert everything.Age.isna().sum() == zero_age
    assert (everything.loc[everything.Age.isna(), 'FROMAGE'] == 0).all()    # source column untouched
    assert dated.attrs['gprm_age']['field'] == 'Age'


def test_unknown_catalogue_and_reconstruction_are_rejected():
    import pytest
    from gprm.datasets.Rocks import Kimberlites

    with pytest.raises(ValueError):
        Kimberlites(catalogue='NotACatalogue')
    with pytest.raises(ValueError):
        Kimberlites(catalogue='Tappe2018', reconstruction='NotAModel')


@pytest.mark.parametrize('catalogue, rows', [('SIO_good', 39399), ('SIO_all', 43969),
                                             ('SIO_shallow', 1719), ('SIO_short', 477),
                                             ('SIO_tall', 2374)])
def test_sio_seamounts_keep_every_row(catalogue, rows):
    """The files have no header, but the loader used to skip 17 lines as if they did, silently
    dropping the first 17 seamounts. Row counts are the ones given in the archive's README."""
    from gprm.datasets.Seafloor import Seamounts

    gdf = Seamounts(catalogue=catalogue)
    assert len(gdf) == rows
    assert gdf.iloc[0]['Name'] == _first_sio_name(catalogue)


def _first_sio_name(catalogue):
    """The name on the first line of the raw file, read independently of the loader."""
    from gprm.datasets.Seafloor import Seamounts

    with open(Seamounts(catalogue=catalogue, load=False)) as f:
        return f.readline().split()[5]


def test_pacific_seamount_ages_2013_names_are_in_the_name_column():
    """In the source rows the name comes before the two-letter chain code, although the file's
    commented header lists them the other way round; the loader used to follow the header."""
    from gprm.datasets.Seafloor import PacificSeamountAges

    gdf = PacificSeamountAges(catalogue='2013')
    assert len(gdf) == 409
    assert gdf.Tag.str.fullmatch(r'[A-Z]{2}').all()
    assert 'Macdonald' in set(gdf.SeamountName)
    assert (gdf.loc[gdf.SeamountName == 'Macdonald', 'SeamountChain'] == 'Austral').all()


def test_pacific_seamount_ages_2021_columns_follow_the_documented_header():
    """'# Lon Lat Age Error Type Ref Name(Sample) Tag Chain' -- checked so the 2013 fix above is
    not mistaken for a problem common to both files."""
    from gprm.datasets.Seafloor import PacificSeamountAges

    gdf = PacificSeamountAges(catalogue='2021')
    assert len(gdf) == 419
    assert gdf.Tag.str.fullmatch(r'[A-Z]{2}').all()
    assert gdf.Ref.str.fullmatch(r'R\d+').all()


def test_kim_wessel_seamounts_skip_only_the_17_line_header():
    """Unlike the SIO files, this one has a real 17-line '#' header (the last line names the
    columns) plus '>' basin separators; 24643 data rows remain. (The header itself says 24,646.)"""
    from gprm.datasets.Seafloor import Seamounts

    gdf = Seamounts(catalogue='KimWessel')
    assert len(gdf) == 24643
    assert gdf.iloc[0]['ID'] == 'KW-00001'
    assert gdf.Longitude.between(-180, 180).all() and gdf.Latitude.between(-90, 90).all()
    assert (gdf.Long == gdf.Longitude).all()     # old gprm name kept as an alias


def _loader(path):
    import importlib
    module, _, function = path.rpartition('.')
    return getattr(importlib.import_module('gprm.datasets.' + module), function)


@pytest.mark.parametrize('call, kwargs, source_to_alias', [
    ('Rocks.BaseMetalDeposits', {'deposit_type': 'VMS'}, {'Lon': 'Longitude', 'Lat': 'Latitude'}),
    ('Rocks.Metamorphism', {}, {'LONGITUDE (˚E)': 'Longitude', 'LATITUDE (˚N)': 'Latitude'}),
    ('Seafloor.PacificSeamountAges', {'catalogue': '2013'},
     {'Lon': 'Long', 'Average_age(Ma)': 'Average_Age_Ma', 'Average_error(Ma)': 'Average_Age_Error_Ma',
      'Name(island,seamount,plateau_or_sample)': 'SeamountName',
      'Island_or_seamount_chain': 'SeamountChain'}),
    ('Seafloor.PacificSeamountAges', {'catalogue': '2021'},
     {'Lon': 'Long', 'Age': 'Average_Age_Ma', 'Error': 'Average_Age_Error_Ma',
      'Name(Sample)': 'SampleName', 'Chain': 'SeamountChain'}),
    ('Seafloor.Seamounts', {'catalogue': 'SIO_good'},
     {'longitude': 'Long', 'latitude': 'Lat', 'height': 'Height', 'radius': 'Radius',
      'base_depth': 'Base_Depth', 'name': 'Name'}),
    ('Zircons.loadDB', {'version': 2021},
     {'Non-iter          age                     (Ma)': 'Non_Iter_Age_Ma',
      'Est. Depos. Age': 'Est_Depos_Age_Ma'}),
    ('Zircons.get_igneous_samples', {'version': 2019},
     {'GPS Longitude': 'Longitude', '206Pb/238U Age (Ma)': '206Pb_238U_Age_Ma',
      'Uncert. (2σ).2': '207Pb_206Pb_Precis'}),
])
def test_source_columns_are_kept_and_gprm_names_are_added_alongside(call, kwargs, source_to_alias):
    """Loaders return every column under the source's own name; gprm's conventional names are
    extra columns holding the same values (see gprm/datasets/_columns.py)."""
    gdf = _loader(call)(**kwargs)
    for source, alias in source_to_alias.items():
        assert source in gdf.columns, source
        assert alias in gdf.columns, alias
        pd_testing = pytest.importorskip('pandas.testing')
        pd_testing.assert_series_equal(gdf[source], gdf[alias], check_names=False)


def test_geochem_keeps_the_source_longitude_and_wraps_only_the_alias():
    from gprm.datasets.Rocks import Geochem

    gdf = Geochem(usecols=['Longitude', 'Latitude'])
    assert gdf.longitude.max() > 180                      # 0-360 rows as in the file
    assert gdf.Longitude.between(-180, 180, inclusive='left').all()
    assert (gdf.latitude == gdf.Latitude).all()


def test_zircons_2018_source_columns_keep_their_raw_values():
    """Two source columns use decimal commas in places; the fix goes on the alias only."""
    from gprm.datasets.Zircons import loadDB

    samples, data = loadDB(version=2018)
    raw = data['207Pb /\n206Pb\nAge\n(Ma)']
    assert raw.astype(str).str.contains(',').any()
    assert data['207Pb_206Pb_Age_Ma'].dtype == float
    assert 'Est. Depos. Age (Ma)' in samples.columns and 'Est_Depos_Age_Ma' in samples.columns


def _pbdb_if_downloaded():
    """pbdb is not fetched by pooch; the user downloads it. Skip if it is not there."""
    from gprm.datasets import cache_path
    from gprm.datasets.Strat import pbdb
    if not cache_path('pbdb', 'pbdb_data.csv').exists():
        pytest.skip('no local pbdb download')
    return pbdb()


AGE_DESCRIBED_LOADERS = [
    ('Rocks.Geochem', lambda: _loader('Rocks.Geochem')()),
    ('Rocks.BaseMetalDeposits', lambda: _loader('Rocks.BaseMetalDeposits')('VMS')),
    ('Rocks.Kimberlites:Faure2010', lambda: _loader('Rocks.Kimberlites')()),
    ('Rocks.Kimberlites:Tappe2018', lambda: _loader('Rocks.Kimberlites')(catalogue='Tappe2018')),
    ('Rocks.Carbonatites', lambda: _loader('Rocks.Carbonatites')()),
    ('Rocks.Metamorphism', lambda: _loader('Rocks.Metamorphism')()),
    ('Seafloor.MagneticPicks', lambda: _loader('Seafloor.MagneticPicks')()),
    ('Seafloor.PacificSeamountAges:2013', lambda: _loader('Seafloor.PacificSeamountAges')('2013')),
    ('Seafloor.PacificSeamountAges:2021', lambda: _loader('Seafloor.PacificSeamountAges')('2021')),
    ('Seafloor.Seamounts:KimWessel', lambda: _loader('Seafloor.Seamounts')('KimWessel')),
    ('Seafloor.Seamounts:SIO', lambda: _loader('Seafloor.Seamounts')('SIO_good')),
    ('Seafloor.LargeIgneousProvinces:Whittaker', lambda: _loader('Seafloor.LargeIgneousProvinces')('Whittaker')),
    ('Seafloor.LargeIgneousProvinces:Johansson', lambda: _loader('Seafloor.LargeIgneousProvinces')('Johansson')),
    ('Seafloor.LargeIgneousProvinces:Johansson_centroids',
     lambda: _loader('Seafloor.LargeIgneousProvinces')('Johansson_centroids')),
    ('Seafloor.LargeIgneousProvinces:UTIG', lambda: _loader('Seafloor.LargeIgneousProvinces')('UTIG')),
    ('Strat.pbdb', _pbdb_if_downloaded),
    ('Strat.PaleoLithology', lambda: _loader('Strat.PaleoLithology')()),
    ('Strat.PaleoReefs', lambda: _loader('Strat.PaleoReefs')()),
    ('Geology.GlobalTectonicMap', lambda: _loader('Geology.fetch_GlobalTectonicMap')()),
    ('Geology.SurfaceGeology', lambda: _loader('Geology.fetch_SurfaceGeology')()),
    ('Zircons.loadDB:2018:samples', lambda: _loader('Zircons.loadDB')(2018)[0]),
    ('Zircons.loadDB:2018:data', lambda: _loader('Zircons.loadDB')(2018)[1]),
    ('Zircons.loadDB:2019', lambda: _loader('Zircons.loadDB')(2019)),
    ('Zircons.loadDB:2021', lambda: _loader('Zircons.loadDB')(2021)),
    ('Zircons.loadDB:2024', lambda: _loader('Zircons.loadDB')(2024)),
    ('Zircons.loadDB:2026', lambda: _loader('Zircons.loadDB')(2026)),
    ('Zircons.load_Hf', lambda: _loader('Zircons.load_Hf')()),
    ('Zircons.get_igneous_samples:2018', lambda: _loader('Zircons.get_igneous_samples')(version=2018)),
    ('Zircons.get_igneous_samples:2019', lambda: _loader('Zircons.get_igneous_samples')(version=2019)),
    ('Zircons.get_igneous_samples:2026', lambda: _loader('Zircons.get_igneous_samples')(version=2026)),
    ('Zircons.get_mafic_felsic_samples', lambda: _loader('Zircons.get_mafic_felsic_samples')('Felsic')),
    ('Zircons.get_sedimentary_samples:2018', lambda: _loader('Zircons.get_sedimentary_samples')(version=2018)),
]


@pytest.mark.parametrize('dataset, load', AGE_DESCRIBED_LOADERS, ids=[d for d, _ in AGE_DESCRIBED_LOADERS])
def test_loader_stamps_an_age_description_naming_only_real_columns(dataset, load):
    """Each loader attaches its entry from gprm/datasets/_ages.py, and every column that entry
    names is present in the table -- so a consumer that trusts the stamp finds what it expects."""
    df = load()
    description = df.attrs.get('gprm_age')
    assert description is not None, 'no age description stamped'
    assert description['dataset'] == dataset

    named = [description['field'], description['source_field'], description['uncertainty']]
    named += list(description['range'] or ()) + list(description['other_ages'])
    missing = [c for c in named if c is not None and c not in df.columns]
    assert not missing, 'described but not in the table: {}'.format(missing)


def test_geochem_age_is_cleaned_without_touching_the_source_or_dropping_rows():
    """2 ages exceed the age of the Earth and 3 are large negatives: NaN in 'Age'. 119 are
    historical eruptions entered as tiny negative ages: 0 in 'Age'. 'age' is left as it was."""
    import pandas as pd
    from gprm.datasets.Rocks import Geochem

    gdf = Geochem(usecols=['Longitude', 'Latitude', 'age'])
    source = pd.to_numeric(gdf['age'], errors='coerce')

    assert gdf.Age.dropna().between(0, 4567).all()
    assert int(gdf.Age.isna().sum() - source.isna().sum()) == 5
    assert int(((source < 0) & (gdf.Age == 0)).sum()) == 119
    assert source.max() == 127000             # the source column is not cleaned
    assert gdf.attrs['gprm_age']['field'] == 'Age'


def test_valdes2021_fetches_every_run_and_checks_its_checksum():
    """109 BRIDGE runs, one annual-mean file each, all verified against pinned sha256."""
    import xarray as xr
    from gprm.datasets.Paleogeography import fetch_Valdes2021

    runs = fetch_Valdes2021()
    assert len(runs) == 109
    assert min(runs) == 0 and max(runs) == 541
    with xr.open_dataset(runs[0.0], decode_times=False) as ds:
        assert 'temp_mm_1_5m' in ds and 'precip_mm_srf' in ds


def test_lihu2022_is_the_figshare_file():
    """754 MB download on first run. Checked against figshare's own published md5."""
    from gprm.datasets.Paleogeography import fetch_LiHu2022

    ds = fetch_LiHu2022(return_xarray=True)
    assert ds.sizes['simulation'] == 55
    assert {'T', 'P'} <= set(ds.data_vars)


def test_paleoreefs_is_pared_version_1():
    from gprm.datasets.Strat import PaleoReefs

    gdf = PaleoReefs()
    assert len(gdf) == 4363
    assert gdf['r_number'].is_unique
    assert {'latit', 'longit', 'intervall', 'biota_main_t', 'pal_lat_scotese'} <= set(gdf.columns)
    assert gdf.crs.to_epsg() == 4326

