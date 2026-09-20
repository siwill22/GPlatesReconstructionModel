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
