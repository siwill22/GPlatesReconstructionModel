"""Shared fixtures.

These tests deliberately avoid the network. The reconstruction-model fetchers are covered by
a separate, network-marked test, because a failure there means an upstream URL has rotted
rather than that gprm is broken, and the two should not look the same in CI.
"""
import textwrap

import pytest


# A minimal but genuine GPlates rotation file: plate 701 moving relative to the anchor plate
# 000, with a 15 degree rotation accumulated by 100 Ma.
MINIMAL_ROTATION_FILE = textwrap.dedent("""\
    701   0.0   0.0   0.0   0.0  000 ! test
    701 100.0  20.0  30.0  15.0  000 ! test
    801   0.0   0.0   0.0   0.0  000 ! test
    801 100.0 -10.0  40.0  -8.0  000 ! test
    """)


@pytest.fixture
def rotation_file(tmp_path):
    """Path to a minimal valid .rot file defining plate ids 701 and 801."""
    path = tmp_path / 'test_model.rot'
    path.write_text(MINIMAL_ROTATION_FILE)
    return str(path)


@pytest.fixture
def reconstruction_model(rotation_file):
    """A ReconstructionModel with a rotation model but no polygons."""
    from gprm import ReconstructionModel

    model = ReconstructionModel('TestModel')
    model.add_rotation_model(rotation_file)
    return model


@pytest.fixture
def sample_points():
    """A GeoDataFrame of age-coded sample points, of the kind a user would bring."""
    import geopandas as gpd
    import pandas as pd

    df = pd.DataFrame({
        'sample_id': ['a', 'b', 'c'],
        'lon': [20.0, -60.0, 130.0],
        'lat': [-25.0, -5.0, 35.0],
        'FROMAGE': [100.0, 80.0, 60.0],
        'TOAGE': [90.0, 70.0, 50.0],
    })
    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs=4326)
