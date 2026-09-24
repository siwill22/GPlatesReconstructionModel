"""The per-dataset age descriptions in gprm/datasets/_ages.py, checked without a network.

Whether each description matches the columns a loader actually returns is checked in
test_dataset_loaders.py, which needs the data.
"""
import re
from pathlib import Path

import pandas as pd
import pytest

from gprm.datasets import _ages

DATASETS_DIR = Path(__file__).parent.parent / 'gprm' / 'datasets'


@pytest.mark.parametrize('name', _ages.datasets())
def test_every_description_is_well_formed(name):
    d = _ages.age_description(name)

    assert d['dataset'] == name
    assert d['field'] is None or isinstance(d['field'], str)
    assert d['event'] is None or d['event'] in _ages.EVENTS
    assert d['validity_window'] in _ages.VALIDITY_WINDOWS
    assert d['range'] is None or (isinstance(d['range'], tuple) and len(d['range']) == 2)
    assert isinstance(d['other_ages'], dict)
    if d['source_field'] is not None:
        assert d['field'] is not None, 'source_field without the derived field it feeds'
    # A dataset whose age lives in FROMAGE must say what FROMAGE holds
    if d['field'] == 'FROMAGE':
        assert d['validity_window'] == 'age to present'


def test_age_description_returns_a_copy():
    d = _ages.age_description('Rocks.Carbonatites')
    d['field'] = 'something else'
    d['other_ages']['x'] = 'y'
    assert _ages.age_description('Rocks.Carbonatites')['field'] == 'Age'
    assert 'x' not in _ages.age_description('Rocks.Carbonatites')['other_ages']


def test_unknown_dataset_names_the_known_ones():
    with pytest.raises(KeyError, match='Rocks.Carbonatites'):
        _ages.age_description('Rocks.NotADataset')


def test_public_accessor():
    from gprm.datasets import age_description
    assert age_description('Rocks.Carbonatites')['event'] == 'emplacement'


def test_stamp_attaches_the_description_to_attrs():
    df = pd.DataFrame({'Age': [1.0]})
    assert _ages.stamp(df, 'Rocks.Carbonatites') is df
    assert df.attrs['gprm_age']['field'] == 'Age'


def test_every_stamp_in_the_loaders_names_a_registered_dataset():
    """Catches a typo in a loader's _stamp(...) key without downloading anything. Keys built at
    run time ('...:' + catalogue, '...:{}'.format(version)) are checked by their fixed prefix."""
    registered = set(_ages.datasets())
    found = 0
    for path in DATASETS_DIR.glob('*.py'):
        text = path.read_text()
        for key in (_stamp_key(call) for call in _stamp_calls(text)):
            found += 1
            if key.endswith(':') or '{}' in key:
                prefix = key.split('{}')[0]
                assert any(r.startswith(prefix) for r in registered), (path.name, key)
            else:
                assert key in registered, (path.name, key)
    assert found > 20


def _stamp_calls(text):
    """The full text of every _stamp(...) call, however it is wrapped over lines."""
    for match in re.finditer(r'_stamp\(', text):
        depth, i = 1, match.end()
        while depth:
            depth += {'(': 1, ')': -1}.get(text[i], 0)
            i += 1
        yield text[match.end():i - 1]


def _stamp_key(call):
    keys = re.findall(r"'((?:Rocks|Seafloor|Strat|Geology|Zircons)\.[^']*)'", call)
    assert len(keys) == 1, call
    return keys[0]


# --- resolve_age_field: which column a consumer takes the sample age from -------------------

def _carbonatite_like():
    """A table shaped like Rocks.Carbonatites: FROMAGE is age + error, not the age."""
    import geopandas as gpd
    df = pd.DataFrame({'Age_ma': [100.0, 50.0], 'Age': [100.0, 50.0], 'Error_ma': [5.0, 2.0],
                       'FROMAGE': [105.0, 52.0], 'TOAGE': [95.0, 48.0], 'PLATEID1': [701, 701]})
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy([10., 20.], [0., 5.]), crs=4326)
    return _ages.stamp(gdf, 'Rocks.Carbonatites')


def test_an_explicit_age_field_wins_over_the_stamp():
    assert _ages.resolve_age_field(_carbonatite_like(), 'FROMAGE', default='age') == 'FROMAGE'


def test_the_stamp_is_used_when_no_field_is_given(recwarn):
    assert _ages.resolve_age_field(_carbonatite_like(), None, default='FROMAGE') == 'Age'
    assert not [w for w in recwarn if 'No age description' in str(w.message)]


def test_without_a_stamp_the_default_is_used_with_a_warning():
    df = pd.DataFrame({'age': [1.0]})
    with pytest.warns(UserWarning, match="'age' is assumed"):
        assert _ages.resolve_age_field(df, None, default='age') == 'age'


def test_a_stamp_lost_in_a_merge_falls_back_with_a_warning():
    """merge() drops attrs -- the case the warning exists for."""
    gdf = _carbonatite_like().merge(pd.DataFrame({'PLATEID1': [701], 'x': [1]}), on='PLATEID1')
    assert 'gprm_age' not in gdf.attrs
    with pytest.warns(UserWarning, match='No age description'):
        assert _ages.resolve_age_field(gdf, None, default='FROMAGE') == 'FROMAGE'


def test_a_dataset_with_no_single_sample_age_raises_rather_than_guessing():
    df = _ages.stamp(pd.DataFrame({'min_ma': [1.0], 'max_ma': [2.0]}), 'Strat.pbdb')
    with pytest.raises(ValueError, match="no single sample age.*min_ma"):
        _ages.resolve_age_field(df, None, default='FROMAGE')


def test_a_stamp_naming_a_dropped_column_raises():
    df = _ages.stamp(pd.DataFrame({'Error_ma': [1.0]}), 'Rocks.Carbonatites')
    with pytest.raises(ValueError, match="'Age', which is not in this table"):
        _ages.resolve_age_field(df, None)


def test_reconstruct_to_time_of_appearance_uses_the_stamped_age_not_fromage(reconstruction_model):
    """End to end on the case that motivated all this: for carbonatites FROMAGE is age + error,
    so reconstructing 'to the time of appearance' by FROMAGE lands every sample at the wrong
    time. With the stamp it goes to Age; stripping the stamp reverts to FROMAGE, with a warning."""
    gdf = _carbonatite_like()

    stamped = reconstruction_model.reconstruct_to_time_of_appearance(gdf)
    by_age = reconstruction_model.reconstruct_to_time_of_appearance(gdf, ReconstructTime='Age')
    assert stamped.geometry.geom_equals(by_age.geometry).all()

    unstamped = gdf.copy()
    unstamped.attrs = {}
    with pytest.warns(UserWarning, match='No age description'):
        by_fromage = reconstruction_model.reconstruct_to_time_of_appearance(unstamped)
    explicit_fromage = reconstruction_model.reconstruct_to_time_of_appearance(
        gdf, ReconstructTime='FROMAGE')
    assert by_fromage.geometry.geom_equals(explicit_fromage.geometry).all()
    assert not stamped.geometry.geom_equals(by_fromage.geometry).any()
