"""ReconstructionModel: catching wrong inputs, one test per fixed defect.

A full audit of this class (its 21 public methods, cross-checked against every ReconstructionModel
call site in the downstream project Geode) found several of its most-used methods silently
returning None or silently producing plausible-looking wrong output for invalid input, and one
live data-corruption bug in copy(). Every case here was reproduced against the code as it stood
before the fix; the docstrings record what the old behaviour was.
"""
import pygplates
import geopandas as gpd
import pytest
from shapely.geometry import Point

from gprm import ReconstructionModel


def test_copy_shallow_does_not_share_mutable_state(reconstruction_model, rotation_file):
    """copy(deep=False) used to alias the list attributes via copy.copy(self); calling any
    add_* method on the copy silently mutated the original's list in place.
    """
    copy_model = reconstruction_model.copy(deep=False)
    copy_model.add_rotation_model(rotation_file)

    assert len(copy_model.rotation_files) == 2
    assert len(reconstruction_model.rotation_files) == 1


def test_info_does_not_crash_on_default_constructor_state():
    """ReconstructionModel().info() used to raise TypeError from '{:s}'.format(None), since
    name defaults to None -- the most basic possible usage sequence.
    """
    model = ReconstructionModel()
    model.info()  # must not raise
    assert '<unnamed>' in repr(model)


def test_reconstruct_rejects_wrong_type_instead_of_returning_none(reconstruction_model):
    """reconstruct's isinstance dispatch had no final else, so any features type other than
    FeatureCollection/GeoDataFrame (including None) silently returned None.
    """
    with pytest.raises(TypeError, match='FeatureCollection.*GeoDataFrame'):
        reconstruction_model.reconstruct('not a feature collection', 100.0)
    with pytest.raises(TypeError):
        reconstruction_model.reconstruct(None, 100.0)


def test_reconstruct_to_time_of_appearance_rejects_wrong_type(reconstruction_model):
    """Same fallthrough-to-None bug as reconstruct, on the sibling method."""
    with pytest.raises(TypeError, match='FeatureCollection.*GeoDataFrame'):
        reconstruction_model.reconstruct_to_time_of_appearance('nope')


def test_assign_plate_ids_rejects_unknown_polygons_value(reconstruction_model):
    """assign_plate_ids' polygons parameter used to silently fall back to 'static' for any
    unrecognised string -- the identical three-way check on the sibling method
    polygon_snapshot was already fixed to raise; this one wasn't.
    """
    gdf = gpd.GeoDataFrame(geometry=[Point(0., 0.)])
    with pytest.raises(ValueError, match="'static', 'coastlines', 'continents'"):
        reconstruction_model.assign_plate_ids(gdf, polygons='statc')


def test_assign_plate_ids_names_the_missing_polygon_type(reconstruction_model):
    """'No polygons found for partitioning' used to not say which polygon type was empty --
    which matters once the typo above no longer silently redirects to 'static'.
    """
    gdf = gpd.GeoDataFrame(geometry=[Point(0., 0.)])
    with pytest.raises(ValueError, match='coastlines'):
        reconstruction_model.assign_plate_ids(gdf, polygons='coastlines')


def test_assign_plate_ids_rejects_wrong_type(reconstruction_model):
    """The final else used to raise a vague ValueError with no mention of what was received."""
    with pytest.raises(TypeError, match='FeatureCollection.*GeoDataFrame'):
        reconstruction_model.assign_plate_ids('nope')


@pytest.mark.parametrize('add_method,argname', [
    ('add_rotation_model', 'rotation_file'),
    ('add_static_polygons', 'static_polygons_file'),
    ('add_dynamic_polygons', 'dynamic_polygons_file'),
    ('add_coastlines', 'coastlines_file'),
    ('add_continent_polygons', 'continent_polygons_file'),
])
def test_add_methods_reject_non_string_path_with_a_named_error(add_method, argname):
    """All five add_* methods built their friendly 'Unable to find file' message with
    '{:s}'.format(x) -- which itself raised a different, unrelated ValueError/TypeError
    whenever x wasn't a string, so the intended message never reached the caller.
    """
    model = ReconstructionModel()
    method = getattr(model, add_method)

    with pytest.raises(TypeError, match=argname):
        method(None)
    with pytest.raises(TypeError, match=argname):
        method(123)


def test_add_rotation_model_still_names_a_missing_file(reconstruction_model):
    """The one case the old check handled correctly must keep working."""
    with pytest.raises(ValueError, match='nonexistent_file.rot'):
        reconstruction_model.add_rotation_model('nonexistent_file.rot')


def test_reconstruct_topological_raises_rather_than_returning_none(reconstruction_model):
    """topological=True used to '#TODO perform a topological reconstruction' then bare
    `return`, silently giving back None -- an un-migrated instance of a bug pattern this
    project has already fixed everywhere else it appeared.
    """
    with pytest.raises(NotImplementedError):
        reconstruction_model.reconstruct(pygplates.FeatureCollection(), 100.0, topological=True)


def test_reconstruct_use_tempfile_branches_agree_on_empty_result(reconstruction_model):
    """The use_tempfile=False branch returned bare None for 'nothing reconstructed', while
    use_tempfile=True returned an empty, correctly-typed GeoDataFrame -- so the two branches
    of the same method, selected only by a flag, disagreed on the result type.
    """
    gdf = gpd.GeoDataFrame({'PLATEID1': [701]}, geometry=[Point(0., 0.)])
    gdf.attrs['_gprm_plate_id_provenance'] = reconstruction_model.name
    # FROMAGE/TOAGE window that excludes reconstruction_time -> nothing reconstructed
    gdf['FROMAGE'] = 10.
    gdf['TOAGE'] = 0.

    result = reconstruction_model.reconstruct(gdf, 500.0, use_tempfile=False)

    assert isinstance(result, gpd.GeoDataFrame)
    assert len(result) == 0


def test_plate_snapshot_rejects_non_numeric_reconstruction_time(reconstruction_model):
    reconstruction_model.dynamic_polygons = []
    with pytest.raises(TypeError, match='reconstruction_time'):
        reconstruction_model.plate_snapshot('100')


def test_plate_snapshot_rejects_unknown_anchor_plate_id(reconstruction_model):
    """An anchor_plate_id not defined anywhere in the rotation model used to silently produce
    an identity-rotated (i.e. meaningless) snapshot instead of raising.
    """
    reconstruction_model.dynamic_polygons = []
    with pytest.raises(ValueError, match='999999'):
        reconstruction_model.plate_snapshot(50.0, anchor_plate_id=999999)


def test_plate_snapshot_accepts_the_default_anchor_plate_id(reconstruction_model):
    """anchor_plate_id=0 (Geode's only real usage pattern) must still work -- 0 is the
    always-valid anchor/unpartitioned-feature id, not necessarily a rotation-file entry.
    """
    reconstruction_model.dynamic_polygons = []
    snapshot = reconstruction_model.plate_snapshot(50.0)
    assert snapshot.plate_count == 0


def test_polygon_snapshot_rejects_unknown_anchor_plate_id(reconstruction_model):
    reconstruction_model.static_polygons = []
    with pytest.raises(ValueError, match='999999'):
        reconstruction_model.polygon_snapshot('static_polygons', 50.0, anchor_plate_id=999999)


def test_reconstruct_rejects_unknown_plate_id_in_feature_collection(reconstruction_model):
    """The FeatureCollection branch of reconstruct never called any plate-id guard at all
    (only the GeoDataFrame branch did) -- an unknown reconstruction plate id silently
    round-tripped through pygplates' identity-rotation fallback as a 'successful' result.
    """
    feature = pygplates.Feature()
    feature.set_geometry(pygplates.PointOnSphere(10., 10.))
    feature.set_reconstruction_plate_id(88888)
    fc = pygplates.FeatureCollection([feature])

    with pytest.raises(ValueError, match='88888'):
        reconstruction_model.reconstruct(fc, 50.0)


def test_reconstruct_to_time_of_appearance_rejects_unknown_plate_id_in_feature_collection(reconstruction_model):
    """Same gap as above, on reconstruct_to_time_of_appearance's FeatureCollection branch."""
    feature = pygplates.Feature()
    feature.set_geometry(pygplates.PointOnSphere(10., 10.))
    feature.set_reconstruction_plate_id(88888)
    feature.set_valid_time(50., 0.)
    fc = pygplates.FeatureCollection([feature])

    with pytest.raises(ValueError, match='88888'):
        reconstruction_model.reconstruct_to_time_of_appearance(fc)


def test_reconstruct_rejects_a_projected_crs(reconstruction_model):
    """A GeoDataFrame in a projected CRS was silently treated as lon/lat degrees -- large
    magnitudes happened to raise deep inside pygplates with no gprm-level context, smaller
    ones would silently reconstruct nonsense geometry.
    """
    gdf = gpd.GeoDataFrame({'PLATEID1': [701]}, geometry=[Point(500000., 4000000.)],
                           crs='EPSG:32633')
    gdf.attrs['_gprm_plate_id_provenance'] = reconstruction_model.name

    with pytest.raises(ValueError, match='projected CRS'):
        reconstruction_model.reconstruct(gdf, 50.0)


def test_reconstruct_accepts_a_geographic_crs(reconstruction_model):
    """The CRS check must not reject ordinary lon/lat data."""
    gdf = gpd.GeoDataFrame({'PLATEID1': [701]}, geometry=[Point(10., 10.)], crs='EPSG:4326')
    gdf.attrs['_gprm_plate_id_provenance'] = reconstruction_model.name

    result = reconstruction_model.reconstruct(gdf, 50.0)
    assert len(result) == 1


def test_reconstruct_accepts_no_crs_at_all(reconstruction_model):
    """Plenty of legitimate GeoDataFrames have no CRS set; the check must only fire when a
    CRS is actually present.
    """
    gdf = gpd.GeoDataFrame({'PLATEID1': [701]}, geometry=[Point(10., 10.)])
    gdf.attrs['_gprm_plate_id_provenance'] = reconstruction_model.name

    result = reconstruction_model.reconstruct(gdf, 50.0)
    assert len(result) == 1


def test_from_agegrid_config_names_the_missing_key(tmp_path):
    """A malformed YAML config used to let a raw, unattributed KeyError through."""
    config_path = tmp_path / 'bad_config.yaml'
    config_path.write_text('InputFiles:\n  MODELDIR: some_model\n')

    with pytest.raises(ValueError, match='input_rotation_filenames'):
        ReconstructionModel.from_agegrid_config(str(config_path))


def test_from_web_service_names_the_missing_dependency():
    """gwsFeatureCollection is not on PyPI at all, so this always raised ModuleNotFoundError
    with no indication that the dependency is unobtainable via pip -- now an ImportError
    that says so.
    """
    model = ReconstructionModel()
    with pytest.raises(ImportError, match='not on PyPI'):
        model.from_web_service()


def test_construct_topological_model_exists_and_is_callable(reconstruction_model):
    """Previously only defined on the class at all if pygplates >= 32 at import time, so
    calling it on an older install raised a plain AttributeError with no indication that
    the real cause was the pygplates version.
    """
    reconstruction_model.dynamic_polygons = []
    # Either succeeds (pygplates >= 32, as in this environment) or raises a clear,
    # version-naming RuntimeError (pygplates < 32) -- never AttributeError.
    try:
        reconstruction_model.construct_topological_model()
    except RuntimeError as error:
        assert 'pygplates >= 32' in str(error)


def test_force_polygon_geometries_warns_about_dropped_features():
    """Features with a reversed (FROMAGE < TOAGE) valid time were silently dropped with no
    count or warning -- the exact real-world case COBTerranes data hits, per the function's
    own comment. pygplates' own set_valid_time refuses to construct a reversed-time feature
    directly, so this proxies get_valid_time() on a real feature to force the case gprm's own
    code has to handle when reading such data from a file.
    """
    from gprm.utils.spatial import force_polygon_geometries

    class _ReversedValidTime:
        def __init__(self, feature):
            self._feature = feature
        def get_all_geometries(self):
            return self._feature.get_all_geometries()
        def get_feature_type(self):
            return self._feature.get_feature_type()
        def get_reconstruction_plate_id(self):
            return self._feature.get_reconstruction_plate_id()
        def get_valid_time(self):
            return (0., 10.)  # reversed: FROMAGE (0.) < TOAGE (10.)

    good = pygplates.Feature()
    good.set_geometry(pygplates.PolylineOnSphere([(0., 0.), (0., 10.), (10., 10.), (10., 0.)]))
    good.set_valid_time(10., 0.)

    bad_source = pygplates.Feature()
    bad_source.set_geometry(pygplates.PolylineOnSphere([(20., 20.), (20., 30.), (30., 30.), (30., 20.)]))
    bad_source.set_valid_time(10., 0.)  # valid at construction; get_valid_time is overridden below

    with pytest.warns(UserWarning, match='1 feature'):
        result = force_polygon_geometries([good, _ReversedValidTime(bad_source)])

    assert len(list(result)) == 1


@pytest.fixture
def model_with_dateline_polygon(rotation_file, tmp_path):
    """A model whose single static polygon spans the antimeridian: 20 degrees wide, from
    lon 170 to lon -170. On the sphere it contains lon 180 and not lon 0. Read as planar
    lon/lat it is instead 340 degrees wide and contains exactly the opposite points, which
    is what makes it separate a spherical containment test from a shapely one.
    """
    from gprm import ReconstructionModel

    feature = pygplates.Feature()
    feature.set_geometry(pygplates.PolygonOnSphere(
        [(-10., 170.), (10., 170.), (10., -170.), (-10., -170.)]))
    feature.set_reconstruction_plate_id(701)
    feature.set_valid_time(600., -999.)

    path = tmp_path / 'dateline_polygon.gpml'
    pygplates.FeatureCollection([feature]).write(str(path))

    model = ReconstructionModel('DatelineTest')
    model.add_rotation_model(rotation_file)
    model.add_static_polygons(str(path))
    return model


def test_assign_plate_ids_tests_containment_on_the_sphere(model_with_dateline_polygon):
    """The GeoDataFrame path used to test containment with geopandas' .overlay(), which is
    planar, so a polygon crossing the antimeridian captured the points on the far side of
    the globe and missed the ones it actually contains.
    """
    gdf = gpd.GeoDataFrame(geometry=[Point(180., 0.), Point(0., 0.)], crs='EPSG:4326')

    result = model_with_dateline_polygon.assign_plate_ids(gdf, polygons='static')

    assert list(result['PLATEID1']) == [701, 0]


def test_assign_plate_ids_preserves_row_count_and_order(model_with_dateline_polygon):
    """.overlay(how='intersection') silently dropped unpartitioned rows whatever
    keep_unpartitioned_features said, and emitted one row per match, so a point falling in
    two polygons was duplicated. Neither can happen now: one row in, one row out.
    """
    gdf = gpd.GeoDataFrame(geometry=[Point(0., 0.), Point(180., 0.), Point(0., 0.)],
                           crs='EPSG:4326')

    result = model_with_dateline_polygon.assign_plate_ids(gdf, polygons='static')
    assert len(result) == 3
    assert list(result['PLATEID1']) == [0, 701, 0]

    dropped = model_with_dateline_polygon.assign_plate_ids(
        gdf, polygons='static', keep_unpartitioned_features=False)
    assert len(dropped) == 1
    assert dropped['PLATEID1'].iloc[0] == 701


def test_assign_plate_ids_overlay_method_still_works_but_warns(model_with_dateline_polygon):
    """The planar path is kept so earlier results can be reproduced, but it has to announce
    itself -- and it still gets the antimeridian exactly backwards, which is the point.
    """
    gdf = gpd.GeoDataFrame(geometry=[Point(180., 0.), Point(0., 0.)], crs='EPSG:4326')

    with pytest.warns(FutureWarning, match='plane'):
        result = model_with_dateline_polygon.assign_plate_ids(
            gdf, polygons='static', method='overlay')

    # Planar: it keeps the point the polygon does not contain, and drops the one it does.
    assert len(result) == 1
    assert result.geometry.x.iloc[0] == 0.


def test_assign_plate_ids_rejects_unknown_method(reconstruction_model):
    """A mistyped method must not silently fall back to either implementation."""
    gdf = gpd.GeoDataFrame(geometry=[Point(0., 0.)])
    with pytest.raises(ValueError, match="'spatial_tree', 'overlay'"):
        reconstruction_model.assign_plate_ids(gdf, method='tree')


def test_assign_plate_ids_uses_every_geometry_of_a_multi_geometry_feature(rotation_file, tmp_path):
    """Feature.get_geometry() returns None when a feature holds more than one geometry, so
    selecting polygons through it dropped every multi-geometry feature -- 279 of Torsvik &
    Cocks 2017's 600 continent polygons, which left more than half its points unpartitioned.
    """
    from gprm import ReconstructionModel

    feature = pygplates.Feature()
    feature.set_geometry([
        pygplates.PolygonOnSphere([(-10., 0.), (10., 0.), (10., 20.), (-10., 20.)]),
        pygplates.PolygonOnSphere([(-10., 100.), (10., 100.), (10., 120.), (-10., 120.)]),
    ])
    feature.set_reconstruction_plate_id(701)
    feature.set_valid_time(600., -999.)
    assert feature.get_geometry() is None  # the trap this guards against

    path = tmp_path / 'multi_geometry.gpml'
    pygplates.FeatureCollection([feature]).write(str(path))

    model = ReconstructionModel('MultiGeometryTest')
    model.add_rotation_model(rotation_file)
    model.add_static_polygons(str(path))

    # One point in each of the feature's two polygons, and one outside both.
    gdf = gpd.GeoDataFrame(geometry=[Point(10., 0.), Point(110., 0.), Point(60., 0.)],
                           crs='EPSG:4326')

    result = model.assign_plate_ids(gdf, polygons='static')

    assert list(result['PLATEID1']) == [701, 701, 0]
