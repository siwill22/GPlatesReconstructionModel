"""One test per defect fixed, so that each stays fixed.

Every case here was reproduced against the code as it stood before the fix; the docstrings
record what the old behaviour was, since that is the part that is hard to recover later.
"""
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, MultiPoint, Point, Polygon


# ---------------------------------------------------------------- dateline wrapping

def test_dateline_split_polyline_keeps_every_piece():
    """Previously returned wrapped[0] only, silently discarding the rest of the line."""
    from gprm.utils.geometry import wrap_polyline_feature

    crossing = pd.Series({'geometry': LineString([(170, -10), (-170, 10)])})

    result = wrap_polyline_feature(crossing)

    assert result.geom_type == 'MultiLineString'
    assert len(result.geoms) == 2


def test_polyline_not_crossing_dateline_is_unchanged():
    from gprm.utils.geometry import wrap_polyline_feature

    ordinary = pd.Series({'geometry': LineString([(10, -10), (20, 10)])})

    assert wrap_polyline_feature(ordinary).geom_type == 'LineString'


def test_dateline_split_polygon_keeps_every_part():
    """Previously printed a warning and returned only the first part."""
    from gprm.utils.geometry import wrap_polygon_feature

    straddling = pd.Series({
        'geometry': Polygon([(170, -10), (-170, -10), (-170, 10), (170, 10)])})

    result = wrap_polygon_feature(straddling)

    assert result.geom_type == 'MultiPolygon'
    assert len(result.geoms) == 2


# ---------------------------------------------------------------- distance to features

def test_empty_feature_list_raises_rather_than_crashing():
    """Previously: TypeError: unsupported operand type(s) for *: 'NoneType' and 'float'.

    nearest_feature returns None for an empty feature set, and the result was multiplied by
    the Earth's radius without being checked.
    """
    import pygplates

    from gprm.utils.geometry import distance_between_reconstructed_points_and_features

    class FakeReconstructedGeometry:
        def get_reconstructed_geometry(self):
            return pygplates.PointOnSphere(10.0, 20.0)

    with pytest.raises(ValueError, match='distance'):
        distance_between_reconstructed_points_and_features([FakeReconstructedGeometry()], [])


# ---------------------------------------------------------------- multipart geometries

def test_apply_reconstruction_does_not_silently_return_none_for_multipart():
    """Previously: an if/elif over Point/LineString/Polygon with no else, so every multipart
    geometry came back as None with no warning. Multipart geometries are common in real
    shapefiles, so this quietly emptied the geometry column."""
    import pygplates

    from gprm.utils.geometry import apply_reconstruction

    rotation_model = pygplates.RotationModel(
        [pygplates.Feature(pygplates.FeatureType.gpml_total_reconstruction_sequence)])
    row = pd.Series({'geometry': MultiPoint([(0, 0), (1, 1)]),
                     'reconstruction_time': 10.0, 'PLATEID1': 0})

    result = apply_reconstruction(row, rotation_model)

    assert result is not None, 'multipart geometry silently became None'


# ---------------------------------------------------------------- plate id guards

def test_reconstruct_without_plate_ids_names_the_fix(reconstruction_model, sample_points):
    """Previously a KeyError from inside a DataFrame.apply, which did not say what to do."""
    with pytest.raises(ValueError, match='assign_plate_ids'):
        reconstruction_model.reconstruct(sample_points, 100.0)


def test_plate_ids_from_another_model_are_refused(reconstruction_model, sample_points):
    gdf = sample_points.copy()
    gdf['PLATEID1'] = 701
    gdf.attrs['gprm_reconstruction_model'] = 'SomeOtherModel'

    with pytest.raises(ValueError, match='SomeOtherModel'):
        reconstruction_model.reconstruct(gdf, 100.0)


def test_unknown_plate_ids_are_refused_when_provenance_was_lost(reconstruction_model,
                                                                sample_points):
    """The provenance stamp lives in GeoDataFrame.attrs, which pandas does not carry through
    every operation, so the plate ids themselves are checked as a fallback."""
    gdf = sample_points.copy()
    gdf['PLATEID1'] = 999901          # not defined in the test rotation file
    gdf.attrs.clear()

    with pytest.raises(ValueError, match='999901'):
        reconstruction_model.reconstruct(gdf, 100.0)


def test_genuine_plate_ids_pass_without_provenance(reconstruction_model, sample_points):
    """The fallback must not produce false positives on legitimate data."""
    gdf = sample_points.copy()
    gdf['PLATEID1'] = 701             # defined in the test rotation file
    gdf.attrs.clear()

    result = reconstruction_model.reconstruct(gdf, 100.0)

    assert result is not None and len(result) > 0


def test_known_plate_ids_reads_the_rotation_files(reconstruction_model):
    assert reconstruction_model.known_plate_ids() == {0, 701, 801}


def test_known_plate_ids_is_none_without_rotation_files():
    """Provenance cannot be checked for a model built without files, and must not pretend to."""
    from gprm import ReconstructionModel

    assert ReconstructionModel('empty').known_plate_ids() is None


# ---------------------------------------------------------------- age handling

def test_midtime_refuses_the_distant_future_sentinel(reconstruction_model, sample_points):
    """Previously: TOAGE of -999 (the GPlates distant-future sentinel) was averaged into
    MidTime, giving e.g. (100 + -999)/2 = -449.5 Ma, a reconstruction age in the future."""
    gdf = sample_points.copy()
    gdf['PLATEID1'] = 701
    gdf['TOAGE'] = -999.0

    with pytest.raises(ValueError, match='-999'):
        reconstruction_model.reconstruct_to_time_of_appearance(gdf, ReconstructTime='MidTime')


def test_midtime_works_with_real_end_ages(reconstruction_model, sample_points):
    gdf = sample_points.copy()
    gdf['PLATEID1'] = 701

    result = reconstruction_model.reconstruct_to_time_of_appearance(gdf,
                                                                    ReconstructTime='MidTime')

    assert list(result['reconstruction_time']) == [95.0, 75.0, 55.0]


def test_missing_age_column_is_reported_clearly(reconstruction_model, sample_points):
    gdf = sample_points.drop(columns=['FROMAGE', 'TOAGE'])
    gdf['PLATEID1'] = 701

    with pytest.raises(ValueError, match='FROMAGE'):
        reconstruction_model.reconstruct_to_time_of_appearance(gdf)


def test_attributes_survive_reconstruction(reconstruction_model, sample_points):
    """The point of working in GeoDataFrames: the sample's own columns must come back."""
    gdf = sample_points.copy()
    gdf['PLATEID1'] = 701

    result = reconstruction_model.reconstruct_to_time_of_appearance(gdf,
                                                                    ReconstructTime='FROMAGE')

    assert list(result['sample_id']) == ['a', 'b', 'c']
    assert 'geometry' in result.columns


# ---------------------------------------------------------------- spherical distance cutoff

@pytest.mark.parametrize('cutoff_degrees,should_find', [(5.0, True), (3.0, False)])
def test_distance_upper_bound_is_interpreted_as_an_angle(cutoff_degrees, should_find):
    """Previously the cutoff was converted from degrees to radians and handed to a KD-tree
    that measures chord length on the unit sphere, so the threshold was too generous by
    ~1.2% at 30 degrees and ~9% at 90."""
    from gprm.utils.sphere import sampleOnSphere

    source_lons = np.arange(-180.0, 180.0, 10.0)
    source_lats = np.zeros_like(source_lons)
    values = np.arange(source_lons.size, dtype=float)

    # the query sits 4 degrees from the nearest source
    distances, _ = sampleOnSphere(source_lons, source_lats, values,
                                  np.array([4.0]), np.array([0.0]),
                                  k=1, distance_upper_bound=cutoff_degrees)

    assert bool(np.isfinite(distances[0])) == should_find


def test_reconstruction_actually_moves_the_points(reconstruction_model, sample_points):
    gdf = sample_points.copy()
    gdf['PLATEID1'] = 701

    result = reconstruction_model.reconstruct_to_time_of_appearance(gdf,
                                                                    ReconstructTime='FROMAGE')

    moved = [not point.equals(original)
             for point, original in zip(result.geometry, gdf.geometry)]
    assert any(moved)


# ---------------------------------------------- non-exhaustive if/elif chains

def test_unknown_method_names_the_valid_options():
    """Previously fell off the end of an if/elif chain and raised UnboundLocalError: mask.
    This one is on Geode's path (prep_oldmap.py calls get_merged_cob_terrane_raster)."""
    from gprm.utils.spatial import get_merged_cob_terrane_raster

    with pytest.raises(ValueError, match="pygplates"):
        get_merged_cob_terrane_raster('unused.gpml', None, 0, sampling=1, method='raterio')


def test_unknown_sampling_method_names_the_valid_options():
    """Previously UnboundLocalError: point_raster_values, raised well after the expensive
    reconstruction had already run."""
    pytest.importorskip('pygmt')
    from gprm.utils.raster import reconstruct_raster

    with pytest.raises(ValueError, match="'scipy', 'gmt', 'stripy'"):
        reconstruct_raster(None, None, None, 0, 10, sampling_method='gtm')


# ---------------------------------------------------------- dataset fetch errors

def test_fetch_failure_names_the_dataset_and_the_url():
    """Previously a dead URL surfaced as a bare requests traceback, with no indication of
    which loader was responsible or where its cache lives."""
    from gprm.datasets import DatasetFetchError
    from gprm.datasets._fetch import retrieve

    def Seamounts():           # stands in for a real loader, to exercise the stack walk
        return retrieve(url='https://gprm.invalid/nope.nc', known_hash=None,
                        path='/tmp/gprm-test-cache')

    with pytest.raises(DatasetFetchError) as excinfo:
        Seamounts()

    message = str(excinfo.value)
    assert 'Seamounts' in message
    assert 'https://gprm.invalid/nope.nc' in message
    assert '/tmp/gprm-test-cache' in message
    assert excinfo.value.__cause__ is not None


def test_dataset_fetch_error_is_a_runtime_error():
    """Subclassing RuntimeError keeps existing broad excepts working."""
    from gprm.datasets import DatasetFetchError

    assert issubclass(DatasetFetchError, RuntimeError)


def test_http_status_is_turned_into_advice():
    from gprm.datasets._fetch import _diagnose

    class FakeResponse:
        status_code = 404

    class FakeError(Exception):
        response = FakeResponse()

    assert 'no longer exists' in _diagnose(FakeError(), 'https://example.invalid/x')
    assert _diagnose(ValueError('hash of downloaded file is different'), 'x') is not None


# --------------------------------------------------------------- dead code / portability

def test_proximity_does_not_import_xrspatial_or_datashader():
    """generate_shadows was the last importer of xrspatial, and had no return statement."""
    import gprm.utils.proximity as proximity

    assert not hasattr(proximity, 'generate_shadows')

    # Parsed rather than grepped: the module docstring and the commented-out previous
    # implementations both still mention xrspatial by name, deliberately.
    import ast

    tree = ast.parse(open(proximity.__file__).read())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split('.')[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split('.')[0])

    assert 'xrspatial' not in imported
    assert 'datashader' not in imported


def test_to_gplates_does_not_hardcode_a_version_or_reject_linux():
    """Previously raised NotImplementedError on Linux, and pinned GPlates 2.3.0 elsewhere."""
    import inspect
    from gprm.GPlatesReconstructionModel import ReconstructionModel

    source = inspect.getsource(ReconstructionModel.to_GPlates)

    assert '2.3.0' not in source
    assert 'NotImplementedError' not in source


def test_to_gplates_rejects_a_path_that_does_not_exist(reconstruction_model):
    from gprm.GPlatesReconstructionModel import ReconstructionModel

    with pytest.raises(FileNotFoundError):
        reconstruction_model.to_GPlates(path_to_gplates='/no/such/gplates')
