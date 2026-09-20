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
def two_plate_subduction_topology(tmp_path):
    """A minimal but genuine topological model: two adjacent closed plate boundaries sharing
    a subduction-zone trench, with plate 200 (east, overriding) rotating about the geographic
    north pole and plate 100 (west, subducting) held fixed. The rotation is pure east-west
    motion at the trench (which runs north-south along the lon=0 meridian), so the sign and
    magnitude of convergence and migration are unambiguous.

    :returns: (rotation_model, topological_features) -- both valid arguments to
        ``pygplates.TopologicalSnapshot`` and to ``ptt.subduction_convergence.subduction_convergence``.
    """
    import pygplates

    # The trench has several vertices so it tessellates into multiple boundary segments,
    # rather than being treated as a single, un-subdivided great circle arc.
    trench_line = pygplates.PolylineOnSphere(
        [(30., 0.), (15., 0.), (0., 0.), (-15., 0.), (-30., 0.)])
    trench_feature = pygplates.Feature(pygplates.FeatureType.gpml_subduction_zone)
    trench_feature.set_geometry(trench_line)
    trench_feature.set_enumeration(pygplates.PropertyName.gpml_subduction_polarity, 'Left')
    trench_feature.set_reconstruction_plate_id(200)

    def edge_feature(p1, p2, plate_id):
        feature = pygplates.Feature()
        feature.set_geometry(pygplates.PolylineOnSphere([p1, p2]))
        feature.set_reconstruction_plate_id(plate_id)
        return feature

    # West plate (100): a simple quadrilateral west of the trench.
    W1, W2, W3, W4 = (-30., -40.), (30., -40.), (30., 0.), (-30., 0.)
    west_edge = edge_feature(W1, W2, 100)
    north_edge_w = edge_feature(W2, W3, 100)
    south_edge_w = edge_feature(W4, W1, 100)

    # East plate (200): a simple quadrilateral east of the trench.
    E1, E2, E3, E4 = (-30., 0.), (30., 0.), (30., 40.), (-30., 40.)
    north_edge_e = edge_feature(E2, E3, 200)
    east_edge = edge_feature(E3, E4, 200)
    south_edge_e = edge_feature(E4, E1, 200)

    def section(feature, reverse):
        # A feature's geometry property name differs by feature type ('gpml:centerLineOf'
        # for a subduction zone, an unclassified name for a plain line); take whichever
        # property the feature actually has, rather than hardcoding one.
        geometry_property_name = next(iter(feature)).get_name()
        return pygplates.GpmlTopologicalLineSection(
            pygplates.GpmlPropertyDelegate(
                feature.get_feature_id(), geometry_property_name, pygplates.GmlLineString),
            reverse)

    west_polygon = pygplates.Feature(pygplates.FeatureType.gpml_topological_closed_plate_boundary)
    west_polygon.set_topological_geometry(pygplates.GpmlTopologicalPolygon([
        section(west_edge, False), section(north_edge_w, False),
        section(trench_feature, False), section(south_edge_w, False)]))
    west_polygon.set_reconstruction_plate_id(100)

    east_polygon = pygplates.Feature(pygplates.FeatureType.gpml_topological_closed_plate_boundary)
    east_polygon.set_topological_geometry(pygplates.GpmlTopologicalPolygon([
        section(trench_feature, True), section(north_edge_e, False),
        section(east_edge, False), section(south_edge_e, False)]))
    east_polygon.set_reconstruction_plate_id(200)

    topological_features = [
        trench_feature, west_edge, north_edge_w, south_edge_w,
        north_edge_e, east_edge, south_edge_e, west_polygon, east_polygon,
    ]

    rotation_path = tmp_path / 'two_plate_subduction.rot'
    rotation_path.write_text(textwrap.dedent("""\
        100   0.0   0.0   0.0   0.0  000 ! west, fixed (subducting)
        100 100.0   0.0   0.0   0.0  000 ! west, fixed (subducting)
        200   0.0   0.0   0.0   0.0  000 ! east, overriding
        200  10.0  90.0   0.0   2.0  000 ! east: rotates about the geographic north pole
        """))
    rotation_model = pygplates.RotationModel(str(rotation_path))

    return rotation_model, topological_features


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
