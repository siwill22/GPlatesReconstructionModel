"""Proximity accuracy.

The proximity functions were previously built on a raster distance transform that did not
wrap the antimeridian. On a global 0.25 degree grid it was wrong by up to 413 km, with the
worst cell sitting exactly at lon -180. These tests check the result against brute-force
great-circle distance, and specifically at the seam, since that is where the old
implementation looked fine everywhere except.
"""
import numpy as np
import pytest
import xarray as xr

from gprm.utils.proximity import (EARTH_RADIUS_M, boundary_proximity, points_proximity,
                                  polyline_proximity)


def brute_force_distance(query_lons, query_lats, source_lons, source_lats):
    """Great-circle distance in metres, computed independently of the code under test."""
    query_lons = np.radians(np.asarray(query_lons, dtype=float))
    query_lats = np.radians(np.asarray(query_lats, dtype=float))
    source_lons = np.radians(np.asarray(source_lons, dtype=float))
    source_lats = np.radians(np.asarray(source_lats, dtype=float))

    # haversine from every query point to every source, then take the minimum
    dlon = query_lons[:, None] - source_lons[None, :]
    dlat = query_lats[:, None] - source_lats[None, :]
    a = (np.sin(dlat / 2) ** 2
         + np.cos(query_lats)[:, None] * np.cos(source_lats)[None, :] * np.sin(dlon / 2) ** 2)
    return (2 * EARTH_RADIUS_M * np.arcsin(np.sqrt(np.clip(a, 0, 1)))).min(axis=1)


@pytest.fixture
def scattered_mask():
    """A 1 degree global binary grid with scattered source cells, plus their coordinates."""
    sampling = 1.0
    lons = np.arange(-180, 180 + sampling, sampling)
    lats = np.arange(-90, 90 + sampling, sampling)

    rng = np.random.default_rng(0)
    rows = rng.integers(0, lats.size, 200)
    cols = rng.integers(0, lons.size, 200)

    mask = np.zeros((lats.size, lons.size))
    mask[rows, cols] = 1

    da = xr.DataArray(mask, coords=[('y', lats), ('x', lons)])
    return da, lons[cols], lats[rows]


def test_boundary_proximity_matches_brute_force(scattered_mask):
    da, source_lons, source_lats = scattered_mask

    result = boundary_proximity(da, inside=False)

    grid_lons, grid_lats = np.meshgrid(da['x'].values, da['y'].values)
    expected = brute_force_distance(grid_lons.ravel(), grid_lats.ravel(),
                                    source_lons, source_lats).reshape(grid_lons.shape)

    np.testing.assert_allclose(result.values, expected, atol=1e-6)


def test_no_error_at_the_antimeridian(scattered_mask):
    """The seam is where the previous implementation was 413 km wrong."""
    da, source_lons, source_lats = scattered_mask

    result = boundary_proximity(da, inside=False)
    grid_lons, grid_lats = np.meshgrid(da['x'].values, da['y'].values)
    expected = brute_force_distance(grid_lons.ravel(), grid_lats.ravel(),
                                    source_lons, source_lats).reshape(grid_lons.shape)

    for column, label in [(0, 'lon -180'), (-1, 'lon +180')]:
        np.testing.assert_allclose(result.values[:, column], expected[:, column], atol=1e-6,
                                   err_msg='error at {}'.format(label))


def test_distance_across_the_seam_is_short():
    """Two points either side of the dateline are close together, not a world apart."""
    result = points_proximity([-179.5], [0.0], spacing=0.5)

    across_the_seam = float(result.sel(x=179.5, y=0.0))
    expected = float(brute_force_distance([179.5], [0.0], [-179.5], [0.0])[0])

    assert across_the_seam == pytest.approx(expected, rel=1e-9)
    assert across_the_seam < 120_000  # one degree of longitude at the equator, not half the globe


def test_boundary_option_is_not_swallowed_by_truthiness(scattered_mask):
    """'both' and 'boundary' are truthy strings, and used to fall into the `if inside:` branch.

    The consequence was that contour_proximity, whose default is inside='boundary', silently
    returned only the inside distance.
    """
    da, _, _ = scattered_mask

    inside = boundary_proximity(da, inside=True)
    boundary = boundary_proximity(da, inside='boundary')
    both = boundary_proximity(da, inside='both')

    assert isinstance(both, tuple) and len(both) == 2
    assert not np.allclose(boundary.values, inside.values)
    np.testing.assert_allclose(boundary.values, both[0].values + both[1].values)


def test_source_cells_are_at_zero_distance(scattered_mask):
    da, _, _ = scattered_mask

    result = boundary_proximity(da, inside=False)

    assert result.values[da.values == 1].max() == 0.0


def test_empty_source_set_gives_infinity():
    da = xr.DataArray(np.zeros((19, 37)),
                      coords=[('y', np.arange(-90, 91, 10.0)), ('x', np.arange(-180, 181, 10.0))])

    result = boundary_proximity(da, inside=False)

    assert np.isinf(result.values).all()


def test_polyline_proximity_is_zero_on_the_line():
    import geopandas as gpd
    from shapely.geometry import LineString

    gdf = gpd.GeoDataFrame(geometry=[LineString([(0, -20), (0, 20)])], crs=4326)

    result = polyline_proximity(gdf, spacing=1.0)

    assert float(result.sel(x=0.0, y=0.0)) == pytest.approx(0.0, abs=1.0)


def test_output_grid_is_gridline_registered():
    result = points_proximity([0.0], [0.0], spacing=1.0)

    assert result.shape == (181, 361)
    assert float(result['x'][0]) == -180.0 and float(result['x'][-1]) == 180.0
    assert float(result['y'][0]) == -90.0 and float(result['y'][-1]) == 90.0
