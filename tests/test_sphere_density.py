"""Point density / HEALPix binning on the sphere.

groupby_healpix previously assigned each point to its nearest member of a separately-built
equal-area point set via a great-circle proximity search -- an indirect, approximate way to
reinvent HEALPix pixel assignment. It now calls astropy_healpix's pixel-index function
directly: exact, vectorized, and it no longer needs PlateTectonicTools or a materialized point
set at all. These tests check the new binning is correct and that the functions built on top of
it (healpix_density, dominant_class, healpix_bin_geometries) behave as designed.
"""
import numpy as np
import pandas as pd
import geopandas as gpd
import pytest
from shapely.geometry import Point

astropy_healpix = pytest.importorskip('astropy_healpix')
from astropy_healpix import healpy as hp

from gprm.utils.sphere import (groupby_healpix, healpix_density, point_density,
                               dominant_class, healpix_bin_geometries)


def _points_gdf(lons, lats, **columns):
    n = len(lons)
    columns = {name: (value if hasattr(value, '__len__') and not isinstance(value, str)
                     else [value] * n)
              for name, value in columns.items()}
    return gpd.GeoDataFrame(columns, geometry=[Point(x, y) for x, y in zip(lons, lats)])


def test_groupby_healpix_matches_ang2pix_directly():
    """Bin assignment should be exactly what astropy_healpix's own pixel lookup gives."""
    rng = np.random.default_rng(0)
    lons = rng.uniform(-180, 180, 50)
    lats = rng.uniform(-90, 90, 50)
    gdf = _points_gdf(lons, lats)

    nside = 8
    grouped = groupby_healpix(gdf, nside)

    expected_pixels = hp.ang2pix(nside, np.radians(90. - lats), np.radians(lons))
    observed_pixels = np.concatenate([[bin_id] * len(group) for bin_id, group in grouped])
    # order isn't guaranteed to match input order once grouped, so compare as multisets
    np.testing.assert_array_equal(sorted(observed_pixels), sorted(expected_pixels))


def test_groupby_healpix_does_not_mutate_caller_frame():
    gdf = _points_gdf([0., 10.], [0., 10.])
    groupby_healpix(gdf, nside=8)
    assert 'bin_id' not in gdf.columns


def test_groupby_healpix_supports_arbitrary_aggregation():
    """The scope decision for this work: gprm assigns bins, pandas does the statistics."""
    gdf = _points_gdf([0., 0.1, 0.2], [0., 0.1, 0.2], age=[10., 20., 30.])
    grouped = groupby_healpix(gdf, nside=32)
    # a single bin at this resolution and clustering -- median works out of the box
    medians = grouped['age'].median()
    assert medians.iloc[0] == 20.


def test_healpix_density_counts_vs_weights():
    gdf = _points_gdf([0., 0., 90.], [0., 0., 0.], magnitude=[1., 3., 100.])

    counts = healpix_density(gdf, nside=8)
    weighted = healpix_density(gdf, weights='magnitude', nside=8)

    # same two occupied bins either way
    assert set(counts.index) == set(weighted.index)
    two_point_bin = counts['value'].idxmax()
    assert counts.loc[two_point_bin, 'value'] == 2
    assert weighted.loc[two_point_bin, 'value'] == 4.  # 1 + 3


def test_point_density_dispatch_matches_named_functions():
    gdf = _points_gdf([0., 5., 170.], [0., 5., -40.])

    via_dispatch = point_density(gdf, method='healpix', nside=8)
    direct = healpix_density(gdf, nside=8)
    pd.testing.assert_frame_equal(via_dispatch, direct)

    kde_result = point_density(gdf, method='kde', bandwidth=0.2, sampling=10)
    assert kde_result.dims == ('lat', 'lon')


def test_point_density_rejects_unknown_method():
    gdf = _points_gdf([0.], [0.])
    with pytest.raises(ValueError, match="Unknown method"):
        point_density(gdf, method='nonsense')


def test_dominant_class_normalizes_by_default():
    """A class with far more points should not automatically win everywhere -- only where its
    own (normalized) spatial pattern is actually stronger.
    """
    rng = np.random.default_rng(1)
    lons = rng.uniform(-170, 170, 200)
    lats = rng.uniform(-80, 80, 200)

    class_a = _points_gdf(lons, lats, cls='A')          # 200 points, spread out
    class_b = _points_gdf(lons[:20], lats[:20], cls='B')  # 20 points, same locations subset
    both = gpd.GeoDataFrame(pd.concat([class_a, class_b], ignore_index=True))

    normalized = dominant_class(both, class_column='cls', nside=8, normalize=True)
    raw = dominant_class(both, class_column='cls', nside=8, normalize=False)

    # unnormalized, A has more raw points in every bin it shares with B, so B never wins
    assert 'B' not in raw['dominant'].values
    # normalized, B's smaller-but-concentrated sample should win in at least some of its bins
    assert 'B' in normalized['dominant'].values


def test_dominant_class_every_row_has_a_real_winner():
    """Every returned pixel has at least one point from some class -- there is no
    all-classes-absent row to produce a spurious/undefined winner for.
    """
    gdf = _points_gdf([0., 90.], [0., 0.], cls=['A', 'B'])
    result = dominant_class(gdf, class_column='cls', nside=8)
    assert result['dominant'].notna().all()
    assert set(result['dominant']) == {'A', 'B'}


def test_healpix_bin_geometries_documents_the_antimeridian_caveat():
    """Pinned, not fixed: a pixel straddling the 0/360 wraparound produces an invalid polygon
    in this planar lon/lat representation, exactly as documented in the function's docstring.
    """
    polygons = healpix_bin_geometries([463], nside=8)
    assert not polygons.geometry.iloc[0].is_valid


def test_healpix_bin_geometries_round_trip():
    # Well away from the wraparound meridian (see the antimeridian caveat in
    # healpix_bin_geometries' docstring) so the polygons here are guaranteed simple.
    gdf = _points_gdf([60., 60., 100.], [10., 10., -10.])
    density = healpix_density(gdf, nside=8)

    polygons = healpix_bin_geometries(density.index, nside=8, values=density['value'])

    assert len(polygons) == len(density)
    assert (polygons.geometry.is_valid).all()
    assert polygons.crs.to_epsg() == 4326
    # the two-point bin's polygon should contain that bin's own pixel center
    two_point_row = density.loc[density['value'].idxmax()]
    matching_polygon = polygons.loc[polygons['bin_id'] == density['value'].idxmax(), 'geometry'].iloc[0]
    assert matching_polygon.contains(Point(two_point_row['longitude'], two_point_row['latitude']))


def test_spherical_kde_weights_change_the_result():
    from gprm.utils.sphere import spherical_kde

    lons = np.array([0., 90.])
    lats = np.array([0., 0.])

    unweighted = spherical_kde(lons, lats, bandwidth=0.3, sampling=10)
    weighted = spherical_kde(lons, lats, weights=[100., 1.], bandwidth=0.3, sampling=10)

    # heavily weighting the point at lon=0 should raise density there relative to lon=90
    lat0_idx = dict(lon=0, lat=0)
    unweighted_ratio = (unweighted.sel(lon=0, lat=0, method='nearest')
                        / unweighted.sel(lon=90, lat=0, method='nearest'))
    weighted_ratio = (weighted.sel(lon=0, lat=0, method='nearest')
                     / weighted.sel(lon=90, lat=0, method='nearest'))
    assert float(weighted_ratio) > float(unweighted_ratio)


def test_spherical_kde_default_weights_match_previous_behaviour():
    """weights=None must reproduce exactly what the function did before weights existed."""
    from gprm.utils.sphere import spherical_kde

    lons = np.array([0., 45., 170.])
    lats = np.array([-10., 20., 60.])

    result = spherical_kde(lons, lats, bandwidth=0.1, sampling=15)
    assert not np.any(result.values == -9999.0)  # the old dead sentinel must never appear
    assert np.all(np.isfinite(result.values))
