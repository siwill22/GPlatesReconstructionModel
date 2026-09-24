"""Tests for the Molchan / alarm-based association analysis.

The Molchan diagram has exact analytic anchors, which is what makes it cheap to test:

* an alarm applied to events scattered uniformly at random has skill 0,
* an alarm that puts every event at distance zero has skill 1/2,
* the skill score is identically ``AUC - 0.5`` of the corresponding ROC curve, where the
  "negative" class is *area* rather than a set of non-events.

Those three pin the whole computation, so they are asserted directly rather than against
recorded numbers from a previous run.
"""
import numpy as np
import pandas as pd
import pytest

from gprm.utils.proximity import points_proximity


# ------------------------------------------------------------------ fixtures

def _equal_area_sample(rng, n):
    """Points uniform on the sphere: uniform in longitude and in sin(latitude)."""
    return (rng.uniform(-180.0, 180.0, n),
            np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, n))))


@pytest.fixture
def rng():
    return np.random.default_rng(20240919)


@pytest.fixture
def target_lons_lats(rng):
    return _equal_area_sample(rng, 150)


@pytest.fixture
def alarm_grid(target_lons_lats):
    """Distance to a scattered set of targets: the 'alarm function'."""
    lons, lats = target_lons_lats
    return points_proximity(x=lons, y=lats, spacing=0.5)


# ------------------------------------------------------- the analytic anchors

def test_random_events_have_no_skill(alarm_grid, rng):
    from gprm.utils.molchan import molchan_test

    lons, lats = _equal_area_sample(rng, 5000)
    points = pd.DataFrame({'Longitude': lons, 'Latitude': lats})

    _, _, skill = molchan_test(alarm_grid, points, interpolater='scipy')

    assert skill == pytest.approx(0.0, abs=0.02)


def test_events_on_the_targets_have_maximum_skill(alarm_grid, target_lons_lats):
    from gprm.utils.molchan import molchan_test

    lons, lats = target_lons_lats
    points = pd.DataFrame({'Longitude': lons, 'Latitude': lats})

    _, _, skill = molchan_test(alarm_grid, points, interpolater='scipy')

    assert skill == pytest.approx(0.5, abs=0.01)


def test_skill_equals_auc_minus_a_half(alarm_grid, target_lons_lats, rng):
    """The Molchan skill score is the ROC AUC of the same alarm function, shifted by 1/2,
    with area playing the role of the negative class. Worth asserting because it is the
    single most useful thing to know about the statistic, and it is not obvious from the
    trapezoidal integration the function performs."""
    from gprm.utils.molchan import molchan_test, scipy_interpolater

    sklearn_metrics = pytest.importorskip('sklearn.metrics')

    # a mixture, so the skill is somewhere strictly between 0 and 1/2
    t_lons, t_lats = target_lons_lats
    idx = rng.integers(0, t_lons.size, 600)
    ev_lons = np.concatenate([t_lons[idx] + rng.normal(0, 4, 600), _equal_area_sample(rng, 900)[0]])
    ev_lats = np.concatenate([np.clip(t_lats[idx] + rng.normal(0, 4, 600), -89, 89),
                              _equal_area_sample(rng, 900)[1]])
    ev_lons = ((ev_lons + 180.0) % 360.0) - 180.0
    points = pd.DataFrame({'Longitude': ev_lons, 'Latitude': ev_lats})

    _, _, skill = molchan_test(alarm_grid, points.copy(), interpolater='scipy')

    # the same quantity as a plain AUC, against an equal-area sample of the study region
    bg_lons, bg_lats = _equal_area_sample(rng, 200000)
    background = scipy_interpolater(alarm_grid, pd.DataFrame({'x': bg_lons, 'y': bg_lats}))
    events = scipy_interpolater(alarm_grid, pd.DataFrame({'x': ev_lons, 'y': ev_lats}))
    background = background[np.isfinite(background)]
    events = events[np.isfinite(events)]

    labels = np.r_[np.ones(events.size), np.zeros(background.size)]
    scores = -np.r_[events, background]          # nearer to a target = stronger alarm
    auc = sklearn_metrics.roc_auc_score(labels, scores)

    assert skill == pytest.approx(auc - 0.5, abs=0.01)


# ------------------------------------------------- truncation must not invent skill

def test_sparse_targets_do_not_manufacture_skill(rng):
    """With only two targets, 45% of the globe lies beyond the default distance_max of
    10,000 km -- half the maximum possible great-circle distance of 20,015 km. Both curves
    used to be truncated there, and the truncated trapezoidal integral reported skill
    +0.36 for events that were uniformly random, i.e. a strong false positive in the
    dangerous direction."""
    from gprm.utils.molchan import molchan_test

    grid = points_proximity(x=np.array([0.0, 10.0]), y=np.array([0.0, 5.0]), spacing=0.5)
    assert float((grid > 1e7).mean()) > 0.3         # the condition that used to break it

    lons, lats = _equal_area_sample(rng, 5000)
    points = pd.DataFrame({'Longitude': lons, 'Latitude': lats})

    _, _, skill = molchan_test(grid, points, interpolater='scipy')

    assert skill == pytest.approx(0.0, abs=0.02)


def test_skill_does_not_depend_on_the_binning_parameters(alarm_grid, rng):
    """distance_max and distance_step describe how the returned curves are sampled. They
    are not supposed to change the statistic, and now do not."""
    from gprm.utils.molchan import molchan_test

    lons, lats = _equal_area_sample(rng, 2000)
    points = pd.DataFrame({'Longitude': lons, 'Latitude': lats})

    _, _, coarse = molchan_test(alarm_grid, points.copy(), interpolater='scipy',
                                distance_step=1e5)
    _, _, fine = molchan_test(alarm_grid, points.copy(), interpolater='scipy',
                              distance_step=1e4)
    _, _, short = molchan_test(alarm_grid, points.copy(), interpolater='scipy',
                               distance_max=5e6)

    assert coarse == pytest.approx(fine, abs=1e-9)
    assert short == pytest.approx(fine, abs=1e-9)


# ------------------------------------------------------------- output contract

def test_molchan_test_output_shape_and_direction(alarm_grid, rng):
    from gprm.utils.molchan import molchan_test

    lons, lats = _equal_area_sample(rng, 500)
    points = pd.DataFrame({'Longitude': lons, 'Latitude': lats})

    grid_fraction, points_fraction, skill = molchan_test(alarm_grid, points,
                                                         interpolater='scipy')

    n = len(np.arange(0.0, 1e7 + 2e4 / 2, 2e4))
    assert len(grid_fraction) == n
    assert len(points_fraction) == n
    assert isinstance(skill, float)

    gf = np.asarray(grid_fraction, dtype=float)
    pf = np.asarray(points_fraction, dtype=float)

    # both are returned in order of decreasing distance, so tau falls 1 -> 0 while the
    # miss rate rises 0 -> 1
    assert gf[0] == pytest.approx(1.0, abs=1e-9)
    assert gf[-1] == pytest.approx(0.0, abs=1e-9)
    assert pf[0] == pytest.approx(0.0, abs=1e-9)
    assert pf[-1] == pytest.approx(1.0, abs=1e-9)
    assert np.all(np.diff(gf) <= 1e-12)
    assert np.all(np.diff(pf) >= -1e-12)
    assert np.all((gf >= -1e-12) & (gf <= 1 + 1e-12))
    assert np.all((pf >= -1e-12) & (pf <= 1 + 1e-12))


def test_molchan_point_agrees_with_the_curve(alarm_grid, rng):
    """molchan_point returns (distance, area fraction better than this point). That area
    fraction is exactly tau evaluated at the point's distance, which is what makes
    Skill = 0.5 - mean(tau) hold."""
    from gprm.utils.molchan import molchan_point, molchan_test

    lons, lats = _equal_area_sample(rng, 400)
    points = pd.DataFrame({'Longitude': lons, 'Latitude': lats})

    taus = []
    for lon, lat in zip(lons, lats):
        one = pd.DataFrame({'Longitude': [lon], 'Latitude': [lat]})
        distance, area_fraction = molchan_point(alarm_grid, one, interpolater='scipy')
        assert np.isfinite(distance)
        taus.append(area_fraction)

    _, _, skill = molchan_test(alarm_grid, points, interpolater='scipy')

    assert 0.5 - np.mean(taus) == pytest.approx(skill, abs=0.01)


def test_molchan_point_can_return_distance_only(alarm_grid):
    from gprm.utils.molchan import molchan_point

    one = pd.DataFrame({'Longitude': [0.0], 'Latitude': [0.0]})
    distance = molchan_point(alarm_grid, one, interpolater='scipy', return_fraction=False)

    assert isinstance(distance, float)
    assert distance >= 0.0


# ---------------------------------------------------- runs without GMT installed

def test_molchan_imports_without_pygmt():
    """molchan used to call require('pygmt') at module scope, so the whole module -- the
    statistics as well as the sampling -- was unavailable on a pip-only install. pygmt is
    not pip-installable because it wraps the GMT C library rather than bundling it."""
    import ast

    import gprm.utils.molchan as molchan

    tree = ast.parse(open(molchan.__file__).read())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name.split('.')[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split('.')[0])

    assert 'pygmt' not in imported

    # require() is still allowed, but only inside a function, for interpolater='pygmt'
    module_level = [n for n in tree.body
                    if isinstance(n, ast.Expr) or isinstance(n, ast.Assign)]
    for node in module_level:
        for call in ast.walk(node):
            if isinstance(call, ast.Call) and getattr(call.func, 'id', None) == 'require':
                pytest.fail('require() is called at module scope, so importing molchan '
                            'still fails when the optional dependency is absent')


def test_default_interpolater_needs_no_gmt(alarm_grid, rng):
    """The default path must work too, not just interpolater='scipy'."""
    from gprm.utils.molchan import molchan_test

    lons, lats = _equal_area_sample(rng, 300)
    points = pd.DataFrame({'Longitude': lons, 'Latitude': lats})

    _, _, skill = molchan_test(alarm_grid, points)

    assert np.isfinite(skill)


# ------------------------------------------- the replaced engines, against GMT

SPHERE_KM2 = 4 * np.pi * 6371.0088 ** 2


def test_total_permissible_area_matches_gmt(alarm_grid):
    """The spherical cell-area sum that replaced pygmt.grdvolume."""
    pygmt = pytest.importorskip('pygmt')
    from gprm.utils.molchan import _area_above

    gmt = pygmt.grdvolume(alarm_grid, contour=[0, 1e7, 2e4], f='g', unit='k')

    ours = float(_area_above(alarm_grid, [-np.inf])[0])

    assert abs(ours - gmt.iloc[0, 1]) / SPHERE_KM2 < 1e-4


def test_area_profile_converges_to_gmt_as_the_grid_refines(target_lons_lats):
    """GMT resolves the contour within each cell; summing whole cells does not, so the two
    differ by a quantisation error rather than by a mistake. Asserting first-order
    convergence pins that down far better than any single tolerance would: if the cell
    areas were wrong -- flat, or missing the polar half-cells -- the error would not shrink
    with the grid."""
    pygmt = pytest.importorskip('pygmt')
    from gprm.utils.molchan import _area_above

    lons, lats = target_lons_lats
    errors = []
    for spacing in (1.0, 0.5, 0.25):
        da = points_proximity(x=lons, y=lats, spacing=spacing)
        gmt = pygmt.grdvolume(da, contour=[0, 1e7, 2e4], f='g', unit='k')
        ours = _area_above(da, gmt.iloc[:, 0].to_numpy())
        errors.append(np.abs(ours - gmt.iloc[:, 1].to_numpy()).max() / SPHERE_KM2)

    assert errors[0] < 0.005
    # each halving of the cell size at least a third off the error
    assert errors[1] < 0.67 * errors[0]
    assert errors[2] < 0.67 * errors[1]


def test_sampling_matches_gmt_grdtrack(alarm_grid, rng):
    """The bilinear sampler that replaced pygmt.grdtrack -nl."""
    pygmt = pytest.importorskip('pygmt')
    from gprm.utils.molchan import scipy_interpolater

    lons, lats = _equal_area_sample(rng, 500)
    frame = pd.DataFrame({'x': lons, 'y': lats})

    gmt = pygmt.grdtrack(grid=alarm_grid, points=frame, no_skip=False,
                         interpolation='l', newcolname='d')['d'].to_numpy()
    ours = scipy_interpolater(alarm_grid, frame)

    both = np.isfinite(gmt) & np.isfinite(ours)
    assert both.sum() > 400
    # a metre, on distances of order 10^6 m
    assert np.abs(ours[both] - gmt[both]).max() < 1.0


def test_sampling_handles_points_outside_the_grid(alarm_grid):
    """grdtrack returns NaN off-grid; scipy's RegularGridInterpolator raises by default."""
    from gprm.utils.molchan import scipy_interpolater

    frame = pd.DataFrame({'x': [0.0, 400.0, -400.0], 'y': [0.0, 0.0, 95.0]})
    result = scipy_interpolater(alarm_grid, frame)

    assert np.isfinite(result[0])
    assert np.isfinite(result[1])        # 400 deg wraps to 40 deg
    assert np.isfinite(result[2])        # -400 wraps to -40; lat 95 clips to the pole


# --------------------------------------------------- unreconstructable points

def test_points_off_the_permissible_region_count_as_misses(rng):
    """A point that cannot be scored (NaN in the alarm grid, e.g. outside a continental
    mask) is deliberately treated as missed at every threshold, so that a model is not
    flattered by discarding the data it cannot explain. Asserted because it is a design
    decision that a reader would otherwise have to infer."""
    from gprm.utils.molchan import molchan_test

    lons, lats = _equal_area_sample(rng, 60)
    grid = points_proximity(x=lons, y=lats, spacing=1.0)
    grid = grid.where(grid['y'] > 0)             # mask the southern hemisphere

    on = pd.DataFrame({'Longitude': lons[lats > 20], 'Latitude': lats[lats > 20]})
    _, _, baseline = molchan_test(grid, on.copy(), interpolater='scipy')

    skills = [baseline]
    for n_invalid in (10, 50, 200):
        padded = pd.concat(
            [on, pd.DataFrame({'Longitude': rng.uniform(-180, 180, n_invalid),
                               'Latitude': -rng.uniform(5, 80, n_invalid)})],
            ignore_index=True)
        _, _, skill = molchan_test(grid, padded, interpolater='scipy')
        skills.append(skill)

    # strictly decreasing: the previous implementation moved these *up* towards +0.5,
    # so a model that could not score most of its data scored better for it
    assert np.all(np.diff(skills) < 0), skills


# ------------------------------------------------ raster lookup by sample age

def test_age_not_on_the_time_steps_is_reported_clearly(alarm_grid):
    """Previously a bare KeyError from a raw dict lookup, with nothing to say which age
    was at fault or what the sequence covered."""
    import geopandas as gpd
    from shapely.geometry import Point

    from gprm.utils.molchan import space_time_distances

    raster_dict = {0.0: alarm_grid, 10.0: alarm_grid, 20.0: alarm_grid}
    gdf = gpd.GeoDataFrame({'age': [14.3]}, geometry=[Point(0.0, 0.0)], crs='EPSG:4326')

    with pytest.raises(ValueError, match='14.3'):
        space_time_distances(raster_dict, gdf)


def test_ages_on_the_time_steps_are_accepted(alarm_grid):
    import geopandas as gpd
    from shapely.geometry import Point

    from gprm.utils.molchan import space_time_distances

    raster_dict = {0.0: alarm_grid, 10.0: alarm_grid, 20.0: alarm_grid}
    gdf = gpd.GeoDataFrame({'age': [10.0, 20.0]},
                           geometry=[Point(0.0, 0.0), Point(30.0, 10.0)], crs='EPSG:4326')

    result = space_time_distances(raster_dict, gdf)

    assert list(result.columns) == ['distance', 'area_fraction']
    assert len(result) == 2
    assert np.all(np.isfinite(result['distance']))


def test_floating_point_noise_in_an_age_still_matches(alarm_grid):
    """An age of 10.0000000001, produced by rounding arithmetic upstream, should not be
    treated as a missing time step."""
    import geopandas as gpd
    from shapely.geometry import Point

    from gprm.utils.molchan import space_time_distances

    raster_dict = {0.0: alarm_grid, 10.0: alarm_grid}
    gdf = gpd.GeoDataFrame({'age': [10.0 + 1e-10]}, geometry=[Point(0.0, 0.0)],
                           crs='EPSG:4326')

    assert len(space_time_distances(raster_dict, gdf)) == 1


# ------------------------------------------- sample_distance_analysis: plate ids and ages

@pytest.fixture
def aged_samples():
    """Samples carrying their own ages in gprm's FROMAGE/TOAGE convention."""
    import geopandas as gpd
    from shapely.geometry import Point

    ages = [10.0, 50.0, 90.0, 130.0]
    gdf = gpd.GeoDataFrame(
        {'sample_id': list('abcd'), 'age': ages, 'FROMAGE': ages, 'TOAGE': [0.0] * 4},
        geometry=[Point(x, 0.0) for x in (-30.0, -10.0, 10.0, 30.0)], crs='EPSG:4326')
    return gdf


@pytest.fixture
def target_lookup():
    """A target at every whole Ma, so the analysis has something to measure against."""
    import pygplates

    def one_target():
        feature = pygplates.Feature()
        feature.set_geometry(pygplates.PointOnSphere(0.0, 0.0))
        return [feature]

    return {float(t): one_target() for t in range(0, 1001)}


def test_supplied_plate_ids_are_not_overwritten(reconstruction_model, aged_samples, target_lookup):
    """It used to call assign_plate_ids unconditionally, replacing ids the caller had chosen
    deliberately -- possibly from a different polygon set."""
    from gprm.utils.molchan import sample_distance_analysis

    gdf = aged_samples.copy()
    gdf['PLATEID1'] = 701
    gdf.attrs['gprm_reconstruction_model'] = reconstruction_model.name

    result = sample_distance_analysis(gdf, reconstruction_model, age_field='age',
                                      time_max=200, targets=target_lookup)

    assert (result['PLATEID1'] == 701).all()


def test_sample_ages_survive_the_analysis(reconstruction_model, aged_samples, target_lookup):
    """assign_plate_ids(copy_valid_times=True) copies the partitioning polygon's valid time
    into FROMAGE and TOAGE, which in gprm's age-coded convention hold the sample's own age.
    The sample ages used to be destroyed, and the row filter then compared each sample
    against its polygon's appearance time while appearing to compare it against its own."""
    from gprm.utils.molchan import sample_distance_analysis

    gdf = aged_samples.copy()
    gdf['PLATEID1'] = 701
    gdf.attrs['gprm_reconstruction_model'] = reconstruction_model.name

    result = sample_distance_analysis(gdf, reconstruction_model, age_field='age',
                                      time_max=200, targets=target_lookup)

    assert list(result['FROMAGE']) == list(aged_samples['age'])
    assert list(result['TOAGE']) == [0.0] * len(aged_samples)
    assert list(result['age']) == list(aged_samples['age'])


def test_plate_ids_from_another_model_are_refused_here_too(reconstruction_model, aged_samples, target_lookup):
    from gprm.utils.molchan import sample_distance_analysis

    gdf = aged_samples.copy()
    gdf['PLATEID1'] = 701
    gdf.attrs['gprm_reconstruction_model'] = 'SomeOtherModel'

    with pytest.raises(ValueError, match='SomeOtherModel'):
        sample_distance_analysis(gdf, reconstruction_model, age_field='age',
                                 time_max=200, targets=target_lookup)


def test_exclusions_are_counted_on_the_result(reconstruction_model, aged_samples, target_lookup):
    """Rows used to disappear with only a print() to say so."""
    from gprm.utils.molchan import sample_distance_analysis

    gdf = aged_samples.copy()
    gdf['PLATEID1'] = 701
    gdf.attrs['gprm_reconstruction_model'] = reconstruction_model.name

    with pytest.warns(UserWarning, match='outside the'):
        result = sample_distance_analysis(gdf, reconstruction_model, age_field='age',
                                          time_max=100, targets=target_lookup)

    counts = result.attrs['sample_distance_analysis_counts']
    assert counts['input'] == 4
    assert counts['outside_time_range'] == 1          # only the 130 Ma sample exceeds 100
    assert counts['analysed'] == len(result) == 3


def test_unknown_targets_names_the_valid_options(reconstruction_model, aged_samples):
    from gprm.utils.molchan import sample_distance_analysis

    gdf = aged_samples.copy()
    gdf['PLATEID1'] = 701
    gdf.attrs['gprm_reconstruction_model'] = reconstruction_model.name

    with pytest.raises(ValueError, match='midoceanridge'):
        sample_distance_analysis(gdf, reconstruction_model, targets='subducton')


def test_topology_lookup_accepts_an_anchor_plate(reconstruction_model):
    """Without this the samples could be anchored but the targets could not, putting the two
    in different frames and making every distance between them wrong."""
    import inspect

    from gprm.utils.spatial import topology_lookup

    assert 'anchor_plate_id' in inspect.signature(topology_lookup).parameters


def test_sample_distance_analysis_takes_the_age_from_the_stamp(reconstruction_model, aged_samples,
                                                               target_lookup):
    """With no age_field, the column named by a gprm loader's age description is used -- here
    'Age', which differs from FROMAGE, as in the Carbonatites table."""
    from gprm.datasets._ages import stamp
    from gprm.utils.molchan import sample_distance_analysis

    gdf = aged_samples.rename(columns={'age': 'Age'}).copy()
    gdf['Age_ma'] = gdf['Age']
    gdf['FROMAGE'] = gdf['Age'] + 5.0
    gdf['Error_ma'] = 5.0
    gdf['PLATEID1'] = 701
    gdf.attrs['gprm_reconstruction_model'] = reconstruction_model.name
    stamp(gdf, 'Rocks.Carbonatites')

    result = sample_distance_analysis(gdf, reconstruction_model, time_max=200,
                                      targets=target_lookup)

    assert list(result['reconstruction_time']) == list(aged_samples['age'])


def test_space_time_distances_warns_when_it_has_to_assume_the_age_column():
    import pandas as pd
    import geopandas as gpd
    from gprm.utils.molchan import space_time_distances

    gdf = gpd.GeoDataFrame({'age': [0.0]}, geometry=gpd.points_from_xy([0.], [0.]), crs=4326)
    with pytest.warns(UserWarning, match="'age' is assumed"):
        with pytest.raises(Exception):      # the empty raster dict; only the warning matters here
            space_time_distances({}, gdf)
