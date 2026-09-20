"""Molchan test and space-time distance analysis for alarm-based association testing.

What this computes
------------------
Given an *alarm function* over the globe -- here always distance to some target, so that
"close to a target" means "high alarm" -- and a set of events, the Molchan diagram plots

    tau = the fraction of the permissible AREA at least as strongly alarmed, against
    mu  = the fraction of events MISSED at that alarm level.

A useless alarm gives mu = 1 - tau, the diagonal. The skill score is the area between the
curve and that diagonal::

    Skill = 0.5 - integral(mu d_tau)

which is 0 for an alarm no better than chance and 0.5 for one that captures every event in
a vanishing area. Equivalently, and rather more usefully::

    Skill = 0.5 - mean(tau evaluated at each event)

This is the same statistic as the ROC AUC of the alarm function, shifted by a half, with
*area* playing the role of the negative class::

    Skill = AUC - 0.5

so it is also the Gini coefficient halved, and it is the "capture-efficiency" or
"prediction-area" curve of mineral prospectivity under another name. Being rank-based, it
is unchanged by any monotone transform of distance.

What it does not tell you
-------------------------
The null hypothesis is that events are placed uniformly at random over the permissible
area. For geological sample data that null is usually false for reasons having nothing to
do with the target: samples are collected where there is outcrop, access and funding, and
that bias is itself spatially structured. A high skill score against a uniform null is
therefore evidence of association, not of a causal or tectonic relationship. Restricting
the permissible region (see ``generate_masked_distance_raster_sequence``) controls the
crudest form of this, and nothing here controls the rest.

The score also carries no confidence interval. The usual one assumes independent events,
which zircon grains from a single sample, or deposits from a single district, are not.

Two further choices are made here rather than by the caller, and neither is announced by
the result:

* The alarm is *distance to the nearest* target, so a sample 500 km from three separate
  trenches scores exactly as one 500 km from a single trench. If the hypothesis concerns
  the amount of subduction rather than the presence of it, a density-weighted alarm would
  be the right one and this is not it.
* Reconstructed sample positions are treated as exact. Palaeo-position uncertainty at
  100+ Ma is of the same order as the distances being tested, and it grows with age, so
  the scores for young and old samples are not equally well determined. The effect is to
  blur association, biasing skill towards zero, by an amount that varies across the
  dataset. Quantifying it would take error propagation through the rotation model, which
  nothing here does.

MIT License

Copyright (c) 2017-2021 Simon Williams
"""
import warnings

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from .proximity import contour_proximity, polyline_proximity, polygons_buffer, points_proximity, boundary_proximity, reconstruct_and_rasterize_polygons
from .create_gpml import gpml2gdf
import xarray as xr
import shapely
import pygplates
from .geometry import apply_nearest_feature, apply_reconstruction, wrap_polygon_feature, wrap_polyline_feature
from .spatial import topology_lookup


from collections import OrderedDict


DEFAULT_DISTANCE_MAX = 1e7
DEFAULT_DISTANCE_STEP = 2e4
DEFAULT_GEOGRAPHIC_EXTENT = [-180.,180.,-90.,90.]
DEFAULT_GEOGRAPHIC_SAMPLING = 0.25

# IUGG mean radius, matching pygplates.Earth.mean_radius_in_kms and utils.proximity
EARTH_RADIUS_KM = 6371.0088

# GPlates uses FROMAGE and TOAGE for two unrelated things: the valid time of a reconstructable
# feature, and -- in gprm's age-coded convention -- the age of the sample itself. Where both are
# in play the plate's valid time gets its own names, so that neither silently overwrites the
# other. See sample_distance_analysis.
PLATE_ID_FIELD = 'PLATEID1'
SAMPLE_APPEARANCE_FIELD = 'FROMAGE'
SAMPLE_DISAPPEARANCE_FIELD = 'TOAGE'
PLATE_APPEARANCE_FIELD = 'plate_FROMAGE'
PLATE_DISAPPEARANCE_FIELD = 'plate_TOAGE'


def scipy_interpolater(da, points):
    """Bilinearly sample a regular lon/lat grid at scattered points.

    Replaces ``pygmt.grdtrack(..., interpolation='l')``, which needed a GMT installation.
    Longitudes are wrapped into the grid's range and latitudes clipped to it, so that a
    point at 190 degrees east is sampled at -170 rather than falling off the edge; points
    that genuinely fall on no-data cells come back NaN, as grdtrack's ``no_skip`` does.

    :param da: xarray DataArray with dimensions (y, x) in degrees.
    :param points: array-like of shape (n, 2), first column longitude, second latitude.
        A DataFrame is accepted and its column *order* is what matters, not the names.
    :returns: 1-D numpy array of sampled values, NaN where the grid has no data.
    """
    coordinates = np.asarray(points, dtype=float)
    lons = coordinates[:, 0]
    lats = coordinates[:, 1]

    x = np.asarray(da['x'], dtype=float)
    y = np.asarray(da['y'], dtype=float)
    values = np.asarray(da.data, dtype=float)

    # RegularGridInterpolator requires ascending coordinates
    if y.size > 1 and y[0] > y[-1]:
        y = y[::-1]
        values = values[::-1, :]
    if x.size > 1 and x[0] > x[-1]:
        x = x[::-1]
        values = values[:, ::-1]

    if x.size > 1 and np.isclose(abs(x[-1] - x[0]), 360.0):
        lons = ((lons - x[0]) % 360.0) + x[0]
    lons = np.clip(lons, x.min(), x.max())
    lats = np.clip(lats, y.min(), y.max())

    interpolate = RegularGridInterpolator((x, y), values.T, method='linear',
                                          bounds_error=False, fill_value=np.nan)
    return interpolate(np.column_stack([lons, lats]))


def _cell_area_weights(da):
    """Area of every cell of a regular lon/lat grid, in square kilometres.

    A cell spanning dlon in longitude and [lat0, lat1] in latitude covers
    ``R^2 * dlon * (sin(lat1) - sin(lat0))`` exactly, so no approximation is involved
    beyond the grid itself. Counting cells instead -- treating the grid as flat -- is wrong
    by up to 10% of the sphere for a target whose alarm region is concentrated in latitude,
    which trenches and ridges are.
    """
    y = np.asarray(da['y'], dtype=float)
    x = np.asarray(da['x'], dtype=float)
    if y.size < 2 or x.size < 2:
        raise ValueError('a grid of at least 2x2 cells is needed to compute cell areas')

    dlon = np.deg2rad(abs(x[1] - x[0]))
    dlat = np.deg2rad(abs(y[1] - y[0]))

    # clipped, so the half-cells at the poles are not counted beyond +/-90
    lower = np.clip(np.deg2rad(y) - dlat / 2.0, -np.pi / 2, np.pi / 2)
    upper = np.clip(np.deg2rad(y) + dlat / 2.0, -np.pi / 2, np.pi / 2)
    row = (EARTH_RADIUS_KM ** 2) * dlon * (np.sin(upper) - np.sin(lower))

    weights = np.repeat(row[:, None], x.size, axis=1)

    # a gridline-registered global grid holds the same meridian at both -180 and +180
    if np.isclose(abs(x[-1] - x[0]), 360.0):
        weights[:, 0] *= 0.5
        weights[:, -1] *= 0.5

    return weights


def _sorted_area_profile(da):
    """Grid values sorted ascending, with the cumulative area at or below each one."""
    values = np.asarray(da.data, dtype=float)
    weights = _cell_area_weights(da)

    finite = np.isfinite(values)
    values = values[finite]
    weights = weights[finite]

    order = np.argsort(values, kind='stable')
    return values[order], np.cumsum(weights[order])


def _area_above(da, contours):
    """Area, in square kilometres, where the grid exceeds each contour.

    Replaces ``pygmt.grdvolume(..., f='g', unit='k')``, whose second column this
    reproduces. Agrees with GMT to better than 0.1% of the sphere.
    """
    values, cumulative = _sorted_area_profile(da)
    contours = np.asarray(contours, dtype=float)
    if values.size == 0:
        return np.zeros(contours.shape)

    total = cumulative[-1]
    index = np.searchsorted(values, contours, side='right')
    at_or_below = np.where(index > 0, cumulative[np.clip(index - 1, 0, None)], 0.0)
    return total - at_or_below


def _alarm_fraction(da, distances):
    """tau: the fraction of permissible area at least as strongly alarmed as each distance.

    A distance of NaN means the event fell somewhere the alarm function does not cover --
    outside a continental mask, say. Such an event is missed at every alarm level, so its
    tau is 1. That is the behaviour the module has always documented; it is only now the
    behaviour it has (see the note in ``molchan_test``).
    """
    values, cumulative = _sorted_area_profile(da)
    distances = np.asarray(distances, dtype=float)
    if values.size == 0:
        return np.full(distances.shape, np.nan)

    index = np.searchsorted(values, distances, side='right')
    at_or_below = np.where(index > 0, cumulative[np.clip(index - 1, 0, None)], 0.0)
    tau = at_or_below / cumulative[-1]
    return np.where(np.isfinite(distances), tau, 1.0)


def _miss_rate(distances, contours):
    """mu: the fraction of events not captured at each alarm level."""
    distances = np.asarray(distances, dtype=float)
    if distances.size == 0:
        return np.ones(np.shape(contours))

    scoreable = np.sort(distances[np.isfinite(distances)])
    captured = np.searchsorted(scoreable, np.asarray(contours, dtype=float), side='left')
    # divided by every event, not only the scoreable ones, so that events the alarm
    # cannot reach count against it
    return 1.0 - captured / distances.size


def _warn_about_failed_times(sequence, what):
    """Say up front which times produced no raster, rather than leaving None in the sequence.

    The workers warn individually, but those warnings come from separate processes and are
    easily missed; a sequence with holes in it otherwise looks complete until something much
    later skips the gaps without comment.
    """
    failed = [time for time, raster in sequence.items() if raster is None]
    if failed:
        shown = ', '.join(str(t) for t in failed[:8])
        warnings.warn(
            'Could not {} at {} of {} times ({}{}). Those entries are None, and any analysis '
            'over this sequence will silently skip them.'.format(
                what, len(failed), len(sequence), shown,
                ', ...' if len(failed) > 8 else ''),
            stacklevel=2)


def _matching_raster_time(raster_dict, age, age_field_name='age'):
    """Find the key of raster_dict matching a sample's age.

    Exact matches are used as they are; a key differing only by floating-point noise is
    accepted. Anything else raises, rather than snapping silently, because an age that
    misses the sequence usually means the ages were never rounded to the time step or are
    in the wrong units, and quietly moving the sample would hide that.
    """
    if age in raster_dict:
        return age

    times = np.asarray(sorted(raster_dict.keys()), dtype=float)
    if times.size == 0:
        raise ValueError('raster_dict is empty, so there is nothing to sample.')

    nearest = times[np.abs(times - float(age)).argmin()]
    if np.isclose(nearest, float(age), rtol=0.0, atol=1e-6):
        return next(k for k in raster_dict if np.isclose(float(k), nearest,
                                                         rtol=0.0, atol=1e-6))

    raise ValueError(
        "No raster in the sequence for {} = {}. The nearest is {}, and the sequence runs "
        "{} to {}. Round the ages onto the sequence's time steps before calling this, the "
        "way sample_distance_analysis does.".format(
            age_field_name, age, nearest, times.min(), times.max()))


def _sample_distances(grid, points, interpolater, buffer_radius=1):
    """Distance of each point from the target, read off the alarm grid."""
    if interpolater not in ('scipy', 'pygmt'):
        raise ValueError(
            "Unknown interpolater {!r}. Choose one of: 'scipy', 'pygmt'.".format(interpolater))

    frame = pd.DataFrame(data=points[['Longitude', 'Latitude']])

    if interpolater == 'pygmt':
        from ._optional import require
        pygmt = require('pygmt', 'sampling grids with grdtrack')
        return pygmt.grdtrack(points=frame, grid=grid, no_skip=False, interpolation='l',
                              radius=buffer_radius, newcolname='dist')['dist'].to_numpy()

    return scipy_interpolater(grid, frame)


def molchan_test(grid,
                 points,
                 distance_max = DEFAULT_DISTANCE_MAX,
                 distance_step = DEFAULT_DISTANCE_STEP,
                 buffer_radius=1,
                 interpolater='scipy',
                 verbose=False):
    """Molchan test for a set of points against one alarm grid.

    :param grid: xarray DataArray of distance to the target, over the permissible region.
        Cells that are NaN are outside the permissible region and are excluded from the
        area; events falling on them count as missed.
    :param points: DataFrame with 'Longitude' and 'Latitude' columns. A 'distance' column
        is added in place, as before.
    :param distance_max: largest distance at which the returned curves are sampled.
    :param distance_step: spacing at which the returned curves are sampled.
    :param buffer_radius: passed to grdtrack; ignored by the default interpolater.
    :param interpolater: 'scipy' (default) or 'pygmt'. These agree to under a metre; the
        default changed because it needs no GMT installation.
    :param verbose: report the permissible area as a fraction of the Earth's surface.
    :returns: (alarm fraction tau, miss rate mu, skill score). The two curves are returned
        in order of *decreasing* distance, as before.

    The skill score no longer depends on ``distance_max`` or ``distance_step``. Those set
    how finely the returned curves are sampled and nothing else; the score is computed
    exactly, as ``0.5 - mean(tau)`` over the events. Two consequences worth knowing if you
    are comparing against numbers from an earlier version:

    * Distances beyond ``distance_max`` used to truncate both curves, and the truncated
      trapezoidal integral reported skill where there was none. With two targets, where
      45% of the globe lies beyond the default 10,000 km, uniformly random events scored
      +0.36 instead of 0.
    * Events the alarm could not score used to move the result *towards* +0.5, the maximum,
      although the code asserted the opposite. They now count as missed, as documented.
    """
    distances = _sample_distances(grid, points, interpolater, buffer_radius)
    points['distance'] = distances

    if verbose:
        total_area = float(_area_above(grid, [-np.inf])[0])
        earth_area = 4 * np.pi * EARTH_RADIUS_KM ** 2
        print('Total permissible area is {:0.1f}% of total Earth surface'.format(
            100 * total_area / earth_area))

    contours = np.arange(0.0, distance_max + distance_step / 2.0, distance_step)

    # tau, the alarm fraction, rising 0 -> 1 with distance
    grid_fraction = pd.Series(1.0 - _area_above(grid, contours) / _area_above(grid, [-np.inf])[0])
    # mu, the miss rate, falling 1 -> 0 with distance
    points_fraction = _miss_rate(distances, contours)

    # Exact, rather than integrated over the sampled curves: for an alarm ranked by
    # distance, 0.5 - integral(mu d_tau) is identically 0.5 - mean(tau at each event).
    Skill = float(0.5 - np.mean(_alarm_fraction(grid, distances)))

    return grid_fraction[::-1], points_fraction[::-1], Skill



def molchan_point(grid, 
                  points, 
                  distance_max = DEFAULT_DISTANCE_MAX, 
                  distance_step = DEFAULT_DISTANCE_STEP, 
                  buffer_radius=1,
                  interpolater='scipy',
                  verbose=False,
                  return_fraction=True):
    """Molchan test for a single point.

    :returns: (distance to the target, tau) where tau is the fraction of the permissible
        area at least as close to a target as this point is -- so a small tau means the
        point sits in a small, strongly alarmed region. ``molchan_test``'s skill score is
        ``0.5 - mean(tau)`` over the events, so the two are directly comparable.

        With ``return_fraction=False``, only the distance.

    ``distance_max`` and ``distance_step`` are accepted for call compatibility and no
    longer affect the result: tau is evaluated exactly rather than interpolated from a
    contour table, so a point beyond ``distance_max`` is no longer clamped to its value
    there.
    """
    distances = _sample_distances(grid, points, interpolater, buffer_radius)
    points['distance'] = distances

    if not return_fraction:
        return float(distances[0])

    if verbose:
        total_area = float(_area_above(grid, [-np.inf])[0])
        earth_area = 4 * np.pi * EARTH_RADIUS_KM ** 2
        print('Total permissible area is {:0.1f}% of total Earth surface'.format(
            100 * total_area / earth_area))

    return float(distances[0]), float(_alarm_fraction(grid, distances)[0])
    
    

def space_time_molchan_test(raster_dict, 
                            point_distances,
                            healpix_resolution=128,
                            distance_max=DEFAULT_DISTANCE_MAX,
                            distance_step=DEFAULT_DISTANCE_STEP,
                            interpolater='scipy'):
    """Molchan test over a sequence of rasters, one per reconstruction time.

    The alarm fraction is taken over the whole space-time volume: an equal-area healpix
    distribution is sampled from every raster in the sequence and the results pooled, so
    each time step contributes equally. That makes the result a statement about
    space-time association. If the events cluster in time and the alarmed area also varies
    with time, a purely temporal coincidence will register as skill here; the statistic
    cannot separate the two.

    :param raster_dict: {reconstruction_time: DataArray of distance to target}.
    :param point_distances: distances already extracted for the events, e.g. from
        ``space_time_distances``. NaN entries count as missed.
    :param healpix_resolution: healpix N for the background sample. The points are equal
        area by construction, so no latitude weighting is needed here.
    :returns: (alarm fraction tau, miss rate mu, skill score), in order of increasing
        distance.

    As in ``molchan_test``, the skill score is computed exactly and does not depend on
    ``distance_max`` or ``distance_step``.
    """
    from gprm import PointDistributionOnSphere
    hp = PointDistributionOnSphere(distribution_type='healpix', N=healpix_resolution)
    hp_dataframe = pd.DataFrame(data={'x':hp.longitude, 'y':hp.latitude})

    space_time_distances = []

    for reconstruction_time in raster_dict.keys():
        raster = raster_dict[reconstruction_time]
        if raster is None:
            continue
        if interpolater == 'scipy':
            smpl = scipy_interpolater(raster, hp_dataframe)
        else:
            from ._optional import require
            pygmt = require('pygmt', 'sampling grids with grdtrack')
            smpl = pygmt.grdtrack(grid=raster, points=hp_dataframe, no_skip=False,
                                  interpolation='l', newcolname='distance')['distance'].to_numpy()
        space_time_distances.extend(np.asarray(smpl)[np.isfinite(smpl)].tolist())

    if not space_time_distances:
        raise ValueError(
            'No raster in the sequence yielded any valid samples, so the alarm fraction '
            'cannot be computed. Check that raster_dict is not empty and that its entries '
            'are not all None or all-NaN.')

    background = np.sort(np.asarray(space_time_distances, dtype=float))
    point_distances = np.asarray(point_distances, dtype=float)

    contours = np.arange(0.0, distance_max + distance_step / 2.0, distance_step)

    # tau, from the equal-area background sample
    grid_fraction = np.searchsorted(background, contours, side='right') / background.size
    # mu, over every event including those the alarm could not score
    point_fraction = _miss_rate(point_distances, contours)

    # tau evaluated at each event, exactly as in molchan_test
    tau_at_events = np.searchsorted(background, point_distances, side='right') / background.size
    tau_at_events = np.where(np.isfinite(point_distances), tau_at_events, 1.0)
    Skill = float(0.5 - np.mean(tau_at_events))

    return grid_fraction, point_fraction, Skill



def combine_raster_sequences(raster_dict1, raster_dict2):
    """
    Given two raster sequences (dictionaries with coincident keys), 
    generate a new raster sequence that multiplies the coincident rasters
    from each sequence
    """
    
    raster_dict3 = OrderedDict()

    for key in raster_dict1.keys():
        raster_dict3[key] = raster_dict1[key] * raster_dict2[key]
        
    return raster_dict3

    

def space_time_distances(raster_dict, gdf, age_field_name='age', 
                         distance_max=DEFAULT_DISTANCE_MAX, 
                         distance_step=DEFAULT_DISTANCE_STEP, 
                         buffer_radius=1,
                         interpolater='scipy'):
    """
    Computes the distances to targets rconstructed to their time of appearance
    from a raster sequence of raster grids
    
    The input gdf is assumed to have reconstructed coordinates in its geometry

    Each sample is looked up in ``raster_dict`` by its age, so the ages must already sit on
    the sequence's time steps. Round them yourself, as ``sample_distance_analysis`` does,
    rather than having it happen silently here; an age matching no raster raises an error
    naming the nearest one available, instead of a bare KeyError.
    """

    results = []

    for i,row in gdf.iterrows():
        reconstruction_time = _matching_raster_time(raster_dict, row[age_field_name],
                                                    age_field_name)
        result = molchan_point(raster_dict[reconstruction_time],
                               pd.DataFrame(data={'Longitude': [row.geometry.x], 
                                                  'Latitude': [row.geometry.y]}),
                               distance_max=distance_max, 
                               distance_step=distance_step, 
                               buffer_radius=buffer_radius, 
                               interpolater=interpolater,
                               )
        results.append(result)

    return pd.DataFrame(data=results, 
                        columns=['distance', 'area_fraction'])

from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from collections import OrderedDict
import numpy as np


def _process_polygon_rasterization(reconstruction_time, features, rotation_model,
                                    sampling, anchor_plate_id, buffer_distance):
    try:
        tmp = reconstruct_and_rasterize_polygons(features, rotation_model,
                                                  reconstruction_time,
                                                  sampling=sampling,
                                                  anchor_plate_id=anchor_plate_id)
        tmp = tmp.where(tmp != 0, np.nan)
        if buffer_distance is not None:
            bn = boundary_proximity(tmp)
            tmp.data[bn.data <= buffer_distance] = 1
        return reconstruction_time, tmp
    except Exception as e:
        # Printed from a worker process, where stdout is easily lost, and the None then sat
        # in the raster sequence to be skipped silently much later. Warned instead, and the
        # caller is told how many times failed.
        warnings.warn('Could not rasterize polygons at {} Ma: {}: {}'.format(
            reconstruction_time, type(e).__name__, e))
        return reconstruction_time, None


def generate_raster_sequence_from_polygons(features,
                                           rotation_model,
                                           reconstruction_times,
                                           sampling=DEFAULT_GEOGRAPHIC_SAMPLING,
                                           buffer_distance=None,
                                           max_workers=None,
                                           anchor_plate_id=0):
    """
    Given some reconstructable polygon features, generates a series of rasterized outputs
    using multiprocessing with progress bar.

    Parameters:
    -----------
    features : object
        Polygon features for reconstruction
    rotation_model : object
        Rotation model for reconstruction
    reconstruction_times : list
        List of reconstruction times
    sampling : float
        Geographic sampling parameter
    buffer_distance : float, optional
        Buffer distance for boundary proximity (default: None)
    max_workers : int, optional
        Maximum number of worker processes (default: None uses ProcessPoolExecutor default)
    anchor_plate_id : int, optional
        Anchor plate ID for reconstruction (default: 0)

    Returns:
    --------
    OrderedDict : Ordered dictionary of rasterized data keyed by reconstruction time
    """
    results = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_time = {
            executor.submit(_process_polygon_rasterization, time, features, rotation_model,
                            sampling, anchor_plate_id, buffer_distance): time
            for time in reconstruction_times
        }
        with tqdm(total=len(reconstruction_times), desc="Processing polygon rasterization") as pbar:
            for future in as_completed(future_to_time):
                reconstruction_time, raster_data = future.result()
                results[reconstruction_time] = raster_data
                pbar.update(1)

    raster_dict = OrderedDict()
    for reconstruction_time in reconstruction_times:
        raster_dict[reconstruction_time] = results[reconstruction_time]

    _warn_about_failed_times(raster_dict, 'rasterize polygons')

    return raster_dict


def _process_distance_raster(reconstruction_time, target_features, reconstruction_model,
                              sampling, region):
    try:
        if isinstance(target_features, dict):
            if not target_features[reconstruction_time]:
                return reconstruction_time, zeros_grid_like(sampling=sampling, region=region) * np.nan
            elif isinstance(target_features[reconstruction_time][0], pygplates.Feature):
                r_target_features = gpml2gdf(target_features[reconstruction_time])
            elif isinstance(target_features[reconstruction_time], pd.GeoDataFrame):
                r_target_features = target_features[reconstruction_time]
        else:
            r_target_features = reconstruction_model.reconstruct(target_features,
                                                                 reconstruction_time,
                                                                 use_tempfile=False)

        if r_target_features is not None:
            if isinstance(r_target_features.geometry.iloc[0], shapely.geometry.point.Point):
                prox_grid = points_proximity(r_target_features.geometry.x,
                                             r_target_features.geometry.y,
                                             spacing=sampling,
                                             region=region)
            elif isinstance(r_target_features.geometry.iloc[0], shapely.geometry.linestring.LineString):
                date_line_wrapper = pygplates.DateLineWrapper(0.0)
                r_target_features['geometry'] = r_target_features.apply(
                    lambda x: wrap_polyline_feature(x, date_line_wrapper), axis=1)
                prox_grid = polyline_proximity(r_target_features, spacing=sampling, region=region)
            elif isinstance(r_target_features.geometry.iloc[0], shapely.geometry.polygon.Polygon):
                date_line_wrapper = pygplates.DateLineWrapper(0.0)
                r_target_features['geometry'] = r_target_features.apply(
                    lambda x: wrap_polygon_feature(x, date_line_wrapper), axis=1)
                prox_grid = polygons_buffer(r_target_features, sampling=sampling, region=region)
            else:
                raise ValueError("Unsupported geometry type in target features.")
        else:
            prox_grid = zeros_grid_like(sampling=sampling, region=region) * np.nan

        return reconstruction_time, prox_grid

    except Exception as e:
        warnings.warn('Could not build a distance raster at {} Ma: {}: {}'.format(
            reconstruction_time, type(e).__name__, e))
        return reconstruction_time, None


def generate_distance_raster_sequence(target_features,
                                      reconstruction_model,
                                      reconstruction_times,
                                      sampling=DEFAULT_GEOGRAPHIC_SAMPLING,
                                      region=DEFAULT_GEOGRAPHIC_EXTENT,
                                      max_workers=None):
    """
    Generate distance raster sequence using multiprocessing with progress bar.

    Parameters:
    -----------
    target_features : dict or other
        Target features for reconstruction
    reconstruction_model : object
        Model for reconstruction
    reconstruction_times : list
        List of reconstruction times
    sampling : float
        Geographic sampling parameter
    region : object
        Geographic region
    max_workers : int, optional
        Maximum number of worker processes (default: None uses ProcessPoolExecutor default)

    Returns:
    --------
    OrderedDict : Ordered dictionary of proximity grids keyed by reconstruction time
    """
    results = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_time = {
            executor.submit(_process_distance_raster, time, target_features,
                            reconstruction_model, sampling, region): time
            for time in reconstruction_times
        }
        with tqdm(total=len(reconstruction_times), desc="Processing reconstruction times") as pbar:
            for future in as_completed(future_to_time):
                reconstruction_time, prox_grid = future.result()
                results[reconstruction_time] = prox_grid
                pbar.update(1)

    prox_grid_sequence = OrderedDict()
    for reconstruction_time in reconstruction_times:
        prox_grid_sequence[reconstruction_time] = results[reconstruction_time]

    _warn_about_failed_times(prox_grid_sequence, 'build a distance raster')

    return prox_grid_sequence


def generate_masked_distance_raster_sequence(
        reconstruction_model, boundary_lookup,
        reconstruction_times, sampling=DEFAULT_GEOGRAPHIC_SAMPLING,
        polygon_buffer_distance=None):
    """Distance-to-target rasters restricted to the continents at each reconstruction time.

    Restricting the permissible region is what keeps the Molchan null honest: without it,
    the alarm is compared against the whole sphere, most of which is ocean where a
    continental sample could never have been collected, and almost any continental dataset
    scores well for that reason alone.

    Two things to be aware of when reporting a result built on this:

    * **The domain moves.** Continents grow, drift and change area through time, so the
      permissible region is different at every step and the alarm fraction tau is pooled
      over a time-varying domain. That is the intended behaviour -- tau is a space-time
      area fraction -- but it means tau is not a fraction of any one fixed map, and two
      runs over different time ranges are not directly comparable.
    * **A continental mask is the crudest possible control for collection bias.** It
      removes the ocean and nothing else. Outcrop, accessibility, national survey coverage
      and publication history are all spatially structured and none of them are addressed
      here, so a positive result remains evidence of association rather than of a tectonic
      relationship.

    :param polygon_buffer_distance: if given, the continental mask is dilated by this
        distance in metres, so that samples slightly off the reconstructed polygon edge
        remain scoreable.
    """
    # Make two raster sequences, where:
    # 1. Mask rasters where the pixels lying within continents (from those lying outside, therefore unreconstructable)
    # 2. Distance rasters from subduction zone geometries

    reconstruction_raster_dict = generate_raster_sequence_from_polygons(
        reconstruction_model.continent_polygons[0],
        reconstruction_model.rotation_model,
        reconstruction_times,
        sampling=sampling,
        buffer_distance=polygon_buffer_distance
    )

    target_distance_dict = generate_distance_raster_sequence(
        boundary_lookup,
        reconstruction_model,
        reconstruction_times,
        sampling=sampling
    )
    
    #Combine the distance rasters with the mask rasters
    target_distance_dict_mask = combine_raster_sequences(target_distance_dict, 
                                                         reconstruction_raster_dict)
    
    return target_distance_dict_mask


def generate_random_distance_sequence(target_distance_dict_mask, healpix_N=64):
    
    from gprm import PointDistributionOnSphere
    hp = PointDistributionOnSphere(distribution_type='healpix', N=healpix_N)
    hp_dataframe = pd.DataFrame(data={'x':hp.longitude, 'y':hp.latitude})

    space_time_distances = []
    for reconstruction_time in target_distance_dict_mask.keys():
        raster = target_distance_dict_mask[reconstruction_time]
        if raster is None:
            continue
        smpl = hp_dataframe.copy()
        smpl['distance'] = scipy_interpolater(raster, hp_dataframe)
        smpl['reconstruction_time'] = reconstruction_time
        space_time_distances.append(smpl.dropna(subset=['distance']))

    #space_time_distances
    return pd.concat(space_time_distances)


def zeros_grid_like(sampling=DEFAULT_GEOGRAPHIC_SAMPLING, 
                    region=DEFAULT_GEOGRAPHIC_EXTENT):
    """
    Helper function to get grid of zeros with expected grid shape and sampling
    """
    # This should return a grid with the expected shape based on sampling and region
    # You'll need to replace this with your actual implementation
    # For example:
    # return create_grid_from_region(region, sampling)
    return xr.DataArray(
        np.full((len(np.arange(region[2], region[3]+sampling, sampling)), len(np.arange(region[0], region[1]+sampling, sampling))), 0.),
        coords={'y': np.arange(region[2], region[3]+sampling, sampling), 'x': np.arange(region[0], region[1]+sampling, sampling)},
        dims=['y', 'x']
        )
    #pass



def sample_distance_analysis(data_df, reconstruction_model, 
                             age_field='age', time_min=0, time_max=1000.,
                             reconstruction_time_step=1, targets='subduction',
                             anchor_plate_id=0):
    """Nearest distance between age-coded samples and reconstructed target features.

    Distances are measured from the geometries themselves rather than sampled off a distance
    raster: more accurate, and slower on large datasets.

    :param data_df: GeoDataFrame of samples, with an age column and point geometries at
        present-day coordinates. If it already carries a ``PLATEID1`` column those ids are
        used as they stand, after checking they came from this model; if not, the samples are
        partitioned here.
    :param age_field: column holding each sample's own age in Ma (default 'age').
    :param targets: 'subduction', 'midoceanridge', 'other', or a {time: features} dict.
    :param anchor_plate_id: Plate held fixed (default 0). Applied to both the samples and the
        target boundaries, so the two stay in one reference frame. Note that it cannot be
        applied to a ``targets`` dict you build yourself -- resolve that with the same anchor.
    :returns: the analysed subset of data_df, with 'reconstruction_time', 'rgeometry' and
        'distance_to_target' added. Counts of what was excluded and why are on
        ``result.attrs['sample_distance_analysis_counts']``.

    **On plate ids.** This used to call ``assign_plate_ids(copy_valid_times=True)``
    unconditionally, which overwrote any ids the caller had assigned -- including ids from a
    different and deliberately chosen polygon set -- and, because that call copies the
    partitioning polygon's valid time into FROMAGE and TOAGE, silently replaced the sample
    ages held in those columns with the polygon's. The row filter that followed then compared
    each sample's age against the polygon's appearance time while appearing to compare it
    against its own, discarding roughly 44% of a uniformly-aged dataset without saying so.

    Now: supplied ids are validated, not replaced; absent ids are assigned for convenience;
    the polygon's valid time is written to ``plate_FROMAGE``/``plate_TOAGE``, never over the
    sample's own; and every exclusion is counted and warned about.

    The filter that drops samples older than their plate's appearance can only be applied when
    the partitioning happened here, since that is the only time the polygon's valid time is
    known. Supply ``PLATEID1`` yourself and you are taken to have made that judgement already.
    """
    
    # If not provided, create a lookup table for the target features
    if isinstance(targets, dict):
        target_lookup = targets
    elif targets in ['subduction', 'midoceanridge', 'other']:
        target_lookup = topology_lookup(reconstruction_model,
                                        np.arange(time_min, time_max+reconstruction_time_step, reconstruction_time_step),
                                        boundary_types=[targets],
                                        anchor_plate_id=anchor_plate_id)
    elif isinstance(targets, str):
        raise ValueError(
            "Unknown targets {!r}. Choose one of: 'subduction', 'midoceanridge', 'other', or "
            'pass a {{reconstruction_time: features}} dict.'.format(targets))
    else:
        raise TypeError(
            'targets must be one of the boundary-type names or a '
            '{reconstruction_time: features} dict, not {}.'.format(type(targets).__name__))
        
    counts = {'input': len(data_df)}

    data_select = data_df[(data_df[age_field] <= time_max) & (data_df[age_field] >= time_min)]
    counts['within_time_range'] = len(data_select)

    if PLATE_ID_FIELD in data_select.columns:
        # The caller has already partitioned, deliberately and possibly against a polygon set
        # of their own choosing, so respect it. Check only that the ids belong to this model,
        # which is the same guard ReconstructionModel.reconstruct applies.
        reconstruction_model._check_plate_ids(data_select)
        data_select = data_select.copy().reset_index(drop=True)
        counts['partitioned_here'] = 0
    else:
        # Nothing assigned, so do it -- more convenient than making the caller run one line
        # before every analysis.
        #
        # assign_plate_ids(copy_valid_times=True) writes the *partitioning polygon's* valid
        # time into FROMAGE and TOAGE, having first dropped whatever was in those columns.
        # In gprm's age-coded convention those same columns hold the *sample's* own age, so
        # the two meanings collide and the sample ages would be destroyed. Renaming them out
        # of the way first is what keeps them: they then travel with their rows through the
        # overlay, which reindexes and so rules out putting them back afterwards.
        stash = {name: '__sample_' + name
                 for name in (SAMPLE_APPEARANCE_FIELD, SAMPLE_DISAPPEARANCE_FIELD)
                 if name in data_select.columns}
        prepared = data_select.rename(columns=stash)

        assigned = reconstruction_model.assign_plate_ids(prepared,
                                                         keep_unpartitioned_features=False,
                                                         copy_valid_times=True)
        counts['partitioned_here'] = len(assigned)
        counts['outside_all_polygons'] = counts['within_time_range'] - len(assigned)

        assigned = assigned.rename(columns={SAMPLE_APPEARANCE_FIELD: PLATE_APPEARANCE_FIELD,
                                            SAMPLE_DISAPPEARANCE_FIELD: PLATE_DISAPPEARANCE_FIELD})
        assigned = assigned.rename(columns={v: k for k, v in stash.items()})

        # A sample cannot sit on a plate that had not appeared yet. This was the intent of the
        # original filter; it just used the column the sample ages live in.
        in_existence = assigned[age_field] <= assigned[PLATE_APPEARANCE_FIELD]
        counts['older_than_their_plate'] = int((~in_existence).sum())
        data_select = assigned[in_existence].reset_index(drop=True)

    counts['analysed'] = len(data_select)

    counts['outside_time_range'] = counts['input'] - counts['within_time_range']

    reasons = [(counts.get(reason, 0), description) for reason, description in (
        ('outside_time_range', 'outside the {}-{} Ma range'.format(time_min, time_max)),
        ('outside_all_polygons', 'not inside any polygon of the model'),
        ('older_than_their_plate', 'older than the plate they sit on'))]
    excluded = [(n, description) for n, description in reasons if n]

    if excluded:
        warnings.warn(
            '{} of {} samples were not analysed: {}. Full counts are on the result, as '
            "result.attrs['sample_distance_analysis_counts'].".format(
                counts['input'] - counts['analysed'], counts['input'],
                '; '.join('{} {}'.format(n, description) for n, description in excluded)),
            stacklevel=2)

    # assign a reconstruction time which is the nearest time step to the age associated with the data point
    data_select['reconstruction_time'] = np.round(data_select[age_field]/reconstruction_time_step)*reconstruction_time_step

    # reconstruct the points
    # Same anchor as the targets were resolved in, or the samples and the boundaries they are
    # measured against would sit in different reference frames.
    data_select['rgeometry'] = data_select.apply(
        lambda x: apply_reconstruction(x, reconstruction_model.rotation_model,
                                       anchor_plate_id=anchor_plate_id), axis=1)

    # Determine the shortest distance to the target features at the associated time
    data_select['distance_to_target'] = data_select.apply(
        lambda x: apply_nearest_feature(x,
                                        target_lookup,
                                        geometry_field='rgeometry',
                                        age_field='reconstruction_time'), axis=1)

    data_select.attrs['sample_distance_analysis_counts'] = counts

    return data_select

