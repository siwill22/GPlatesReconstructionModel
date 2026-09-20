"""
Point distributions on sphere, spherical interpolation, and data binning utilities.

MIT License

Copyright (c) 2017-2021 Simon Williams

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import numpy as np
import pygplates
from scipy import spatial
# pygmt is imported per-function rather than at module level: it costs ~1.7 s to
# import, and most of this module does not need it.
import xarray as xr

import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon


def marsaglias_method(N=10000):

    ## Marsaglia's method
    dim = 3

    norm = np.random.normal
    normal_deviates = norm(size=(dim, N))

    radius = np.sqrt((normal_deviates**2).sum(axis=0))
    points = normal_deviates/radius

    return points


def fibonacci_sphere(N=10000):
# https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/44164075
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(N):
        y = 1 - (i / float(N - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points.append((x, y, z))

    return np.vstack(points).T


def golden_spiral(N=10000):
# https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/44164075
    indices = np.arange(0, N, dtype=float) + 0.5

    phi = np.arccos(1 - 2*indices/N)
    theta = np.pi * (1 + 5**0.5) * indices

    x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)

    return np.vstack((x,y,z))


def points_on_sphere(N, distribution_type='marsaglia'):
    '''
    Function to call one of several methods and return Long/
    Lat arrays of points distributed on sphere
    N controls number of points
    distribution_type can be: 'marsaglia' (or 'random')
                              'fibonacci'
                              'spiral'
    '''

    if distribution_type in ['marsaglia','random']:
        points = marsaglias_method(N)
    elif distribution_type=='fibonacci':
        points = fibonacci_sphere(N)
    elif distribution_type=='spiral':
        points = golden_spiral(N)
    else:
        raise ValueError('unrecognised method for point on sphere generation')

    Long=[]
    Lat=[]
    for xyz in points.T:
        LL = pygplates.PointOnSphere((xyz))
        Lat.append(LL.to_lat_lon()[0])
        Long.append(LL.to_lat_lon()[1])

    return np.array(Long), np.array(Lat)


def random_points_feature(N,filename=None):
# function to call Marsaglia's method and return
# feature collection or save to file

    points = marsaglias_method(N)

    #multipoint = pygplates.MultiPointOnSphere((points.T))
    multipoint_feature = pygplates.Feature()
    multipoint_feature.set_geometry(pygplates.MultiPointOnSphere((points.T)))
    multipoint_feature.set_name("Random Points from Marsaglia's method")

    multipoint_feature_collection = pygplates.FeatureCollection(multipoint_feature)

    if filename is not None:
        multipoint_feature_collection.write(filename)
    else:
        return multipoint_feature_collection


def rtp2xyz(r, theta, phi):
    # if only one value, shape will be empty, hence the next if statement
    if r.size==1:
        rdim=1
    else:
        rdim = r.shape[0]
    rst = r * np.sin(theta)
    xout = np.zeros((rdim,3))
    xout[:,0] = rst * np.cos(phi)       # x
    xout[:,1] = rst * np.sin(phi)       # y
    xout[:,2] = r * np.cos(theta)       # z

    return xout


def create_tree_for_spherical_data(inputLons, inputLats, inputVals, n=16):

    ithetas = np.radians(90.-inputLats)
    iphis   = np.radians(inputLons)
    irs     = np.ones(np.shape(ithetas))
    nodes = []

    ixyzs=rtp2xyz(irs.ravel(), ithetas.ravel(), iphis.ravel())
    tree = spatial.cKDTree(ixyzs, n)

    return tree


def sampleOnSphere(inputLons, inputLats, inputVals, sample_points_lons, sample_points_lats, tree=None, n=16, k=1, distance_upper_bound=np.inf):
    """Nearest-neighbour lookup on the sphere. Longitudes and latitudes are in DEGREES.

    :param distance_upper_bound: Maximum separation to consider, in degrees of arc.
    :returns: (d, l) where l are indices into the input arrays and d are the corresponding
        CHORD lengths on the unit sphere. Chord length is monotonic in great-circle distance,
        so the ranking and the indices are exact; convert with 2*arcsin(d/2) if the angle
        itself is wanted.
    """
    # The tree is built on unit-sphere xyz, so query() measures chord length, not arc. A
    # threshold given as an angle therefore has to be converted to the corresponding chord
    # (2*sin(theta/2)); treating the angle as a chord directly overstates the cutoff by ~1.2%
    # at 30 degrees and ~9% at 90.
    if not np.isnan(distance_upper_bound) and np.isfinite(distance_upper_bound):
        distance_upper_bound = 2.0 * np.sin(np.radians(distance_upper_bound) / 2.0)

    if (tree is None):
        tree = create_tree_for_spherical_data(inputLons, inputLats, inputVals, n=n)

    othetas = np.radians(90.-sample_points_lats)
    ophis = np.radians(sample_points_lons)
    oxyzs = rtp2xyz(np.ones(np.shape(othetas)), othetas, ophis)

    d,l = tree.query(oxyzs, k=k, distance_upper_bound=distance_upper_bound)

    return d,l


def healpix_mesh(nSide):
    """
    create a set of healpix points, returned as numpy arrays of the longitudes and latitudes
    """
    #import healpy as hp
    from astropy_healpix import healpy as hp
    othetas,ophis = hp.pix2ang(nSide,np.arange(12*nSide**2))
    othetas = np.pi/2-othetas
    ophis[ophis>np.pi] -= np.pi*2

    # ophis -> longitude, othetas -> latitude
    return np.degrees(ophis), np.degrees(othetas)


## Some functions for data binning
def groupby_healpix(gdf, nside, order='ring'):
    """
    Bin points in a (geo)dataframe into HEALPix pixels and return a groupby object.

    Assigns each row of ``gdf`` (point geometry, lon/lat) directly to its HEALPix pixel at the
    given ``nside`` -- exact and equal-area by construction, no proximity search involved -- and
    returns ``gdf.groupby('bin_id')``. Any pandas reduction works from there: ``grouped.mean()``,
    ``grouped['age'].quantile(0.5)``, or ``grouped.apply(my_stat_fn)`` for something bespoke like
    a combined detrital-zircon age spectrum per bin. gprm's job stops at correct bin assignment;
    generic statistics are what pandas is already the right tool for.

    :param gdf: GeoDataFrame with point geometry.
    :param nside: HEALPix resolution parameter (number of pixels = 12*nside**2).
    :param order: 'ring' or 'nested' HEALPix pixel ordering.
    :returns: gdf.groupby('bin_id') (a copy of gdf is used, the caller's frame is not mutated).
    """
    from ._optional import require
    require('astropy_healpix', 'HEALPix pixel binning')
    from astropy_healpix import healpy as hp

    pixel_ids = hp.ang2pix(nside, np.radians(90. - gdf.geometry.y.to_numpy()),
                          np.radians(gdf.geometry.x.to_numpy()), nest=(order == 'nested'))

    gdf = gdf.copy()
    gdf['bin_id'] = pixel_ids
    return gdf.groupby('bin_id')


def healpix_density(points, weights=None, nside=32, order='ring'):
    """
    Count (or sum weights of) points falling in each HEALPix pixel at the given nside.

    :param points: GeoDataFrame with point geometry.
    :param weights: optional per-point weights -- either an array/Series the same length as
        ``points``, or the name of a column in ``points``. If None, each point counts as 1.
    :param nside: HEALPix resolution parameter (number of pixels = 12*nside**2).
    :param order: 'ring' or 'nested' HEALPix pixel ordering.
    :returns: DataFrame indexed by HEALPix pixel id ('bin_id'), with a 'value' column (counts,
        or the weighted sum) and the pixel-center 'longitude'/'latitude', in degrees. Only
        pixels with at least one point are included.
    """
    if isinstance(weights, str):
        weights = points[weights]

    if weights is None:
        grouped = groupby_healpix(points, nside, order=order)
        result = grouped.size().rename('value').to_frame()
    else:
        points = points.assign(_gprm_weight=np.asarray(weights))
        grouped = groupby_healpix(points, nside, order=order)
        result = grouped['_gprm_weight'].sum().rename('value').to_frame()

    from astropy_healpix import healpy as hp
    lon, lat = hp.pix2ang(nside, result.index.to_numpy(), nest=(order == 'nested'), lonlat=True)
    result['longitude'] = lon
    result['latitude'] = lat
    return result


def point_density(points, method='healpix', weights=None, **kwargs):
    """
    Coherent entry point for point density / binned counts on the sphere.

    :param points: GeoDataFrame with point geometry.
    :param method: 'healpix' for equal-area binned counts (see healpix_density), or 'kde' for a
        continuous spherical kernel density surface (see spherical_kde).
    :param weights: optional per-point weights -- an array/Series the same length as ``points``,
        or the name of a column in ``points``.
    :param kwargs: forwarded to the selected method: ``nside``/``order`` for 'healpix';
        ``bandwidth``/``sampling``/``region`` for 'kde'.
    """
    if isinstance(weights, str):
        weights = points[weights].to_numpy()

    if method == 'healpix':
        return healpix_density(points, weights=weights, **kwargs)
    elif method == 'kde':
        return spherical_kde(points.geometry.x.to_numpy(), points.geometry.y.to_numpy(),
                             weights=weights, **kwargs)
    else:
        raise ValueError("Unknown method {!r} for point_density. Choose 'healpix' or 'kde'.".format(method))


def dominant_class(points, class_column, nside=32, order='ring', normalize=True):
    """
    For each HEALPix pixel, find which class (from ``points[class_column]``) is most represented.

    Bins each class's points into the same HEALPix pixel grid independently (via
    healpix_density), then compares. By default each class's per-pixel counts are normalized to
    that class's own total count first, so a class with more points overall does not win
    everywhere by default regardless of spatial pattern -- pass ``normalize=False`` to compare
    raw counts/weights instead. Ties are broken by column order (first class wins).

    :param points: GeoDataFrame with point geometry.
    :param class_column: name of the column in ``points`` giving each point's class/category.
    :param nside: HEALPix resolution parameter.
    :param order: 'ring' or 'nested' HEALPix pixel ordering.
    :param normalize: if True (default), compare each class's counts as a fraction of that
        class's total before taking the argmax.
    :returns: DataFrame indexed by pixel id, with one column per class (raw counts), a
        'dominant' column (the winning class label), and 'longitude'/'latitude' pixel centers.
        Only pixels with at least one point from any class are included.
    """
    classes = points[class_column].unique()

    per_class = {}
    for cls in classes:
        subset = points[points[class_column] == cls]
        per_class[cls] = healpix_density(subset, nside=nside, order=order)['value']

    # Each class's Series only has entries for pixels where that class has >=1 point, so their
    # union (this concat) can never produce an all-NaN row: every row has a real point from at
    # least one class. Missing entries for other classes are genuine zeros, not unknowns.
    table = pd.DataFrame(per_class).fillna(0.)

    if normalize:
        comparison = table.divide(table.sum(axis=0), axis='columns')
    else:
        comparison = table

    result = table.copy()
    result['dominant'] = comparison.idxmax(axis='columns')

    from astropy_healpix import healpy as hp
    lon, lat = hp.pix2ang(nside, result.index.to_numpy(), nest=(order == 'nested'), lonlat=True)
    result['longitude'] = lon
    result['latitude'] = lat
    return result


def healpix_bin_geometries(pixel_ids, nside, values=None, order='ring', step=4):
    """
    Return a GeoDataFrame of exact HEALPix pixel boundary polygons for the given pixel ids.

    The counterpart to groupby_healpix/healpix_density: turns bin ids (and whatever values you
    computed for them, e.g. via ``grouped.median()`` or a bespoke ``.apply()``) back into
    plottable/exportable spherical geometry, without needing to touch HEALPix boundary math.

    :param pixel_ids: array-like of HEALPix pixel indices (e.g. from healpix_density's index,
        or groupby_healpix's 'bin_id' column) -- not necessarily every pixel on the sphere.
    :param nside: HEALPix resolution parameter the pixel ids belong to.
    :param values: optional array-like, same length as pixel_ids, attached as a 'value' column.
    :param order: 'ring' or 'nested', matching whatever assigned the pixel ids.
    :param step: boundary vertices per edge (4*step corners per pixel).
    :returns: GeoDataFrame with one row per pixel id, polygon geometry, and (if given) 'value'.

    .. warning::
       This builds each pixel's polygon directly from its lon/lat corners. A pixel that straddles
       the antimeridian (lon +/-180 deg) will produce a wrapped/self-intersecting polygon in this
       planar representation -- the same caveat that applies to other planar-geometry helpers in
       this codebase (e.g. the latitude-band clipping in utils.spatial). Not handled here.
    """
    from ._optional import require
    require('astropy_healpix', 'HEALPix pixel geometry')
    from astropy_healpix import healpy as hp

    pixel_ids = np.asarray(pixel_ids)
    lon, lat = hp.boundaries_lonlat(pixel_ids, step, nside, order=order)
    lon = lon.to_value('deg')
    lat = lat.to_value('deg')

    polygons = [Polygon(zip(lon[i], lat[i])) for i in range(len(pixel_ids))]

    result = gpd.GeoDataFrame({'bin_id': pixel_ids}, geometry=polygons, crs='EPSG:4326')
    if values is not None:
        result['value'] = np.asarray(values)
    return result


def plot_groups(pixel_ids, nside, bin_values, fig=None, filename=None, order='ring',
                grid_resolution=0.2, color_range=None, cmap='hot', reverse=True,
                pen='0.1p,gray50', transparency=0, **kwargs):
    """
    Generate a visual representation of HEALPix-binned data (e.g. from healpix_density, or a
    caller's own ``groupby_healpix(...).agg()``). The result can either be added to a pygmt
    figure or saved to a GIS file (format taken from the ``filename`` extension, e.g. shp,
    geojson, gpkg -- whatever geopandas' ``to_file`` supports).

    :param pixel_ids: HEALPix pixel ids, e.g. healpix_density's result index.
    :param nside: HEALPix resolution parameter the pixel ids belong to.
    :param bin_values: one value per pixel id, used for the fill colour.
    :param order: 'ring' or 'nested', matching whatever assigned the pixel ids.
    """
    from ._optional import require
    pygmt = require('pygmt', 'plotting spatially binned data')

    pixel_ids = np.asarray(pixel_ids)
    bin_values = np.asarray(bin_values)
    polygons = healpix_bin_geometries(pixel_ids, nside, values=bin_values, order=order)

    if filename:
        polygons.to_file(filename)

    if fig:
        from astropy_healpix import healpy as hp

        grid_lon, grid_lat = np.meshgrid(np.arange(-180., 180., grid_resolution),
                                        np.arange(-90., 90., grid_resolution))
        grid_pixel = hp.ang2pix(nside, np.radians(90. - grid_lat.ravel()),
                                np.radians(grid_lon.ravel()), nest=(order == 'nested'))

        full_sky = np.full(12 * nside**2, np.nan)
        full_sky[pixel_ids] = bin_values
        grid_z = full_sky[grid_pixel].reshape(grid_lon.shape)

        ds = xr.DataArray(grid_z, coords=[('lat',grid_lat[:,0]), ('lon',grid_lon[0,:])], name='z')

        if not color_range:
            color_range = (np.nanmin(bin_values), np.nanmax(bin_values))
            reverse = True
        pygmt.makecpt(cmap=cmap, series='{:f}/{:f}'.format(color_range[0],color_range[1]),
                      reverse=reverse, background='o')

        fig.grdimage(ds, transparency=transparency, cmap=True, nan_transparent=True)
        fig.plot(data=polygons, pen=pen, transparency=transparency, close=True, **kwargs)


def spherical_kde(lons, lats, weights=None, bandwidth=0.05, sampling=0.5, region='d'):
    """
    Generate a kernel density map in spherical coordinates.

    Uses a haversine-metric KDE, so the density estimate itself is a genuine great-circle
    calculation; the output grid is a plain equal-angle lon/lat mesh, which is not equal-area
    (cells shrink towards the poles), so per-cell values here are not directly comparable to a
    HEALPix bin count without an area correction.

    :param lons, lats: point coordinates, in degrees.
    :param weights: optional per-point weights, same length as lons/lats.
    :param bandwidth: KDE bandwidth, in radians of great-circle distance.
    :param sampling: output grid spacing, in degrees.
    :param region: 'd' for [-180,180], 'g' for [0,360] longitude range.
    """
    from ._optional import require
    require('sklearn', 'spherical kernel density estimation')
    from sklearn.neighbors import KernelDensity

    if region=='d':
        region = [-180,180,-90,90]
    elif region=='g':
        region = [0,360,-90,90]

    xgrid = np.arange(region[0],region[1]+sampling,sampling)
    ygrid = np.arange(region[2],region[3]+sampling,sampling)
    X, Y = np.meshgrid(xgrid, ygrid)
    xy = np.radians(np.vstack([Y.ravel(), X.ravel()]).T)

    latlon = np.vstack([lats, lons]).T

    kde = KernelDensity(bandwidth=bandwidth, metric='haversine')
    kde.fit(np.radians(latlon), sample_weight=weights)

    Z = np.exp(kde.score_samples(xy)).reshape(X.shape)

    ds = xr.DataArray(Z, coords=[('lat',ygrid), ('lon',xgrid)], name='z')

    return ds
