'''
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
'''

import numpy as np
import pygplates
from scipy import spatial
import os
import tempfile
import pygmt
import xarray as xr

import pandas as pd


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
# function to call one of several methods and return Long/
# Lat arrays of points distributed on sphere

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
    multipoint_feature.set_name('Random Points from Marsaglia''s method')

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

    # if distance_upper_bound specified, assume that it is specified in degrees, convert to radians
    if not np.isnan(distance_upper_bound):
        distance_upper_bound = np.radians(distance_upper_bound)

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
def groupby_healpix(gdf, equal_area_points, return_point_indices=True):
    """
    Given a (geo)dataframe with irregularly distributed point in lat,long, 
    and equal area point distribution, this functions will bin the dataframe point
    based on the equal area pixel extents (using the underlying healpix function).
    The function returns a groupby object that can be interrogated to get the 
    statistics of points within each bin, and (optionally) a list of point indices  
    """
    # TODO force input to be dataframe??
    bin_counts, bin_indices = equal_area_points.point_feature_heatmap(
        [pygplates.PointOnSphere(point) for point in zip(gdf.geometry.y,
                                                         gdf.geometry.x)], return_indices=True)

    point_indices = np.unique(bin_indices)
    gdf['bin_id'] = bin_indices

    grouped_points = gdf.groupby(by=['bin_id'])

    #binned_df = pd.DataFrame(columns=df.columns)
    #for i,group in enumerate(grouped_points.groups):
    #    points_selection = grouped_points.get_group(group)
    #    binned_df.loc[i] = points_selection.median()

    if return_point_indices:
        return grouped_points, point_indices
    else:
        return grouped_points


def plot_groups(equal_area_points, bin_values, fig=None, filename=None, grid_resolution=0.2,
                color_range=None, cmap='hot', reverse=True, pen='0.1p,gray50', transparency=0, **kwargs):

    """
    Generate a visual representation of spatially binned data generated by 'groupby_healpix'.
    The result can either be added to a pygmt figure or saved to a GIS file (with the file
    type taken from the filename extension of the optional 'filename' parameter, e.g. shp, gmt, geojson) 
    """

    points = pygplates.MultiPointOnSphere(zip(equal_area_points.latitude,equal_area_points.longitude)).to_xyz_array() 

    radius = 1
    center = np.array([0, 0, 0])
    sv = spatial.SphericalVoronoi(points, radius, center)
    sv.sort_vertices_of_regions()

    polygon_features = []
    for region,zval in zip(sv.regions,bin_values):
        polygon = np.vstack((sv.vertices[region],sv.vertices[region][0,:]))
        polygon_feature = pygplates.Feature()
        polygon_feature.set_geometry(pygplates.PolygonOnSphere(polygon))
        polygon_feature.set_shapefile_attribute('zval', zval)
        polygon_features.append(polygon_feature)

    if filename:
        return_file = True
        pygplates.FeatureCollection(polygon_features).write(filename)

    else:
        return_file = False
        plot_file = tempfile.NamedTemporaryFile(delete=False, suffix='.gmt')
        plot_file.close()
        filename = plot_file.name
        pygplates.FeatureCollection(polygon_features).write(filename)

    if fig:
        grid_lon, grid_lat = np.meshgrid(np.arange(-180.,180.,grid_resolution),np.arange(-90.,90.,grid_resolution))
    
        d,l = sampleOnSphere(np.radians(equal_area_points.longitude),
                            np.radians(equal_area_points.latitude),
                            np.array(bin_values),
                            np.radians(grid_lon).ravel(),
                            np.radians(grid_lat).ravel(),
                            k=1)
        grid_z = np.array(bin_values)[l].reshape(grid_lon.shape)
        
        #spherical_triangulation = stripy.sTriangulation(lons=np.radians(equal_area_points.longitude), lats=np.radians(equal_area_points.latitude))
        #grid_z,_ = spherical_triangulation.interpolate_nearest(np.radians(grid_lon).ravel(), np.radians(grid_lat).ravel(), np.array(bin_values))

        ds = xr.DataArray(grid_z.reshape(grid_lon.shape), coords=[('lat',grid_lat[:,0]), ('lon',grid_lon[0,:])], name='z')

        pygmt.config(COLOR_FOREGROUND='white', COLOR_BACKGROUND='black')
        if not color_range:
            color_range = (np.nanmin(bin_values), np.nanmax(bin_values))
            reverse = True
        pygmt.makecpt(cmap=cmap, series='{:f}/{:f}'.format(color_range[0],color_range[1]), reverse=reverse)

        # This line would allow the polygons to be plotted directly with a colormap, but tends to crash when 
        # healpix of N=32 or greater is input
        #fig.plot(data=filename, pen=pen, color='+z', cmap=True, a='Z=zval', close=True, **kwargs)
        fig.grdimage(ds, transparency=transparency, cmap=True)
        fig.plot(data=filename, pen=pen, transparency=transparency, close=True, **kwargs)


    if not return_file:
        os.unlink(plot_file.name)


def spherical_kde(lons,lats, bandwidth=0.05, sampling=0.5, region='d'):
    """
    generate a kernel density map in spherical coordinates
    """
    from sklearn.neighbors import KernelDensity

    if region=='d':
        region = [-180,180,-90,90]
    elif region=='g':
        region = [0,360,-90,90]

    xgrid = np.arange(region[0],region[1]+sampling,sampling)
    ygrid = np.arange(region[2],region[3]+sampling,sampling)
    X, Y = np.meshgrid(xgrid, ygrid)
    grid = (np.ones(X.shape).ravel()==1)
    xy = np.vstack([Y.ravel(), X.ravel() ]).T
    xy = np.radians(xy[grid])
    
    latlon = np.vstack([lats, lons]).T
    
    kde = KernelDensity(bandwidth=bandwidth, metric='haversine')
    kde.fit(np.radians(latlon))

    # evaluate only on the land: -9999 indicates ocean
    Z = np.full(grid.shape[0], -9999.0)
    Z[grid] = np.exp(kde.score_samples(xy))
    Z = Z.reshape(X.shape)

    #coords = [('lat',ygrid), ('lon',xgrid)]
    ds = xr.DataArray(Z, coords=[('lat',ygrid), ('lon',xgrid)], name='z')
    
    return ds
