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

import os
import math
import numpy as np
import pygplates
from ptt.utils import points_in_polygons
from ptt.utils import points_spatial_tree
from ptt.utils.proximity_query import find_closest_geometries_to_points_using_points_spatial_tree

from .fileio import write_xyz_file

import scipy.interpolate as spi
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from scipy.spatial import cKDTree

import pandas as pd
import xarray as xr
import pygmt


def xyzfile_to_spatial_tree_of_points(xyzfile):

    data = np.loadtxt(xyzfile)

    points = [pygplates.PointOnSphere(lat, lon) for lat, lon in zip(data[:,1],data[:,0])]

    return points,data


def reconstruct_raster_stage(static_polygon_features,
                             rotation_model,
                             time_from,
                             time_to,
                             uniform_recon_points,
                             spatial_tree_of_uniform_recon_points,
                             anchor_plate_id=0):

    print('Reconstruct static polygons...')

    # Reconstruct the multipoint feature.
    recon_static_polygon_features = []
    pygplates.reconstruct(static_polygon_features, rotation_model, recon_static_polygon_features, time_to, anchor_plate_id=anchor_plate_id)

    # Extract the polygons and plate IDs from the reconstructed static polygons.
    recon_static_polygons = []
    recon_static_polygon_plate_ids = []
    for recon_static_polygon_feature in recon_static_polygon_features:
        recon_plate_id = recon_static_polygon_feature.get_feature().get_reconstruction_plate_id()
        recon_polygon = recon_static_polygon_feature.get_reconstructed_geometry()

        recon_static_polygon_plate_ids.append(recon_plate_id)
        recon_static_polygons.append(recon_polygon)

    print('Find static polygons...')

    # Find the reconstructed static polygon (plate IDs) containing the uniform (reconstructed) points.
    #
    # The order (and length) of 'recon_point_plate_ids' matches the order (and length) of 'uniform_recon_points'.
    # Points outside all static polygons return a value of None.
    recon_point_plate_ids = points_in_polygons.find_polygons_using_points_spatial_tree(
            uniform_recon_points, spatial_tree_of_uniform_recon_points, recon_static_polygons, recon_static_polygon_plate_ids)

    print('Group by polygons...')

    # Group recon points with plate IDs so we can later create one multipoint per plate.
    recon_points_grouped_by_plate_id = {}
    for point_index, point_plate_id in enumerate(recon_point_plate_ids):
        # Reject any points outside all reconstructed static polygons.
        if point_plate_id is None:
            continue

        # Add empty list to dict if first time encountering plate ID.
        if point_plate_id not in recon_points_grouped_by_plate_id:
            recon_points_grouped_by_plate_id[point_plate_id] = []

        # Add to list of points associated with plate ID.
        recon_point = uniform_recon_points[point_index]
        recon_points_grouped_by_plate_id[point_plate_id].append(recon_point)

    print('Reverse reconstruct points...')

    # Reconstructed points.
    recon_point_lons = []
    recon_point_lats = []

    # Present day points associated with reconstructed points.
    point_lons = []
    point_lats = []

    # Create a multipoint feature for each plate ID and reverse-reconstruct it to get present-day points.
    #
    # Iterate over key/value pairs in dictionary.
    for plate_id, recon_points_in_plate in recon_points_grouped_by_plate_id.items():
        # Reverse reconstructing a multipoint is much faster than individually reverse-reconstructing points.
        multipoint_feature = pygplates.Feature()
        multipoint_feature.set_geometry(pygplates.MultiPointOnSphere(recon_points_in_plate))
        multipoint_feature.set_reconstruction_plate_id(plate_id)

        # Reverse reconstruct the multipoint feature.
        pygplates.reverse_reconstruct(multipoint_feature, rotation_model, time_to, anchor_plate_id=anchor_plate_id)

        #Forward reconstruct multipoint to
        multipoint_at_from_time = []
        pygplates.reconstruct(multipoint_feature,rotation_model,multipoint_at_from_time,time_from, anchor_plate_id=anchor_plate_id)

        # Extract reverse-reconstructed geometry.
        multipoint = multipoint_at_from_time[0].get_reconstructed_geometry()

        # Collect present day and associated reconstructed points.
        for point_index, point in enumerate(multipoint):
            lat, lon = point.to_lat_lon()
            point_lons.append(lon)
            point_lats.append(lat)

            recon_point = recon_points_in_plate[point_index]
            recon_lat, recon_lon = recon_point.to_lat_lon()
            recon_point_lons.append(recon_lon)
            recon_point_lats.append(recon_lat)


    #print('Sample present-day grid...')
    # Query present-day grid using present-day points.
    #
    # TODO: Note sure what happens in regions where there's no data in grid (need to ignore those points).
    #data = data_grid.ev(point_lons, point_lats)
    #data = [1.0] * len(recon_point_lons)
    #data = sample_grid_using_scipy(point_lons,point_lats,grdfile)

    return recon_point_lons,recon_point_lats,point_lons,point_lats


def run_grid_pip(recon_time, points, polygons, rotation_model, anchor_plate_id=0):

    reconstructed_polygons = []
    pygplates.reconstruct(polygons, rotation_model, reconstructed_polygons, recon_time, anchor_plate_id=anchor_plate_id)
    rpolygons = []
    for polygon in reconstructed_polygons:
        if polygon.get_reconstructed_geometry():
            rpolygons.append(polygon.get_reconstructed_geometry())
    polygons_containing_points = points_in_polygons.find_polygons(points, rpolygons)
    lat = []
    lon = []
    zval = []
    for pcp,point in zip(polygons_containing_points,points):
        lat.append(point.get_latitude())
        lon.append(point.get_longitude())
        if pcp is not None:
            zval.append(1)
        else:
            zval.append(0)
    return zval


def run_grid_pnp(recon_time, points, spatial_tree_of_uniform_recon_points, polygons, rotation_model,
                 anchor_plate_id=0, distance_threshold_radians=2):

    reconstructed_polygons = []
    pygplates.reconstruct(polygons, rotation_model, reconstructed_polygons, recon_time, anchor_plate_id=anchor_plate_id)
    rpolygons = []
    for polygon in reconstructed_polygons:
        if polygon.get_reconstructed_geometry():
            rpolygons.append(polygon.get_reconstructed_geometry())
    res = find_closest_geometries_to_points_using_points_spatial_tree(points,
                                                                     spatial_tree_of_uniform_recon_points,
                                                                     rpolygons,
                                                                     distance_threshold_radians = distance_threshold_radians,
                                                                     geometries_are_solid = True)

    pnp_test = []
    for index in res:
        if index is not None:
            pnp_test.append(1)
        else:
            pnp_test.append(0)

    return pnp_test


###########################################################


def reconstruct_raster(raster, static_polygons, rotation_model, time_from, time_to,
                       grid_sampling=1., anchor_plate_id=0, sampling_method='gmt', return_type='xarray'):

    grid_longitudes, grid_latitudes = np.meshgrid(np.arange(-180.,180.0001,grid_sampling), np.arange(-90.,90.0001,grid_sampling))
    grid_longitudes = grid_longitudes.flatten()
    grid_latitudes = grid_latitudes.flatten()

    points = [pygplates.PointOnSphere(point) for point in zip(grid_latitudes, grid_longitudes)]

    spatial_tree_of_uniform_recon_points = points_spatial_tree.PointsSpatialTree(points)

    (time_to_point_lons,
     time_to_point_lats,
     time_from_point_lons,
     time_from_point_lats) = reconstruct_raster_stage(static_polygons,
                                                      rotation_model,
                                                      time_from,
                                                      time_to,
                                                      points,
                                                      spatial_tree_of_uniform_recon_points,
                                                      anchor_plate_id)

    if isinstance(raster, xr.core.dataarray.DataArray):
        point_raster_values = pygmt.grdtrack(grid=raster, 
                                            points=pd.DataFrame(data={'x':time_from_point_lons,
                                                                    'y':time_from_point_lats}), 
                                            newcolname='z')['z']

    # retained for backward compatibility
    else:
        if sampling_method=='scipy':
            point_raster_values = raster.sample(time_from_point_lons, time_from_point_lats)
        elif sampling_method=='gmt':
            point_raster_values = raster.sample_using_gmt(time_from_point_lons, time_from_point_lats)
        elif sampling_method=='stripy':
            point_raster_values = raster.sample_using_stripy(time_from_point_lons, time_from_point_lats)
    
    if return_type=='ndarray':
        return xyz2grd(time_to_point_lons,time_to_point_lats,point_raster_values,
                       grid_longitudes,grid_latitudes)
    else:
        return time_to_point_lons, time_to_point_lats, point_raster_values


def xyz2grd(point_lons,point_lats,point_zvals,grid_lons,grid_lats):
    """
    Taking as input lists/flat arrays of longitudes, latitudes, and z values - already on a regular
    grid, convert to a 2D array containing the same values (suitable to be saved as netcdf for example)
    """
    # https://stackoverflow.com/questions/30655749/how-to-set-a-maximum-distance-between-points-for-interpolation-when-using-scipy
    
    # TODO change this to the pygmt implementation?

    if grid_lons.ndim == 1:
        grid_lons, grid_lats = np.meshgrid(grid_lons, grid_lats)

    grid_sampling = np.abs(grid_lats[1]-grid_lats[0])
    xy = np.vstack((point_lons, point_lats)).T

    # Construct kd-tree, functionality copied from scipy.interpolate
    tree = cKDTree(xy)
    xi = _ndim_coords_from_arrays((grid_lons, grid_lats), ndim=xy.shape[1])
    dists, indexes = tree.query(xi)

    grid_interp = spi.griddata(xy, point_zvals,
                              (grid_lons, grid_lats),
                              method='nearest')

    # Copy original result but mask missing values with NaNs
    result = grid_interp[:]
    result[dists > grid_sampling/2.] = np.nan

    return result


def to_anchor_plate(grid, reconstruction_model, reconstruction_time,
                    new_anchor_plate_id, old_anchor_plate_id=0,
                    region=None, spacing=1.):
    """
    Given an xarray raster object, derive a 
    new raster with all values rotated into the frame of a different reference plate.
    
    Optionally, specify a new region and grid sampling (default is to take the same region
    and grid sampling as the input)
    """

    # should be a better way to do this (ie pass the gridfile to grdtrack directly)
    if not isinstance(grid, xr.core.dataarray.DataArray):
        if os.path.isfile(grid):
            grid = xr.open_dataarray(grid)
        else:
            raise ValueError('Unrecognised input grid object: {}'.format(grid))
        
    # Sort out the coordinates (move this to a generic function??)
    coord_keys = [key for key in grid.coords.keys()]

    if region:
        coords = [('lat',np.arange(region[2],region[3]+spacing, spacing)), ('lon',np.arange(region[0],region[1]+spacing, spacing))]
        XX,YY = np.meshgrid(np.arange(region[0],region[1]+spacing, spacing),
                            np.arange(region[2],region[3]+spacing, spacing))
    else:
        coords = [('lat',grid.coords[coord_keys[0]].data), ('lon',grid.coords[coord_keys[1]].data)]
        XX,YY = np.meshgrid(grid.coords[coord_keys[1]].data, grid.coords[coord_keys[0]].data)

    # us a multipoint as the basis for the raster coordinates. 
    # multipoints vertices are the coordinates of the adjusted grid - then we unrotate them
    # to the source reference frame, sample the grid and assign values to new raster
    mp = pygplates.MultiPointOnSphere(zip(YY.flatten().tolist(),XX.flatten().tolist()))

    adjustment_rotation = reconstruction_model.rotation_model.get_rotation(reconstruction_time,
                                                                           new_anchor_plate_id,
                                                                           0,
                                                                           old_anchor_plate_id)

    multipoint_at_from_time = adjustment_rotation * mp

    point_raster_values = pygmt.grdtrack(grid=grid, 
                                         points=pd.DataFrame(data={'x':multipoint_at_from_time.to_lat_lon_array()[:,1],
                                                                   'y':multipoint_at_from_time.to_lat_lon_array()[:,0]}), 
                                         newcolname='z',
                                         no_skip=True)

    ds = xr.DataArray(np.array(point_raster_values['z']).reshape(XX.shape),
                      coords=coords,
                      name='z')

    return ds

