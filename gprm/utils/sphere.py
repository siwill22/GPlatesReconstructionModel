import numpy as np
import pygplates
from scipy import spatial
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
    ophis   = np.radians(sample_points_lons)
    oxyzs=rtp2xyz(np.ones(np.shape(othetas)), othetas, ophis)

    d,l = tree.query(oxyzs, k=k, distance_upper_bound=distance_upper_bound)

    return d,l


def healpix_mesh(nSide):
    #import healpy as hp
    from astropy_healpix import healpy as hp
    othetas,ophis = hp.pix2ang(nSide,np.arange(12*nSide**2))
    othetas = np.pi/2-othetas
    ophis[ophis>np.pi] -= np.pi*2

    # ophis -> longitude, othetas -> latitude
    return np.degrees(ophis), np.degrees(othetas)


## Some functions for data binning
def groupby_healpix(df, equal_area_points, return_point_indices=True):
    # Given a (geo)dataframe with irregularly distributed point in lat,long, 
    # and equal area point distribution, this functions will bin the dataframe point
    # based on the equal area pixel extents (using the underlying healpix function).
    # The function returns a groupby object that can be interrogated to get the 
    # statistics of points within each bin, and (optionally) a list of point indices  

    # TODO force input to be dataframe??
    bin_counts, bin_indices = equal_area_points.point_feature_heatmap(
        [pygplates.PointOnSphere(point) for point in zip(df['Latitude'],
                                                         df['Longitude'])], return_indices=True)

    point_indices = np.unique(bin_indices)
    df['bin_id'] = bin_indices

    grouped_points = df.groupby(by=['bin_id'])

    #binned_df = pd.DataFrame(columns=df.columns)
    #for i,group in enumerate(grouped_points.groups):
    #    points_selection = grouped_points.get_group(group)
    #    binned_df.loc[i] = points_selection.median()

    if return_point_indices:
        return grouped_points, point_indices
    else:
        return grouped_points


def create_circles(centre_points, circle_radius_degrees=0.5, polygon=None):
    # create small circles centred on the provide points, with user-specified radii
    # optionally, limit the points to those which lie within a provided polygon geometry
    # requires pyshtools
    import pyshtools as pysh

    if polygon:
        pts_in_polygon = []
        for pt in centre_points.multipoint.get_points():
            if geom.is_point_in_polygon(pt):
                pts_in_polygon.append(pt.to_lat_lon())
        #print(pts_in_polygon[0])
    
        circle_list = []
        for lat,lon in pts_in_polygon:
            circle_list.append(pysh.utils.MakeCircleCoord(lat,lon,circle_radius_degrees))
            
    else:
        circle_list = []
        for lon,lat in zip(centre_points.longitude, centre_points.latitude):
            circle_list.append(pysh.utils.MakeCircleCoord(lat,lon,circle_radius_degrees))

    return circle_list

    
def create_graticule(graticule_spacing=1, xlims=[-180,180], ylims=[-90,90], along_line_point_density=0.1, polygon=None):
    # create a latitude-longitude graticule, in the form of line geometries along each meridian and parallel 
    # optionally specify the graticule spacing, extent on the globe, sampling of the vertices along each line
    # optionally clip the graticule to a user-provided polygon geometry
    graticule_lines = []
    for xlevel in np.arange(xlims[0], xlims[1]+graticule_spacing, graticule_spacing):
        graticule_lines.append(np.array([(ylevel,xlevel) for ylevel in np.arange(ylims[0], ylims[1]+along_line_point_density, along_line_point_density)]))
        
    for ylevel in np.arange(ylims[0], ylims[1]+graticule_spacing, graticule_spacing):
        graticule_lines.append(np.array([(ylevel,xlevel) for xlevel in np.arange(xlims[0], xlims[1]+along_line_point_density, along_line_point_density)]))

    if polygon:
        graticule_lines_gdf = geometries_to_geodataframe(graticule_lines, geometry_type='PolyLine')
        polygon_gdf = geometries_to_geodataframe([polygon.to_lat_lon_array()])
        graticule_lines_in_polygon_gdf = gpd.overlay(graticule_lines_gdf, polygon_gdf, how='intersection', keep_geom_type=False)
        graticule_lines = dataframe_to_geometries(graticule_lines_in_polygon_gdf)

    return graticule_lines