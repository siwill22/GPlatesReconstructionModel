import pygplates
import points_in_polygons
import points_spatial_tree
from proximity_query import *
from create_gpml import create_gpml_regular_long_lat_mesh
import numpy as np
from skimage import measure


def merge_polygons(polygons,rotation_model,
                   time=0,sampling=1.,area_threshold=None,filename=None,
                   return_raster=False):

    multipoints = create_gpml_regular_long_lat_mesh(sampling)
    grid_dims = (int(180/sampling)+1,int(360/sampling)+1)

    for multipoint in multipoints:
        for mp in multipoint.get_all_geometries():
            points = mp.to_lat_lon_point_list()

    bi = run_grid_pip(time,points,polygons,rotation_model,grid_dims)
    
    if return_raster:
        return bi
    
    else:
        # To handle edge effects, pad grid before making contour polygons  
        ## --- start
        pad_hor = np.zeros((1,bi.shape[1]))
        pad_ver = np.zeros((bi.shape[0]+2,1))
        pad1 = np.vstack((pad_hor,bi,pad_hor))      # add row of zeros to top and bottom
        pad2 = np.hstack((pad_ver,pad1,pad_ver))    # add row of zeros to left and right
        #pad3 = np.hstack((pad2,pad_ver))
        contours = measure.find_contours(pad2, 0.5, fully_connected='low')
        ## --- end
    
        contour_polygons = []
        contour_features = []
    
        for n,cp in enumerate(contours):
        
            # To handle edge effects again - strip off parts of polygon
            # due to padding, and adjust from image coordinates to long/lat
            # --- start
            cp[:,1] = (cp[:,1]*sampling)-sampling
            cp[:,0] = (cp[:,0]*sampling)-sampling
            cp[np.where(cp[:,0]<0.),0] = 0
            cp[np.where(cp[:,0]>180.),0] = 180
            cp[np.where(cp[:,1]<0.),1] = 0
            cp[np.where(cp[:,1]>360.),1] = 360
            ## --- end
        
            cpf = pygplates.PolygonOnSphere(zip(cp[:,0]-90,cp[:,1]-180))
            contour_polygons.append(cpf)
        
            feature = pygplates.Feature()
            feature.set_geometry(cpf)
            contour_features.append(feature)

        if filename is not None:
            pygplates.FeatureCollection(contour_features).write(filename)

        else:
            return contour_features


def force_polygon_geometries(input_features):
# given any pygplates feature collection, creates an output feature collection
# where all geometries are polygons based on the input geometries
# intended for use in forcing features that are strictly polylines to close

    polygons = []
    for feature in input_features: 
        for geom in feature.get_all_geometries():
            polygon = pygplates.Feature()
            polygon.set_geometry(pygplates.PolygonOnSphere(geom))
            polygon.set_reconstruction_plate_id(feature.get_reconstruction_plate_id())
            # some features in COBTerranes had invalid time ranges - these with throw an error if 
            # we try to create a new feature with same times
            if feature.get_valid_time()[0]>=feature.get_valid_time()[1]:
                polygon.set_valid_time(feature.get_valid_time()[0],feature.get_valid_time()[1])
                polygons.append(polygon)
    polygon_features = pygplates.FeatureCollection(polygons)

    return polygon_features


def polygon_area_threshold(polygons,area_threshold):
    
    polygons_larger_than_threshold = []
    for polygon in polygons:
        if polygon.get_geometry() is not None:
            if polygon.get_geometry().get_area()>area_threshold:
                polygons_larger_than_threshold.append(polygon)

    return polygons_larger_than_threshold


#This is a function to do fast point in polygon text
def run_grid_pip(time,points,polygons,rotation_model,grid_dims):

    reconstructed_polygons = []
    pygplates.reconstruct(polygons,rotation_model,reconstructed_polygons,time)

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

    bi = np.array(zval).reshape(grid_dims[0],grid_dims[1])

    return bi


# Function to run efficient point in/near polygons
# returns two numbers - one is distance to polygon edge,
# other is distance to polygon where distance is zero if inside
def run_grid_pnp(recon_time, 
                 points, 
                 spatial_tree_of_uniform_recon_points, 
                 polygons, 
                 rotation_model, 
                 distance_threshold_radians=2):

    reconstructed_polygons = []
    pygplates.reconstruct(polygons, rotation_model, reconstructed_polygons, recon_time)
    rpolygons = []
    for polygon in reconstructed_polygons:
        if polygon.get_reconstructed_geometry():
            rpolygons.append(polygon.get_reconstructed_geometry())
                
    res1 = find_closest_geometries_to_points_using_points_spatial_tree(points,
                                                                    spatial_tree_of_uniform_recon_points,
                                                                    rpolygons,
                                                                    distance_threshold_radians = distance_threshold_radians,
                                                                    geometries_are_solid = False)

    distance_to_polygon_boundary = np.array(zip(*res1)[0])

    # Make a copy of list of distances.
    distance_to_polygon = list(distance_to_polygon_boundary)

    # Set distance to zero for any points inside a polygon (leave other points unchanged).
    res2 = points_in_polygons.find_polygons_using_points_spatial_tree(points,
                                                                     spatial_tree_of_uniform_recon_points,
                                                                     rpolygons)
    for point_index, rpolygon in enumerate(res2):
        # If not inside any polygons then result will be None.
        if rpolygon:
            distance_to_polygon[point_index] = 0.0

    distance_to_polygon = np.array(distance_to_polygon)

    return distance_to_polygon,distance_to_polygon_boundary

