import math
from call_system_command import call_system_command
import pygplates
import numpy as np



def get_rotation_to_equator(point):
    """Get rotation from 'point' to nearest position on equator."""
    
    rotate_from_north_to_point = pygplates.FiniteRotation(pygplates.PointOnSphere.north_pole, point)
    
    if rotate_from_north_to_point.represents_identity_rotation():
        # Point coincides with North Pole, so just choose any rotation pole on the equator to rotate point with.
        # The point could end up anywhere along equator.
        rotate_from_north_to_point_pole = pygplates.PointOnSphere(0, 0)
        rotate_from_north_to_point_angle = 0.5 * math.pi # Rotate 90 degrees from North Pole to equator.
    else:
        rotate_from_north_to_point_pole, rotate_from_north_to_point_angle = rotate_from_north_to_point.get_euler_pole_and_angle()
    
    return pygplates.FiniteRotation(
            rotate_from_north_to_point_pole,
            0.5 * math.pi - rotate_from_north_to_point_angle)


def get_major_axis_orientation_angle(lat_lon_centroid, lat_lon_points):
    """
    Find the major axis using principle component analysis of a sequence of lat/lon points and
    return angle of major axis (in radians) with respect to the equator (along positive longitude).
    
    Using lat/lon space is probably not the best - might be better to do an equal-area projection instead
    (or even project onto a tangent plane, but that might be problematic for very long polygons).
    It works pretty well here because the input data is near the equator (ie, not too close to poles hopefully).
    
    TODO: Need to handle longitude wraparound of input lat/lon points somehow.
    """
    
    #
    # The following code was copied from:
    #
    #    Principal Components Analysis in 2D
    #            Miguel A. Lerma
    #            January 16, 2014
    #
    #    www.math.northwestern.edu/~mlerma/papers/princcomp2d.pdf
    #
    
    sumx, sumy, sumxx, sumyy, sumxy = 0.0, 0.0, 0.0, 0.0, 0.0
    
    for lat, lon in lat_lon_points:
        x, y = lon, lat
        sumx += x
        sumy += y
        sumxx += x*x
        sumyy += y*y
        sumxy += x*y
    
    num_points = len(lat_lon_points)
    
    # baricenter
    lat_centroid, lon_centroid = lat_lon_centroid
    xbar, ybar = lon_centroid, lat_centroid
    
    # variances and covariance
    varx = sumxx / num_points - xbar * xbar
    vary = sumyy / num_points - ybar * ybar
    covxy = sumxy / num_points - xbar * ybar
    sumvars = varx + vary
    diffvars = varx - vary
    discriminant = diffvars * diffvars + 4.0 * covxy * covxy
    sqrtdiscr = math.sqrt(discriminant)
    
    # eigenvalues
    lambdaplus = (sumvars + sqrtdiscr) / 2.0
    lambdaminus = (sumvars - sqrtdiscr) / 2.0
    # eigenvectors - these are the components of the two vectors
    aplus = varx + covxy - lambdaminus
    bplus = vary + covxy - lambdaminus
    aminus = varx + covxy - lambdaplus
    bminus = vary + covxy - lambdaplus
    # normalizing the vectors
    denomPlus = math.sqrt(aplus*aplus + bplus*bplus)
    denomMinus= math.sqrt(aminus*aminus + bminus*bminus)
    aParallel = aplus/denomPlus
    bParallel = bplus/denomPlus
    aNormal = aminus/denomMinus
    bNormal = bminus/denomMinus
    # semi axes
    k = 2 # scale factor
    majoraxis = k*math.sqrt(np.abs(lambdaplus))
    minoraxis = k*math.sqrt(np.abs(lambdaminus))
    
    
    # We just take x and y components of the 'major' eigenvector and calculate its angle of rotation.
    major_axis_angle = math.atan2(bplus, aplus)
    
    return major_axis_angle


def orientation_polygon_along_equator(polygon):
    """Rotate the polygon such that its major axis is aligned along the equator."""
    
    polygon_centroid = polygon.get_boundary_centroid()
    
    # Rotate polygon so its centroid in on equator (nearest point on equator to original centroid).
    rotation_to_equator = get_rotation_to_equator(polygon_centroid)
    equator_polygon = rotation_to_equator * polygon
    equator_polygon_centroid = rotation_to_equator * polygon_centroid
    
    # Tessellate polygon on equator so we have enough points to do PCA analysis of polygon's boundary.
    tesselated_equator_polygon = equator_polygon.to_tessellated(math.radians(1))
    major_axis_orientation_angle_radians = get_major_axis_orientation_angle(
            equator_polygon_centroid.to_lat_lon(),
            tesselated_equator_polygon.to_lat_lon_list())
    
    # Rotate polygon such that major axis is aligned with equator.
    reorient_major_axis = pygplates.FiniteRotation(
            equator_polygon_centroid,
            -major_axis_orientation_angle_radians)
    oriented_polygon = reorient_major_axis * equator_polygon
    
    #print 'Rotation to equator %s' % rotation_to_equator
    #print 'Rotation to reorient major axis %s' % reorient_major_axis
    
    composed_rotation =  reorient_major_axis * rotation_to_equator
    
    return oriented_polygon,composed_rotation


def orientation_polygon_along_meridian(polygon):

    median_longitude = np.median(polygon.to_lat_lon_array()[:,1])

    reorient_major_axis = pygplates.FiniteRotation(
            pygplates.PointOnSphere(90,0),
            -np.radians(median_longitude))
    oriented_polygon = reorient_major_axis * polygon

    return oriented_polygon,reorient_major_axis


def find_rotations_from_pca(polygon_features,num_iterations=20):

    total_applied_rotations = []
    oriented_polygon_features = []

    for polygon_feature in polygon_features:
    
        polygon = polygon_feature.get_geometry()
    
        # Re-orient polygon such that its major axis is aligned along the equator.
        oriented_polygon,composed_rotation1 = orientation_polygon_along_meridian(polygon)

        #oriented_polygon,composed_rotation1 = orientation_polygon_along_equator(polygon)
        for i in range(num_iterations):   
            oriented_polygon,composed_rotation2 = orientation_polygon_along_equator(oriented_polygon)
            total_applied_rotation = composed_rotation2 * composed_rotation1
            composed_rotation1 = total_applied_rotation

        oriented_polygon,reorient_major_axis = orientation_polygon_along_meridian(oriented_polygon)
        total_applied_rotation = reorient_major_axis * total_applied_rotation


        total_applied_rotations.append(total_applied_rotation)

        # Clone original feature and replace cloned polygon with reoriented polygon.
        oriented_polygon_feature = polygon_feature.clone()
        oriented_polygon_feature.set_geometry(oriented_polygon)

        oriented_polygon_features.append(oriented_polygon_feature)

    return total_applied_rotations,oriented_polygon_features


def grdcontour2feature(grdfile,clevel):

    # call GMT to get a single contour at the specified value of clevel
    call_system_command(['gmt',
                         'grdcontour',
                         grdfile,
                         '-C+%0.8f' % clevel,
                         '-S4',
                         '-Dcontour_%c.txt',
                         '-V'])

    # read in the GMT delimited xyz ascii file, 
    # create a list of lists with polygon coordinates
    f = open('./contour_C.txt', 'r')

    polygons = []
    contourlist = []
    for line in f:
        if line[0] == '>':
            if len(contourlist)>0:
                polygons.append(contourlist)
            contourlist = []
        else:
            line = line.split()
            contourlist.append([float(j) for j in line])
            #break

    # create gplates-format features
    polyline_features = []
    for p in polygons:
        pf = pygplates.PolylineOnSphere(zip(zip(*p)[1],zip(*p)[0]))
        polyline_features.append(pf)

    # use join to handle polylines split across dateline
    joined_polyline_features = pygplates.PolylineOnSphere.join(polyline_features)

    # force polylines to be polygons
    joined_polygon_features = []
    for geom in joined_polyline_features:
        polygon = pygplates.Feature()
        polygon.set_geometry(pygplates.PolygonOnSphere(geom))
        joined_polygon_features.append(polygon)
    
    return joined_polygon_features
