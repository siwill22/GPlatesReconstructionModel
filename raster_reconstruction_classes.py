from reconstruction_classes import *

#import pygplates
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import paleogeography as pg
from pigplates import sphere_tools as pigsph


class GplatesRaster(object):

    def __init__(self, filename, reconstruction_time=0.):

        self.gridX,self.gridY,self.gridZ = pg.load_netcdf(filename)
        self.reconstruction_time = reconstruction_time
        self.source_filename = filename

    def plot(self):
        plt.figure(figsize=(16,6))
        plt.contourf(self.gridX, self.gridY, self.gridZ,
                     20, cmap=plt.cm.BrBG_r)
        plt.gca().set_aspect('equal')
        plt.colorbar()
        plt.show()

    def sample(self, point_lons, point_lats):

        Xg, Yg = np.meshgrid(self.gridX,self.gridY)
        d,l = pigsph.sampleOnSphere(Xg.flatten(),
                                    Yg.flatten(), 
                                    self.gridZ.flatten(),
                                    np.array(point_lons),
                                    np.array(point_lats),
                                    k=4)

        #print d,l
        # based on http://earthpy.org/interpolation_between_grids_with_ckdtree.html
        # note also that where d is zero, we get a divide by zero error - hence, these 
        # values are (currently) set to one 
        w = np.divide(1.,d**2, out=np.ones_like(d), where=d!=0)
        point_z = np.sum(w * self.gridZ.flatten().ravel()[l],axis=1) / np.sum(w,axis=1)

        return point_z

    def cross_section(self, PtLons, PtLats):
        return CrossSection(self, PtLons, PtLats)


class CrossSection(object):

    def __init__(self, target_raster, PtLons, PtLats):

        self.GreatCirclePoints,self.ProfilePoints,arc_distance = pg.create_profile_points(PtLons,PtLats)
        # create an array of distances along profile in km, starting at zero
        
        self.profileX_kms = np.arange(0,self.ProfilePoints.shape[0])*arc_distance

        # extract the values from the (smoothed) topography grid along the profile
        self.grid_values = pg.create_slice(target_raster.gridX,
                                       target_raster.gridY,
                                       target_raster.gridZ,
                                       self.GreatCirclePoints, self.ProfilePoints)

        self.cross_section_geometry = pygplates.PolylineOnSphere(self.GreatCirclePoints)
        self.source_filename = target_raster.source_filename

    def plate_boundary_intersections(self, shared_boundary_sections):

        (self.subduction_intersections,
         self.ridge_intersections,
         self.other_intersections) = pg.plate_boundary_intersections(self.cross_section_geometry,
                                                                     shared_boundary_sections,
                                                                     self.profileX_kms)



class PresentDayAgeGrid(object):
    
    def __init__(self):
        GplatesRaster.__init__(self)


class TimeDependentRasterSequence(object):

    def __init__(self, name):

        self.name = name
