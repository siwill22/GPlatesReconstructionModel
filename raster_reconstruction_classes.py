from reconstruction_classes import *

#import pygplates
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import tempfile
import xarray as xr
from ptt.utils.call_system_command import call_system_command
import paleogeography as pg
#import stripy
import sphere_tools


class GplatesRaster(object):

    def __init__(self, filename, reconstruction_time=0., z_field_name='z'):

        self.gridX,self.gridY,self.gridZ = pg.load_netcdf(filename, z_field_name)
        self.reconstruction_time = reconstruction_time
        self.source_filename = filename

    def plot(self, show=False):
        plt.figure(figsize=(16,6))
        plt.contourf(self.gridX, self.gridY, self.gridZ,
                     20, cmap=plt.cm.BrBG_r)
        plt.gca().set_aspect('equal')
        plt.colorbar()
        if show:
            plt.show()

    def sample(self, point_lons, point_lats, order=0):

        LonGrid, LatGrid = np.meshgrid(self.gridX,self.gridY)
        d,l = sphere_tools.sampleOnSphere(LonGrid.flatten(),
                                    LatGrid.flatten(),
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

    def sample_using_gmt(self, point_lons, point_lats, extrapolate=False):

        dataout = np.vstack((np.asarray(point_lons),np.asarray(point_lats))).T
        xyzfile = tempfile.NamedTemporaryFile()
        grdtrack_file = tempfile.NamedTemporaryFile()

        np.savetxt(xyzfile.name,dataout)
        # Note the a -T option would find the nearest valid grid value,
        # if the point falls on a NaN grid node
        # adding -T+e returns the distance to the node
        if extrapolate:
            call_system_command(['gmt','grdtrack',xyzfile.name,'-G%s' % self.source_filename, '-T', '-nl','-V','>', grdtrack_file.name])
        else:
            call_system_command(['gmt','grdtrack',xyzfile.name,'-G%s' % self.source_filename, '-nl','-V','>', grdtrack_file.name])
        G=[]
        with open(grdtrack_file.name) as f:
            for line in f:
                if line[0] == '>':
                    continue
                else:
                    tmp = line.split()
                    G.append(float(tmp[2]))

        f.close()
        return np.array(G)

    def sample_using_stripy(self, point_lons, point_lats, order=0):
        import stripy

        LonGrid, LatGrid = np.meshgrid(self.gridX,self.gridY)
        tri = stripy.sTriangulation(lons=np.radians(LonGrid.flatten()),
                                    lats=np.radians(LatGrid.flatten()))

        point_z = tri.interpolate(np.radians(point_lons),np.radians(point_lats),
                                  zdata=self.gridZ.flatten(),
                                  order=order)

        return point_z[0]


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
