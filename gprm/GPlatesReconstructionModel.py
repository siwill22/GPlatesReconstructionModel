'''
MIT License

Copyright (c) 2017-2020 Simon Williams

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

#from xarray.core.utils import to_0d_object_array
import pygplates

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import os
from io import StringIO
from pprint import pprint
import tempfile
import copy
import xarray as xr

import gprm.utils as utils

from ptt.utils.proximity_query import find_closest_geometries_to_points
from gprm.utils.geometry import distance_between_reconstructed_points_and_features

import ptt.subduction_convergence as sc
from ptt.utils.call_system_command import call_system_command
from ptt.resolve_topologies import resolve_topologies as topology2gmt

import warnings


DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'Data')


class ReconstructionModel(object):
    """
    Class to contain the various elements of a GPlates-format reconstruction model, 
    including some or all of plate rotation models, topologies, and reconstructable 
    geometries

    Attributes
    ----------
    name : str
        a string with a name for the reconstruction model
    rotation_model : pygplates.RotationModel
        the rotation model in memory
    rotation_files : list of str
        list of filename(s) of rotation files used to build the
        rotation model
    static_polygons : list of str
        list of static polygon feature collections
    static_polygon_files : list of str
        list of static polygon filename(s)
    dynamic_polygons : list of pygplates.FeatureCollection
        dynamic polygon features loaded into memory
    dynamic_polygon_files : list of str
        list of dynamic polygon filename(s)
    coastlines : list of str
        list of coastline feature collections
    coastlines_files : list of str
        list of coastline filename(s)
    continent_polygons : list of str
        list of continent polygons feature collections
    continent_polygons_files : list of str
        list of continent polygons filename(s)
    """

    def __init__(self, name=None):
        self.name = name
        self.rotation_model = []    # creates a new empty list for each
        self.rotation_files = []
        self.static_polygons = []
        self.static_polygon_files = []
        self.dynamic_polygons = []
        self.dynamic_polygon_files = []
        self.coastlines = []
        self.coastlines_files = []
        self.continent_polygons = []
        self.continent_polygons_files = []

    def info(self, show_full_paths=False):

        print('Name: {:s}'.format(self.name))

        for item in [('Rotation', self.rotation_files),
                     ('Static Polygon', self.static_polygon_files), 
                     ('Coastlines', self.coastlines_files),
                     ('Continent Polygon', self.continent_polygons), 
                     ('Dynamic Polygon', self.dynamic_polygon_files)]:
            print('{:s} Files:'.format(item[0]))
            for f in item[1]:
                if show_full_paths:
                    print('   - {:s}'.format(f))
                else:
                    print('   - {:s}'.format(os.path.split(f)[1]))


    def add_rotation_model(self, rotation_file, replace=False):
        """
        Add a rotation model by specifying the path and filename to a .rot file_extension.
        Can be called multiple times to add a series of file into a single object instance.

        :param rotation_file: (str) name of rotation file
        :param replace: (bool, optional) A flag to specify whether to add to existing rotation model (if present), or 
            replace the current contents (default is False)
        """
        if not os.path.isfile(rotation_file):
            raise ValueError('Unable to find file {:s}'.format(rotation_file))

        if replace:
            self.rotation_model = []
            self.rotation_files = []

        self.rotation_files.append(rotation_file)
        self.rotation_model = pygplates.RotationModel(self.rotation_files)

    def add_static_polygons(self, static_polygons_file, replace=False):
        """
        Add a set of static polygons to the reconstruction model object by specifying
        path and file to a GPlates compatible file format (gpml, gpmlz, shp, gmt)
        """
        if not os.path.isfile(static_polygons_file):
            raise ValueError('Unable to find file {:s}'.format(static_polygons_file))

        if replace:
            self.static_polygons = []
            self.static_polygon_files = []

        self.static_polygon_files.append(static_polygons_file)
        self.static_polygons.append(static_polygons_file)  # Should this be loaded into memory like dynamic polygons??

    def add_dynamic_polygons(self, dynamic_polygons_file, replace=False):
        """
        Add topology files to be used in resolving topological polygons.
        Can be called multiple times to add a series of file into a single object instance.
        """
        if not os.path.isfile(dynamic_polygons_file):
            raise ValueError('Unable to find file {:s}'.format(dynamic_polygons_file))

        if replace:
            self.dynamic_polygons = []
            self.dynamic_polygon_files = []

        self.dynamic_polygon_files.append(dynamic_polygons_file)
        self.dynamic_polygons = [pygplates.FeatureCollection(dpfile) for dpfile in self.dynamic_polygon_files]

    def add_coastlines(self, coastlines_file, replace=False):
        """
        Add a set of coastline polygons to the reconstruction model object by specifying
        path and file to a GPlates compatible file format (gpml, gpmlz, shp, gmt)
        """
        if not os.path.isfile(coastlines_file):
            raise ValueError('Unable to find file {:s}'.format(coastlines_file))

        if replace:
            self.coastlines = []
            self.coastlines_files = []

        self.coastlines_files.append(coastlines_file)
        self.coastlines.append(coastlines_file)

    def add_continent_polygons(self, continent_polygons_file, replace=False):
        """
        Add a set of continent polygons to the reconstruction model object by specifying
        path and file to a GPlates compatible file format (gpml, gpmlz, shp, gmt)
        """
        if not os.path.isfile(continent_polygons_file):
            raise ValueError('Unable to find file {:s}'.format(continent_polygons_file))
        
        if replace:
            self.continent_polygons = []
            self.continent_polygons_files = []

        self.continent_polygons_files.append(continent_polygons_file)
        self.continent_polygons.append(continent_polygons_file)

    def from_web_service(self, model='MULLER2016', url='https://gws.gplates.org'):
        """
        Add a reconstruction model directly from the GPlates web service.
        """
        import gwsFeatureCollection
        self.rotation_model = gwsFeatureCollection.FeatureCollection(model=model, layer='rotations', url=url)
        self.static_polygons = gwsFeatureCollection.FeatureCollection(model=model, layer='static_polygons', url=url)
        self.dynamic_polygons = gwsFeatureCollection.FeatureCollection(model=model, layer='plate_polygons', url=url)

    def copy(self, deep=False):
        """
        Make a copy of an existing reconstruction_model
        """
        if deep:
            return copy.deepcopy(self)
        else:
            return copy.copy(self)

    def plate_snapshot(self, reconstruction_time, anchor_plate_id=0):
        """
        Generate a snapshot of a topological reconstruction model
        Returns an object of the PlateSnapshot class, containing the resolved plate polygons
        """
        resolved_topologies = []
        resolved_topological_sections = []
        pygplates.resolve_topologies(self.dynamic_polygons,
                                     self.rotation_model,
                                     resolved_topologies,
                                     reconstruction_time,
                                     resolved_topological_sections,
                                     anchor_plate_id=anchor_plate_id)

        return PlateSnapshot(resolved_topologies,
                             resolved_topological_sections,
                             self.rotation_model,
                             reconstruction_time,
                             anchor_plate_id)

    def polygon_snapshot(self, polygon_type, reconstruction_time, anchor_plate_id=0):
        """
        Create a set of reconstructed polygons for a specific reconstruction time 
        """

        if polygon_type == 'coastlines':
            polygons_to_reconstruct = self.coastlines
        elif polygon_type == 'continents':
            polygons_to_reconstruct = self.continent_polygons
        elif polygon_type == 'static_polygons':
            polygons_to_reconstruct = self.static_polygons
        else:
            print('some error msg')
        
        reconstructed_polygons = []
        pygplates.reconstruct(polygons_to_reconstruct,
                              self.rotation_model,
                              reconstructed_polygons,
                              reconstruction_time, anchor_plate_id=anchor_plate_id)

        # TODO detect if reconstructed_polygons is empty, return None if so?

        return ReconstructedPolygonSnapshot(reconstructed_polygons,
                                            self.rotation_model,
                                            reconstruction_time,
                                            anchor_plate_id)

    def rotation_table(self, plate_id_list=None, asdataframe=False):
        """
        Return an object containing the rotation parameters for some or all of the 
        current rotation model, either as a pandas dataframe or a list of lists (one 
        list per finite rotation entry)

        :param plate_id_list: (list of ints, optional) list of plate_ids to include in the output. 
            Default is None, in which case the full rotation table is returned
        :param asdataframe: (bool, optional) whether to return a pandas data frame (default is False) 
        """
        rotation_features = utils.rotation.generate_rotation_feature(self.rotation_files)

        return utils.rotation.get_rotation_table(rotation_features, plate_id_list=plate_id_list, asdataframe=asdataframe)


    if pygplates.Version.get_imported_version() >= pygplates.Version(32):
        def construct_topological_model(self, anchor_plate_id=0,
                            default_resolve_topology_parameters=pygplates.ResolveTopologyParameters(enable_strain_rate_clamping=True)):

            # tell the object to generate a TopologicalModel object based on the already
            # assigned rotations and topologies
            self.topological_model = pygplates.TopologicalModel(
                self.dynamic_polygons,
                self.rotation_model,
                anchor_plate_id=anchor_plate_id,
                # Enable strain rate clamping to better control crustal stretching factors...
                default_resolve_topology_parameters=default_resolve_topology_parameters)


    def reconstruct(self, features, reconstruction_time, anchor_plate_id=0, topological=False):
        """
        Reconstruct feature collection or a geopandas dataframe using the 
        """

        print(isinstance(features, gpd.GeoDataFrame))

        if not topological:
            if isinstance(features, pygplates.FeatureCollection):

                # TODO assign plate ids if not available already (or option selected)

                reconstructed_features = []
                pygplates.reconstruct(features, self.rotation_model, 
                                      reconstructed_features, reconstruction_time, anchor_plate_id=anchor_plate_id)
                return reconstructed_features

            elif isinstance(features, gpd.GeoDataFrame):

                # select only the features that are valid at reconstruction time?
                # convert geometries to gpml (features?)
                # reconstruct
                # somehow map reconstructed features back to original attribute table

                temporary_file = tempfile.NamedTemporaryFile(delete=True, suffix='.geojson')
                temporary_file.close()

                temporary_file_r = tempfile.NamedTemporaryFile(delete=True, suffix='.geojson')
                temporary_file_r.close()

                features.to_file(temporary_file.name, driver='GeoJSON')  #GeoJSON

                pygplates.reconstruct(temporary_file.name, self.rotation_model, 
                                      temporary_file_r.name, reconstruction_time,
                                      anchor_plate_id=anchor_plate_id)

                reconstructed_gdf = gpd.read_file(temporary_file_r.name)

                # The reconstructed file will have various extra columns, of which the name
                # of the temporary file is definitelt not useful so we delete it
                # (checking for the unlikely event of the column 'FILE1' already existing)
                if not 'FILE1' in features.columns:
                    reconstructed_gdf.drop(columns=['FILE1'], inplace=True)

                os.unlink(temporary_file.name)

                return reconstructed_gdf


        else:
             
             #TODO perform a topological reconstruction

             return



    def reconstruct_vector():

        return



    def assign_plate_ids(self, features, polygons='static', copy_valid_times=False):
        """
        assign plate ids to a data set using polygons from the ReconstructionModel 
        """
        if polygons=='continents':
            partitioning_polygon_features = self.continent_polygons
        elif polygons=='coastlines':
            partitioning_polygon_features = self.coastlines
        else:
            partitioning_polygon_features = self.static_polygons
        if not partitioning_polygon_features:
            raise ValueError('No polygons found for partitioning')

        if isinstance(features, pygplates.FeatureCollection):
            if copy_valid_times:
                properties_to_copy = [pygplates.PartitionProperty.reconstruction_plate_id,
                                      pygplates.PartitionProperty.valid_time_period]
            else:
                properties_to_copy = [pygplates.PartitionProperty.reconstruction_plate_id]
            return pygplates.FeatureCollection(
                pygplates.partition_into_plates(partitioning_polygon_features,
                                                   self.rotation_model,
                                                   features,
                                                   properties_to_copy=properties_to_copy))

        elif isinstance(features, gpd.GeoDataFrame):

            # TODO handle cases where static polygons are spread across multiple feature collections
            polygon_gdf = utils.create_gpml.gpml2gdf(pygplates.FeatureCollection(partitioning_polygon_features[0]))
            # TODO handle the FROMAGE and TOAGE
            if copy_valid_times:
                polygon_gdf = polygon_gdf[['geometry', 'PLATEID1', 'FROMAGE', 'TOAGE']]
            else:
                polygon_gdf = polygon_gdf[['geometry', 'PLATEID1']]

            return gpd.overlay(features, polygon_gdf, how='intersection', keep_geom_type=False)





class ReconstructedPolygonSnapshot(object):
    """
    Class to contain a set of reconstructed polygon features
    """
    def __init__(self,
                 reconstructed_polygons,
                 rotation_model,
                 reconstruction_time,
                 anchor_plate_id):

        self.reconstructed_polygons = reconstructed_polygons
        self.rotation_model = rotation_model
        self.reconstruction_time = reconstruction_time
        self.anchor_plate_id = anchor_plate_id

    def _plot(self, fig, pen='black', color='wheat', **kwargs):
        """
        DEPRECATED
        plot the reconstructed polygons into a pygmt map figure

        :param fig: (pygmt.Figure) pygmt figure object to plot to
        :param kwargs: (optional) set of additonal keyword arguments to pass to the pygmt 'plot' command
        """
        for polygon in self.reconstructed_polygons:
            data = polygon.get_reconstructed_geometry().to_lat_lon_array()
            fig.plot(x=data[:,1],y=data[:,0], 
                     pen=pen, color=color, **kwargs)

    def plot(self, fig, pen='black', color='wheat', **kwargs):

        if not self.reconstructed_polygons:
            print('No polygons to plot')
            return

        # plotting is generally faster if saved to temporary file
        plot_file = tempfile.NamedTemporaryFile(delete=False, suffix='.xy')
        plot_file.close()

        features = []
        for feature in self.reconstructed_polygons:
            f = pygplates.Feature()
            f.set_geometry(feature.get_reconstructed_geometry())
            features.append(f)

        pygplates.FeatureCollection(features).write(plot_file.name)
        fig.plot(data = plot_file.name, pen=pen, color=color, **kwargs)

        os.unlink(plot_file.name)

    #TODO add a 'to_file' method? But that should be an option when creating the snapshot




## TODO
     # topology_checks - area, segment lengths... [should go in plate snapshot??]
     # continents --> continent area, masking between continents and oceans
     # age grid(s) associated with reconstruction model



class PlateSnapshot(object):
    """
    Class containing a snapshot of a topological plate reconstruction.
    An instance of this class is initialized from a ReconstructionModel instance.
    """

    def __init__(self,
                 resolved_topologies,
                 resolved_topological_sections,
                 rotation_model,
                 reconstruction_time,
                 anchor_plate_id):

        self.reconstruction_time = reconstruction_time
        self.anchor_plate = anchor_plate_id
        self.rotation_model = rotation_model
        self.resolved_topologies = resolved_topologies
        self.resolved_topological_sections = resolved_topological_sections
        self.plate_count = len(resolved_topologies)
        self.plate_ids = [resolved_topology.get_resolved_feature().get_reconstruction_plate_id() \
                          for resolved_topology in resolved_topologies]
        self.plate_areas = [resolved_topology.get_resolved_geometry().get_area() * pygplates.Earth.mean_radius_in_kms**2\
                            for resolved_topology in resolved_topologies]
        self.plate_perimeters = [resolved_topology.get_resolved_geometry().get_arc_length() * pygplates.Earth.mean_radius_in_kms\
                            for resolved_topology in resolved_topologies]
        self.plate_centroids = [resolved_topology.get_resolved_geometry().get_interior_centroid()\
                                for resolved_topology in resolved_topologies]

    def get_boundary_features(self, boundary_types=['subduction','midoceanridge','other']):
        """
        Get the boundary features from a topological reconstruction plot_snapshot.
        Optionally specify to only return boundaries of type 'SubductionZone', 'MidOceanRidge',
        or other boundaries not of these two types.
        """
        # return a list of boundary features, optionally matching a certain boundary type
        resolved_boundary_segments = []
        for resolved_topological_section in self.resolved_topological_sections:
            for shared_sub_segment in resolved_topological_section.get_shared_sub_segments():
                if shared_sub_segment.get_feature().get_feature_type() == pygplates.FeatureType.gpml_subduction_zone:
                    if 'subduction' in boundary_types:
                        resolved_boundary_segments.append(shared_sub_segment.get_resolved_feature())
                elif shared_sub_segment.get_feature().get_feature_type() == pygplates.FeatureType.gpml_mid_ocean_ridge:
                    if 'midoceanridge' in boundary_types:
                        resolved_boundary_segments.append(shared_sub_segment.get_resolved_feature())
                else:
                    if 'other' in boundary_types:
                        resolved_boundary_segments.append(shared_sub_segment.get_resolved_feature())

        return resolved_boundary_segments

    def proximity_to_boundaries(self,reconstructed_point_features,
                                boundary_types=['subduction','midoceanridge','other']):
        """
        Given a set of reconstructed point features, returns the distance of each
        point to the nearest boundary of the specified type(s)
        """

        boundary_features = self.get_boundary_features(boundary_types=boundary_types)

        (reconstructed_lon,
         reconstructed_lat,
         distances) = distance_between_reconstructed_points_and_features(reconstructed_point_features,
                                                                         boundary_features)

        return reconstructed_lon, reconstructed_lat, distances

    def velocity_field(self, velocity_domain_features=None, velocity_type='both', delta_time=1.):
        """
        From the plate snapshot, generates a set of velocity vectors, which are return_closest_index
        as a instance of the VelocityField class.
        """
        if velocity_domain_features is None:
            velocity_domain_features = [PointDistributionOnSphere(distribution_type='healpix',N=32).meshnode_feature]

        # All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
        all_domain_points = []
        all_velocities_en = []
        all_velocities_magaz = []
        plate_ids = []

        rotation_model = pygplates.RotationModel(self.rotation_model)

        # Partition our velocity domain features into our topological plate polygons at the current 'time'.
        plate_partitioner = pygplates.PlatePartitioner(self.resolved_topologies,
                                                       rotation_model)

        for velocity_domain_feature in velocity_domain_features:

            # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
            # Iterate over them all.
            for velocity_domain_geometry in velocity_domain_feature.get_geometries():

                for velocity_domain_point in velocity_domain_geometry.get_points():

                    all_domain_points.append(velocity_domain_point)

                    partitioning_plate = plate_partitioner.partition_point(velocity_domain_point)
                    if partitioning_plate:

                        # We need the newly assigned plate ID to get the equivalent stage rotation of that tectonic plate.
                        partitioning_plate_id = partitioning_plate.get_feature().get_reconstruction_plate_id()

                        # Get the stage rotation of partitioning plate from 'time + delta_time' to 'time'.
                        equivalent_stage_rotation = rotation_model.get_rotation(self.reconstruction_time,
                                                                                partitioning_plate_id,
                                                                                self.reconstruction_time + delta_time)

                        # Calculate velocity at the velocity domain point.
                        # This is from 'time + delta_time' to 'time' on the partitioning plate.
                        velocity_vectors = pygplates.calculate_velocities(
                            [velocity_domain_point],
                            equivalent_stage_rotation,
                            delta_time)

                        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
                        velocities = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                                [velocity_domain_point],
                                velocity_vectors)
                        all_velocities_en.append(velocities[0])

                        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
                        velocities = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                                [velocity_domain_point],
                                velocity_vectors)
                        all_velocities_magaz.append(velocities[0])

                        plate_ids.append(partitioning_plate_id)

                    else:
                        # If point is not within a polygon, set velocity and plate_id to zero
                        all_velocities_en.append((0,0,0))
                        all_velocities_magaz.append((0,0,0))
                        plate_ids.append(0)

        vel_east=[]
        vel_north=[]
        vel_mag = []
        vel_azim = []
        for velocity_vector in all_velocities_en:
            if getattr(velocity_vector,'get_x',None) is not None:
                vel_east.append(velocity_vector.get_y())
                vel_north.append(velocity_vector.get_x())
            else:
                vel_east.append(0.)
                vel_north.append(0.)

        for velocity_vector in all_velocities_magaz:
            vel_mag.append(velocity_vector[0])
            vel_azim.append(velocity_vector[1])

        pt_lon = []
        pt_lat = []
        for pt in all_domain_points:
            pt_lon.append(pt.to_lat_lon()[1])
            pt_lat.append(pt.to_lat_lon()[0])

        return VelocityField(pt_lat,pt_lon,vel_east,vel_north,vel_mag,vel_azim,plate_ids)

    def plot(self, ax=None, projection=None,
             linewidth=2):
        """
        Plot topological plate boundaries into a cartopy map.
        Optionally specify the projection (default = Mollweide) and linewidth (default = 2).
        """
        #ccrs.Mollweide()
        import cartopy
        import cartopy.crs as ccrs
        if not projection:
            projection=ccrs.Mollweide()
        if not ax:
            fig = plt.figure(figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1, projection=projection)
            ax.set_global()

        for boundary in self.resolved_topological_sections:
            shared_sub_segments = boundary.get_shared_sub_segments()
            for shared_sub_segment in shared_sub_segments:
                if shared_sub_segment.get_resolved_feature().get_geometry() is not None:
                    if shared_sub_segment.get_resolved_feature().get_feature_type() == pygplates.FeatureType.gpml_subduction_zone:
                        boundary_color = 'black'
                    elif shared_sub_segment.get_resolved_feature().get_feature_type() == pygplates.FeatureType.gpml_mid_ocean_ridge:
                        boundary_color = 'red'
                    else:
                        boundary_color = 'lightblue'
                    ax.plot(shared_sub_segment.get_resolved_feature().get_geometry().to_lat_lon_array()[:,1],
                            shared_sub_segment.get_resolved_feature().get_geometry().to_lat_lon_array()[:,0],
                            transform=ccrs.Geodetic(), color=boundary_color, linewidth=linewidth)

        return ax

    # pygmt functions
    def plot_subduction_zones(self, fig, color='black', gap=10, size=4, **kwargs):
        """
        plot subduction zones into a pygmt map figure

        :param fig: (pygmt.Figure) pygmt figure object to plot to
        :param color: (string, optional) color of the subduction lines and triangle (default is 'black') 
        :param gap: (float, optional) gap between trianges (default is 10p)
        :param size: (float, optional) size of triangles (default is 4p)
        :param kwargs: (optional) additional arguments to pygmt plot command
        """

        resolved_boundary_segments = self.get_boundary_features(['subduction'])

        plot_file = tempfile.NamedTemporaryFile(delete=False, suffix='.xy')
        plot_file.close()

        features = []
        for resolved_boundary_segment in resolved_boundary_segments:
            if resolved_boundary_segment.get_geometry() is not None:
                if resolved_boundary_segment.get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)=='Left':
                    resolved_boundary_segment.set_geometry(pygplates.PolylineOnSphere(resolved_boundary_segment.get_geometry().to_lat_lon_list()[::-1]))
                features.append(resolved_boundary_segment)
        
        if not features:
            print('No subduction zones to plot')
            return

        pygplates.FeatureCollection(features).write(plot_file.name)
        fig.plot(data = plot_file.name, color=color, style='f{:f}p/{:f}p+r+t'.format(float(gap), float(size)), **kwargs)

        os.unlink(plot_file.name)


    def plot_mid_ocean_ridges(self, fig, pen='0.5p,red', **kwargs):
        """
        plot mid ocean ridges into a pygmt map figure

        :param fig: (pygmt.Figure) pygmt figure object to plot to
        :param color: (string, optional) color of the subduction lines and triangle (default is 'red') 
        :param kwargs: (optional) additional arguments to pygmt plot command
        """
        resolved_boundary_segments = self.get_boundary_features(['midoceanridge'])

        plot_file = tempfile.NamedTemporaryFile(delete=False, suffix='.xy')
        plot_file.close()

        features = []
        for resolved_boundary_segment in resolved_boundary_segments:
            if resolved_boundary_segment.get_geometry() is not None:
                features.append(resolved_boundary_segment)

        if not features:
            print('No mid-ocean ridges to plot')
            return

        pygplates.FeatureCollection(features).write(plot_file.name)
        fig.plot(data = plot_file.name, pen=pen, **kwargs)

        os.unlink(plot_file.name)


    def plot_other_boundaries(self, fig, pen='0.5p,gray70', **kwargs):
        """
        plot plate boundaries of 'other' types (ie not subduction zone 
        or mid ocean ridge) into a pygmt map figure

        :param fig: (pygmt.Figure) pygmt figure object to plot to
        :param color: (string, optional) color of the subduction lines and triangle (default is 'gray70') 
        :param kwargs: (optional) additional arguments to pygmt plot command
        """
        resolved_boundary_segments = self.get_boundary_features(['other'])

        plot_file = tempfile.NamedTemporaryFile(delete=False, suffix='.xy')
        plot_file.close()

        features = []
        for resolved_boundary_segment in resolved_boundary_segments:
            if resolved_boundary_segment.get_geometry() is not None:
                features.append(resolved_boundary_segment)
        
        if not features:
            print('No plate boundaries to plot')
            return

        pygplates.FeatureCollection(features).write(plot_file.name)
        fig.plot(data = plot_file.name, pen=pen, **kwargs)

        os.unlink(plot_file.name)

    #TODO plot polygons
    def plot_deformation_zones(self, fig, pen='0.5p,gray10', color='p14+b+r300', **kwargs):

        plot_file = tempfile.NamedTemporaryFile(delete=False, suffix='.xy')
        plot_file.close()

        features = []
        for topology in self.resolved_topologies:
            if isinstance(topology, pygplates.ResolvedTopologicalNetwork):
                features.append(topology.get_resolved_feature())

        if not features:
            print('No deformation zones to plot')
            return

        pygplates.FeatureCollection(features).write(plot_file.name)
        fig.plot(data = plot_file.name, pen=pen, color=color, **kwargs)

        os.unlink(plot_file.name)


class MotionPathFeature:
    """
    Class to define a motion path feature.
    """

    def __init__(self, path_times=np.arange(0.,201.,10.), 
                 reconstruction_plate_id=0, 
                 seed_points=None, lats=None, longs=None,
                 relative_plate_id=0):
        """
        create a motion path feature 

        :param seed_points: (tuple or list of tuples) lat,lon coordinates of the seed point(s)
        :param lats:
        :param longs:
        :param path_times: (array)
        :param reconstruction_plate_id: (int)
        :param relative_plate_id: (int, optional)
        """
        if seed_points:
            if type(seed_points) is tuple:
                seed_points = [seed_points]
        elif lats and longs:
            seed_points = []
            for x,y in zip(lats,longs):
                seed_points.append((x,y))
        else:
            raise ValueError('Unrecognised format for seed point coordinates')


        seed_points_at_digitisation_time = pygplates.MultiPointOnSphere(seed_points)
        motion_path_feature = pygplates.Feature.create_motion_path(seed_points_at_digitisation_time,
                                                                   path_times,
                                                                   valid_time=(pygplates.GeoTimeInstant.create_distant_past(), pygplates.GeoTimeInstant.create_distant_future()),
                                                                   relative_plate = relative_plate_id,
                                                                   reconstruction_plate_id = reconstruction_plate_id)

        self.seed_points = seed_points
        self.path_times = path_times
        self.motion_path_feature = motion_path_feature

    def reconstruct_motion_path(self, reconstruction_model, reconstruction_time=0, anchor_plate_id=0):
        """
        generate reconstructed trails from the motion path feature according to a specified reconstruction model
        """

        reconstructed_motion_paths = []
        pygplates.reconstruct(self.motion_path_feature, reconstruction_model.rotation_model,
                              reconstructed_motion_paths, reconstruction_time,
                              anchor_plate_id=anchor_plate_id,
                              reconstruct_type=pygplates.ReconstructType.motion_path)

        trails = []
        for reconstructed_motion_path in reconstructed_motion_paths:
            trails.append(reconstructed_motion_path.get_motion_path().to_lat_lon_array())

        return trails

    def rate(self, reconstruction_model, reconstruction_time=0):

        reconstructed_motion_paths = []
        pygplates.reconstruct(self.motion_path_feature, reconstruction_model.rotation_model,
                              reconstructed_motion_paths, reconstruction_time,
                              reconstruct_type=pygplates.ReconstructType.motion_path)

        #dists = []
        rates = []
        for reconstructed_motion_path in reconstructed_motion_paths:
            dist = []
            for segment in reconstructed_motion_path.get_motion_path().get_segments():
                dist.append(segment.get_arc_length()*pygplates.Earth.mean_radius_in_kms)
            #dists.append(dist[::-1])
            
            # Get rate of motion as distance per Myr
            rates.append(np.asarray(dist[::-1])/np.diff(self.path_times))

        return rates

    def step_plot(self, reconstruction_model, reconstruction_time=0, show=False):

        rates = self.rate(reconstruction_model, reconstruction_time=0)
        
        step_rates = []
        for rate in rates:
            step_rate = np.zeros(len(rate)*2)
            step_rate[::2] = rate
            step_rate[1::2] = rate
            step_rates.append(step_rate)

        step_time = np.zeros(len(rates[0])*2)
        step_time[::2] = self.path_times[:-1]
        step_time[1::2] = self.path_times[1:]
        
        if show:
            fig = plt.figure(figsize=(10,4))
            plt.plot(step_time,np.array(step_rates).T)
            plt.xlabel('Reconstruction Time (Myr)')
            plt.ylabel('Rate of motion (mm/yr)')
            plt.gca().invert_xaxis()
            plt.show()
        else:
            return np.array(step_time), np.array(step_rates).squeeze()


class FlowlineFeature:
    
    def __init__(self, path_times=np.arange(0.,151.,5.), 
                 seed_points=None, lats=None, longs=None,
                 left_plate=None, right_plate=None):
        """
        Create a reconstructable flowline feature

        :param seed_points: (tuple or list of tuples) lat,lon coordinates of the seed point(s)
        :param lats: 
        :param longs:
        :param path_times: (array)
        :param reconstruction_plate_id: (int)
        :param relative_plate_id: (int, optional)
        """
        
        if seed_points:
            if type(seed_points) is tuple:
                seed_points = [seed_points]
        elif lats and longs:
            seed_points = []
            for x,y in zip(lats,longs):
                seed_points.append((x,y))
        else:
            raise ValueError('Unrecognised format for seed point coordinates')

        # CREATE FLOWLINE
        # POINTS ON THE FLOWLINE
        multi_point = pygplates.MultiPointOnSphere(seed_points)

        #reverse_reconstruct=(rotation_model, 0, 1)

        flowline_feature = pygplates.Feature(pygplates.FeatureType.create_gpml('Flowline'))
        flowline_feature.set_geometry(multi_point)
        flowline_feature.set_times(path_times)
        flowline_feature.set_valid_time(np.max(path_times), np.min(path_times))
        flowline_feature.set_left_plate(left_plate)
        flowline_feature.set_right_plate(right_plate)
        #flowline_feature.set_geometry(multi_point, reverse_reconstruct=(rotation_model,0))
        
        self.seed_point = seed_points
        self.path_times = path_times
        self.flowline_feature = flowline_feature
        
        
    def reconstruct_flowline(self, reconstruction_model, reconstruction_time=0, anchor_plate_id=0):
        """
        reconstruct the flowline and return a list of reconstructed flowline geometries

        :param reconstruction_model: object of type ReconstructionModel containing rotation model
        :param reconstruction_time: time at which the reconstructed flowline is defined (typically 
                                    0 for overlaying on present day bathymetry, but can be any age value)
                                    [default is 0]
        :param anchor_plate_id: plate_id of fixed plate, relevant if reconstruction_time is not 0
                                [default is 0]
        """

        reconstructed_flowlines = []
        pygplates.reconstruct(self.flowline_feature, reconstruction_model.rotation_model, 
                              reconstructed_flowlines, reconstruction_time,
                              anchor_plate_id=anchor_plate_id, 
                              reconstruct_type=pygplates.ReconstructType.flowline)

        flowlines = []
        for reconstructed_flowline in reconstructed_flowlines:
            flowlines.append(reconstructed_flowline.get_left_flowline().to_lat_lon_array())
            flowlines.append(reconstructed_flowline.get_right_flowline().to_lat_lon_array())

        return flowlines

    def rate(self, reconstruction_model, reconstruction_time=0):
        """
        calculate half=spreading rate along flowline assumung symmetric spreading

        :param reconstruction_model: object of type ReconstructionModel containing rotation model
        :param reconstruction_time: time at which the reconstructed flowline is defined (typically 
                                    0 for overlaying on present day bathymetry, but can be any age value)
                                    [default is 0]
        """

        reconstructed_flowlines = []
        pygplates.reconstruct(self.flowline_feature, reconstruction_model.rotation_model, 
                              reconstructed_flowlines, reconstruction_time,
                              anchor_plate_id=0, 
                              reconstruct_type=pygplates.ReconstructType.flowline)

        rates = []
        for reconstructed_flowline in reconstructed_flowlines:
            # Iterate over the left flowline points
            flowlinearray_left = np.empty([0,0])
            for left_point in reconstructed_flowline.get_left_flowline():
                flowlinearray_left = np.append(flowlinearray_left, left_point.to_lat_lon_array())
            # Iterate over the right flowline points
            flowlinearray_right = np.empty([0,0])
            for right_point in reconstructed_flowline.get_right_flowline():
                flowlinearray_right = np.append(flowlinearray_right, right_point.to_lat_lon_array())
            
            tmp = reconstructed_flowline.get_left_flowline()
            dist = []
            for segment in tmp.get_segments():
                dist.append(segment.get_arc_length()*pygplates.Earth.mean_radius_in_kms)

            rate = 2*np.asarray(dist)/np.diff(self.path_times)  # *2 for full spreading rate
            
            rates.append(rate)
            
        return rates
            
    def step_plot(self, reconstruction_model, reconstruction_time=0, show=False):
        
        rates = self.rate(reconstruction_model, reconstruction_time=0)
        
        step_rates = []
        for rate in rates:
            step_rate = np.zeros(len(rate)*2)
            step_rate[::2] = rate
            step_rate[1::2] = rate
            step_rates.append(step_rate)

        step_time = np.zeros(len(rate)*2)
        step_time[::2] = self.path_times[:-1]
        step_time[1::2] = self.path_times[1:]

        if show:
            fig = plt.figure(figsize=(10,4))
            plt.plot(step_time,np.array(step_rates).T)
            plt.xlabel('Reconstruction Time (Myr)')
            plt.ylabel('Full Spreading Rate (mm/yr)')   ## IS this 
            plt.gca().invert_xaxis()
            plt.show()
        else:
            return np.array(step_time), np.array(step_rates).squeeze()



class PlateTree(object):
    """
    Class describing a geographic representation of a plate hierarchy,
    where the hierarchy joins the centroid points of plate polygons.
    """
    #TODO handle dynamic polygons as well as static

    def __init__(self, reconstruction_model):
        self.reconstruction_model = reconstruction_model

    def plot_snapshot(self, reconstruction_time, figsize=(14,9), show=True):
        """
        simple snapshot of a platetree snapshot at a specified reconstruction time, rendered in matplotlib
        """

        utils.platetree.plot_snapshot(self.reconstruction_model.static_polygons,
                                      self.reconstruction_model.rotation_model,
                                      reconstruction_time, 
                                      figsize=figsize,
                                      show=show)

    def to_gpml(self, reconstruction_times, filename):
        """
        Save a platetree object to a vector file with a GPlates-compatible file type
        """

        if isinstance(reconstruction_times, (float,int)):
            reconstruction_times = [reconstruction_times]

        for reconstruction_time in reconstruction_times:

            (uniq_plates_from_polygons,
             chains,
             reconstruction_tree,
             reconstructed_polygons) = utils.platetree.tree_snapshot(self.reconstruction_model.static_polygons,
                                                                     self.reconstruction_model.rotation_model,
                                                                     reconstruction_time)

            # for now, the valid time is set to be plus/minus 0.5 Myr
            tree_features = utils.platetree.create_hierarchy_features(chains,reconstructed_polygons,tree_features,
                                                                      valid_time=(reconstruction_time+0.5,reconstruction_time-0.5))

        pygplates.FeatureCollection(tree_features).write(filename)



class VelocityField(object):
    """
    Class containing velocity field vectors for a reconstruction snapshot.
    """
    def __init__(self, pt_lat,pt_lon,vel_east,vel_north,vel_mag,vel_azim,plate_ids):
        self.longitude = pt_lon
        self.latitude = pt_lat
        self.plate_id = plate_ids
        self.velocity_east = vel_east
        self.velocity_north = vel_north
        self.velocity_magnitude = vel_mag
        self.velocity_azimuth = vel_azim


    def rms_velocity(self, plate_id_selection=None):
        """
        compute the rms velocity for a specific plate, or list of plates, or all plates within
        the current plate snapshot
        """

        if plate_id_selection is None:
            return np.sqrt(np.mean(np.square(np.asarray(self.velocity_magnitude))))

        elif type(plate_id_selection) is list:
            index = np.array([plate_id in plate_id_selection for plate_id in self.plate_id])

        else:
            index = np.array([plate_id == plate_id_selection for plate_id in self.plate_id])

        #print index
        return np.sqrt(np.mean(np.square(np.asarray(self.velocity_magnitude)[index])))


    def plot(self, fig, scaling=500., style="v0.1c+e", pen="0.1p,black", color="black", **kwargs):
        """
        Plot velocity vectors to a pygmt figure
        """

        tmp = np.vstack((self.longitude, self.latitude,
                        np.degrees(self.velocity_azimuth),
                        np.array(self.velocity_magnitude)/scaling)).T

        fig.plot(data=tmp, style=style, pen=pen, color=color, **kwargs)




class SubductionConvergence(object):
    """
    Class for holding the analysis of subduction zone kinematics for a single time
    or range of time, contained within a pandas data frame
    """
    def __init__(self, reconstruction_model,
                 reconstruction_times,
                 threshold_sampling_distance_radians,
                 velocity_delta_time=1,
                 anchor_plate_id=0):
        """
        Initiate a subduction convergence object, which is essentially a pandas dataframe
        with the results of subduction kinematic analysis for a single or multiple time snapshots

        :param reconstruction_model: The ReconstructionModel object containing the rotation model and dynamic polygons
        :param reconstruction_times: (int, or float, or list or array of ints or floats)
        :param threshold_sampling_distance_radians:
        :param velocity_delta_time: (default=1)
        :param anchor_plate_id: (default=0)
        """

        # Data frame template defining the column names
        DataFrameTemplate = ('lon','lat','conv_rate','conv_obliq','migr_rate',
                             'migr_obliq','arc_length','arc_azimuth',
                             'subducting_plate','overriding_plate','time')

        # Create an empty dataframe to concatenate results to
        df_AllTimes = pd.DataFrame(columns=DataFrameTemplate)

        if isinstance(reconstruction_times, (float,int)):
            reconstruction_times = [reconstruction_times]

        for reconstruction_time in reconstruction_times:

            result = sc.subduction_convergence(
                reconstruction_model.rotation_model,
                reconstruction_model.dynamic_polygons,
                threshold_sampling_distance_radians,
                reconstruction_time,
                velocity_delta_time,
                anchor_plate_id)

            # Make a flat list of subduction stats to input into the proximity test
            subduction_data = []
            for data in result:
                subduction_data.append(data+(reconstruction_time,))

            df = pd.DataFrame(subduction_data, columns = DataFrameTemplate)

            # append dataframe
            df_AllTimes = df_AllTimes.append(df)

        # convert list array to dataframe
        self.df = df_AllTimes
        self.reconstruction_model = reconstruction_model


    def plot(self, variable='convergence rate'):
        """
        simple plot of subduction kinematics with matplotlib
        """
        plt.figure(figsize=(12,5))
        if variable in ['convergence rate','cr']:
            self.df.plot(kind='scatter', x='lon', y='lat', c='conv_rate', 
                         colormap='magma', edgecolor=None)
            plt.title('convergence rate')
        if variable in ['convergence obliquity','co']:
            self.df.plot(kind='scatter', x='lon', y='lat', c='conv_obliq', 
                         colormap='magma', edgecolor=None)
            plt.title('migration rate')
        if variable in ['migration rate','mr']:
            self.df.plot(kind='scatter', x='lon', y='lat', c='migr_rate', 
                         colormap='magma', edgecolor=None)
            plt.title('migration rate')
        if variable in ['migration obliquity','mo']:
            self.df.plot(kind='scatter', x='lon', y='lat', c='migr_obliq', 
                         colormap='magma', edgecolor=None)
            plt.title('migration rate')
        plt.show()

    # TODO make the histograms weighted by segment length
    def hist(self, variable='convergence rate',bins=50):
        """
        plot histogram of subduction parameters
        """

        plt.figure(figsize=(12,5))
        if variable in ['convergence rate','cr']:
            plt.hist(self.df.conv_rate,bins=bins)
            plt.title('convergence rate')
        if variable in ['convergence obliquity','co']:
            plt.hist(self.df.conv_obliq,bins=bins)
            plt.title('migration rate')
        if variable in ['migration rate','mr']:
            plt.hist(self.df.migr_rate,bins=bins)
            plt.title('migration rate')
        if variable in ['migration obliquity','mo']:
            plt.hist(self.df.migr_obliq,bins=bins)
            plt.title('migration rate')
        plt.show()



class AgeCodedPointDataset(object):

    def __init__(self, source, field_mapping = None):
        """
        Initiate an AgeCodedPointDataset class

        This can be:
        1. Any feature collection file readable by GPlates
        2. Any csv (TODO add support for any pandas-readable file)
        3. The paleobiology database web service (TODO add more generic support for web services)
        """

        try:
            filename, file_extension = os.path.splitext(source)
        except:
            filename = None; file_extension = None

        if file_extension in ['.shp','.gpml','.gpmlz','.gmt']:
            feature_collection = pygplates.FeatureCollection(source)

            self._point_features = feature_collection

            DataFrameTemplate = ['lon','lat','name','description','reconstruction_plate_id','from_age','to_age']

            # Get attribute (other than coordinate) names from first feature
            for feature in feature_collection:
                if feature.get_shapefile_attributes():
                    for attribute in feature.get_shapefile_attributes():
                        DataFrameTemplate.append(attribute)
                break

            result = []
            for feature in feature_collection:
                tmp = []
                tmp.append(feature.get_geometry().to_lat_lon()[1])
                tmp.append(feature.get_geometry().to_lat_lon()[0])
                tmp.append(feature.get_name())
                tmp.append(feature.get_description())
                tmp.append(feature.get_reconstruction_plate_id())
                tmp.append(feature.get_valid_time()[0])
                tmp.append(feature.get_valid_time()[1])
                if feature.get_shapefile_attributes():
                    for attribute in feature.get_shapefile_attributes():
                        tmp.append(feature.get_shapefile_attribute(attribute))
                result.append(tmp)

            self._df = pd.DataFrame(result,columns=DataFrameTemplate)
            self._field_mapping = {'latitude_field':'lat', 'longitude_field':'lon',
                                   'max_age_field':'from_age', 'min_age_field':'to_age'}

        else:
            if file_extension == '.csv':
                self._df = pd.read_csv(source)
            elif "http://" in source or "https://" in source:
                import requests
                r = requests.get(source)
                self._df = pd.read_csv(StringIO(r.text))
                field_mapping = {'latitude_field':'lat', 'longitude_field':'lng',
                                 'max_age_field':'max_ma', 'min_age_field':'min_ma'}
            elif isinstance(source,pd.DataFrame):
                self._df = source

            self._field_mapping = field_mapping

            self._point_features = []
            for index,row in self._df.iterrows():
                point = pygplates.PointOnSphere(float(row[field_mapping['latitude_field']]),
                                                float(row[field_mapping['longitude_field']]))
                point_feature = pygplates.Feature()
                point_feature.set_geometry(point)
                point_feature.set_reconstruction_plate_id(0)
                try:
                    point_feature.set_valid_time(row[field_mapping['max_age_field']],-999.)
                except:
                    warnings.warn('Unable to set valid time for row %d' % index)
                    point_feature.set_valid_time(-998,-999.)
                self._point_features.append(point_feature)


    def assign_reconstruction_model(self, reconstruction_model, polygons='static'):
        """
        assign plate ids to a point data set using an existing ReconstructionModel class
        """
        if polygons=='continents':
            partitioning_polygon_features = reconstruction_model.continent_polygons
        elif polygons=='coastlines':
            partitioning_polygon_features = reconstruction_model.coastlines
        else:
            partitioning_polygon_features = reconstruction_model.static_polygons

        if not partitioning_polygon_features:
            raise ValueError('No polygons found for partitioning')
        partitioned_point_features = pygplates.partition_into_plates(partitioning_polygon_features,
                                                                     reconstruction_model.rotation_model,
                                                                     self._point_features)
        self._point_features = partitioned_point_features
        self.reconstruction_model = reconstruction_model


    def reconstruct(self,reconstruction_time,anchor_plate_id=0):
        """
        reconstruct point data to specified time (and optionally with specified anchor
        plate id)
        """

        reconstructed_features = []
        pygplates.reconstruct(self._point_features,
                              self.reconstruction_model.rotation_model,
                              reconstructed_features,
                              reconstruction_time,
                              anchor_plate_id=anchor_plate_id)

        return reconstructed_features


    def plot_reconstructed(self,reconstruction_time,anchor_plate_id=0):
        """
        Quick plot of points reconstructed to specified time (and optionally with
        specified anchor plate id)
        """

        reconstructed_features = []
        pygplates.reconstruct(self._point_features,
                              self.reconstruction_model.rotation_model,
                              reconstructed_features,
                              reconstruction_time,
                              anchor_plate_id=anchor_plate_id)

        plt.figure()
        for reconstructed_feature in reconstructed_features:
            plt.plot(reconstructed_feature.get_reconstructed_geometry().to_lat_lon()[1],
                     reconstructed_feature.get_reconstructed_geometry().to_lat_lon()[0],'ro')
            plt.axis([-180,180,-90,90])
        plt.title('%0.2f Ma' % reconstruction_time)
        plt.show()


    def reconstruct_to_time_of_appearance(self,ReconstructTime='BirthTime',anchor_plate_id=0):
        """
        Reconstruct points to time of appearance corresponding to each point feature
        """

        rotation_model = pygplates.RotationModel(self.reconstruction_model.rotation_model)
        recon_points = []
        for point_feature in self._point_features:
            if ReconstructTime == 'MidTime':
                reconstruction_time = (point_feature.get_valid_time()[0]+point_feature.get_valid_time()[1])/2.
            else:
                reconstruction_time = point_feature.get_valid_time()[0]
            if point_feature.get_reconstruction_plate_id()!=0:
                point_rotation = rotation_model.get_rotation(reconstruction_time,
                                                             point_feature.get_reconstruction_plate_id(),
                                                             anchor_plate_id=anchor_plate_id)
                reconstructed_point = point_rotation * point_feature.get_geometry()
                recon_points.append([reconstructed_point.to_lat_lon()[1],
                                     reconstructed_point.to_lat_lon()[0],
                                     reconstruction_time])

        return recon_points


    def spatial_binning(self, reconstruction_time=None, anchor_plate_id=0, binsize=10., axis=None):
        """
        spatial binning within regular long-lat boxes,
        [cf Zeigler++ 2003 Lethaia; Cao++ 2018 Geol.Mag]

        """

        bin_edges=(np.arange(-180,180+binsize,binsize),
                   np.arange(-90,90+binsize,binsize))

        if reconstruction_time is None:
            result = np.histogram2d(self._df[self._field_mapping['longitude_field']],
                                    self._df[self._field_mapping['latitude_field']],
                                    bin_edges)
        else:
            reconstructed_features = []
            pygplates.reconstruct(self._point_features,
                                  self.reconstruction_model.rotation_model,
                                  reconstructed_features,
                                  reconstruction_time,
                                  anchor_plate_id=anchor_plate_id)

            result = np.histogram2d([feature.get_reconstructed_geometry().to_lat_lon()[1] for feature in reconstructed_features],
                                    [feature.get_reconstructed_geometry().to_lat_lon()[0] for feature in reconstructed_features],
                                    bin_edges)


        if axis == None:
            return result,bin_edges
        elif axis in ['latitude',0]:
            return np.nansum(result[0]/result[0],axis=0),bin_edges[1]
        elif axis in ['longitude',1]:
            return np.nansum(result[0]/result[0],axis=1),bin_edges[0]


class PointDistributionOnSphere(object):

    """
    Class to handle point distributions on the sphere
    """
    def __init__(self, distribution_type='random', N=10000, nest=False):
        """
        Initiate a point distribution on the sphere

        distribution_type: 'healpix' or 'random' [default='random']

        N: number controlling point density. If 'distribution_type' is 'healpix',
            N must be a factor of 2 and controls the healpix density. If
            'distribution_type' is 'random', N is the number of points returned.

        nest: logical, only applies for healpix distributions, and controls
              the ordering of the healpix points [default='False']

        """

        if distribution_type=='healpix':
            try:
                #import healpy as hp
                from astropy_healpix import healpy as hp
                othetas,ophis = hp.pix2ang(N,np.arange(12*N**2),nest=nest)
                othetas = np.pi/2-othetas
                ophis[ophis>np.pi] -= np.pi*2

                # ophis -> longitude, othetas -> latitude
                self.longitude = np.degrees(ophis)
                self.latitude = np.degrees(othetas)
            except:
                warnings.warn('unable to import module for healpix generation, trying pregenerated point files')
                features = pygplates.FeatureCollection('{:s}/healpix_mesh_{:d}.gpmlz'.format(DATA_DIR,N))
                for feature in features:
                    for geometry in feature.get_all_geometries():
                        self.latitude = geometry.to_lat_lon_array()[:,0]
                        self.longitude = geometry.to_lat_lon_array()[:,1]

        elif distribution_type in ['random','marsaglia','fibonacci','spiral']:
            # function to call Marsaglia's method and return Long/
            # Lat arrays

            Long, Lat = utils.sphere.points_on_sphere(N=N, distribution_type=distribution_type)

            self.longitude = np.array(Long)
            self.latitude = np.array(Lat)

        else:
            raise ValueError('unrecognized distribution type ''{:s}'' for PointDistributionOnSphere generation'.format(distribution_type))

        self.multipoint = pygplates.MultiPointOnSphere(zip(self.latitude,self.longitude))
        self.meshnode_feature = pygplates.Feature(pygplates.FeatureType.create_from_qualified_string('gpml:MeshNode'))
        self.meshnode_feature.set_geometry(self.multipoint)

    def to_file(self, filename):
        """
        Save to a file. File format is dictated from extension (options are gpml, gpmlz, gmt, shp)
        """

        pygplates.FeatureCollection(self.meshnode_feature).write(filename)

    def mask(self, reconstructable_polygons, rotation_model, reconstruction_time=0., 
             masking='outside', preserve_polygon_attributes=False):
        """
        create a masked version of the point distribution using reconstructed polygons 
        [NB currently returns multipoint feature(s), not a PointDistributionOnSphere object]
        """

        masked_pts = utils.spatial.rasterise_polygons(reconstructable_polygons, 
                                                      rotation_model, 
                                                      reconstruction_time, 
                                                      raster_domain_points=self.meshnode_feature,
                                                      masking=masking)

        if preserve_polygon_attributes:
            return masked_pts
        else:
            merge_points = []
            for points in masked_pts:
                merge_points.extend(points.get_geometry().to_lat_lon_list())
            return pygplates.MultiPointOnSphere(merge_points)


    #TODO - move this to be a method of AgeCodedPointDataset
    def point_feature_heatmap(self, target_features, return_indices=False):
        """
        Given a AgeCodedPointDataset class object, returns a heatmap showing the number
        of points for which each point in the point distribution is the closest.
        Most useful where the point distribution is equal area.
        """

        res = find_closest_geometries_to_points(target_features,
                                                [self.multipoint],
                                                return_closest_position = True,
                                                return_closest_index = True)

        bin_indices = list(zip(*res))[2]
        bin_counts = []
        for j in range(len(self.multipoint.get_points())):
            bin_counts.append(bin_indices.count(j))

        if return_indices:
            return bin_counts, bin_indices
        else:
            return bin_counts


# From raster_reconstruction_classes
class GPlatesRaster(object):
    """
    Class for holding raster data for conveient use with reconstructed data 
    """
    # TODO create a way to initialise from an array
    # TODO save to netcdf
    # TODO grdimage method

    def __init__(self, filename, reconstruction_time=0., z_field_name='z'):
        """
        initialise a GPlatesRaster object from a netcdf file defined in geographic coordinates
        """
        self.gridX,self.gridY,self.gridZ = utils.fileio.load_netcdf(filename, z_field_name)
        self.reconstruction_time = reconstruction_time
        self.source_filename = filename

    def plot(self, show=False, levels=20, extend='both', cmap=plt.cm.BrBG_r, **kwargs):
        """
        generate a quick plot of the raster using matplotlib 
        """
        plt.figure(figsize=(16,6))
        plt.contourf(self.gridX, self.gridY, self.gridZ,
                     levels, extend=extend, cmap=cmap, **kwargs)
        plt.gca().set_aspect('equal')
        plt.colorbar()
        if show:
            plt.show()

    def sample(self, point_lons, point_lats, order=0):
        """
        sample the raster at specififed lon,lat points and return the z values
        """
        LonGrid, LatGrid = np.meshgrid(self.gridX,self.gridY)
        d,l = utils.sphere.sampleOnSphere(LonGrid.flatten(),
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

    
    def reconstruct(self, reconstruction_model, to_time, from_time=0, grid_sampling=1.,
                    anchor_plate_id=0, sampling_method='scipy', return_points=False):
        """
        TODO make default sampling be same as input
        """

        from_time = 0.
        to_time = np.float(to_time)

        (reconstructed_point_lons,
         reconstructed_point_lats,
         point_zvals) = utils.raster.reconstruct_raster(self, 
                                                        reconstruction_model.static_polygons, 
                                                        reconstruction_model.rotation_model,
                                                        from_time, to_time, 
                                                        grid_sampling=grid_sampling)

        if return_points:
            return reconstructed_point_lons, reconstructed_point_lats, point_zvals
        else:
            return utils.raster.xyz2grd(reconstructed_point_lons,reconstructed_point_lats,
                                        point_zvals,self.gridX,self.gridY)

    #def cross_section(self, PtLons, PtLats):
    #    return CrossSection(self, PtLons, PtLats)


class CrossSection(object):

    def __init__(self, target_raster, PtLons, PtLats):

        self.GreatCirclePoints,self.ProfilePoints,arc_distance = utils.paleogeography.create_profile_points(PtLons,PtLats)
        # create an array of distances along profile in km, starting at zero

        self.profileX_kms = np.arange(0,self.ProfilePoints.shape[0])*arc_distance

        # extract the values from the (smoothed) topography grid along the profile
        self.grid_values = utils.paleogeography.create_slice(target_raster.gridX,
                                                             target_raster.gridY,
                                                             target_raster.gridZ,
                                                             self.GreatCirclePoints, self.ProfilePoints)

        self.cross_section_geometry = pygplates.PolylineOnSphere(self.GreatCirclePoints)
        self.source_filename = target_raster.source_filename

    def plate_boundary_intersections(self, shared_boundary_sections):

        (self.subduction_intersections,
         self.ridge_intersections,
         self.other_intersections) = utils.paleogeography.plate_boundary_intersections(self.cross_section_geometry,
                                                                                       shared_boundary_sections,
                                                                                       self.profileX_kms)




# reconstructable datasets
class litho1_scalar_coverage(object):

    def __init__(self, distribution_type='healpix', N=32):
        import litho1pt0 as litho
        self.points = PointDistributionOnSphere(distribution_type,N)
        self.litho = litho
        self.layer_keys = litho.l1_layer_decode.items()
        self.value_keys = litho.l1_data_decode.items()

    def write_layer_depth_to_scalar_coverage(self, filename, layer_names='All'):

        if layer_names == 'All':
            layer_names = [name[0] for name in self.layer_keys]

        scalar_coverages = {}
        for layer_name in layer_names:

            layerZ = litho.layer_depth(self.points.latitude, self.points.longitude, layer_name)
            scalar_coverages[pygplates.ScalarType.create_gpml(layer_name)] = layerZ

        ct_feature = pygplates.Feature()
        ct_feature.set_geometry((self.points.multipoint,scalar_coverages))
        ct_feature.set_name('litho1.0 layers')

        pygplates.FeatureCollection(ct_feature).write(filename)

    def write_layer_thickness_to_scalar_coverage(self, top_layer_name='CRUST1-TOP', bottom_layer_name='CRUST3-BOTTOM'):

        #if layer_names == 'All':
        #    layer_names = [name[0] for name in self.layer_keys]
        top_layer_depth = litho.layer_depth(self.points.latitude, self.points.longitude, top_layer_name)
        bottom_layer_depth = litho.layer_depth(self.points.latitude, self.points.longitude, bottom_layer_name)
        
        scalar_coverage = {}
        #layer_name = 'THICKNESS.{:s}--{:s}'.format(top_layer_name, bottom_layer_name)
        layer_name = 'CrustalThickness'
        layerZ = bottom_layer_depth-top_layer_depth
        scalar_coverage[pygplates.ScalarType.create_gpml(layer_name)] = layerZ

        ct_feature = pygplates.Feature()
        ct_feature.set_geometry((self.points.multipoint,scalar_coverage))
        ct_feature.set_name('litho1.0 layer thickness, {:s} to {:s}'.format(top_layer_name, bottom_layer_name))
        
        return ct_feature
