import pygplates
import numpy as np
import healpy as hp
import pandas as pd
import matplotlib.pyplot as plt
import os
import requests
from io import StringIO

from proximity_query import find_closest_geometries_to_points
import ptt.subduction_convergence as sc
import platetree as ptree

import gwsFeatureCollection


class ReconstructionModel(object):
    
    def __init__(self, name):
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

    def add_rotation_model(self, rotation_file):
        self.rotation_files.append(rotation_file)
        self.rotation_model = pygplates.RotationModel(self.rotation_files)
        
    def add_static_polygons(self, static_polygons_file):
        self.static_polygon_files.append(static_polygons_file)
        self.static_polygons.append(static_polygons_file)  # Should this be loaded into memory like dynamic polygons??
        
    def add_dynamic_polygons(self, dynamic_polygons_file):
        self.dynamic_polygon_files.append(dynamic_polygons_file)
        self.dynamic_polygons = [pygplates.FeatureCollection(dpfile) for dpfile in self.dynamic_polygon_files]

    def add_coastlines(self, coastlines_file):
        self.coastlines_files.append(coastlines_file)
        self.coastlines.append(coastlines_file)

    def add_continent_polygons(self, continent_polygons_file):
        self.continent_polygons_files.append(continent_polygons_file)
        self.continent_polygons.append(continent_polygons_file)

    def from_web_service(self, model='MULLER2016', url='https://gws.gplates.org'):
        self.rotation_model = gwsFeatureCollection.FeatureCollection(model=model, layer='rotations', url=url)
        self.static_polygons = gwsFeatureCollection.FeatureCollection(model=model, layer='static_polygons', url=url)
        self.dynamic_polygons = gwsFeatureCollection.FeatureCollection(model=model, layer='plate_polygons', url=url)
    
    def plate_snapshot(self, reconstruction_time, anchor_plate_id=0):
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


## TODO
     # topology_checks - area, segment lengths... [should go in plate snapshot??]
     # continents --> continent area, masking between continents and oceans
     # age grid(s) associated with reconstruction model



class PlateSnapshot(object):

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
        self.plate_areas = [resolved_topology.get_resolved_geometry().get_area() * pygplates.Earth.mean_radius_in_kms\
                            for resolved_topology in resolved_topologies]
        self.plate_perimeters = [resolved_topology.get_resolved_geometry().get_arc_length() * pygplates.Earth.mean_radius_in_kms\
                            for resolved_topology in resolved_topologies] 
        self.plate_centroids = [resolved_topology.get_resolved_geometry().get_interior_centroid()\
                                for resolved_topology in resolved_topologies]

    def velocity_field(self, velocity_domain_features=None, velocity_type='both', delta_time=1.):

        if velocity_domain_features is None:
            velocity_domain_features = PointDistributionOnSphere(distribution_type='healpix',N=32).meshnode_feature

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


class PlateTree(object):
    #TODO handle dynamic polygons as well as static

    def __init__(self, reconstruction_model):
        self.reconstruction_model = reconstruction_model

    def plot_snapshot(self, reconstruction_time):
        
        ptree.plot_snapshot(self.reconstruction_model.static_polygons,
                            self.reconstruction_model.rotation_model,
                            reconstruction_time)

    def to_gpml(self, reconstruction_times, filename):

        if isinstance(reconstruction_times, (float,int)):
            reconstruction_times = [reconstruction_times]

        for reconstruction_time in reconstruction_times:

            (uniq_plates_from_polygons, 
             chains, 
             reconstruction_tree,
             reconstructed_polygons) = tree_snapshot(polygons,
                                                     rotation_model,
                                                     reconstruction_time)

            # for now, the valid time is set to be plus/minus 0.5 Myr
            tree_features = pt.create_hierarchy_features(chains,reconstructed_polygons,tree_features,
                                                         valid_time=(recon_time+0.5,recon_time-0.5))

        pygplates.FeatureCollection(tree_features).write(filename)



class VelocityField(object):
    
    def __init__(self, pt_lat,pt_lon,vel_east,vel_north,vel_mag,vel_azim,plate_ids):
        self.longitude = pt_lon
        self.latitude = pt_lat
        self.plate_id = plate_ids
        self.velocity_east = vel_east
        self.velocity_north = vel_north
        self.velocity_magnitude = vel_mag
        self.velocity_azimuth = vel_azim

    def rms_velocity(self, plate_id_selection=None):

        if plate_id_selection is None:
            return np.sqrt(np.mean(np.square(np.asarray(self.velocity_magnitude))))

        elif type(plate_id_selection) is list:
            index = [plate_id in plate_id_selection for plate_id in self.plate_id]

        else:
            index = [plate_id == plate_id_selection for plate_id in self.plate_id]
        
        #print index
        return np.sqrt(np.mean(np.square(np.asarray(self.velocity_magnitude)[index])))



class SubductionConvergence(object):
    
    def __init__(self, reconstruction_model,
                 reconstruction_times,
                 threshold_sampling_distance_radians,
                 velocity_delta_time=1,
                 anchor_plate_id=0,
                 time_step=1):

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
        plt.figure(figsize=(12,5))
        if variable in ['convergence rate','cr']:
            plt.scatter(self.df.lon, self.df.lat,c=self.df.conv_rate,edgecolors='')
            plt.title('convergence rate')
        if variable in ['convergence obliquity','co']:
            plt.scatter(self.df.lon, self.df.lat,c=self.df.migr_obliq,edgecolors='')
            plt.title('migration rate')
        if variable in ['migration rate','mr']:
            plt.scatter(self.df.lon, self.df.lat,c=self.df.migr_rate,edgecolors='')
            plt.title('migration rate')
        if variable in ['migration obliquity','mo']:
            plt.scatter(self.df.lon, self.df.lat,c=self.df.migr_obliq,edgecolors='')
            plt.title('migration rate')
        plt.colorbar()
        plt.show()

    # TODO make the histograms weighted by segment length
    def hist(self, variable='convergence rate',bins=50):

        plt.figure(figsize=(12,5))
        if variable in ['convergence rate','cr']:
            plt.hist(self.df.conv_rate,bins=bins)
            plt.title('convergence rate')
        if variable in ['convergence obliquity','co']:
            plt.scatter(self.df.migr_obliq,bins=bins)
            plt.title('migration rate')
        if variable in ['migration rate','mr']:
            plt.scatter(self.df.migr_rate,bins=bins)
            plt.title('migration rate')
        if variable in ['migration obliquity','mo']:
            plt.scatter(self.df.migr_obliq,bins=bins)
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

        filename, file_extension = os.path.splitext(source)

        if file_extension in ['.shp','.gpml','.gpmlz','.gmt']:
            feature_collection = pygplates.FeatureCollection(source)

            self._point_features = feature_collection

            DataFrameTemplate = ['lon','lat','name','description','reconstruction_plate_id','from_age','to_age']

            # Get attribute (other than coordinate) names from first feature
            for feature in feature_collection:
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
                point_feature.set_valid_time(row[field_mapping['max_age_field']],-999.)
                self._point_features.append(point_feature)


    def assign_reconstruction_model(self,reconstruction_model):
        """
        assign plate ids to a point data set using an existing ReconstructionModel class
        """

        partitioned_point_features = pygplates.partition_into_plates(reconstruction_model.static_polygons,
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
            if ReconstructTime is 'MidTime':
                time = (point_feature.get_valid_time()[0]+point_feature.get_valid_time()[1])/2.
            else:
                time = point_feature.get_valid_time()[0]
            if point_feature.get_reconstruction_plate_id()!=0:
                point_rotation = rotation_model.get_rotation(time,
                                                             point_feature.get_reconstruction_plate_id(),
                                                             anchor_plate_id=anchor_plate_id)
                reconstructed_point = point_rotation * point_feature.get_geometry()
                recon_points.append([reconstructed_point.to_lat_lon()[1],
                                     reconstructed_point.to_lat_lon()[0],
                                     time])
            
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
    def __init__(self, distribution_type='random', N=10000):
        """
        Initiate a point distribution on the sphere
        
        distribution_type: 'healpix' or 'random' [default='random']
        
        N: number controlling point density. If 'distribution_type' is 'healpix',
            N must be a factor of 2 and controls the healpix density. If 
            'distribution_type' is 'random', N is the number of points returned.
        
        """

        if distribution_type=='healpix':               
            othetas,ophis = hp.pix2ang(N,np.arange(12*N**2))
            othetas = np.pi/2-othetas
            ophis[ophis>np.pi] -= np.pi*2

            # ophis -> longitude, othetas -> latitude
            self.longitude = np.degrees(ophis)
            self.latitude = np.degrees(othetas)

        elif distribution_type=='random':
            # function to call Marsaglia's method and return Long/
            # Lat arrays

            ## Marsaglia's method
            dim = 3
            norm = np.random.normal
            normal_deviates = norm(size=(dim, N))

            radius = np.sqrt((normal_deviates**2).sum(axis=0))
            points = normal_deviates/radius

            Long=[], Lat=[]
            for xyz in points.T:
                LL = pygplates.PointOnSphere((xyz))
                Lat.append(LL.to_lat_lon()[0])
                Long.append(LL.to_lat_lon()[1])

            self.longitude = np.array(Long)
            self.latitude = np.array(Lat)

        self.multipoint = pygplates.MultiPointOnSphere(zip(self.latitude,self.longitude))
        self.meshnode_feature = pygplates.Feature(pygplates.FeatureType.create_from_qualified_string('gpml:MeshNode'))
        self.meshnode_feature.set_geometry(self.multipoint)

    def to_file(self, filename):
        """
        Save to a file. File format is dictated from extension (options are gpml, gpmlz, gmt, shp)
        """

        pygplates.FeatureCollection(self.meshnode_feature).write(filename)

#TODO - move this to be a method of AgeCodedPointDataset
    def point_feature_heatmap(self, target_features):
        """
        Given a AgeCodedPointDataset class object, returns a heatmap showing the number
        of points for which each point in the point distribution is the closest. 
        Most useful where the point distribution is equal area.
        """

        res = find_closest_geometries_to_points(target_features,
                                                [self.multipoint],
                                                return_closest_position = True,
                                                return_closest_index = True)

        bin_indices = zip(*res)[2]
        bin_counts = []
        for j in range(len(self.multipoint.get_points())):
            bin_counts.append(bin_indices.count(j))

        return bin_counts


