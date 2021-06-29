import pygplates
import numpy as np
import ptt

from . import tools as _tools



class PlateReconstruction(object):

    def __init__(self, rotation_model=None, topology_features=None, static_polygons=None):


        rotation_model = pygplates.RotationModel(rotation_model)

        default_topology_features = pygplates.FeatureCollection()
        for topology in topology_features:
            default_topology_features.add( pygplates.FeatureCollection(topology) )

        self.rotation_model = rotation_model
        self.topology_features = default_topology_features
        self.static_polygons = static_polygons


    def tesselate_subduction_zones(self, time, tessellation_threshold_radians=0.001, anchor_plate_id=0):
        subduction_data = ptt.subduction_convergence.subduction_convergence(
            self.rotation_model,
            self.topology_features,
            tessellation_threshold_radians,
            float(time),
            anchor_plate_id=anchor_plate_id)
        subduction_data = np.vstack(subduction_data)
        return subduction_data


    def tesselate_mid_ocean_ridges(self, time, tessellation_threshold_radians=0.001, anchor_plate_id=0):
        ridge_data = ptt.ridge_spreading_rate.spreading_rates(
            self.rotation_model,
            self.topology_features,
            float(time),
            tessellation_threshold_radians,
            anchor_plate_id=anchor_plate_id)
        ridge_data = np.vstack(ridge_data)
        return ridge_data

    def reconstruct(self, feature, to_time, from_time=0, anchor_plate_id=0, **kwargs):
        from_time, to_time = float(from_time), float(to_time)
        if from_time != 0.0:
            raise NotImplementedError("Soon...")

        reconstructed_features = []
        pygplates.reconstruct(feature, self.rotation_model, reconstructed_features, to_time,\
            anchor_plate_id=anchor_plate_id, **kwargs)
        return reconstructed_features


    def get_point_velocities(self, lons, lats, time, delta_time=1.0):
        """ function to make a velocity mesh nodes at an arbitrary set of points defined in Lat
        Lon and Lat are assumed to be 1d arrays. """
        # Add points to a multipoint geometry

        reconstruction_time = float(reconstruction_time)

        multi_point = pygplates.MultiPointOnSphere([(float(lat),float(lon)) for lat, lon in zip(lats,lons)])

        # Create a feature containing the multipoint feature, and defined as MeshNode type
        meshnode_feature = pygplates.Feature(pygplates.FeatureType.create_from_qualified_string('gpml:MeshNode'))
        meshnode_feature.set_geometry(multi_point)
        meshnode_feature.set_name('Velocity Mesh Nodes from pygplates')

        velocity_domain_features = pygplates.FeatureCollection(meshnode_feature)
        
        # NB: at this point, the feature could be written to a file using
        # output_feature_collection.write('myfilename.gpmlz')
        
        
        # All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
        all_domain_points = []
        all_velocities = []

        # Partition our velocity domain features into our topological plate polygons at the current 'time'.
        plate_partitioner = pygplates.PlatePartitioner(topology_features, rotation_model, reconstruction_time)

        for velocity_domain_feature in velocity_domain_features:
            # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
            # Iterate over them all.
            for velocity_domain_geometry in velocity_domain_feature.get_geometries():

                for velocity_domain_point in velocity_domain_geometry.get_points():

                    all_domain_points.append(velocity_domain_point)

                    partitioning_plate = plate_partitioner.partition_point(velocity_domain_point)
                    if partitioning_plate:

                        # We need the newly assigned plate ID
                        # to get the equivalent stage rotation of that tectonic plate.
                        partitioning_plate_id = partitioning_plate.get_feature().get_reconstruction_plate_id()

                        # Get the stage rotation of partitioning plate from 'time + delta_time' to 'time'.
                        equivalent_stage_rotation = rotation_model.get_rotation(reconstruction_time,
                                                                                partitioning_plate_id,
                                                                                reconstruction_time + delta_time)

                        # Calculate velocity at the velocity domain point.
                        # This is from 'time + delta_time' to 'time' on the partitioning plate.
                        velocity_vectors = pygplates.calculate_velocities(
                            [velocity_domain_point],
                            equivalent_stage_rotation,
                            delta_time)

                        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples
                        # (one tuple per point).
                        velocities =pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                                [velocity_domain_point],
                                velocity_vectors)
                        all_velocities.append((velocities[0].get_x(), velocities[0].get_y()))

                    else:
                        all_velocities.append((0,0))
                        
        return np.array(all_velocities)



class Points(object):

    def __init__(self, PlateReconstruction_object, lons, lats, time=0, plate_id=None):
        self.lons = lons
        self.lats = lats
        self.time = time

        # get Cartesian coordinates
        self.x, self.y, self.z = _tools.lonlat2xyz(lons, lats)

        # scale by average radius of the Earth
        self.x *= pygplates.Earth.mean_radius_in_kms
        self.y *= pygplates.Earth.mean_radius_in_kms
        self.z *= pygplates.Earth.mean_radius_in_kms

        # store concatenated arrays
        self.lonlat = np.c_[self.lons, self.lats]
        self.xyz = np.c_[self.x, self.y, self.z]


        rotation_model = PlateReconstruction_object.rotation_model
        static_polygons = PlateReconstruction_object.static_polygons
        self.PlateReconstruction_object = PlateReconstruction_object

        features = points_to_features(lons, lats, plate_id)

        if plate_id:
            self.features = features
        else:
            # partition using static polygons
            # being careful to observe 'from time'
            partitioned_features = pygplates.partition_into_plates(static_polygons, rotation_model, features)
            self.features = partitioned_features

        self.FeatureCollection = pygplates.FeatureCollection(self.features)


    def reconstruct(self, time, anchor_plate_id=0, **kwargs):
        from_time = self.time
        to_time = time
        reconstructed_features = self.PlateReconstruction_object.reconstruct(
            self.features, to_time, from_time, anchor_plate_id=anchor_plate_id, **kwargs)

        rlons, rlats = _tools.extract_feature_lonlat(reconstructed_features)
        return rlons, rlats

    def plate_velocity(self, time, delta_time=1):
        return self.PlateReconstruction_object.get_point_velocities(self.lons, self.lats, time, delta_time)


    def save(self, filename):
        self.FeatureCollection.save(filename)