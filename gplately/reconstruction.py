import pygplates
import numpy as np
import ptt


_rotation_filenames= _data_dir + '/Alps_Mesh_Rotations_2019_v2.rot',\
                     _data_dir + '/Andes_Flat_Slabs_Rotations_2019_v2.rot',\
                     _data_dir + '/Andes_Rotations_2019_v2.rot',\
                     _data_dir + '/Australia_Antarctica_Mesh_Rotations_2019_v2.rot',\
                     _data_dir + '/Australia_North_Zealandia_Rotations_2019_v2.rot',\
                     _data_dir + '/Eurasia_Arabia_Mesh_Rotations_2019_v2.rot',\
                     _data_dir + '/Global_250-0Ma_Rotations_2019_v2.rot',\
                     _data_dir + '/North_America_Flat_Slabs_Rotations_2019_v2.rot',\
                     _data_dir + '/North_America_Mesh_Rotations_2019_v2.rot',\
                     _data_dir + '/North_China_Mesh_Rotations_2019_v2.rot',\
                     _data_dir + '/South_Atlantic_Rotations_2019_v2.rot',\
                     _data_dir + '/Southeast_Asia_Rotations_2019_v2.rot'


_topology_filenames = input_directory + '/Global_Mesozoic-Cenozoic_PlateBoundaries_2019_v2.gpml',\
                      input_directory + '/Alps_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/America_Anyui_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Andes_Flat_Slabs_Topologies_2019_v2.gpml',\
                      input_directory + '/Andes_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Arctic_Eurasia_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Australia_Antarctica_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Australia_North_Zealandia_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Coral_Sea_Topologies_2019_v2.gpml',\
                      input_directory + '/East_African_Rift_Deforming_Mesh_and_Topologies_2019_v2.gpml',\
                      input_directory + '/East-West_Gondwana_Deforming_Mesh_and_Topologies_2019_v2.gpml',\
                      input_directory + '/Eurasia_Arabia_Deforming_Mesh_and_Topologies_2019_v2.gpml',\
                      input_directory + '/Greater_India_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Inactive_Meshes_and_Topologies_2019_v2.gpml',\
                      input_directory + '/North_America_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/North_Atlantic_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/North_China_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Northern_Andes_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Papua_New_Guinea_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Scotia_Deforming_Mesh_and_Topologies_2019_v2.gpml',\
                      input_directory + '/Siberia_Eurasia_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/South_Atlantic_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/South_China_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/South_Zealandia_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/Southeast_Asia_Mesh_Topologies_2019_v2.gpml',\
                      input_directory + '/West_Antarctica_Zealandia_Mesh_Topologies_2019_v2.gpml', \
                      input_directory + '/Western_Tethys_Tectonic_Boundary_Topologies_2019_v2.gpml'



class PlateReconstruction(object):

    def __init__(self, rotation_model=None, topology_features=None):


        if rotation_model is None:
            rotation_model = _rotation_filenames
        if topology_features is None:
            topology_features = _topology_filenames


        rotation_model = pygplates.RotationModel(_rotation_filenames)

        default_topology_features = pygplates.FeatureCollection()
        for topology in topology_features:
            default_topology_features.add( pygplates.FeatureCollection(_topology_filename) )

        self.rotation_model = rotation_model
        self.topology_features = default_topology_features


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

    def reconstruct(self, feature, to_time, from_time=0, anchor_plate_id=0):
        pass




class Points(object):

    def __init__(self, lons, lats, PlateReconstruction_object, time=0, plate_id=None):
        self.lons = lons
        self.lats = lats
        self.time = time

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


    def save(self, filename):
        pygplates
