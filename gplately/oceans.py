import pygplates
import numpy as np
from pydoc import plain
import stripy 
import ptt
#import pandas as pd
#import os
import multiprocessing
import glob

from . import reconstruction
from . import plot
from ptt import separate_ridge_transform_segments
#
#
#
#
#
#
#
# TO DO:
# - URGENT: Do we need the ocean basin seed points to be plate-partitioned?
# - Generalise the gridding extent 
# - Remove need for user-input re: gplot and model; make it internal 
# (although this is not urgent, should still work nonetheless but users will need to supply a time that is consistent with the time attributed to gplot ALWAYS)
#   - Will require user input for the file_collection string? Depends on whether this will be an isolated .py script.
#
#
#
#
# How time is iterated through. will we call an object once per time? or will a list iterate through an attributed min, max time?
# min_time: 0.
#    max_time: 170.
#    mor_time_step: 1.
#    gridding_time_step: 1.
#
#
#
class SeafloorGrid(object):

    """A class with tools to track static and dynamic data on global ocean basins 
    through geological time.
    """

    def __init__(
        self, 
        PlateReconstruction_object=None, 
        PlotTopologies_object=None,
        refinement_levels=5, 
        ridge_sampling=0.1,
        ridge_time_step=1.,
        time=0):

        # Provides a rotation model, topology features and reconstruction time for 
        # the SeafloorGrid
        self.PlateReconstruction_object = PlateReconstruction_object
        self._PlotTopologies_object = PlotTopologies_object
        self.refinement_levels = refinement_levels
        self.ridge_sampling = ridge_sampling
        self._time = time
        # self.PlotTopologies_object.time = self._time

        ## TO DO LATER:
        # The following lines will have to deal with how the input data looks like - 
        # is it netCDF? CSV?
        
        """
        # Ensure the data to plot on the SeafloorGrid is either a numpy array or filename. 
        if filename is None and array is None:
            raise ValueError("Supply either a filename or numpy array of data to plot on the seafloor")
        elif filename and array:
            raise ValueError("Supply either a filename or numpy array of data to plot on the seafloor")

        elif filename is not None:
        self.data, lons, lats = read_netcdf_grid(filename, return_grids=True, resample=resample)
        self.extent = [lons.min(), lons.max(), lats.min(), lats.max()]
        self.lons = lons
        self.lats = lats

        elif array is not None:
        if extent is None:
            extent = [-180,180,-90,90]
        self.data = array
        self.extent = extent
        self.lons = np.linspace(extent[0], extent[1], self.data.shape[1])
        self.lats = np.linspace(extent[2], extent[3], self.data.shape[0])
        """

    # Allow SeafloorGrid time to be updated, and to update the internally-used 
    # PlotTopologies' time attribute too. If PlotTopologies is used outside the
    # object, its `time` attribute is not updated. 
    @property
    def time(self):
        """ The reconstruction time."""
        return self._time


    @property
    def PlotTopologiesTime(self):
        return self._PlotTopologies_object.time


    @time.setter
    def time(self, var):
        if var >= 0:
            self.update_time(var)
        else:
            raise ValueError("Enter a valid time >= 0")


    def update_time(self, time):
        self._time = float(time)
        self._PlotTopologies_object.time = float(time)



    def get_initial_ocean_seed_points(self, save_directory=None):
        """ Create an ocean basin seed point domain for the specified extent using a
        pre-defined mesh for an icosahedral triangulation. 

        Notes
        -----
        Accesses continental polygons from the continent shapefile or GPML file 
        attributed to the `PlotTopologies_object`. The object automatically resolves
        the continental polygons to the `time` set in `SeafloorGrid.time` attribute.
        See the `PlotTopologies` object for more information. Once continents are
        resolved to `time`, Plate Tectonic Tools' point-in-polygon spatial tree 
        identifies ocean basin points that lie outside them.

        Outputs the ocean basin seed point mesh as a GPML file with the filename:
        "ocean_basin_seed_points_{}Ma.gpml" if a `save_directory` is passed.
        Otherwise, the mesh is returned as a pyGPlates FeatureCollection object.
        """

        # Create the ocean basin mesh using a fine icosahedral spherical mesh
        icosahedral_ocean_basin_mesh = stripy.spherical_meshes.icosahedral_mesh(
            self.refinement_levels, 
            include_face_points=False, 
            trisection=False, 
            tree=False
        )
        # Get lons and lats of mesh, and turn them into a MultiPointOnSphere
        ocean_lats = np.rad2deg(icosahedral_ocean_basin_mesh.lats)
        ocean_lons = np.rad2deg(icosahedral_ocean_basin_mesh.lons)
        multi_point = pygplates.MultiPointOnSphere(zip(ocean_lats,ocean_lons))

        # Collect continental polygon features and reconstructed geometries
        polygons = []
        polygon_features = []
        for reconstructed_continental_geometry in self._PlotTopologies_object.continents:
            polygons.append(
                reconstructed_continental_geometry.get_reconstructed_geometry()
            )
            polygon_features.append(
                reconstructed_continental_geometry.get_feature()
            )
        # Determine which continental polygons contain points from the isocahedral mesh
        continental_polygon_features_containing_points = ptt.utils.points_in_polygons.find_polygons(
            multi_point, polygons, polygon_features, all_polygons=True
        )
        # Ocean basin points aren't in the continental polygons
        points_in_ocean = []
        for point_index, polygon_feature_list in enumerate(continental_polygon_features_containing_points):
            if not polygon_feature_list:
                points_in_ocean.append(multi_point[point_index])

        # Turn ocean basin points into a FeatureCollection of MultiPointOnSphere geometries
        multi_point_features = []
        multi_point_feature = pygplates.Feature()
        multi_point_feature.set_geometry(
            pygplates.MultiPointOnSphere(points_in_ocean)
        )
        multi_point_features.append(multi_point_feature)
        ocean_basin_point_mesh = pygplates.FeatureCollection(multi_point_features)

        # Determine whether to save to GPML or return as a FeatureCollection
        if save_directory:
            full_directory = save_directory+"/ocean_basin_seed_points_{}Ma.gpml".format(self.time)
            ocean_basin_point_mesh.write(filename)
        else:
            return ocean_basin_point_mesh


        """
        FOR LATER = PLATE PARTITIONING THE OCEAN BASIN POINTS?
        plate_partitioner = pygplates.PlatePartitioner(pg_features, rotation_model, reconstruction_time=time)

        if masking is not None:
        pg_points = plate_partitioner.partition_features(raster_domain,
         partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned,
         properties_to_copy=[pygplates.PropertyName.gpml_shapefile_attributes])
        if masking == 'Outside':
        pg_points = pg_points[0]
        elif masking == 'Inside':
        pg_points = pg_points[1]
        """


    def get_mid_ocean_ridge_seedpoints(self):
        """ Resolve mid-ocean ridges to the current `time`, and shift their shared sub-segments
        them slightly off the ridge using their stage rotation. 

        Adapted from an age gridding workflow by Simon Williams, John Cannon and Nicky Wright.

        Returns
        -----
        shifted_mor_points : list of pygplates.PointOnSphere objects
            A list containing ridge seed points that have been slightly rotated away from
            ridge locations at the current timestep.
        """

        # Get the rotation model ascribed to the plate model. Topology features
        # are already resolved to `time`. 
        rotation_model = self.PlateReconstruction_object.rotation_model
        topology_features = self.PlateReconstruction_object.topology_features
        topology_features_extracted = pygplates.FeaturesFunctionArgument(topology_features)

        # Resolve topologies to the current time.
        resolved_topologies = []
        shared_boundary_sections = []
        pygplates.resolve_topologies(
            topology_features_extracted.get_features(), 
            rotation_model, 
            resolved_topologies, 
            self.time, 
            shared_boundary_sections
        )
        shifted_mor_points = []

        # pygplates.ResolvedTopologicalSection objects.
        for shared_boundary_section in shared_boundary_sections:
            if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('MidOceanRidge'):
                
                spreading_feature = shared_boundary_section.get_feature()

                # Find the stage rotation of the spreading feature in the 
                # frame of reference of its geometry at the current 
                # reconstruction time (the MOR is currently actively spreading).
                # The stage pole can then be directly geometrically compared 
                # to the *reconstructed* spreading geometry.
                stage_rotation = separate_ridge_transform_segments.get_stage_rotation_for_reconstructed_geometry(
                    spreading_feature, rotation_model, self.time
                )
                if not stage_rotation:
                    # Skip current feature - it's not a spreading feature.
                    continue

                # Get the stage pole of the stage rotation.
                # Note that the stage rotation is already in frame of 
                # reference of the *reconstructed* geometry at the spreading time.
                stage_pole, _ = stage_rotation.get_euler_pole_and_angle()

                # One way rotates left and the other right, but don't know 
                # which - doesn't matter in our example though.
                rotate_slightly_off_mor_one_way = pygplates.FiniteRotation(
                    stage_pole, 
                    np.radians(0.01)
                )
                rotate_slightly_off_mor_opposite_way = rotate_slightly_off_mor_one_way.get_inverse()
                
                # Iterate over the shared sub-segments.
                for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():

                    # Tessellate MOR section.
                    mor_points = pygplates.MultiPointOnSphere(
                        shared_sub_segment.get_resolved_geometry().to_tessellated(np.radians(self.ridge_sampling))
                    )
                    for point in mor_points.get_points()[1:-1]:
                        # Append shifted geometries (one with points rotated one way 
                        # and the other rotated the opposite way).
                        shifted_mor_points.append(
                            rotate_slightly_off_mor_one_way * point
                        )
                        shifted_mor_points.append(
                            rotate_slightly_off_mor_opposite_way * point    
                        )
        return shifted_mor_points




