import pygplates
import numpy as np
from pydoc import plain
import stripy 
import ptt
import multiprocessing
import glob
import os
from skimage import measure

from . import reconstruction
from . import plot
from . import grids
from . import tools
from ptt import separate_ridge_transform_segments


def create_icosahedral_mesh(refinement_levels):
    # Create the ocean basin mesh using Stripy's icosahedral spherical mesh
    icosahedral_global_mesh = stripy.spherical_meshes.icosahedral_mesh(
        refinement_levels, 
        include_face_points=False, 
        trisection=False, 
        tree=False
    )
    # Get lons and lats of mesh, and turn them into a MultiPointOnSphere
    lats_arr = np.rad2deg(icosahedral_global_mesh.lats)
    lons_arr = np.rad2deg(icosahedral_global_mesh.lons)
    multi_point = pygplates.MultiPointOnSphere(zip(lats_arr,lons_arr))

    return multi_point, icosahedral_global_mesh


def ensure_polygon_geometry(reconstructed_polygons, rotation_model, time):
    """ Ensure COB terrane geometries are polygons with reconstruction plate 
    IDs and valid times before they are used to:
    - identify ocean basin points, and 
    - form a continental mask for `ReconstructByTopologies`.
    """
    continent_FeatCol = []
    # self._PlotTopologies_object.continents
    for n in reconstructed_polygons:
        continent_FeatCol.append(n.get_feature())

    polygon_feats = pygplates.FeatureCollection(continent_FeatCol)

    # From GPRM's force_polygon_geometries(); ascribe feature attributes
    # like valid times and plate IDs
    polygons = []
    for feature in polygon_feats: 
        for geom in feature.get_all_geometries():
            polygon = pygplates.Feature(
                feature.get_feature_type()
            )
            polygon.set_geometry(
                pygplates.PolygonOnSphere(geom)
            )
            polygon.set_reconstruction_plate_id(
                feature.get_reconstruction_plate_id()
            )
            # Avoid features in COBTerranes with invalid time
            if feature.get_valid_time()[0]>=feature.get_valid_time()[1]:
                polygon.set_valid_time(
                    feature.get_valid_time()[0],feature.get_valid_time()[1]
                )
                polygons.append(polygon)
    cobter_polygon_features = pygplates.FeatureCollection(polygons)

    # Turn the feature collection back into ReconstructedFeatureGeometry 
    # objects otherwise it will not work with PIP
    reconstructed_cobter_polygons = []
    pygplates.reconstruct(
        cobter_polygon_features, 
        rotation_model, 
        reconstructed_cobter_polygons, 
        time
    )
    return reconstructed_cobter_polygons


def point_in_polygon_routine(multi_point, COB_polygons):
    """ Perform Plate Tectonic Tools' point in polygon routine to determine
    which points from a MultiPointOnSphere feature are in or out of polygons 
    from `COB_polygons`.

    Returns
    -------
    point_mesh_in : FeatureCollection of MultiPointOnSphere objects
        Point features that are within COB terrane polygons.
    point_mesh_out : FeatureCollection of MultiPointOnSphere objects
        Point features that are outside COB terrane polygons.
    zvals : list
        A binary list. If an entry is == 0, its corresponing point in the 
        MultiPointOnSphere object is on the ocean. If == 1, the point is 
        in the COB terrane polygon.
    """

    # Collect continental polygon features and reconstructed geometries
    # By default, looks for points in polygons
    polygons = []
    polygon_features = []
    for reconstructed_continental_geometry in COB_polygons:
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
    # Look for points inside polygons, i.e. points in COB terrane polygons
    points_in_arr = []
    points_out_arr = []
    zvals = []
    for point_index, polygon_feature_list in enumerate(continental_polygon_features_containing_points):
        if polygon_feature_list:
            points_in_arr.append(multi_point[point_index])
            zvals.append(1)
        else:
            points_out_arr.append(multi_point[point_index])
            zvals.append(0)

    return pygplates.MultiPointOnSphere(points_in_arr), pygplates.MultiPointOnSphere(points_out_arr), zvals


def _extract_point_feature_attributes_for_rbt(ocean_basin, mor_all_times, time_arr):
    """ Used to extract feature attributes from all point features over time, like ocean basin
    seed points or MOR segment points.

    TO DO: Expand this for general shapefile attribute allocation.
    """
    # Attributes to extract
    active_points = []  
    appearance_time = []  
    birth_lat = []  # latitude_of_crust_formation
    prev_lat = []
    prev_lon = []
    spreading_rates = []

    # ocean_basin is a PointOnSphere object already
    for feature in ocean_basin:
        active_points.append(feature.get_geometry())
        appearance_time.append(feature.get_valid_time()[0])
        birth_lat.append(
            feature.get_geometry().to_lat_lon_list()[0][0]
        )
        prev_lat.append(
            feature.get_geometry().to_lat_lon_list()[0][0]
        )
        prev_lon.append(
            feature.get_geometry().to_lat_lon_list()[0][1]
        )
        try:
            spreading_rate = float(feature.get_shapefile_attribute("SPREADING_RATE"))
        except Exception:
            spreading_rate = np.nan
        spreading_rates.append(spreading_rate)
    del ocean_basin

    for i, time in enumerate(time_arr):
        features = pygplates.FeatureCollection(
            mor_all_times[i]
        )
        for feature in features:
            # seeds_from_topologies.append(feature)
            if feature.get_valid_time()[0]<time_arr[0]:
                active_points.append(
                    feature.get_geometry()
                )
                appearance_time.append(
                    feature.get_valid_time()[0]
                )
                birth_lat.append(
                    feature.get_geometry().to_lat_lon_list()[0][0]
                )
                prev_lat.append(
                    feature.get_geometry().to_lat_lon_list()[0][0]
                )
                prev_lon.append(
                    feature.get_geometry().to_lat_lon_list()[0][1]
                )
                try:
                    spreading_rate = float(feature.get_shapefile_attribute("SPREADING_RATE"))
                except Exception:
                    spreading_rate = np.nan
                spreading_rates.append(spreading_rate)

    return active_points, appearance_time, birth_lat, prev_lat, prev_lon, spreading_rates


# TO-DO: extend this to >1 attribute, and consider time-dependence of attributes. Improve speed?
def set_shapefile_attribute(
    pygplates_feature, 
    shapefile_attribute, 
    shapefile_attribute_name, 
    overwrite_attribute=False
):
    new_attribute = {shapefile_attribute_name : shapefile_attribute}
    
    # If the feature has no shapefile attribute dictionary, make one
    if pygplates_feature.get_shapefile_attributes() is None:
        pygplates_feature.set_shapefile_attributes()
    
    # If the shapefile attribute does not already exist in the dictionary, add it
    if shapefile_attribute_name not in pygplates_feature.get_shapefile_attributes().keys():
        pygplates_feature.set_shapefile_attributes(
            {**pygplates_feature.get_shapefile_attributes(), **new_attribute}
        )
    # If the shapefile attribute exists already but needs to be overwritten
    elif (shapefile_attribute_name in pygplates_feature.get_shapefile_attributes().keys()
         and
         overwrite_attribute):
        pygplates_feature.set_shapefile_attributes(
            {**pygplates_feature.get_shapefile_attributes(), **new_attribute}
        )
    return pygplates_feature


class SeafloorGrid(object):

    """A class to generate grids that track data atop global ocean basin points 
    (which emerge from mid ocean ridges) through geological time.
    """

    def __init__(
        self, 
        PlateReconstruction_object, 
        PlotTopologies_object,
        max_time,
        min_time,
        ridge_time_step,
        save_directory=None,
        file_collection=None,
        refinement_levels=5, 
        ridge_sampling=0.1,
        subduction_collision_parameters = (5.0, 10.0),
        initial_ocean_mean_spreading_rate = 75.
        ):

        # Provides a rotation model, topology features and reconstruction time for 
        # the SeafloorGrid
        self.PlateReconstruction_object = PlateReconstruction_object
        self.rotation_model = self.PlateReconstruction_object.rotation_model
        self.topology_features = self.PlateReconstruction_object.topology_features
        self._PlotTopologies_object = PlotTopologies_object
        self.save_directory = save_directory
        self.file_collection = file_collection

        # Topological parameters
        self.refinement_levels = refinement_levels
        self.ridge_sampling = ridge_sampling
        self.subduction_collision_parameters = subduction_collision_parameters
        self.initial_ocean_mean_spreading_rate = initial_ocean_mean_spreading_rate

        # Temporal parameters
        self._max_time = max_time
        self.min_time = min_time
        self.ridge_time_step = ridge_time_step
        self.time_array = np.arange(self._max_time, self.min_time-0.1, -self.ridge_time_step)

        # If PlotTopologies' time attribute is not equal to the maximum time in the 
        # seafloor grid reconstruction tree, make it equal. This will ensure the time
        # for continental masking is consistent.
        if self._PlotTopologies_object.time != self._max_time:
            self._PlotTopologies_object.time = self._max_time

        # Essential features and meshes for the SeafloorGrid
        self._PlotTopologies_object.continents = PlotTopologies_object.continents
        self.icosahedral_multi_point, self.icosahedral_global_mesh = create_icosahedral_mesh(self.refinement_levels)


    # Allow SeafloorGrid time to be updated, and to update the internally-used 
    # PlotTopologies' time attribute too. If PlotTopologies is used outside the
    # object, its `time` attribute is not updated. 
    @property
    def max_time(self):
        """ The reconstruction time."""
        return self._max_time


    @property
    def PlotTopologiesTime(self):
        return self._PlotTopologies_object.time


    @max_time.setter
    def max_time(self, var):
        if var >= 0:
            self.update_time(var)
        else:
            raise ValueError("Enter a valid time >= 0")


    def update_time(self, max_time):
        self._max_time = float(max_time)
        self._PlotTopologies_object.time = float(max_time)


    def create_initial_ocean_seed_points(self):
        """ Create the initial ocean basin seed point domain (at `max_time` only) 
        using Stripy's icosahedral triangulation with the specified 
        `self.refinement_levels`. 

        Notes
        ----- 
        This point mesh represents ocean basin seafloor that was produced
        before `SeafloorGrid.max_time`, and thus has unknown properties like valid
        time and spreading rate. As time passes, the plate reconstruction model sees 
        points emerging from MORs. These new points spread to occupy the ocean basins, 
        moving the initial filler points closer to subduction zones and continental 
        polygons with which they can collide. If a collision is detected by 
        `PlateReconstruction`s ReconstructByTopologies object, these points are deleted. 
        Optimally, if a reconstruction tree spans a large time range, these initial mesh 
        points completely disappear, leaving behind a mesh of well-defined MOR-emerged ocean
        basin points that data can be attributed to.

        `create_initial_ocean_seed_points` accesses continental polygons from the continent 
        shapefile or GPML file  attributed to the `PlotTopologies_object`. Ideally this should 
        be a COB terrane file. 
        Automatically resolves continental polygons to the `SeafloorGrid.time` attribute
        using PlotTopologies.
        Plate Tectonic Tools' point-in-polygon routine identifies ocean basin points 
        that lie:
        * outside the polygons (for the ocean basin point domain)
        * inside the polygons (for the continental mask)

        Outputs the ocean basin seed point mesh as a GPML file with the filename:
        "ocean_basin_seed_points_{}Ma.gpml" if a `save_directory` is passed.
        Otherwise, the mesh is returned as a pyGPlates FeatureCollection object.

        Returns
        -------
        ocean_basin_point_mesh : pygplates.FeatureCollection of pygplates.MultiPointOnSphere
            A feature collection of point objects on the ocean basin.
        """
        print("Generating global point mesh...")

        # Ensure COB terranes have reconstruction IDs and valid times
        COB_polygons = ensure_polygon_geometry(
            self._PlotTopologies_object.continents, 
            self.rotation_model,
            self._max_time
        )
        # zval is a binary array encoding whether a point 
        # coordinate is within a COB terrane polygon or not.
        # Use the icosahedral mesh MultiPointOnSphere attribute
        _, ocean_basin_point_mesh, zvals = point_in_polygon_routine(
            self.icosahedral_multi_point, 
            COB_polygons
        )
        
        # Plates to partition with
        plate_partitioner = pygplates.PlatePartitioner(
            COB_polygons,
            self.rotation_model, 
        )
        print("Partitioning global mesh by COB terrane plate polygons...")

        # Plate partition the ocean basin points
        meshnode_feature = pygplates.Feature(
            pygplates.FeatureType.create_from_qualified_string('gpml:MeshNode')
        )
        meshnode_feature.set_geometry(
            ocean_basin_point_mesh
            #multi_point
        )
        ocean_basin_meshnode = pygplates.FeatureCollection(meshnode_feature)

        paleogeography = plate_partitioner.partition_features(
            ocean_basin_meshnode, partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned,
            properties_to_copy=[pygplates.PropertyName.gpml_shapefile_attributes]
        )
        ocean_points = paleogeography[1]  # Separate those inside polygons
        continent_points = paleogeography[0] # Separate those outside polygons

        # Determine age of ocean basin points using their proximity to MOR features
        # and an assumed globally-uniform ocean basin mean spreading rate
        print("Determining initial ocean basin point ages as of {} Ma...".format(self._max_time))

        resolved_topologies = []
        shared_boundary_sections = []
        pygplates.resolve_topologies(
            self.topology_features, 
            self.rotation_model, 
            resolved_topologies, 
            self._max_time, 
            shared_boundary_sections
        )
        pX,pY,pZ = tools.find_distance_to_nearest_ridge(
            resolved_topologies, 
            shared_boundary_sections, 
            #global_points
            ocean_points
        )

        # Divide spreading rate by 2 to use half the mean spreading rate
        pAge = np.array(pZ) / (self.initial_ocean_mean_spreading_rate/2.)

        initial_ocean_point_features = []
        initial_ocean_multipoints = []

        for point in zip(pX,pY,pAge):

            point_feature = pygplates.Feature()
            #point_feature.set_name(str(initial_ocean_mean_spreading_rate))
            point_feature.set_geometry(pygplates.PointOnSphere(point[1], point[0]))

            # note that we add 'time' to the age at the time of computation
            # to get the valid time in Ma
            point_feature.set_valid_time(point[2]+self._max_time, -1)
            point_feature = set_shapefile_attribute(point_feature, self.initial_ocean_mean_spreading_rate, "SPREADING_RATE")  # Seems like static data
            initial_ocean_point_features.append(point_feature)
            initial_ocean_multipoints.append(point_feature.get_geometry())

        # print(initial_ocean_point_features)
        multi_point_feature = pygplates.MultiPointOnSphere(initial_ocean_multipoints)

        if self.save_directory:
            if self.file_collection is not None:
                full_directory = "{}/{}_ocean_basin_seed_points_{}Ma.gpml".format(
                    self.save_directory,
                    self.file_collection, 
                    self._max_time
                )
            else:
                full_directory = "{}/ocean_basin_seed_points_{}Ma.gpml".format(
                    self.save_directory, 
                    self._max_time
                )
            pygplates.FeatureCollection(initial_ocean_point_features).write(full_directory)

        return pygplates.FeatureCollection(initial_ocean_point_features), multi_point_feature


    def get_mid_ocean_ridge_seedpoints(self, time):
        """ Resolve mid-ocean ridges to the current `time`, and shift their shared sub-segments
        them slightly off the ridge to the left and to the right using their stage rotation. 

        Adapted from an age gridding workflow by Simon Williams, John Cannon and Nicky Wright.

        Returns
        -----
        mor_point_features : FeatureCollection
            Ridge seed points that have been slightly rotated away from
            ridge locations at the current timestep.
        """

        # Topology features are already resolved to `time`. 
        topology_features_extracted = pygplates.FeaturesFunctionArgument(self.topology_features)

        # Resolve topologies to the current time.
        resolved_topologies = []
        shared_boundary_sections = []
        pygplates.resolve_topologies(
            topology_features_extracted.get_features(), 
            self.rotation_model, 
            resolved_topologies, 
            time, 
            shared_boundary_sections
        )
        shifted_mor_points = []
        point_spreading_rates = []

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
                    spreading_feature, self.rotation_model, time
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

                    coords = mor_points.to_lat_lon_list()
                    lats = [i[0] for i in coords]
                    lons = [i[1] for i in coords]
                    left_plate = shared_boundary_section.get_feature().get_left_plate(None)
                    right_plate = shared_boundary_section.get_feature().get_right_plate(None)
                    if left_plate is not None and right_plate is not None:
                        spreading_rates = tools.calculate_spreading_rates(
                            time=time,
                            lons=lons,
                            lats=lats,
                            left_plates=[left_plate] * len(lons),
                            right_plates=[right_plate] * len(lons),
                            rotation_model=self.rotation_model,
                            delta_time=self.ridge_time_step,
                        )

                    else:
                        spreading_rates = [np.nan] * len(lons)
                    for point, rate in zip(
                        mor_points.get_points()[1:-1],
                        spreading_rates[1:-1],
                    ):
                        shifted_mor_points.append(rotate_slightly_off_mor_one_way * point)
                        shifted_mor_points.append(rotate_slightly_off_mor_opposite_way * point)
                        point_spreading_rates.extend([rate] * 2)

        # Summarising get_isochrons_for_ridge_snapshot;
        # Write out the ridge point born at 'ridge_time' but their position at 'ridge_time - time_step'.
        mor_point_features = []
        for curr_point, spreading_rate in zip(shifted_mor_points, point_spreading_rates):
            feature = pygplates.Feature()
            feature.set_geometry(curr_point)
            feature.set_valid_time(time, -999)  # delete - time_step
            #feature.set_name(str(spreading_rate))
            feature = set_shapefile_attribute(feature, spreading_rate, "SPREADING_RATE")  # make spreading rate a shapefile attribute
            mor_point_features.append(feature)

        mor_points = pygplates.FeatureCollection(mor_point_features)

        # The following file is size-intensive
        #mor_points.write('{}/MOR_plus_one_points_{:0.2f}.gpml'.format(
        #    self.save_directory, 
        #    time)
        #)
        return mor_points



    def create_continental_mask(self, time):
        """ Create a continental mask for use as a specified collision type in 
        `ReconstructByTopologies`. 

        If `ReconstructByTopologies` identifies a continental collision 
        between oceanic points and the boundaries of this continental
        mask at `time`, those points are deleted at `time`. 

        The continental mask is also saved to "/continent_mask_{}Ma.nc"
        if a `save_directory` is passed. Otherwise, the final grid is returned 
        as a NumPy ndarray object.

        Returns
        -------
        final_grid : ndarray
            A masked grid with 1=continental point, and 0=ocean point, for all points
            on the full global icosahedral mesh.
        """

        # Ensure COB terranes have reconstruction IDs and valid times
        self._PlotTopologies_object.time = time
        COB_polygons = ensure_polygon_geometry(
            self._PlotTopologies_object.continents,
            self.rotation_model,
            time) 

        # zval is a binary array encoding whether a point 
        # coordinate is within a COB terrane polygon or not
        _, _, zvals = point_in_polygon_routine(
            self.icosahedral_multi_point, 
            COB_polygons
        )

        # Interpolate the zval binaries onto the icosahedral global mesh
        boolean_identifier, _ = self.icosahedral_global_mesh.interpolate(
            self.icosahedral_global_mesh.lons, 
            self.icosahedral_global_mesh.lats, order=3, 
            zdata=np.array(zvals)
        )

        # A regular grid to interpolate the spherical mesh onto
        # TO-DO: resX, resY and extent should be user-defined?
        resX = 200
        resY = 100
        extent_globe = np.radians([-180,180,-90,90])
        grid_lon = np.linspace(extent_globe[0], extent_globe[1], resX)
        grid_lat = np.linspace(extent_globe[2], extent_globe[3], resY)

        # Interpolate the icosahedral-meshed zval binaries onto the regular global extent grid
        grid_z1 = self.icosahedral_global_mesh.interpolate_to_grid(
            grid_lon, grid_lat, boolean_identifier
        )
        # Ensure the PIP binaries are integers
        final_grid = np.rint(grid_z1)

        if self.save_directory is not None:
            if self.file_collection is not None:
                full_directory = "{}/{}_continent_mask_{}Ma.nc".format(
                    self.save_directory, 
                    self.file_collection, 
                    time
                )
            else:
                full_directory = "{}/continent_mask_{}Ma.nc".format(
                    self.save_directory, 
                    time
                )
            grids.write_netcdf_grid(
                full_directory, 
                final_grid, 
                extent=[-180,180,-90,90]
            )
        #print("Continental mask for {} Ma done!".format(time))
        return final_grid


    def prepare_for_reconstruction_by_topologies(self):
        """ Prepare three main auxiliary files for seafloor data gridding:
        * Initial ocean seed points (at `max_time`)
        * MOR points (from `max_time` to `min_time`)
        * Continental masks (from `max_time` to `min_time`)
        """

        initial_ocean_seed_points, initial_ocean_seed_points_mp = self.create_initial_ocean_seed_points()

        print("Finished building initial_ocean_seed_points!")
        time_array = np.arange(self._max_time, self.min_time-1, -self.ridge_time_step)
        all_mor_features = []
        all_continental_masks = []

        for time in self.time_array:
            all_mor_features.append(
                self.get_mid_ocean_ridge_seedpoints(time)
            )
            print("Finished building MOR seedpoints at {} Ma!".format(time))
            
            all_continental_masks.append(
                self.create_continental_mask(time)
            )
            print("Finished building a continental mask at {} Ma!".format(time))

        active_points, appearance_time, \
        birth_lat, prev_lat, prev_lon, spreading_rates  = _extract_point_feature_attributes_for_rbt(
            initial_ocean_seed_points,
            all_mor_features,
            self.time_array
        )
        return active_points, appearance_time, birth_lat, prev_lat, prev_lon, spreading_rates


    def reconstruct_by_topologies(self):
        """ Obtain all active ocean seed points at `time` - these are 
        points that have not been consumed at subduction zones or have not
        collided with continental polygons.
        """

        print("Preparing all initial files...")

        # Obtain all info from the ocean seed points and MOR, store in
        # one array.
        active_points, appearance_time, birth_lat, \
        prev_lat, prev_lon, spreading_rates = self.prepare_for_reconstruction_by_topologies()

        ####  Begin reconstruction by topology process:
        # Indices for all active points at the current time step
        point_id = range(len(active_points))

        # Specify the collision detection
        default_collision = reconstruction.DefaultCollision(
            feature_specific_collision_parameters = [
            (pygplates.FeatureType.gpml_subduction_zone, self.subduction_collision_parameters)
            ]
        )
        # In addition to the default subduction detection, also detect continental collisions
        if self.file_collection is not None:
            collision_spec = reconstruction.ContinentCollision(
                self.save_directory+"/"+self.file_collection+"_continent_mask_{}Ma.nc", 
                default_collision
            )
        else:
            collision_spec = reconstruction.ContinentCollision(
                self.save_directory+"/continent_mask_{}Ma.nc", 
                default_collision
            )

        # Call the reconstruct by topologies object
        topology_reconstruction = reconstruction.ReconstructByTopologies(
            self.rotation_model, 
            self.topology_features,
            self._max_time, 
            self.min_time, 
            self.ridge_time_step,
            active_points, 
            point_begin_times=appearance_time,
            detect_collisions = collision_spec
        )
        # Initialise the reconstruction.
        topology_reconstruction.begin_reconstruction()

        # Loop over the reconstruction times until reached end of the reconstruction time span, or
        # all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        reconstruction_data = []
        while True:
            print('Reconstruct by topologies: working on time {:0.2f} Ma'.format(
                topology_reconstruction.get_current_time())
            )

            # Collect latitudes and longitudes of currently active points in the ocean basin mesh
            curr_points = topology_reconstruction.get_active_current_points()
            #nan_points, nan_point_indices = topology_reconstruction.get_in_continent_indices()

            curr_lat_lon_points = [point.to_lat_lon() for point in curr_points]

            if curr_lat_lon_points:
                curr_latitudes, curr_longitudes = zip(*curr_lat_lon_points)

                # TO BE REPLACED W/ USER-INPUT DATA LATER
                seafloor_age = []
                spreading_rate_snapshot = []

                # Time-dependent point attributes
                birth_lat_snapshot = []
                point_id_snapshot = []
                prev_lat_snapshot = []
                prev_lon_snapshot = []

                for point_index,current_point in enumerate(topology_reconstruction.get_all_current_points()):

                    # Look at all active points (these have not collided with a continent or trench)
                    if current_point is not None:
                        seafloor_age.append(
                                appearance_time[point_index] - topology_reconstruction.get_current_time()
                            )
                        birth_lat_snapshot.append(birth_lat[point_index])
                        point_id_snapshot.append(point_id[point_index])
                        prev_lat_snapshot.append(prev_lat[point_index])
                        prev_lon_snapshot.append(prev_lon[point_index])
                        
                        # TO-DO: add a general method to extract shapefile attributes with user-input data
                        spreading_rate_snapshot.append(spreading_rates[point_index])
                        
                        prev_lat[point_index] = current_point.to_lat_lon()[0]
                        prev_lon[point_index] = current_point.to_lat_lon()[1]

                # TO-DO: Look for storage-efficient alternative(s) to .csv intermediate/aux files to store 
                # these arrays (unique to each reconstruction time).
                header = 'CURRENT_LONGITUDES,'\
                'CURRENT_LATITUDES,'\
                'SEAFLOOR_AGE,'\
                'SPREADING_RATE_SNAPSHOT,'\
                'BIRTH_LAT_SNAPSHOT,'\
                'POINT_ID_SNAPSHOT'

                zippeddata = list(zip(curr_longitudes, curr_latitudes, seafloor_age, spreading_rate_snapshot, birth_lat_snapshot, point_id_snapshot, 
                #any extra user-input data here), 
                ))
                if self.file_collection is not None:
                    np.savetxt('{:s}/{}_gridding_input_{:0.1f}Ma.csv'.format(self.save_directory, self.file_collection, topology_reconstruction.get_current_time()), 
                        zippeddata, header=header, delimiter=',', comments='', fmt='%s'
                    )
                else:
                    np.savetxt('{:s}/gridding_input_{:0.1f}Ma.csv'.format(self.save_directory, topology_reconstruction.get_current_time()), 
                        zippeddata, header=header, delimiter=',', comments='', fmt='%s'
                    )

            if not topology_reconstruction.reconstruct_to_next_time():
                break

        print('Reconstruction done for {}!'.format(topology_reconstruction.get_current_time()))
        # return reconstruction_data
