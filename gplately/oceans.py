import pygplates
import numpy as np
from pydoc import plain
import stripy 
import ptt
import multiprocessing
import glob
from skimage import measure

from . import reconstruction
from . import plot
from . import grids
from ptt import separate_ridge_transform_segments
from gplately.utils import reconstruct_by_topologies as rbt


def _create_icosahedral_mesh(refinement_levels):
        # Create the ocean basin mesh using a fine icosahedral spherical mesh
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


def _ensure_polygon_geometry(reconstructed_polygons, rotation_model, time):
    """ Ensure COB terrane geometries have reconstruction plate IDs and valid times
    before they are used to:
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


def _point_in_polygon_routine(multi_point, COB_polygons):
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

    # Turn ocean basin points into a FeatureCollection of MultiPointOnSphere geometries
    multi_point_features_in = []
    multi_point_features_out = []
    multi_point_feature_in = pygplates.Feature()
    multi_point_feature_out = pygplates.Feature()

    multi_point_feature_in.set_geometry(
        pygplates.MultiPointOnSphere(points_in_arr)
    )
    multi_point_feature_out.set_geometry(
        pygplates.MultiPointOnSphere(points_out_arr)
    )
    multi_point_features_in.append(multi_point_feature_in)
    multi_point_features_out.append(multi_point_feature_out)
    point_mesh_in = pygplates.FeatureCollection(multi_point_features_in)
    point_mesh_out = pygplates.FeatureCollection(multi_point_features_out)

    return point_mesh_in, point_mesh_out, zvals


def _extract_point_feature_attributes(point_feature_collections):
    """ Used to extract feature attributes from point features, like ocean basin
    seed points or MOR segment points.
    """
    # Attributes to extract
    active_points = []  # current_point # TODO more verbose names for cp and at
    appearance_time = []  # appearance_time
    birth_lat = []  # latitude_of_crust_formation
    prev_lat = []
    prev_lon = []

    for point_feature_collection in point_feature_collections:
        for feature in point_feature_collection:
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
    return active_points, appearance_time, birth_lat, prev_lat, prev_lon





class SeafloorGrid(object):

    """A class with tools to track static and dynamic data on global ocean basins 
    through geological time.
    """

    def __init__(
        self, 
        PlateReconstruction_object, 
        PlotTopologies_object,
        max_time,
        min_time,
        ridge_time_step,
        refinement_levels=5, 
        ridge_sampling=0.1,
        subduction_collision_parameters = (5.0, 10.0), #Should raise an error if not supplied.
        ):

        # Provides a rotation model, topology features and reconstruction time for 
        # the SeafloorGrid
        self.PlateReconstruction_object = PlateReconstruction_object
        self.rotation_model = self.PlateReconstruction_object.rotation_model
        self.topology_features = self.PlateReconstruction_object.topology_features
        self._PlotTopologies_object = PlotTopologies_object
        # self._PlotTopologies_object.continents = PlotTopologies_object.continents


        # Topological parameters
        self.refinement_levels = refinement_levels
        self.ridge_sampling = ridge_sampling
        self.subduction_collision_parameters = subduction_collision_parameters

        # Temporal parameters
        self._max_time = max_time
        self.min_time = min_time
        self.ridge_time_step = ridge_time_step
        self.time_array = np.arange(self._max_time, self.min_time-1, -self.ridge_time_step)


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


    def _create_initial_ocean_seed_points(self, save_directory=None):
        """ Create the initial ocean basin seed point domain (at `max_time` only) 
        using an icosahedral triangulation with the specified 
        `self.refinement_levels`, 

        Notes
        ----- 
        Accesses continental polygons from the continent shapefile or GPML file 
        attributed to the `PlotTopologies_object`. The object automatically resolves
        the continental polygons to the `time` set in `SeafloorGrid.time` attribute
        (see the `PlotTopologies` object for more information). 
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

        # Create MultiPointOnSphere from stripy's icosahedral spherical mesh
        multi_point, icosahedral_global_mesh = _create_icosahedral_mesh(self.refinement_levels)

        # Ensure COB terranes have reconstruction IDs and valid times
        COB_polygons = _ensure_polygon_geometry(
            self._PlotTopologies_object.continents, 
            self.rotation_model,
            self._max_time
        )
        # zval is a binary array encoding whether a point 
        # coordinate is within a COB terrane polygon or not
        ocean_basin_point_mesh, _, _ = _point_in_polygon_routine(
            multi_point, 
            COB_polygons
        )
        ### PREPARE THE OCEAN BASIN SEED POINTS
        if save_directory:
            full_directory = save_directory+"/ocean_basin_seed_points_{}Ma.gpml".format(self._max_time)
            ocean_basin_point_mesh.write(filename)
        return ocean_basin_point_mesh


    ### FUNCTIONS TO REPEAT IN TIME LOOPS/POOL? 
    def _get_mid_ocean_ridge_seedpoints(self, time):
        """ Resolve mid-ocean ridges to the current `time`, and shift their shared sub-segments
        them slightly off the ridge using their stage rotation. 

        Adapted from an age gridding workflow by Simon Williams, John Cannon and Nicky Wright.

        Returns
        -----
        mor_point_features : FeatureCollection
            Ridge seed points that have been slightly rotated away from
            ridge locations at the current timestep.
        """

        # Get the rotation model ascribed to the plate model. Topology features
        # are already resolved to `time`. 
        rotation_model = self.rotation_model
        topology_features = self.topology_features
        topology_features_extracted = pygplates.FeaturesFunctionArgument(topology_features)

        # Resolve topologies to the current time.
        resolved_topologies = []
        shared_boundary_sections = []
        pygplates.resolve_topologies(
            topology_features_extracted.get_features(), 
            rotation_model, 
            resolved_topologies, 
            time, 
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
                    spreading_feature, rotation_model, time
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

        # Summarising get_isochrons_for_ridge_snapshot;
        # Write out the ridge point born at 'ridge_time' but their position at 'ridge_time - time_step'.
        mor_point_features = []
        for curr_point in shifted_mor_points:
            feature = pygplates.Feature()
            feature.set_geometry(curr_point)
            feature.set_valid_time(time, -999)  # delete - time_step
            #feature.set_name(str(spreading_rate))
            mor_point_features.append(feature)
        pygplates.FeatureCollection(mor_point_features)

        return mor_point_features


    def _create_continental_mask(self, time, save_directory=None):
        """ Create a continental mask for use as a specified collision type in 
        `ReconstructByTopologies`. 

        If `ReconstructByTopologies` identifies a continental collision 
        between points on the ocean basin seed point domain and continental
        mask at `time`, those points are not valid at that `time`. 

        The continental mask is also saved to "/continent_mask_{}Ma.nc"
        if a `save_directory` is passed. Otherwise, the final grid is returned 
        as a NumPy ndarray object.

        Returns
        -------
        final_grid : ndarray
            A masked grid with 1=continental point, and 0=ocean point, for all points
            on the full global icosahedral mesh.
        """
        # Create MultiPointOnSphere from stripy's icosahedral spherical mesh
        multi_point, icosahedral_global_mesh = _create_icosahedral_mesh(self.refinement_levels)
        
        # Ensure COB terranes have reconstruction IDs and valid times
        self._PlotTopologies_object.time = time
        COB_polygons = _ensure_polygon_geometry(
            self._PlotTopologies_object.continents,
            self.rotation_model,
            time) 

        # zval is a binary array encoding whether a point 
        # coordinate is within a COB terrane polygon or not
        _, _, zvals = _point_in_polygon_routine(
            multi_point, 
            COB_polygons
        )

        # Interpolate the zval binaries onto the icosahedral global mesh
        boolean_identifier, _ = icosahedral_global_mesh.interpolate(
            icosahedral_global_mesh.lons, 
            icosahedral_global_mesh.lats, order=3, 
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
        grid_z1 = icosahedral_global_mesh.interpolate_to_grid(
            grid_lon, grid_lat, boolean_identifier
        )
        # Ensure the PIP binaries are integers
        final_grid = np.rint(grid_z1)

        if save_directory:
            full_directory = save_directory+"/continent_mask_{}Ma.nc".format(self._max_time)
            gplately.grids.write_netcdf_grid(
                full_directory, final_grid, extent=[-180,180,-90,90]
            )
        return final_grid


    def prepare_for_reconstruction_by_topologies(self, save_directory=None):
        """ Prepare three main auxiliary files for seafloor data gridding:
        * Initial ocean seed points (at `max_time`)
        * MOR points (from `max_time` to `min_time`)
        * Continental masks (from `max_time` to `min_time`)
        """

        initial_ocean_seed_points = self._create_initial_ocean_seed_points(save_directory)
        print("Finished building initial_ocean_seed_points!")
        time_array = np.arange(self._max_time, self.min_time-1, -self.ridge_time_step)
        all_mor_features = []
        all_continental_masks = []

        for time in self.time_array:
            all_mor_features.append(
                self._get_mid_ocean_ridge_seedpoints(time)
            )
            print("Finished building MOR seedpoints at {} Ma!".format(time))
            all_continental_masks.append(
                self._create_continental_mask(time, save_directory)
            )
            print("Finished building a continental mask at {} Ma!".format(time))
        return initial_ocean_seed_points, all_mor_features, all_continental_masks



    def reconstruct_by_topologies(self):
        """ Obtain all active ocean seed points at `time` - these are 
        points that have not been consumed at subduction zones or have not
        collided with continental polygons.
        """
        rotation_model = self.rotation_model

        # Obtain all info from the ocean seed points and MOR, store in
        # one array.
        seeds_at_start_time, continent_mask_file_pattern = self._get_initial_ocean_seed_points_and_continental_mask()
        mor_at_start_time = self.get_mid_ocean_ridge_seedpoints()

        active_points, appearance_time, \
        birth_lat, prev_lat, prev_lon = _extract_point_feature_attributes(
            [seeds_at_start_time, mor_at_start_time]
        )
        # Conserve memory by deleting ocean seeds
        del seeds_at_start_time

        ####  Begin reconstruction by topology process:
        # Indices for all active points at the current time step
        point_id = range(len(active_points))

        # Specify the collision detection
        default_collision = rbt.DefaultCollision(
            feature_specific_collision_parameters = [
            (pygplates.FeatureType.gpml_subduction_zone, self.subduction_collision_parameters)
            ]
        )
        # In addition to the default subduction detection, also detect
        # continental collisions
        collision_spec = rbt.ContinentCollision(
            continent_mask_file_pattern, default_collision
        )

        # Call the reconstruct by topologies object
        topology_reconstruction = rbt.ReconstructByTopologies(
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
        while True:
            print('reconstruct by topologies: working on time {:0.2f} Ma'.format(
                topology_reconstruction.get_current_time())
            )

            # Collect latitudes and longitudes of currently active points in the ocean basin mesh
            curr_points = topology_reconstruction.get_active_current_points()
            curr_lat_lon_points = [point.to_lat_lon() for point in curr_points]
            if curr_lat_lon_points:
                curr_latitudes, curr_longitudes = zip(*curr_lat_lon_points)

                # TO BE REPLACED W/ USER-INPUT DATA LATER
                seafloor_age = []
                #spreading_rate_snapshot = []

                # Time-dependent point attributes
                birth_lat_snapshot = []
                point_id_snapshot = []
                prev_lat_snapshot = []
                prev_lon_snapshot = []
                for point_index,current_point in enumerate(topology_reconstruction.get_all_current_points()):
                    if current_point is not None:
                        #all_birth_ages.append(at[point_index])
                        seafloor_age.append(
                            appearance_time[point_index] - topology_reconstruction.get_current_time()
                        )
                        birth_lat_snapshot.append(birth_lat[point_index])
                        point_id_snapshot.append(point_id[point_index])
                        prev_lat_snapshot.append(prev_lat[point_index])
                        prev_lon_snapshot.append(prev_lon[point_index])
                        
                        # TO BE REPLACED W/ USER-INPUT DATA LATER
                        #spreading_rate_snapshot.append(spreading_rates[point_index])
                        
                        prev_lat[point_index] = current_point.to_lat_lon()[0]
                        prev_lon[point_index] = current_point.to_lat_lon()[1]


                # TO-DO: Remove .xy dependency, move straight to gplately grids object, or if more convenient
                # save as file for now since each array is renewed per reconstruction time. 
                write_xyz_file('{:s}/gridding_input/gridding_input_{:0.1f}Ma.xy'.format("/Users/laurenilano/Downloads",
                                                                         topology_reconstruction.get_current_time()),
                               zip(curr_longitudes,
                               curr_latitudes,
                               seafloor_age,
                               birth_lat_snapshot,
                               point_id_snapshot,
                               spreading_rate_snapshot))

            if not topology_reconstruction.reconstruct_to_next_time():
                break

        print('done')
        return curr_longitudes, curr_latitudes, seafloor_age, birth_lat_snapshot, point_id_snapshot, spreading_rate_snapshot



