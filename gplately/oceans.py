""" A module to generate grids of seafloor age, seafloor spreading rate 
and other oceanic data from the `gplately.PlateReconstruction` and 
`gplately.PlotToplogies` objects. 

Gridding methods in this script have been adapted from Simon Williams' 
development repository for an 
(auto-age-gridding workflow)[https://github.com/siwill22/agegrid-0.1].

Methodology
-----------
There are two main steps in grid generation:
1) Preparation for reconstruction by topologies
2) Reconstruction by topologies

The preparation step involves building a:
* global domain of initial points, 
* continental mask per timestep,
* set of points that emerge to the left and right of mid-ocean 
ridge segments per timestep, as well as the z-value to allocate to these
points.

First, the global domain of initial points is created using Stripy's 
icosahedral triangulated mesh. The number of points in this mesh can be
controlled using a `refinement_levels` integer (the larger this integer,
the more resolved the continent masks will be). 

These points are spatially partitioned by plate ID so they can be passed
into a point-in-polygon routine. This identifies points that lie within
continental polygon boundaries and those that are in the ocean. From this,
continental masks are built per timestep, and the initial seed points are
allocated ages at the first reconstruction timestep `max_time`. Each point's 
initial age is calculated by dividing its proximity to the nearest
MOR segment by half its assumed spreading rate. This spreading rate 
(`initial_ocean_mean_spreading_rate`) is assumed to be uniform for all points.

These initial points momentarily fill the global ocean basin. 

Ridge topologies are resolved at each reconstruction time step and partitioned
into segments with a valid stage rotation. Each segment is divided into points 
at a specified ridge sampling spacing (`ridge_sampling`). Each point is 
ascribed a latitude, longitude, spreading rate and age (from plate reconstruction 
model files, as opposed to ages of the initial ocean mesh points), a point index 
and the general z-value that will be gridded onto it. 


Reconstruction by topologies involves determining which points are active and 
inactive (collided with a continent or subducted at a trench) for each reconstruction 
time step. This is done using `PlateReconstruction`s `ReconstructByTopologies` object.
If an ocean point with a certain velocity on one plate ID transitions into another 
rigid plate ID at another timestep (with another velocity), the velocity difference 
between both plates is calculated. The point may have subducted/collided with a continent 
if this velocity difference is higher than a specified velocity threshold (which can be
controlled with ``). To ascertain whether the point should be deactivated, a displacement
test is conducted. If the proximity of the point's previous time position to the polygon 
boundary it is approaching is higher than a set distance threshold, then the point is 
far enough away from the boundary that it cannot be subducted or consumed by it, and 
hence the point is still active. Otherwise, it is deemed inactive and deleted from the 
ocean basin mesh. 

With each reconstruction time step, points from mid-ocean ridges (which have more 
accurate spreading rates and attributed valid times) will spread across the ocean
floor. Eventually, points will be pushed into continental boundaries or subduction 
zones, where they are deleted. Ideally, all initial ocean points (from the Stripy 
icosahedral mesh) should be deleted over time. However, not all will be deleted - 
such points typically reside near continental boundaries. This happens if the 
emerged ridge points do not spread far enough to "phase out" these points at 
collision regions - likely due to insufficient reconstruction detail. These 
undeleted points form artefacts of anomalously high seafloor age that append 
over the reconstruction time range. 

Once reconstruction by topologies determines the ocean basin snapshot per timestep,
a data frame of all longitudes, latitudes, seafloor ages, spreading rates and any other
attributed z values will be written to a gridding input file per timestep.

Each active longitude, latitude and chosen z value (identified by a gridding input file
column index integer, i.e. `2` is seafloor age) is gridded using nearest-neighbour
interpolation and written to a netCDF4 format.

Classes
-------
* SeafloorGrid

"""
import pygplates
import numpy as np
from pydoc import plain
import stripy 
import ptt
import multiprocessing
import glob
import os
from skimage import measure
from scipy.interpolate import griddata
import pandas as pd

from . import reconstruction
from . import plot
from . import grids
from . import tools
from ptt import separate_ridge_transform_segments

# -------------------------------------------------------------------------
# Auxiliary functions for SeafloorGrid

def create_icosahedral_mesh(refinement_levels):
    """ Define a global point mesh with Stripy's
    [icosahedral triangulated mesh][https://github.com/underworldcode/stripy/blob/294354c00dd72e085a018e69c345d9353c6fafef/stripy/spherical_meshes.py#L27] 
    and turn all mesh domains into pyGPlates MultiPointOnSphere types.

    This global mesh will be masked with a set of continental or COB terrane
    polygons to define the ocean basin at a given reconstruction time.
    The `refinement_levels` integer is proportional to the resolution of the
    mesh and the ocean/continent boundary.

    Parameters
    ----------
    refinement_levels : int
        Refine the number of points in the triangulation. The larger the 
        refinement level, the sharper the ocean basin resolution.

    Returns
    -------
    multi_point : instance of <pygplates.MultiPointOnSphere>
        The longitues and latitudes that make up the icosahedral ocean mesh
        collated into a MultiPointOnSphere object. 
    icosahedral_global_mesh : instance of <stripy.spherical_meshes.icosahedral_mesh>
        The original global icosahedral triangulated mesh. 
    """

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
    """ Ensure COB terrane/continental polygon geometries are polygons 
    with reconstruction plate IDs and valid times.

    Notes
    -----
    This step must be done so that the initial set of ocean basin points 
    (the Stripy icosahedral mesh) can be partitioned into plates using
    each reconstruction plate ID for the given plate `model`. 

    This allows for an oceanic point-in-continental 
    polygon query for every identified plate ID. See documentation for
    `point_in_polygon_routine` for more details.

    `ensure_polygon_geometry` works as follows:
    COB terrane/continental polygons are assumed to have been reconstructed 
    already in `reconstructed_polygons` (a list of 
    type <pygplates.ReconstructedFeatureGeometry>). The list contents are 
    turned into a <pygplates.FeatureCollection> to be ascribed a 
    `PolygonOnSphere` geometry, a reconstruction plate ID, and a valid time. 
    Once finished, this feature collection is turned back into a list of 
    instance <pygplates.ReconstructedFeatureGeometry> and returned. 

    This revert must be completed for compatibility with the subsequent 
    point-in-polygon routine. 

    Parameters
    ----------
    reconstructed_polygons : list of instance <pygplates.ReconstructedFeatureGeometry>
        If used in `SeafloorGrid`, these are automatically obtained from the 
        `PlotTopologies.continents` attribute (the reconstructed continental 
        polygons at the current reconstruction time).

    rotation_model : instance of <pygplates.RotationModel>
        A parameter for turning the <pygplates.FeatureCollection> back into a
        list of instance <pygplates.ReconstructedFeatureGeometry> for 
        compatibility with the point-in-polygon routine. 

    """
    continent_FeatCol = []
    # self._PlotTopologies_object.continents
    for n in reconstructed_polygons:
        continent_FeatCol.append(n.get_feature())

    polygon_feats = pygplates.FeatureCollection(continent_FeatCol)

    # From GPRM's force_polygon_geometries(); set feature attributes
    # like valid times and plate IDs to each masking polygon
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
    """ Perform Plate Tectonic Tools' point in polygon routine to partition
    points in a `multi_point` MultiPointOnSphere feature based on whether
    they are inside or outside the polygons in `COB_polygons`.

    Notes
    -----
    Assuming the `COB_polygons` have passed through `ensure_polygon_geometry`,
    each polygon should have a plate ID assigned to it.

    This PIP routine serves two purposes for `SeafloorGrid`:

    1) It identifies continental regions in the icosahedral global mesh 
    MultiPointOnSphere feature and 'erases' in-continent oceanic points
    for the construction of a continental mask at each timestep;

    2) It identifies oceanic points in the icosahedral global mesh. 
    These points will be passed to a function that calculates each point's
    proximity to its nearest MOR segment (if any) within the polygonal domain 
    of its allocated plate ID. Each distance is divided by half the 
    `initial_ocean_mean_spreading_rate` (an attribute of `SeafloorGrids`) to 
    determine a simplified seafloor age for each point. 

    Number 2) only happens once at the start of the gridding process to 
    momentarily fill the gridding region with initial ocean points that have
    set ages (albeit not from a plate model file). After multiple time steps 
    of reconstruction, the ocean basin will be filled with new points (with 
    plate-model prescribed ages) that emerge from ridge topologies. 


    Returns
    -------
    pygplates.MultiPointOnSphere(points_in_arr) : instance <pygplates.MultiPointOnSphere>
        Point features that are within COB terrane polygons.
    pygplates.MultiPointOnSphere(points_out_arr) : instance <pygplates.MultiPointOnSphere>
        Point features that are outside COB terrane polygons.
    zvals : list
        A binary list. If an entry is == 0, its corresponing point in the 
        MultiPointOnSphere object is on the ocean. If == 1, the point is 
        in the COB terrane polygon.
    """

    # Collect reconstructed geometries and features of continental polygons
    # Both are needed for PTT's PIP routine.
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
        # If inside, set zval of masking grid to True
        if polygon_feature_list:
            points_in_arr.append(multi_point[point_index])
            zvals.append(1)
        # If outside, set zval of masking grid to False
        else:
            points_out_arr.append(multi_point[point_index])
            zvals.append(0)

    return pygplates.MultiPointOnSphere(points_in_arr), pygplates.MultiPointOnSphere(points_out_arr), zvals


def _extract_point_feature_attributes_for_rbt(ocean_basin, mor_all_times, time_arr):
    """ Collect all properties (like latitude, longitude, 
    geometry, shapefile attributes, and valid times) of features in the
    pygplates.MultiPointOnSphere object containing the 
    initial ocean basin seed points, and a list of all tessellated 
    points from mid-ocean ridges resolved for all times in `time_arr`.

    Both the initial points and all-time MOR points form the entirety
    of ocean basin points through `time_arr` to reconstruct. 

    TO DO: Expand this for general shapefile attribute allocation, or
    use pandas data frames instead?
    """
    # Attributes to extract
    active_points = []  
    appearance_time = []  
    birth_lat = []  # latitude_of_crust_formation
    prev_lat = []
    prev_lon = []
    spreading_rates = []
    #p_indices = []

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
            #p_index = int(feature.get_shapefile_attribute("INDEX"))
        except Exception:
            spreading_rate = np.nan
        spreading_rates.append(spreading_rate)
        #p_indices.append(p_index)
    del ocean_basin

    # Extract all MOR properties per timestep
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


# TO-DO: Replace this with a pandas dataframe per ocean point? Use unique point ID to track
# points through time. 
def set_shapefile_attribute(
    pygplates_feature, 
    shapefile_attribute, 
    shapefile_attribute_name, 
    overwrite_attribute=False
):
    """ Allocate a z-value to a point using its pygplates shapefile attribute
    dictionary.
    """
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

    Parameters
    ----------
    PlateReconstruction_object : instance of <gplately.PlateReconstruction>
        A GPlately PlateReconstruction object with a <pygplates.RotationModel> and
        a <pygplates.FeatureCollection> containing topology features. 
    PlotTopologies_object : instance of <gplately.PlotTopologies>
        A GPlately PlotTopologies object with a continental polygon or COB terrane
        polygon file to mask grids with.
    max_time : float
        The maximum time for age gridding.
    min_time : float
        The minimum time for age gridding.
    ridge_time_step : float
        The delta time for resolving ridges (and thus age gridding).
    save_directory : str, default None'
        The top-level directory to save all outputs to.
    file_collection : str, default None
        A string to identify the plate model used (will be automated later).
    refinement_levels : int, default 5
        Control the number of points in the icosahedral mesh (higher integer
        means higher resolution of continent masks).
    ridge_sampling : float, default 0.5
        Spatial resolution for points that emerge from ridges.
    resX : int, default 2000
        X resolution of the regular grid to interpolate the icosahedral meshed points onto.
    resY : int, default 1000
        Y resolution of the regular grid to interpolate the icosahedral meshed points onto.
    subduction_collision_parameters : len-2 tuple of float, default (5.0, 10.0)
        A 2-tuple of (threshold velocity delta in kms/my, threshold distance to boundary 
        per My in kms/my)
    initial_ocean_mean_spreading_rate : float, default 75.
        A spreading rate to uniformly allocate to points that define the initial ocean 
        basin. These points will have inaccurate ages, but most of them will be phased
        out after points with plate-model prescribed ages emerge from ridges and spread 
        to push them towards collision boundaries (where they are deleted).
    checkpoint_MORS_delta_time : float, default 5.
        A delta time at which to checkpoint MOR seedpoints by writing them into gpmlz
        format. If ever gridding preparation is interrupted, it will resume seed pointing
        and continental mask building at checkpoint_MORS_delta_time + 1. Ma.  
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
        ridge_sampling=0.5,
        resX = 2000,
        resY = 1000,
        subduction_collision_parameters = (5.0, 10.0),
        initial_ocean_mean_spreading_rate = 75.,
        checkpoint_MORS_delta_time = 5.
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

        # Gridding parameters
        self.resX = resX
        self.resY = resY

        # Temporal parameters
        self._max_time = max_time
        self.min_time = min_time
        self.ridge_time_step = ridge_time_step
        self.time_array = np.arange(self._max_time, self.min_time-0.1, -self.ridge_time_step)
        self.checkpoint_MORS_delta_time = checkpoint_MORS_delta_time

        # If PlotTopologies' time attribute is not equal to the maximum time in the 
        # seafloor grid reconstruction tree, make it equal. This will ensure the time
        # for continental masking is consistent.
        if self._PlotTopologies_object.time != self._max_time:
            self._PlotTopologies_object.time = self._max_time

        # Essential features and meshes for the SeafloorGrid
        self.continental_polygons = ensure_polygon_geometry(
            self._PlotTopologies_object.continents,
            self.rotation_model,
            self._max_time
        )
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

        The ocean mesh starts off as a global-spanning Stripy icosahedral mesh. 
        `create_initial_ocean_seed_points` passes the automatically-resolved-to-current-time
        continental polygons from the `PlotTopologies_object`'s `continents` attribute 
        (which can be from a COB terrane file or a continental polygon file) into 
        Plate Tectonic Tools' point-in-polygon routine. It identifies ocean basin points 
        that lie:
        * outside the polygons (for the ocean basin point domain)
        * inside the polygons (for the continental mask)

        Points from the mesh outside the continental polygons make up the ocean basin seed 
        point mesh. The masked mesh is outputted as a compressed GPML (GPMLZ) file with 
        the filename: "ocean_basin_seed_points_{}Ma.gpmlz" if a `save_directory` is passed.
        Otherwise, the mesh is returned as a pyGPlates FeatureCollection object.

        Notes
        ----- 
        This point mesh represents ocean basin seafloor that was produced
        before `SeafloorGrid.max_time`, and thus has unknown properties like valid
        time and spreading rate. As time passes, the plate reconstruction model sees 
        points emerging from MORs. These new points spread to occupy the ocean basins, 
        moving the initial filler points closer to subduction zones and continental 
        polygons with which they can collide. If a collision is detected by 
        `PlateReconstruction`s `ReconstructByTopologies` object, these points are deleted. 

        Ideally, if a reconstruction tree spans a large time range, **all** initial mesh 
        points would collide with a continent or be subducted, leaving behind a mesh of 
        well-defined MOR-emerged ocean basin points that data can be attributed to. 
        However, some of these initial points situated close to contiental boundaries are 
        retained through time - these form point artefacts with anomalously high ages. Even 
        deep-time plate models (e.g. 1 Ga) will have these artefacts - removing them would 
        require more detail to be added to the reconstruction model.

        Returns
        -------
        ocean_basin_point_mesh : pygplates.FeatureCollection of pygplates.MultiPointOnSphere
            A feature collection of point objects on the ocean basin.
        """
        print("Generating global point mesh...")

        # Ensure COB terranes at max time have reconstruction IDs and valid times
        COB_polygons = ensure_polygon_geometry(
            self._PlotTopologies_object.continents,
            self.rotation_model,
            self._max_time) 

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
        print("Partitioning global mesh by PlotTopologies' continental plate polygons...")

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

        # We need resolved topologies at the `max_time` to pass into the proximity
        # function
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
            ocean_points,
        )

        # Divide spreading rate by 2 to use half the mean spreading rate
        pAge = np.array(pZ) / (self.initial_ocean_mean_spreading_rate/2.)

        initial_ocean_point_features = []
        initial_ocean_multipoints = []

        for point in zip(pX,pY,pAge):

            point_feature = pygplates.Feature()
            point_feature.set_geometry(pygplates.PointOnSphere(point[1], point[0]))

            # Add 'time' to the age at the time of computation, to get the valid time in Ma
            point_feature.set_valid_time(point[2]+self._max_time, -1)

            # For now: custom zvals are added as shapefile attributes - will attempt pandas data frames
            point_feature = set_shapefile_attribute(point_feature, self.initial_ocean_mean_spreading_rate, "SPREADING_RATE")  # Seems like static data
            initial_ocean_point_features.append(point_feature)
            initial_ocean_multipoints.append(point_feature.get_geometry())

        # print(initial_ocean_point_features)
        multi_point_feature = pygplates.MultiPointOnSphere(initial_ocean_multipoints)

        if self.save_directory:
            if self.file_collection is not None:
                full_directory = "{}/{}_ocean_basin_seed_points_{}Ma.gpmlz".format(
                    self.save_directory,
                    self.file_collection, 
                    self._max_time
                )
            else:
                full_directory = "{}/ocean_basin_seed_points_{}Ma.gpmlz".format(
                    self.save_directory, 
                    self._max_time
                )
            pygplates.FeatureCollection(initial_ocean_point_features).write(full_directory)

        return pygplates.FeatureCollection(initial_ocean_point_features), multi_point_feature


    def _checkpoint_MOR(self, mor_points, time):
        """Write MOR seed points to GPMLZ format every `delta_time` to failsafe gridding 
        preparation if it is interrupted. 
        """
        if (time % self.checkpoint_MORS_delta_time == 0):
            mor_points.write('{}/MOR_plus_one_points_{:0.2f}.gpmlz'.format(
                self.save_directory, 
                time
                )
            )
            print("Checkpointed MOR seedpoints at {} Ma!".format(time))


    def get_mid_ocean_ridge_seedpoints(self, time):
        """ Resolve mid-ocean ridges to the current `time`, divide them into their shared 
        sub-segments, and further divide them into points. Rotate these points to the left
        and right of the ridge using their stage rotation so that they spread from the ridge.

        Z-value allocation to each point is done here. In future, a function (like
        the spreading rate function) to calculate general z-data will be an input parameter.

        Notes
        -----
        This assumes that points spread from ridges symmetrically, with the exception of 
        large ridge jumps at successive timesteps. Therefore, z-values allocated to ridge-emerging
        points will appear symmetrical until changes in spreading ridge geometries create 
        asymmetries.

        In future, this will have a checkpoint save feature so that execution
        (which occurs during preparation for ReconstructByTopologies and can take several hours) 
        can be safeguarded against run interruptions. 

        Parameters
        ----------
        time : float
            The time at which to resolve ridge points and stage-rotate them off the ridge.

        Returns
        -------
        mor_point_features : FeatureCollection
            All ridge seed points that have emerged from all ridge topologies at `time`. 
            These points have spread by being slightly rotated away from
            ridge locations at `time`.

        References
        ----------
        get_mid_ocean_ridge_seedpoints() has been adapted from 
        https://github.com/siwill22/agegrid-0.1/blob/master/automatic_age_grid_seeding.py#L117.
        """

        # Topology features from `PlotTopologies`. 
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
        #point_indices = []

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
                
                subsegment_index = []
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
                        # Get the spreading rates for all points in this sub segment
                        spreading_rates, subsegment_index = tools.calculate_spreading_rates(
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

                    # Loop through all but the 1st and last points in the current sub segment
                    for point, rate in zip(
                        mor_points.get_points()[1:-1],
                        spreading_rates[1:-1],
                    ):
                        # Add the point "twice" to the main shifted_mor_points list; once for a L-side
                        # spread, another for a R-side spread. Then add the same spreading rate twice
                        # to the list - this therefore assumes spreading rate is symmetric. 
                        shifted_mor_points.append(rotate_slightly_off_mor_one_way * point)
                        shifted_mor_points.append(rotate_slightly_off_mor_opposite_way * point)
                        point_spreading_rates.extend([rate] * 2)
                        #point_indices.extend(subsegment_index)

        # Summarising get_isochrons_for_ridge_snapshot;
        # Write out the ridge point born at 'ridge_time' but their position at 'ridge_time - time_step'.
        mor_point_features = []
        for i, (curr_point, spreading_rate) in enumerate(zip(shifted_mor_points, point_spreading_rates)):
            feature = pygplates.Feature()
            feature.set_geometry(curr_point)
            feature.set_valid_time(time, -999)  # delete - time_step
            #feature.set_name(str(spreading_rate))
            feature = set_shapefile_attribute(feature, spreading_rate, "SPREADING_RATE")  # make spreading rate a shapefile attribute
            #feature = set_shapefile_attribute(feature, point_indices[i], "INDEX")
            mor_point_features.append(feature)

        mor_points = pygplates.FeatureCollection(mor_point_features)

        # Checkpoint MORs if `time` is a multiple of `checkpoint_MORS`
        self._checkpoint_MOR(mor_points, time)

        return mor_points



    def create_continental_mask(self, time):
        """ Create a continental mask to define the ocean basin, as well as
        use as a set of continental collision boundaries in 
        `ReconstructByTopologies`. 

        If `ReconstructByTopologies` identifies a continental collision 
        between oceanic points and the boundaries of this continental
        mask at `time`, those points are deleted at `time`. 

        The continental mask is also saved to "/continent_mask_{}Ma.nc" as a 
        compressed netCDF4 file if a `save_directory` is passed. Otherwise, 
        the final grid is returned as a NumPy ndarray object.

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

        # A regular grid to interpolate the icosahedral mesh onto
        extent_globe = np.radians([-180,180,-90,90])
        grid_lon = np.linspace(extent_globe[0], extent_globe[1], self.resX)
        grid_lat = np.linspace(extent_globe[2], extent_globe[3], self.resY)

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

        Returns lists of all attributes for the initial ocean point mesh and
        all ridge points for all times in the reconstruction time array.
        """

        # INITIAL OCEAN SEED POINT MESH ----------------------------------------------------
        initial_ocean_seed_points, initial_ocean_seed_points_mp = self.create_initial_ocean_seed_points()
        print("Finished building initial_ocean_seed_points!")

        # MOR SEED POINTS AND CONTINENTAL MASKS --------------------------------------------
        time_array = np.arange(self._max_time, self.min_time-1, -self.ridge_time_step)
        all_mor_features = []
        all_continental_masks = []

        # Build all continental masks and spreading ridge points (with z values)
        for time in self.time_array:
            all_mor_features.append(
                self.get_mid_ocean_ridge_seedpoints(time)
            )
            print("Finished building MOR seedpoints at {} Ma!".format(time))
            
            all_continental_masks.append(
                self.create_continental_mask(time)
            )
            print("Finished building a continental mask at {} Ma!".format(time))

        # ALL-TIME POINTS AND FEATURES -----------------------------------------------------
        # Extract all feature attributes for all reconstruction times into lists
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

        All active points' latitudes, longitues, seafloor ages, spreading rates and all 
        other general z-values are saved to a gridding input file (.npz).
        """
        print("Preparing all initial files...")

        # Obtain all info from the ocean seed points and all MOR points through time, store in
        # arrays
        active_points, appearance_time, birth_lat, \
        prev_lat, prev_lon, spreading_rates = self.prepare_for_reconstruction_by_topologies()

        ####  Begin reconstruction by topology process:
        # Indices for all points (`active_points`) that have existed from `max_time` to `min_time`.
        point_id = range(len(active_points))

        # Specify the default collision detection region as subduction zones
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

        # Loop over the reconstruction times until the end of the reconstruction time span, or until
        # all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        reconstruction_data = []
        while True:
            print('Reconstruct by topologies: working on time {:0.2f} Ma'.format(
                topology_reconstruction.get_current_time())
            )

            # NOTE:
            # topology_reconstruction.get_active_current_points() and topology_reconstruction.get_all_current_points()
            # are different. The former is a subset of the latter, and it represents all points at the timestep that
            # have not collided with a continental or subduction boundary. The remainders in the latter are inactive 
            # (NoneType) points, which represent the collided points.

            # We need to access active point data from topology_reconstruction.get_all_current_points() because it has 
            # the same length as the list of all initial ocean points and MOR seed points that have ever emerged from 
            # spreading ridge topologies through `max_time` to `min_time`. Therefore, it protects the time and space 
            # order in which all MOR points through time were seeded by pyGPlates. At any given timestep, not all these 
            # points will be active, but their indices are retained. Thus, z value allocation, point latitudes and 
            # longitudes of active points will be correctly indexed if taking it from 
            # topology_reconstruction.get_all_current_points().
            curr_points = topology_reconstruction.get_active_current_points()
            curr_points_including_inactive = topology_reconstruction.get_all_current_points()


            # Collect latitudes and longitudes of currently ACTIVE points in the ocean basin
            curr_lat_lon_points = [point.to_lat_lon() for point in curr_points]

            if curr_lat_lon_points:

                # Get the number of active points at this timestep. 
                num_current_points = len(curr_points)

                # ndarray to fill with active point lats, lons and zvalues
                # FOR NOW, the number of gridding input columns is 6:
                # 0 = longitude
                # 1 = latitude
                # 2 = seafloor age
                # 3 = spreading rate
                # 4 = birth latitude snapshot
                # 5 = point id
                # In future, we will add to 6 by len(input params of zvals to SeafloorGrid)
                gridding_input_data = np.empty((num_current_points, 6))

                # Lons and lats are first and second columns of the ndarray respectively
                gridding_input_data[:, 1], gridding_input_data[:,0] = zip(*curr_lat_lon_points)

                # NOTE: We need a single index to access data from curr_points_including_inactive AND allocate 
                # this data to an ndarray with a number of rows equal to num_current_points. This index will 
                # append +1 after each loop through curr_points_including_inactive. 
                i = 0

                # Get indices and points of all points at `time`, both active and inactive (which are NoneType points that
                # have undergone continental collision or subduction at `time`). 
                for point_index,current_point in enumerate(curr_points_including_inactive):

                    # Look at all active points (these have not collided with a continent or trench)
                    if current_point is not None:

                        # Seafloor age
                        gridding_input_data[i, 2] = appearance_time[point_index] - topology_reconstruction.get_current_time()
                        # Birth latitude (snapshot)
                        gridding_input_data[i, 4] = birth_lat[point_index]
                        # Point ID (snapshot)
                        gridding_input_data[i, 5] = point_id[point_index] # The ID of a corresponding point from the original list of all MOR-resolved points
                        # Spreading rate
                        gridding_input_data[i, 3] = spreading_rates[point_index]

                        i += 1

                if self.file_collection is not None:
                    np.savez_compressed('{:s}/{}_gridding_input_{:0.1f}Ma'.format(self.save_directory, self.file_collection, topology_reconstruction.get_current_time()), 
                        CURRENT_LONGITUDES = gridding_input_data[:,0],
                        CURRENT_LATITUDES = gridding_input_data[:,1],
                        SEAFLOOR_AGE = gridding_input_data[:,2],
                        SPREADING_RATE_SNAPSHOT = gridding_input_data[:,3],
                        BIRTH_LAT_SNAPSHOT = gridding_input_data[:,4],
                        POINT_ID_SNAPSHOT = gridding_input_data[:,5]
                    )
                else:
                    np.savez_compressed('{:s}/gridding_input_{:0.1f}Ma.csv'.format(self.save_directory, topology_reconstruction.get_current_time()), 
                        CURRENT_LONGITUDES = gridding_input_data[:,0],
                        CURRENT_LATITUDES = gridding_input_data[:,1],
                        SEAFLOOR_AGE = gridding_input_data[:,2],
                        SPREADING_RATE_SNAPSHOT = gridding_input_data[:,3],
                        BIRTH_LAT_SNAPSHOT = gridding_input_data[:,4],
                        POINT_ID_SNAPSHOT = gridding_input_data[:,5]
                    )

            if not topology_reconstruction.reconstruct_to_next_time():
                break

        print('Reconstruction done for {}!'.format(topology_reconstruction.get_current_time()))
        # return reconstruction_data


    def lat_lon_z_to_netCDF(self, gridding_data_index=2, time_arr=None):
        """ Write ndarray grids containing z values from index `gridding_data_index`
        for a given time range in `time_arr`.

        By default, `gridding_data_index` is 2, which corresponds to seafloor age.
        `gridding_data_index` = 3 corresponds to spreading rate.

        Saves all grids to compressed netCDF format in the attributed directory. Grids
        can be read into ndarray format using `gplately.grids.read_netcdf_grid()`.

        Parameters
        ----------
        gridding_data_index : int, default 2
            An index for a column in the ReconstructByTopologies gridding input files.
            * 2 = seafloor age
            * 3 = spreading rate
            * 6 and onwards will be the general data.
        time_arr : list of float, default None
            A time range to turn lons, lats and z-values into netCDF4 grids. If not provided,
            `time_arr` defaults to the full `time_array` provided to `SeafloorGrids`. 
        """

        # User can put any time array within SeafloorGrid bounds, but if none
        # is provided, it defaults to the attributed time array
        if time_arr is None:
            time_arr = self.time_array
        for time in time_arr:
            # Read the gridding input made by ReconstructByTopologies:
            if self.file_collection is not None:
                gridding_input = '{:s}/{}_gridding_input_{:0.1f}Ma.npz'.format(
                    self.save_directory, 
                    self.file_collection, 
                    time
                )
            else:
                gridding_input = '{:s}/gridding_input_{:0.1f}Ma.npz'.format(self.save_directory, time)

            # Use pandas to load in lons, lats and z values from npz files
            npz = np.load(gridding_input)
            curr_data = pd.DataFrame.from_dict({item: npz[item] for item in npz.files}, orient='index').transpose()

            # Drop duplicate latitudes and longitudes
            unique_data = curr_data.drop_duplicates(subset=["CURRENT_LONGITUDES", "CURRENT_LATITUDES"])

            # Acquire lons, lats and zvalues for each time
            lons = unique_data["CURRENT_LONGITUDES"].to_list()
            lats = unique_data["CURRENT_LATITUDES"].to_list()
            zdata = np.array(unique_data[unique_data.columns[gridding_data_index]].to_list())
            #zdata = np.where(zdata > 375, float("nan"), zdata), to deal with vmax in the future
            zdata = np.nan_to_num(zdata)
            
            # Create a regular grid on which to interpolate lats, lons and zdata
            extent_globe = [-180,180,-90,90]
            grid_lon = np.linspace(extent_globe[0], extent_globe[1], self.resX)
            grid_lat = np.linspace(extent_globe[2], extent_globe[3], self.resY)
            X, Y = np.meshgrid(grid_lon, grid_lat)

            # Interpolate lons, lats and zvals over a regular grid using nearest
            # neighbour interpolation
            Z = griddata((lons, lats), zdata, (X, Y), method='nearest')
            
            # Access continental grids from the save directory
            if self.save_directory is not None:
                if self.file_collection is not None:
                    full_directory = "{}/{}_continent_mask_{}Ma.nc".format(
                        self.save_directory, 
                        self.file_collection, 
                        time
                    )
                    grid_output_dir = "{}/{}_{}_grid_{}Ma.nc".format(
                        self.save_directory, 
                        self.file_collection, 
                        str(unique_data.columns[gridding_data_index]),
                        time
                    )
                else:
                    full_directory = "{}/continent_mask_{}Ma.nc".format(
                        self.save_directory, 
                        time
                    )
                    grid_output_dir = "{}/grid{}Ma.nc".format(
                        self.save_directory, 
                        time
                    )

            # Identify regions in the grid in the continental mask
            cont_mask = grids.Raster(
                    filename=str(full_directory))
            grd = cont_mask.interpolate(X, Y) > 0.5
            Z[grd] = np.nan

            grids.write_netcdf_grid(
                    grid_output_dir, 
                    Z, 
                    extent=[-180,180,-90,90]
                )
            print("netCDF grid for {} Ma complete!".format(time))
        return


