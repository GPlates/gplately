#
#    Copyright (C) 2024 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
Auxiliary functions for SeafloorGrid
"""

import math
import os

import numpy as np
import pygplates

from .. import ptt as _ptt


def create_icosahedral_mesh(refinement_levels):
    """Define a global point mesh with Stripy's
    [icosahedral triangulated mesh](https://github.com/underworldcode/stripy/blob/294354c00dd72e085a018e69c345d9353c6fafef/stripy/spherical_meshes.py#L27)
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
    import stripy

    # Create the ocean basin mesh using Stripy's icosahedral spherical mesh
    icosahedral_global_mesh = stripy.spherical_meshes.icosahedral_mesh(
        refinement_levels, include_face_points=False, trisection=False, tree=False
    )
    # Get lons and lats of mesh, and turn them into a MultiPointOnSphere
    lats_arr = np.rad2deg(icosahedral_global_mesh.lats)
    lons_arr = np.rad2deg(icosahedral_global_mesh.lons)
    multi_point = pygplates.MultiPointOnSphere(zip(lats_arr, lons_arr))

    return multi_point, icosahedral_global_mesh


def ensure_polygon_geometry(reconstructed_polygons, rotation_model, time):
    """Ensure COB terrane/continental polygon geometries are polygons
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
            polygon = pygplates.Feature(feature.get_feature_type())
            polygon.set_geometry(pygplates.PolygonOnSphere(geom))
            polygon.set_reconstruction_plate_id(feature.get_reconstruction_plate_id())
            # Avoid features in COBTerranes with invalid time
            if feature.get_valid_time()[0] >= feature.get_valid_time()[1]:
                polygon.set_valid_time(
                    feature.get_valid_time()[0], feature.get_valid_time()[1]
                )
                polygons.append(polygon)
    cobter_polygon_features = pygplates.FeatureCollection(polygons)

    # Turn the feature collection back into ReconstructedFeatureGeometry
    # objects otherwise it will not work with PIP
    reconstructed_cobter_polygons = []
    pygplates.reconstruct(
        cobter_polygon_features, rotation_model, reconstructed_cobter_polygons, time
    )
    return reconstructed_cobter_polygons


def point_in_polygon_routine(multi_point, COB_polygons):
    """Perform Plate Tectonic Tools' point in polygon routine to partition
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
    # Convert MultiPointOnSphere to array of PointOnSphere
    multi_point = np.array(multi_point.get_points(), dtype="object")

    # Collect reconstructed geometries of continental polygons
    polygons = np.empty(len(COB_polygons), dtype="object")
    for ind, i in enumerate(COB_polygons):
        if isinstance(i, pygplates.ReconstructedFeatureGeometry):
            geom = i.get_reconstructed_geometry()
        elif isinstance(i, pygplates.GeometryOnSphere):
            geom = i
        else:  # e.g. ndarray of coordinates
            geom = pygplates.PolygonOnSphere(i)
        polygons[ind] = geom
    proxies = np.ones(polygons.size)

    pip_result = _ptt.utils.points_in_polygons.find_polygons(
        multi_point, polygons, proxies, all_polygons=False
    )  # 1 for points in polygons, None for points outside
    zvals = np.array(
        pip_result,
        dtype="float",
    ).ravel()
    zvals[np.isnan(zvals)] = 0.0
    zvals = zvals.astype("int")
    points_in_arr = multi_point[zvals == 1]
    points_out_arr = multi_point[zvals != 1]

    return (
        pygplates.MultiPointOnSphere(points_in_arr),
        pygplates.MultiPointOnSphere(points_out_arr),
        zvals,
    )

# FROM RECONSTRUCT_BY_TOPOLOGIES.PY
class DefaultCollision(object):
    """
    Default collision detection function class (the function is the '__call__' method).
    """

    DEFAULT_GLOBAL_COLLISION_PARAMETERS = (7.0, 10.0)
    """
    Default collision parameters for all feature types.

    This is a 2-tuple of (threshold velocity delta in kms/my, threshold distance to boundary per My in kms/my):
    Here we default to the same constants used internally in GPlates 2.0 (ie, 7.0 and 10.0).
    """

    def __init__(
        self,
        global_collision_parameters=DEFAULT_GLOBAL_COLLISION_PARAMETERS,
        feature_specific_collision_parameters=None,
    ):
        """
        global_collision_parameters: The collision parameters to use for any feature type not specified in 'feature_specific_collision_parameters'.
                                     Should be a 2-tuple of (threshold velocity delta in kms/my, threshold distance to boundary per My in kms/my).
                                     The first threshold parameter means:
                                        A point that transitions from one plate to another can disappear if the change in velocity exceeds this threshold.
                                     The second threshold parameter means:
                                        Only those transitioning points exceeding the threshold velocity delta and that are close enough to a plate boundary can disappear.
                                        The distance is proportional to the relative velocity (change in velocity), plus a constant offset based on the threshold distance to boundary
                                        to account for plate boundaries that change shape significantly from one time step to the next
                                        (note that some boundaries are meant to do this and others are a result of digitisation).
                                        The actual distance threshold used is (threshold_distance_to_boundary + relative_velocity) * time_interval
                                     Defaults to parameters used in GPlates 2.0, if not specified.

        feature_specific_collision_parameters: Optional sequence of collision parameters specific to feature types.
                                               If specified then should be a sequence of 2-tuples, with each 2-tuple specifying (feature_type, collision_parameters).
                                               And where each 'collision_parameters' is a 2-tuple of (threshold velocity delta in kms/my, threshold distance to boundary per My in kms/my).
                                                   See 'global_collision_parameters' for details on these thresholds.
                                               Any feature type not specified here defaults to using 'global_collision_parameters'.
        """

        # Convert list of (feature_type, collision_parameters) tuples to a dictionary.
        if feature_specific_collision_parameters:
            self.feature_specific_collision_parameters = dict(
                feature_specific_collision_parameters
            )
        else:
            self.feature_specific_collision_parameters = dict()
        # Fallback for any feature type not specified in the optional feature-specific list.
        self.global_collision_parameters = global_collision_parameters

        # Used to improve performance by caching velocity stage rotations in a dict (for a specific reconstruction time).
        self.velocity_stage_rotation_dict = {}
        self.velocity_stage_rotation_time = None

    def __call__(
        self,
        rotation_model,
        time,
        reconstruction_time_interval,
        prev_point,
        curr_point,
        prev_topology_plate_id,
        prev_resolved_plate_boundary,
        curr_topology_plate_id,
        curr_resolved_plate_boundary,
    ):
        """
        Returns True if a collision was detected.

        If transitioning from a rigid plate to another rigid plate with a different plate ID then
        calculate the difference in velocities and continue testing as follows
        (otherwise, if there's no transition, then the point is still active)...

        If the velocity difference is below a threshold then we assume the previous plate was split,
        or two plates joined. In this case the point has not subducted (forward in time) or
        been consumed by a mid-ocean (backward in time) and hence is still active.

        If the velocity difference is large enough then we see if the distance of the *previous* position
        to the polygon boundary (of rigid plate containing it) exceeds a threshold.
        If the distance exceeds the threshold then the point is far enough away from the boundary that it
        cannot be subducted or consumed by it and hence the point is still active.
        However if the point is close enough then we assume the point was subducted/consumed
        (remember that the point switched plate IDs).
        Also note that the threshold distance increases according to the velocity difference to account for fast
        moving points (that would otherwise tunnel through the boundary and accrete onto the other plate).
        The reason for testing the distance from the *previous* point, and not from the *current* point, is:

          (i)  A topological boundary may *appear* near the current point (such as a plate split at the current time)
               and we don't want that split to consume the current point regardless of the velocity difference.
               It won't get consumed because the *previous* point was not near a boundary (because before split happened).
               If the velocity difference is large enough then it might cause the current point to transition to the
               adjacent split plate in the *next* time step (and that's when it should get consumed, not in the current time step).
               An example of this is a mid-ocean ridge suddenly appearing (going forward in time).

          (ii) A topological boundary may *disappear* near the current point (such as a plate merge at the current time)
               and we want that merge to consume the current point if the velocity difference is large enough.
               In this case the *previous* point is near a boundary (because before plate merged) and hence can be
               consumed (provided velocity difference is large enough). And since the boundary existed in the previous
               time step, it will affect position of the current point (and whether it gets consumed or not).
               An example of this is a mid-ocean ridge suddenly disappearing (going backward in time).

        ...note that items (i) and (ii) above apply both going forward and backward in time.
        """

        # See if a collision occurred.
        if (
            curr_topology_plate_id != prev_topology_plate_id
            and prev_topology_plate_id is not None
            and curr_topology_plate_id is not None
        ):
            #
            # Speed up by caching velocity stage rotations in a dict.
            #
            if time != self.velocity_stage_rotation_time:
                # We've just switched to a new time so clear the cache.
                #
                # We only cache stage rotations for a specific time.
                # We only really need to cache different plate IDs at the same 'time', so this avoids caching for all times
                # (which would also require including 'time' in the key) and using memory unnecessarily.
                self.velocity_stage_rotation_dict.clear()
                self.velocity_stage_rotation_time = time
            prev_location_velocity_stage_rotation = (
                self.velocity_stage_rotation_dict.get(prev_topology_plate_id)
            )
            if not prev_location_velocity_stage_rotation:
                prev_location_velocity_stage_rotation = rotation_model.get_rotation(
                    time + 1, prev_topology_plate_id, time
                )
                self.velocity_stage_rotation_dict[prev_topology_plate_id] = (
                    prev_location_velocity_stage_rotation
                )
            curr_location_velocity_stage_rotation = (
                self.velocity_stage_rotation_dict.get(curr_topology_plate_id)
            )
            if not curr_location_velocity_stage_rotation:
                curr_location_velocity_stage_rotation = rotation_model.get_rotation(
                    time + 1, curr_topology_plate_id, time
                )
                self.velocity_stage_rotation_dict[curr_topology_plate_id] = (
                    curr_location_velocity_stage_rotation
                )

            # Note that even though the current point is not inside the previous boundary (because different plate ID), we can still
            # calculate a velocity using its plate ID (because we really should use the same point in our velocity comparison).
            prev_location_velocity = pygplates.calculate_velocities(
                (curr_point,),
                prev_location_velocity_stage_rotation,
                1,
                pygplates.VelocityUnits.kms_per_my,
            )[0]
            curr_location_velocity = pygplates.calculate_velocities(
                (curr_point,),
                curr_location_velocity_stage_rotation,
                1,
                pygplates.VelocityUnits.kms_per_my,
            )[0]

            delta_velocity = curr_location_velocity - prev_location_velocity
            delta_velocity_magnitude = delta_velocity.get_magnitude()

            # If we have feature-specific collision parameters then iterate over the boundary sub-segments of the *previous* topological boundary
            # and test proximity to each sub-segment individually (with sub-segment feature type specific collision parameters).
            # Otherwise just test proximity to the entire boundary polygon using the global collision parameters.
            if self.feature_specific_collision_parameters:
                for (
                    prev_boundary_sub_segment
                ) in prev_resolved_plate_boundary.get_boundary_sub_segments():
                    # Use feature-specific collision parameters if found (falling back to global collision parameters).
                    (
                        threshold_velocity_delta,
                        threshold_distance_to_boundary_per_my,
                    ) = self.feature_specific_collision_parameters.get(
                        prev_boundary_sub_segment.get_feature().get_feature_type(),
                        # Default to global collision parameters if no collision parameters specified for sub-segment's feature type...
                        self.global_collision_parameters,
                    )

                    # Since each feature type could use different collision parameters we must use the current boundary sub-segment instead of the boundary polygon.
                    if self._detect_collision_using_collision_parameters(
                        reconstruction_time_interval,
                        delta_velocity_magnitude,
                        prev_point,
                        prev_boundary_sub_segment.get_resolved_geometry(),
                        threshold_velocity_delta,
                        threshold_distance_to_boundary_per_my,
                    ):
                        # Detected a collision.
                        return True
            else:
                # No feature-specific collision parameters so use global fallback.
                (
                    threshold_velocity_delta,
                    threshold_distance_to_boundary_per_my,
                ) = self.global_collision_parameters

                # Since all feature types use the same collision parameters we can use the boundary polygon instead of iterating over its sub-segments.
                if self._detect_collision_using_collision_parameters(
                    reconstruction_time_interval,
                    delta_velocity_magnitude,
                    prev_point,
                    prev_resolved_plate_boundary.get_resolved_boundary(),
                    threshold_velocity_delta,
                    threshold_distance_to_boundary_per_my,
                ):
                    # Detected a collision.
                    return True

        return False

    def _detect_collision_using_collision_parameters(
        self,
        reconstruction_time_interval,
        delta_velocity_magnitude,
        prev_point,
        prev_boundary_geometry,
        threshold_velocity_delta,
        threshold_distance_to_boundary_per_my,
    ):
        if delta_velocity_magnitude > threshold_velocity_delta:
            # Add the minimum distance threshold to the delta velocity threshold.
            # The delta velocity threshold only allows those points that are close enough to the boundary to reach
            # it given their current relative velocity.
            # The minimum distance threshold accounts for sudden changes in the shape of a plate boundary
            # which are no supposed to represent a new or shifted boundary but are just a result of the topology
            # builder/user digitising a new boundary line that differs noticeably from that of the previous time period.
            distance_threshold_radians = (
                (threshold_distance_to_boundary_per_my + delta_velocity_magnitude)
                * reconstruction_time_interval
                / pygplates.Earth.equatorial_radius_in_kms
            )
            distance_threshold_radians = min(distance_threshold_radians, math.pi)
            distance_threshold_radians = max(distance_threshold_radians, 0.0)

            distance = pygplates.GeometryOnSphere.distance(
                prev_point,
                prev_boundary_geometry,
                distance_threshold_radians=float(distance_threshold_radians),
            )
            if distance is not None:
                # Detected a collision.
                return True

        return False


DEFAULT_COLLISION = DefaultCollision()


class ContinentCollision(object):
    """
    Continental collision detection function class (the function is the '__call__' method).
    """

    def __init__(
        self,
        grd_output_dir,
        chain_collision_detection=DEFAULT_COLLISION,
        verbose=False,
    ):
        """
        grd_output_dir: The directory containing the continental grids.

        chain_collision_detection: Another collision detection class/function to reference if we find no collision.
                                   If None then no collision detection is chained. Defaults to the default collision detection.

        verbose: Print progress messages
        """

        self.grd_output_dir = grd_output_dir
        self.chain_collision_detection = chain_collision_detection
        self.verbose = verbose

        # Load a new grid each time the reconstruction time changes.
        self.grid_time = None

    @property
    def grid_time(self):
        return self._grid_time

    @grid_time.setter
    def grid_time(self, time):
        from ..grids import read_netcdf_grid

        if time is None:
            self._grid_time = time
        else:
            filename = "{:s}".format(self.grd_output_dir.format(time))
            if self.verbose:
                print(
                    "Points masked against grid: {0}".format(os.path.basename(filename))
                )
            gridZ, gridX, gridY = read_netcdf_grid(filename, return_grids=True)
            self.gridZ = gridZ
            self.ni, self.nj = gridZ.shape
            self.xmin = np.nanmin(gridX)
            self.xmax = np.nanmax(gridX)
            self.ymin = np.nanmin(gridY)
            self.ymax = np.nanmax(gridY)
            self._grid_time = float(time)

    def __call__(
        self,
        rotation_model,
        time,
        reconstruction_time_interval,
        prev_point,
        curr_point,
        prev_topology_plate_id,
        prev_resolved_plate_boundary,
        curr_topology_plate_id,
        curr_resolved_plate_boundary,
    ):
        """
        Returns True if a collision with a continent was detected, or returns result of
        chained collision detection if 'self.chain_collision_detection' is not None.
        """
        # Load the grid for the current time if encountering a new time.
        if time != self.grid_time:
            self.grid_time = time
            self.continent_deletion_count = 0

        # Sample mask grid, which is one over continents and zero over oceans.
        point_lat, point_lon = curr_point.to_lat_lon()
        point_i = (self.ni - 1) * ((point_lat - self.ymin) / (self.ymax - self.ymin))
        point_j = (self.nj - 1) * ((point_lon - self.xmin) / (self.xmax - self.xmin))
        point_i_uint = np.rint(point_i).astype(np.uint)
        point_j_uint = np.rint(point_j).astype(np.uint)
        try:
            mask_value = self.gridZ[point_i_uint, point_j_uint]
        except IndexError:
            point_i = np.clip(np.rint(point_i), 0, self.ni - 1).astype(np.int_)
            point_j = np.clip(np.rint(point_j), 0, self.nj - 1).astype(np.int_)
            mask_value = self.gridZ[point_i, point_j]
        if mask_value >= 0.5:
            # Detected a collision.
            self.continent_deletion_count += 1
            return True

        # We didn't find a collision, so ask the chained collision detection if it did (if we have anything chained).
        if self.chain_collision_detection:
            return self.chain_collision_detection(
                rotation_model,
                time,
                reconstruction_time_interval,
                prev_point,
                curr_point,
                prev_topology_plate_id,
                prev_resolved_plate_boundary,
                curr_topology_plate_id,
                curr_resolved_plate_boundary,
            )

        return False


def reconstruct_points(
    rotation_features_or_model,
    topology_features,
    reconstruction_begin_time,
    reconstruction_end_time,
    reconstruction_time_interval,
    points,
    point_begin_times=None,
    point_end_times=None,
    point_plate_ids=None,
    detect_collisions=DEFAULT_COLLISION,
):
    """
    Function to reconstruct points using the ReconstructByTopologies class below.

    For description of parameters see the ReconstructByTopologies class below.
    """

    topology_reconstruction = ReconstructByTopologies(
        rotation_features_or_model,
        topology_features,
        reconstruction_begin_time,
        reconstruction_end_time,
        reconstruction_time_interval,
        points,
        point_begin_times,
        point_end_times,
        point_plate_ids,
        detect_collisions,
    )

    return topology_reconstruction.reconstruct()

class ReconstructByTopologies(object):
    """
    Class to reconstruct geometries using topologies.

    Currently only points are supported.

    use_plate_partitioner: If True then use pygplates.PlatePartitioner to partition points,
                           otherwise use faster points_in_polygons.find_polygons().
    """

    use_plate_partitioner = False

    def __init__(
        self,
        rotation_features_or_model,
        topology_features,
        reconstruction_begin_time,
        reconstruction_end_time,
        reconstruction_time_interval,
        points,
        point_begin_times=None,
        point_end_times=None,
        point_plate_ids=None,
        detect_collisions=DEFAULT_COLLISION,
    ):
        """
        rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).

        topology_features: Topology feature collection(s), or list of features, or filename(s) or any combination of those.

        detect_collisions: Collision detection function, or None. Defaults to DEFAULT_COLLISION.
        """

        # Turn rotation data into a RotationModel (if not already).
        if not isinstance(rotation_features_or_model, pygplates.RotationModel):
            rotation_model = pygplates.RotationModel(rotation_features_or_model)
        else:
            rotation_model = rotation_features_or_model
        self.rotation_model = rotation_model

        # Turn topology data into a list of features (if not already).
        self.topology_features = pygplates.FeaturesFunctionArgument(
            topology_features
        ).get_features()

        # Set up an array of reconstruction times covering the reconstruction time span.
        self.reconstruction_begin_time = reconstruction_begin_time
        self.reconstruction_end_time = reconstruction_end_time
        if reconstruction_time_interval <= 0.0:
            raise ValueError("'reconstruction_time_interval' must be positive.")
        # Reconstruction can go forward or backward in time.
        if self.reconstruction_begin_time > self.reconstruction_end_time:
            self.reconstruction_time_step = -reconstruction_time_interval
        else:
            self.reconstruction_time_step = reconstruction_time_interval
        # Get number of times including end time if time span is a multiple of time step.
        # The '1' is because, for example, 2 time intervals is 3 times.
        # The '1e-6' deals with limited floating-point precision, eg, we want (3.0 - 0.0) / 1.0 to be 3.0 and not 2.999999 (which gets truncated to 2).
        self.num_times = 1 + int(
            math.floor(
                1e-6
                + float(self.reconstruction_end_time - self.reconstruction_begin_time)
                / self.reconstruction_time_step
            )
        )
        # It's possible the time step is larger than the time span, in which case we change it to equal the time span.
        # This guarantees there'll be at least one time step (which has two times; one at either end of interval).
        if self.num_times == 1:
            self.num_times = 2
            self.reconstruction_time_step = (
                self.reconstruction_end_time - self.reconstruction_begin_time
            )
        self.reconstruction_time_interval = math.fabs(self.reconstruction_time_step)

        self.last_time_index = self.num_times - 1

        self.points = points
        self.num_points = len(points)

        # Use the specified point begin times if provided (otherwise use 'inf').
        self.point_begin_times = point_begin_times
        if self.point_begin_times is None:
            self.point_begin_times = [float("inf")] * self.num_points
        elif len(self.point_begin_times) != self.num_points:
            raise ValueError(
                "Length of 'point_begin_times' must match length of 'points'."
            )

        # Use the specified point end times if provided (otherwise use '-inf').
        self.point_end_times = point_end_times
        if self.point_end_times is None:
            self.point_end_times = [float("-inf")] * self.num_points
        elif len(self.point_end_times) != self.num_points:
            raise ValueError(
                "Length of 'point_end_times' must match length of 'points'."
            )

        # Use the specified point plate IDs if provided (otherwise use '0').
        # These plate IDs are only used when a point falls outside all resolved topologies during a time step.
        self.point_plate_ids = point_plate_ids
        if self.point_plate_ids is None:
            self.point_plate_ids = [0] * self.num_points
        elif len(self.point_plate_ids) != self.num_points:
            raise ValueError(
                "Length of 'point_plate_ids' must match length of 'points'."
            )

        self.detect_collisions = detect_collisions

    def reconstruct(self):
        # Initialise the reconstruction.
        self.begin_reconstruction()

        # Loop over the reconstruction times until reached end of the reconstruction time span, or
        # all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        while self.reconstruct_to_next_time():
            pass

        return self.get_active_current_points()

    def begin_reconstruction(self):
        self.current_time_index = 0

        # Set up point arrays.
        # Store active and inactive points here (inactive points have None in corresponding entries).
        self.prev_points = [None] * self.num_points
        self.curr_points = [None] * self.num_points
        self.next_points = [None] * self.num_points

        # Each point can only get activated once (after deactivation it cannot be reactivated).
        self.point_has_been_activated = [False] * self.num_points
        self.num_activated_points = 0

        # Set up topology arrays (corresponding to active/inactive points at same indices).
        self.prev_topology_plate_ids = [None] * self.num_points
        self.curr_topology_plate_ids = [None] * self.num_points
        self.prev_resolved_plate_boundaries = [None] * self.num_points
        self.curr_resolved_plate_boundaries = [None] * self.num_points

        # Array to store indices of points found in continents
        self.in_continent_indices = [None] * self.num_points
        self.in_continent_points = [None] * self.num_points

        self.deletedpoints = []

        self._activate_deactivate_points()
        self._find_resolved_topologies_containing_points()

    def get_current_time(self):
        return (
            self.reconstruction_begin_time
            + self.current_time_index * self.reconstruction_time_step
        )

    def get_all_current_points(self):
        return self.curr_points

    def get_active_current_points(self):
        # Return only the active points (the ones that are not None).
        return [point for point in self.get_all_current_points() if point is not None]

    def get_in_continent_indices(self):
        return self.in_continent_points, self.in_continent_indices

    def reconstruct_to_next_time(self):
        # If we're at the last time then there is no next time to reconstruct to.
        if self.current_time_index == self.last_time_index:
            return False

        # If all points have been previously activated, but none are currently active then we're finished.
        # This means all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        if self.num_activated_points == self.num_points and not any(self.curr_points):
            return False

        # Cache stage rotations by plate ID.
        # Use different dicts since using different rotation models and time steps, etc.
        reconstruct_stage_rotation_dict = {}

        current_time = self.get_current_time()

        # Iterate over all points to reconstruct them to the next time step.
        for point_index in range(self.num_points):
            curr_point = self.curr_points[point_index]
            if curr_point is None:
                # Current point is not currently active.
                # So we cannot reconstruct to next time.
                self.next_points[point_index] = None
                continue

            # Get plate ID of resolved topology containing current point
            # (this was determined in last call to '_find_resolved_topologies_containing_points()').
            curr_plate_id = self.curr_topology_plate_ids[point_index]
            if curr_plate_id is None:
                # Current point is currently active but it fell outside all resolved polygons.
                # So instead we just reconstruct using its plate ID (that was manually assigned by the user/caller).
                curr_plate_id = self.point_plate_ids[point_index]

            # Get the stage rotation that will move the point from where it is at the current time to its
            # location at the next time step, based on the plate id that contains the point at the current time.

            # Speed up by caching stage rotations in a dict.
            stage_rotation = reconstruct_stage_rotation_dict.get(curr_plate_id)
            if not stage_rotation:
                stage_rotation = self.rotation_model.get_rotation(
                    # Positive/negative time step means reconstructing backward/forward in time.
                    current_time + self.reconstruction_time_step,
                    curr_plate_id,
                    current_time,
                )
                reconstruct_stage_rotation_dict[curr_plate_id] = stage_rotation

            # Use the stage rotation to reconstruct the tracked point from position at current time
            # to position at the next time step.
            self.next_points[point_index] = stage_rotation * curr_point

        #
        # Set up for next loop iteration.
        #
        # Rotate previous, current and next point arrays.
        # The new previous will be the old current.
        # The new current will be the old next.
        # The new next will be the old previous (but values are ignored and overridden in next time step; just re-using its memory).
        self.prev_points, self.curr_points, self.next_points = (
            self.curr_points,
            self.next_points,
            self.prev_points,
        )
        # Swap previous and current topology arrays.
        # The new previous will be the old current.
        # The new current will be the old previous (but values are ignored and overridden in next time step; just re-using its memory).
        self.prev_topology_plate_ids, self.curr_topology_plate_ids = (
            self.curr_topology_plate_ids,
            self.prev_topology_plate_ids,
        )
        self.prev_resolved_plate_boundaries, self.curr_resolved_plate_boundaries = (
            self.curr_resolved_plate_boundaries,
            self.prev_resolved_plate_boundaries,
        )

        # Move the current time to the next time.
        self.current_time_index += 1
        current_time = self.get_current_time()

        self._activate_deactivate_points()
        self._find_resolved_topologies_containing_points()

        # Iterate over all points to detect collisions.
        if self.detect_collisions:
            for point_index in range(self.num_points):
                prev_point = self.prev_points[point_index]
                curr_point = self.curr_points[point_index]
                if prev_point is None or curr_point is None:
                    # If current point is not currently active then no need to detect a collision for it (to deactivate it).
                    # Also previous point might just have been activated now, at end of current time step, and hence
                    # not active at beginning of time step.
                    continue

                # Get plate IDs of resolved topology containing previous and current point
                # (this was determined in last call to '_find_resolved_topologies_containing_points()').
                #
                # Note that could be None, so the collision detection needs to handle that.
                prev_plate_id = self.prev_topology_plate_ids[point_index]
                curr_plate_id = self.curr_topology_plate_ids[point_index]

                # Detect collisions at the end of the current time step since we need previous, and current, points and topologies.
                # De-activate point (in 'curr_points') if subducted (forward in time) or consumed back into MOR (backward in time).
                if self.detect_collisions(
                    self.rotation_model,
                    current_time,
                    self.reconstruction_time_interval,
                    prev_point,
                    curr_point,
                    prev_plate_id,
                    self.prev_resolved_plate_boundaries[point_index],
                    curr_plate_id,
                    self.curr_resolved_plate_boundaries[point_index],
                ):
                    # An inactive point in 'curr_points' becomes None.
                    # It may have been reconstructed from the previous time step to a valid position
                    # but now we override that result as inactive.
                    self.curr_points[point_index] = None
                    # self.curr_points.remove(self.curr_points[point_index])
                    self.deletedpoints.append(point_index)

        # We successfully reconstructed to the next time.
        return True

    def _activate_deactivate_points(self):
        current_time = self.get_current_time()

        # Iterate over all points and activate/deactivate as necessary depending on each point's valid time range.
        for point_index in range(self.num_points):
            if self.curr_points[point_index] is None:
                if not self.point_has_been_activated[point_index]:
                    # Point is not active and has never been activated, so see if can activate it.
                    if (
                        current_time <= self.point_begin_times[point_index]
                        and current_time >= self.point_end_times[point_index]
                    ):
                        # The initial point is assumed to be the position at the current time
                        # which is typically the point's begin time (approximately).
                        # But it could be the beginning of the reconstruction time span (specified in constructor)
                        # if that falls in the middle of the point's valid time range - in this case the
                        # initial point position is assumed to be in a position that is some time *after*
                        # it appeared (at its begin time) - and this can happen, for example, if you have a
                        # uniform grids of points at some intermediate time and want to see how they
                        # reconstruct to either a younger or older time (remembering that points can
                        # be subducted forward in time and consumed back into a mid-ocean ridge going
                        # backward in time).
                        self.curr_points[point_index] = self.points[point_index]
                        self.point_has_been_activated[point_index] = True
                        self.num_activated_points += 1
            else:
                # Point is active, so see if can deactivate it.
                if not (
                    current_time <= self.point_begin_times[point_index]
                    and current_time >= self.point_end_times[point_index]
                ):
                    self.curr_points[point_index] = None

    def _find_resolved_topologies_containing_points(self):

        current_time = self.get_current_time()

        # Resolve the plate polygons for the current time.
        resolved_topologies = []
        pygplates.resolve_topologies(
            self.topology_features,
            self.rotation_model,
            resolved_topologies,
            current_time,
        )

        if ReconstructByTopologies.use_plate_partitioner:
            # Create a plate partitioner from the resolved polygons.
            plate_partitioner = pygplates.PlatePartitioner(
                resolved_topologies, self.rotation_model
            )
        else:
            # Some of 'curr_points' will be None so 'curr_valid_points' contains only the valid (not None)
            # points, and 'curr_valid_points_indices' is the same length as 'curr_points' but indexes into
            # 'curr_valid_points' so we can quickly find which point (and hence which resolved topology)
            # in 'curr_valid_points' is associated with the a particular point in 'curr_points'.
            curr_valid_points = []
            curr_valid_points_indices = [None] * self.num_points
            for point_index, curr_point in enumerate(self.curr_points):
                if curr_point is not None:
                    curr_valid_points_indices[point_index] = len(curr_valid_points)
                    curr_valid_points.append(curr_point)
            # For each valid current point find the resolved topology containing it.
            resolved_topologies_containing_curr_valid_points = (
                _ptt.utils.points_in_polygons.find_polygons(
                    curr_valid_points,
                    [
                        resolved_topology.get_resolved_boundary()
                        for resolved_topology in resolved_topologies
                    ],
                    resolved_topologies,
                )
            )

        # Iterate over all points.
        for point_index, curr_point in enumerate(self.curr_points):
            if curr_point is None:
                # Current point is not currently active - so skip it.
                self.curr_topology_plate_ids[point_index] = None
                self.curr_resolved_plate_boundaries[point_index] = None
                continue

            # Find the plate id of the polygon that contains 'curr_point'.
            if ReconstructByTopologies.use_plate_partitioner:
                curr_polygon = plate_partitioner.partition_point(curr_point)
            else:
                curr_polygon = resolved_topologies_containing_curr_valid_points[
                    # Index back into 'curr_valid_points' and hence also into
                    # 'resolved_topologies_containing_curr_valid_points'.
                    curr_valid_points_indices[point_index]
                ]
            self.curr_resolved_plate_boundaries[point_index] = curr_polygon

            # If the polygon is None, that means (presumably) that it fell into a crack between
            # topologies. So it will be skipped and thrown away from future iterations.
            if curr_polygon is None:
                self.curr_topology_plate_ids[point_index] = None
                continue

            # Set the plate ID of resolved topology containing current point.
            self.curr_topology_plate_ids[point_index] = (
                curr_polygon.get_feature().get_reconstruction_plate_id()
            )