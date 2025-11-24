#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
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
from abc import ABC, abstractmethod
import logging
import math
import os
from typing import Union

import numpy as np
import pygplates

from .. import ptt as _ptt
from ..grids import read_netcdf_grid

logger = logging.getLogger("gplately")


class _DefaultCollision(object):
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
            prev_location_velocity = pygplates.calculate_velocities(  # type: ignore
                (curr_point,),
                prev_location_velocity_stage_rotation,
                1,
                pygplates.VelocityUnits.kms_per_my,
            )[0]
            curr_location_velocity = pygplates.calculate_velocities(  # type: ignore
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


_DEFAULT_COLLISION = _DefaultCollision()


class _ContinentCollision(object):
    """
    Continental collision detection function class (the function is the '__call__' method).
    """

    def __init__(
        self,
        grd_output_dir,
        chain_collision_detection=_DEFAULT_COLLISION,
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


class _ReconstructByTopologies(ABC):
    """Reconstruct geometries using topologies. Currently only points are supported."""

    @abstractmethod
    def _begin_reconstruction_impl(self):
        """Called from 'begin_reconstruction()' and implemented in derived class.

        Derived class should call 'self._activate_deactivate_points_using_their_begin_end_times()' to introduce points
        whose begin time is equal to (or older than) the initial time.

        Derived class can optionally do anything else needed before starting reconstruction over the time intervals.
        """
        pass  # This is an abstract method, no implementation here.

    @abstractmethod
    def _reconstruct_to_next_time_impl(self):
        """Called from 'reconstruct_to_next_time()' and implemented in derived class.

        Derived class should reconstruct the current active points in 'self.curr_points' to the next points 'self.next_points' from
        the current time 'self.get_current_time()' to the next time 'self.get_current_time() + self.reconstruction_time_step'.
        And collision detected for points should deactivate those points (ie, set None in 'self.next_points').

        Derived class then should cycle the three point arrays in preparation for the next time interval.
        Essentially do this:
            self.prev_points, self.curr_points, self.next_points = self.curr_points, self.next_points, self.prev_points

        Derived class then should increment the current time index (to update the current time to the next time).
        Essentially do this:
            self.current_time_index += 1

        Finally, derived class needs to call 'self._activate_deactivate_points_using_their_begin_end_times()' to introduce points
        whose begin time is the new current time, and remove points whose end time is the new current time.
        """
        pass  # This is an abstract method, no implementation here.

    def __init__(
        self,
        reconstruction_begin_time,
        reconstruction_end_time,
        reconstruction_time_interval,
        points,
        point_begin_times: Union[np.ndarray, list, None] = None,
        point_end_times: Union[np.ndarray, list, None] = None,
    ):

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
        # It's possible the time step is larger than the time span, in which case we change it to equal the time span
        # unless the reconstruction begin and end times are equal (in which case there'll be only one reconstruction snapshot).
        # This guarantees there'll be at least one time step (which has two times; one at either end of interval).
        if (
            self.num_times == 1
            and self.reconstruction_end_time != self.reconstruction_begin_time
        ):
            self.num_times = 2
            self.reconstruction_time_step = (
                self.reconstruction_end_time - self.reconstruction_begin_time
            )
        self.reconstruction_time_interval = math.fabs(self.reconstruction_time_step)

        self.last_time_index = self.num_times - 1

        self.points = points
        self.num_points = len(points)

        # Use the specified point begin times if provided (otherwise use 'inf').
        if point_begin_times is None:
            self.point_begin_times = np.full(self.num_points, np.inf)
        else:
            # Make sure numpy array (if not already).
            self.point_begin_times = np.asarray(point_begin_times)
            if len(self.point_begin_times) != self.num_points:
                raise ValueError(
                    "Length of 'point_begin_times' must match length of 'points'."
                )

        # Use the specified point end times if provided (otherwise use '-inf').
        if point_end_times is None:
            self.point_end_times = np.full(self.num_points, -np.inf)
        else:
            # Make sure numpy array (if not already).
            self.point_end_times = np.asarray(point_end_times)
            if len(self.point_end_times) != self.num_points:
                raise ValueError(
                    "Length of 'point_end_times' must match length of 'points'."
                )

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
        self.point_has_been_activated = np.zeros(self.num_points, dtype=bool)
        self.num_activated_points = 0

        # Call derived class implementation.
        self._begin_reconstruction_impl()

    def get_current_time_index(self):
        return self.current_time_index

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

    def reconstruct_to_next_time(self):
        # If we're at the last time then there is no next time to reconstruct to.
        if self.current_time_index == self.last_time_index:
            return False

        # If all points have been previously activated, but none are currently active then we're finished.
        # This means all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        if self.num_activated_points == self.num_points and not any(self.curr_points):
            return False

        # Call derived class implementation.
        self._reconstruct_to_next_time_impl()

        # We successfully reconstructed to the next time.
        return True

    def _activate_deactivate_points_using_their_begin_end_times(self):
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


class _ReconstructByTopologiesImpl(_ReconstructByTopologies):
    """Reconstruct geometries using topologies. Currently only points are supported."""

    use_plate_partitioner = False
    """If the use_plate_partitioner is True then use pygplates.PlatePartitioner to partition points,
        otherwise use faster points_in_polygons.find_polygons()."""

    def __init__(
        self,
        rotation_features_or_model,
        topology_features,
        reconstruction_begin_time,
        reconstruction_end_time,
        reconstruction_time_interval,
        points,
        *,
        point_begin_times: Union[np.ndarray, list, None] = None,
        point_end_times: Union[np.ndarray, list, None] = None,
        point_plate_ids: Union[np.ndarray, list, None] = None,
        detect_collisions=_DEFAULT_COLLISION,
    ):
        """
        rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).

        topology_features: Topology feature collection(s), or list of features, or filename(s) or any combination of those.

        detect_collisions: Collision detection function, or None. Defaults to _DEFAULT_COLLISION.
        """

        super().__init__(
            reconstruction_begin_time,
            reconstruction_end_time,
            reconstruction_time_interval,
            points,
            point_begin_times,
            point_end_times,
        )

        self.rotation_model = pygplates.RotationModel(rotation_features_or_model)

        # Turn topology data into a list of features (if not already).
        self.topology_features = pygplates.FeaturesFunctionArgument(
            topology_features
        ).get_features()

        # Use the specified point plate IDs if provided (otherwise use '0').
        # These plate IDs are only used when a point falls outside all resolved topologies during a time step.
        if point_plate_ids is None:
            self.point_plate_ids = np.zeros(self.num_points, dtype=int)
        else:
            # Make sure numpy array (if not already).
            self.point_plate_ids = np.asarray(point_plate_ids)
            if len(self.point_plate_ids) != self.num_points:
                raise ValueError(
                    "Length of 'point_plate_ids' must match length of 'points'."
                )

        self.detect_collisions = detect_collisions

    def _begin_reconstruction_impl(self):
        # Set up topology arrays (corresponding to active/inactive points at same indices as the point arrays).
        self.prev_topology_plate_ids = [None] * self.num_points
        self.curr_topology_plate_ids = [None] * self.num_points
        self.prev_resolved_plate_boundaries = [None] * self.num_points
        self.curr_resolved_plate_boundaries = [None] * self.num_points

        # Activate any original points whose begin time is equal to (or older than) the newly updated current time.
        # Deactivate any activated points whose end time is equal to (or younger than) the newly updated current time.
        #
        # Note: This should be done before calling 'self._find_resolved_topologies_containing_points()' so that new
        #       points just appearing at their begin time will have associated topologies (that contain them).
        self._activate_deactivate_points_using_their_begin_end_times()

        self._find_resolved_topologies_containing_points()

    def _reconstruct_to_next_time_impl(self):
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
        #
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
        #
        # Move the current time to the next time.
        self.current_time_index += 1
        current_time = self.get_current_time()

        # Activate any original points whose begin time is equal to (or older than) the newly updated current time.
        # Deactivate any activated points whose end time is equal to (or younger than) the newly updated current time.
        #
        # Note: This should be done before calling 'self._find_resolved_topologies_containing_points()' so that new
        #       points just appearing at their begin time will have associated topologies (that contain them).
        self._activate_deactivate_points_using_their_begin_end_times()

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

    def _find_resolved_topologies_containing_points(self):
        current_time = self.get_current_time()

        # Resolve the plate polygons for the current time.
        resolved_topologies = []
        pygplates.resolve_topologies(  # type: ignore
            self.topology_features,
            self.rotation_model,
            resolved_topologies,
            current_time,
        )

        if _ReconstructByTopologiesImpl.use_plate_partitioner:
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
                    curr_valid_points_indices[point_index] = len(curr_valid_points)  # type: ignore
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
            if _ReconstructByTopologiesImpl.use_plate_partitioner:
                curr_polygon = plate_partitioner.partition_point(curr_point)  # type: ignore
            else:
                curr_polygon = resolved_topologies_containing_curr_valid_points[  # type: ignore
                    # Index back into 'curr_valid_points' and hence also into
                    # 'resolved_topologies_containing_curr_valid_points'.
                    curr_valid_points_indices[point_index]  # type: ignore
                ]  # type: ignore
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


class _ReconstructByTopologiesImplV2(_ReconstructByTopologies):
    """An improved version of _ReconstructByTopologiesImpl.

    This uses a different approach to collision detection of seed points
    (to determine if they get subducted going forward in time and hence deactivated).
    Previously we detected if the velocity delta of a seed point exceeded a threshold (when that point
    transitioned from one plate to another) - and the point had to be close enough to the boundary.
    Instead, we now take the plate containing a seed point at the current time and resolve it at the next time
    (ie, current time plus time step) and see if its resolved shape still contains that seed point
    (reconstructed to the next time). If it does then the seed point remains active (otherwise it's deactivated).
    Note that the *same* plate is resolved at the current and next times. Usually when you resolve topologies
    at a different time you will get different plates (eg, there could be plate merges or splits).
    So, to resolve the *same* plate, some fenagling is required to ensure that a resolved plate boundary at the
    current time can also be resolved at the next time. The end result is that any boundary around a plate that
    converges will naturally deactivate seed points near that boundary (without having to use thresholds, etc).
    And this doesn't require plate boundaries to be labelled as converging (ie, subduction zones).
    """

    def __init__(
        self,
        rotation_features_or_model,
        topology_features,
        reconstruction_begin_time,
        reconstruction_end_time,
        reconstruction_time_interval,
        continent_mask_filepath_format,
        points,
        *,
        point_begin_times: Union[np.ndarray, list, None] = None,
        point_end_times: Union[np.ndarray, list, None] = None,
    ):
        """
        continent_mask_filepath_format: str
            The format of the file path of the continental mask grids that is converted to a
            file path using ``continent_mask_filepath_format.format(time)``.
        """

        super().__init__(
            reconstruction_begin_time,
            reconstruction_end_time,
            reconstruction_time_interval,
            points,
            point_begin_times,
            point_end_times,
        )

        self.rotation_model = pygplates.RotationModel(rotation_features_or_model)  # type: ignore

        # Turn topology data into a list of features (if not already).
        self.topology_features = pygplates.FeaturesFunctionArgument(  # type: ignore
            topology_features
        ).get_features()

        self.continent_mask_filepath_format = continent_mask_filepath_format

    def _begin_reconstruction_impl(self):
        # Activate any original points whose begin time is equal to (or older than) the newly updated current time.
        # Deactivate any activated points whose end time is equal to (or younger than) the newly updated current time.
        self._activate_deactivate_points_using_their_begin_end_times()

    def _reconstruct_to_next_time_impl(self):

        current_time = self.get_current_time()
        # Positive/negative time step means reconstructing backward/forward in time.
        next_time = current_time + self.reconstruction_time_step

        curr_topological_snapshot = pygplates.TopologicalSnapshot(  # type: ignore
            self.topology_features, self.rotation_model, current_time
        )
        curr_resolved_topologies = curr_topological_snapshot.get_resolved_topologies()

        # Some of 'curr_points' will be None, so 'curr_active_points' contains only the active (not None)
        # points, and 'curr_active_points_indices' is the same length as 'curr_active_points' but indexes
        # into 'curr_points' so we can quickly find which point in 'curr_points' is associated with a
        # particular point in 'curr_active_points'.
        curr_active_points = []
        curr_active_points_indices = []
        for point_index, curr_point in enumerate(self.curr_points):
            if curr_point is None:
                # Current point is not currently active, so we cannot reconstruct it to the next point.
                self.next_points[point_index] = None
                continue
            curr_active_points.append(curr_point)
            curr_active_points_indices.append(point_index)

        # For each active current point find the resolved topology containing it.
        curr_active_point_resolved_topologies = (
            _ptt.utils.points_in_polygons.find_polygons(
                curr_active_points,
                [
                    resolved_topology.get_resolved_boundary()
                    for resolved_topology in curr_resolved_topologies
                ],
                curr_resolved_topologies,
            )
        )

        # Map of current resolved topologies to the active points contained within them (their active point indices).
        map_curr_resolved_topology_to_active_point_indices = {}

        # Iterate over the resolved topologies containing all currently active points and get a list of the unique resolved topologies.
        #
        # Note: A resolved topology for a currently active point can be None if it fell outside all resolved topologies.
        #       In this case we'll just deactivate the point. Previously we would keep it active but just not reconstruct it
        #       (ie, the current and next positions would be the same). However the point might end up in weird locations.
        #       So it's best to just remove it. And this shouldn't happen very often at all (for topologies with global coverage).
        curr_active_point_index = 0
        while curr_active_point_index < len(curr_active_points):

            curr_active_point_resolved_topology = curr_active_point_resolved_topologies[
                curr_active_point_index
            ]
            # See if current active point fell outside all current resolved topologies.
            if curr_active_point_resolved_topology is None:
                # Index into the active and inactive points 'self.curr_points'.
                curr_point_index = curr_active_points_indices[curr_active_point_index]

                # Current point is currently active but it fell outside all resolved topologies.
                # So we deactivate it.
                self.next_points[curr_point_index] = None

                # Also remove evidence that the current point is active.
                del curr_active_points[curr_active_point_index]
                del curr_active_points_indices[curr_active_point_index]
                del curr_active_point_resolved_topologies[curr_active_point_index]

                # Continue to next active point.
                #
                # Note: We don't increment the active point index.
                #       We just removed an active point which essentially does the same thing.
                continue

            # If resolved topology has not been encountered yet then give it an empty list of active point indices.
            if (
                curr_active_point_resolved_topology
                not in map_curr_resolved_topology_to_active_point_indices
            ):
                map_curr_resolved_topology_to_active_point_indices[
                    curr_active_point_resolved_topology
                ] = []

            # Add the index of the currently active point to the resolved topology containing it.
            map_curr_resolved_topology_to_active_point_indices[
                curr_active_point_resolved_topology
            ].append(curr_active_point_index)

            # Increment to the next active point.
            curr_active_point_index += 1

        # The current resolved topologies that contain active points.
        curr_active_resolved_topologies = list(
            map_curr_resolved_topology_to_active_point_indices.keys()
        )

        # The *next* topology features and their topological section features.
        next_topological_features = self._NextTopologicalFeatures(next_time)
        # For each current resolved topology, create a corresponding *next* topological boundary feature and
        # all the topological section features referenced by them.
        for curr_resolved_topology in curr_active_resolved_topologies:
            next_topological_features.add_current_topology(curr_resolved_topology)

        # Resolve the topologies to the *next* time step.
        next_topological_snapshot = pygplates.TopologicalSnapshot(  # type: ignore
            next_topological_features.next_topological_section_features
            + next_topological_features.next_topological_boundary_features,
            self.rotation_model,
            next_time,
        )

        # First, map the *next* *boundary* features to their indices into 'next_topological_boundary_features'.
        #
        # This will also be the indices into 'curr_active_resolved_topologies' (due to order of insertion).
        map_next_topology_boundary_feature_to_index = {
            next_topology_feature: index
            for index, next_topology_feature in enumerate(
                next_topological_features.next_topological_boundary_features
            )
        }
        #
        # Then, map each current resolved topology to its next resolved topology.
        map_curr_to_next_resolved_topologies = {}
        for (
            next_resolved_topology
        ) in next_topological_snapshot.get_resolved_topologies():
            index = map_next_topology_boundary_feature_to_index.get(
                next_resolved_topology.get_feature()
            )
            # It's possible that the next topology could not get resolved at the next time step.
            # This generally shouldn't happen though.
            if index is not None:
                curr_resolved_topology = curr_active_resolved_topologies[index]
                map_curr_to_next_resolved_topologies[curr_resolved_topology] = (
                    next_resolved_topology
                )

        # Iterate over all currently active points to reconstruct them to the next time step and
        # see if they've been consumed by a convergent plate boundary.
        #
        # Note: Active points are grouped by the current resolved topology containing them.
        for (
            curr_resolved_topology,
            curr_resolved_topology_active_point_indices,
        ) in map_curr_resolved_topology_to_active_point_indices.items():

            # The next resolved topology associated with the current resolved topology.
            next_resolved_topology = map_curr_to_next_resolved_topologies.get(
                curr_resolved_topology
            )
            if next_resolved_topology is None:
                # The current topology could not be resolved to the next time step.
                # So we deactivate all points that it contains. This generally shouldn't happen though.
                for (
                    curr_active_point_index
                ) in curr_resolved_topology_active_point_indices:
                    curr_point_index = curr_active_points_indices[
                        curr_active_point_index
                    ]
                    self.next_points[curr_point_index] = None
                continue

            # Get the boundary polygon of the next resolved topology.
            #
            # If the next resolved topology is consistent with the current resolved topology then we're fine,
            # otherwise we need to replace the boundary of the next resolved topology with something more consistent.
            next_resolved_topology_boundary = self._NextTopologicalBoundary(
                curr_resolved_topology,
                next_resolved_topology,
                self.rotation_model,
            ).get_next_resolved_topology_boundary()

            # Iterate over the currently active points contained by the current resolved topology.
            for curr_active_point_index in curr_resolved_topology_active_point_indices:

                # Reconstruct the currently active point from its position at current time to its position at the next time step.
                curr_active_point = curr_active_points[curr_active_point_index]
                next_active_point = curr_resolved_topology.reconstruct_point(
                    curr_active_point, next_time
                )

                # Index into the active and inactive points 'self.curr_points'.
                curr_point_index = curr_active_points_indices[curr_active_point_index]

                # See if the location (of the current active point) reconstructed to the *next* time step
                # is inside the current topology resolved to the *next* time step.
                #
                # If it is outside then it is subducting going forward in time (or consumed by a mid-ocean ridge going backward in time).
                # It doesn't necessarily have to happen at a subduction zone (or mid-ocean ridge). It can be any part of a plate boundary
                # that is convergent forward in time (or divergent forward in time; which is convergent backward in time).
                if not next_resolved_topology_boundary.is_point_in_polygon(
                    next_active_point
                ):
                    self.next_points[curr_point_index] = None
                    continue

                # The current active point was not consumed by a plate boundary.
                self.next_points[curr_point_index] = next_active_point

        # See if any 'next' points collide with the continents and deactivate those that do.
        self._detect_continent_collisions(next_time)

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
        #
        # Move the current time to the next time.
        self.current_time_index += 1
        current_time = self.get_current_time()

        # Activate any original points whose begin time is equal to (or older than) the newly updated current time.
        # Deactivate any activated points whose end time is equal to (or younger than) the newly updated current time.
        self._activate_deactivate_points_using_their_begin_end_times()

    def _detect_continent_collisions(self, next_time):
        # Get those 'next' points that are active.
        next_active_points = []
        next_active_points_indices = []
        for point_index, next_point in enumerate(self.next_points):
            if next_point:
                next_active_points.append(next_point)
                next_active_points_indices.append(point_index)

        # Convert points to lat/lon.
        points_lat = np.empty(len(next_active_points))
        points_lon = np.empty(len(next_active_points))
        for active_point_index, point in enumerate(next_active_points):
            point_lat, point_lon = point.to_lat_lon()
            points_lat[active_point_index] = point_lat
            points_lon[active_point_index] = point_lon

        # Read the continent mask at 'next_time'.
        continent_mask_filepath = self.continent_mask_filepath_format.format(next_time)
        gridZ, gridX, gridY = read_netcdf_grid(
            continent_mask_filepath, return_grids=True
        )
        ni, nj = gridZ.shape
        xmin = np.nanmin(gridX)
        xmax = np.nanmax(gridX)
        ymin = np.nanmin(gridY)
        ymax = np.nanmax(gridY)

        # Sample continent mask grid, which is one over continents and zero over oceans.
        points_i = (ni - 1) * ((points_lat - ymin) / (ymax - ymin))
        points_j = (nj - 1) * ((points_lon - xmin) / (xmax - xmin))
        points_i_uint = np.rint(points_i).astype(np.uint)
        points_j_uint = np.rint(points_j).astype(np.uint)
        try:
            mask_values = gridZ[points_i_uint, points_j_uint]
        except IndexError:
            points_i = np.clip(np.rint(points_i), 0, ni - 1).astype(np.int_)
            points_j = np.clip(np.rint(points_j), 0, nj - 1).astype(np.int_)
            mask_values = gridZ[points_i, points_j]

        # Deactivate any points that sampled inside the continent mask.
        for active_point_index in np.where(mask_values >= 0.5)[0]:
            self.next_points[next_active_points_indices[active_point_index]] = None

    class _NextTopologicalFeatures(object):
        """The *next* topology boundary features and their topological section features.

        These are essentially the same as the *current* topology features (and their topological section features) but
        with their valid time periods modified (if needed) to include the *next* time.
        """

        def __init__(self, next_time):

            self.next_time = next_time
            # All topology *boundary* features and their referenced topological *section* features ready to be resolved at the next time.
            #
            # Note: Topological *section* features can include topological lines.
            self.next_topological_boundary_features = []
            self.next_topological_section_features = []

            # Map of each current topological section feature to its corresponding *next* topological section feature.
            self._map_curr_to_next_topological_section_feature = {}

        def add_current_topology(self, curr_resolved_topology):

            # Iterate over the boundary sub-segments of the current resolved topology and create a new topological
            # section feature for each one (if its topological section feature hasn't been encountered yet).
            next_boundary_section_property_values = []
            for (
                curr_boundary_sub_segment
            ) in curr_resolved_topology.get_boundary_sub_segments():  # type: ignore

                curr_boundary_section_feature = curr_boundary_sub_segment.get_feature()
                # If the current boundary section feature hasn't been encountered yet then create an associated
                # *next* boundary section feature that is either the same as the current boundary section feature
                # (if it's still valid at the *next* time) or a clone of the current boundary section feature with
                # the valid time range modified to be valid at the *next* time.
                if (
                    curr_boundary_section_feature
                    not in self._map_curr_to_next_topological_section_feature
                ):
                    curr_boundary_section_recon_geom = (
                        curr_boundary_sub_segment.get_topological_section()
                    )
                    # Boundary section is either a ReconstructedFeatureGeometry or a ResolvedTopologicalLine.
                    #
                    # If it's a ResolvedTopologicalLine then we need to iterate through its sub-segments and create a
                    # new line section for each one, and finally create a new topological line feature from them.
                    if isinstance(
                        curr_boundary_section_recon_geom, pygplates.ResolvedTopologicalLine  # type: ignore
                    ):

                        # Iterate over the sub-segments of the current resolved topological line and create a new line
                        # section feature for each one (if its topological section feature hasn't been encountered yet).
                        next_line_section_property_values = []
                        for (
                            curr_line_sub_segment
                        ) in curr_boundary_section_recon_geom.get_line_sub_segments():  # type: ignore
                            curr_line_section_feature = (
                                curr_line_sub_segment.get_feature()
                            )
                            # If the current line section feature hasn't been encountered yet then create an associated
                            # *next* line section feature that is either the same as the current line section feature
                            # (if it's still valid at the *next* time) or a clone of the current line section feature with
                            # the valid time range modified to be valid at the *next* time.
                            if (
                                curr_line_section_feature
                                not in self._map_curr_to_next_topological_section_feature
                            ):
                                # Line section must be a ReconstructedFeatureGeometry (ie, can't be ResolvedTopologicalLine)
                                # since a topological line cannot have sections that are also topological lines.

                                # If feature is not valid at 'next_time' then clone the feature and set a valid time that includes 'next_time'.
                                if not curr_line_section_feature.is_valid_at_time(
                                    self.next_time
                                ):
                                    next_line_section_feature = (
                                        curr_line_section_feature.clone()
                                    )
                                    next_line_section_feature.set_valid_time(
                                        pygplates.GeoTimeInstant.create_distant_past(), 0.0  # type: ignore
                                    )
                                else:
                                    next_line_section_feature = (
                                        curr_line_section_feature
                                    )

                                self.next_topological_section_features.append(
                                    next_line_section_feature
                                )

                                # Map the current to the next (line section feature).
                                self._map_curr_to_next_topological_section_feature[
                                    curr_line_section_feature
                                ] = next_line_section_feature

                            # Get the next line section feature associated with the current one.
                            next_line_section_feature = (
                                self._map_curr_to_next_topological_section_feature[
                                    curr_line_section_feature
                                ]
                            )

                            # The property name of the geometry referenced by the next line section.
                            next_line_section_geometry_property_name = (
                                curr_line_sub_segment.get_topological_section()
                                .get_property()
                                .get_name()
                            )
                            # Whether to reverse the geometry referenced by the next line section.
                            next_line_section_reverse_order = (
                                curr_line_sub_segment.was_geometry_reversed_in_topology()
                            )

                            # Create the *next* line section associated with the current one.
                            next_line_section_property_value = pygplates.GpmlTopologicalSection.create(  # type: ignore
                                next_line_section_feature,
                                geometry_property_name=next_line_section_geometry_property_name,
                                reverse_order=next_line_section_reverse_order,
                                topological_geometry_type=pygplates.GpmlTopologicalLine,  # type: ignore
                            )
                            if not next_line_section_property_value:
                                raise RuntimeError(
                                    "Unable to create topological section for a topological line."
                                )
                            next_line_section_property_values.append(
                                next_line_section_property_value
                            )

                        # Create a topological line feature that enables the current topological line
                        # to be resolved to the *next* time step.
                        #
                        # Note: 'valid_time' argument is not specified which means valid for all time.
                        #       This works for us since we only need it to be valid at the *next* time step.
                        next_line_feature = pygplates.Feature(  # type: ignore
                            pygplates.FeatureType.gpml_unclassified_feature  # type: ignore
                        )
                        next_line_feature.set_topological_geometry(
                            pygplates.GpmlTopologicalLine(  # type: ignore
                                next_line_section_property_values
                            ),
                            # Use the same geometry property name as the current boundary sub-segment so
                            # that it can be found by the topological polygon(s) that will reference it...
                            curr_boundary_sub_segment.get_topological_section()
                            .get_property()
                            .get_name(),
                        )

                        next_boundary_section_feature = next_line_feature

                    else:  # boundary section is a ReconstructedFeatureGeometry ...

                        # If feature is not valid at 'next_time' then clone the feature and set a valid time that includes 'next_time'.
                        if not curr_boundary_section_feature.is_valid_at_time(
                            self.next_time
                        ):
                            next_boundary_section_feature = (
                                curr_boundary_section_feature.clone()
                            )
                            next_boundary_section_feature.set_valid_time(
                                pygplates.GeoTimeInstant.create_distant_past(), 0.0  # type: ignore
                            )
                        else:
                            next_boundary_section_feature = (
                                curr_boundary_section_feature
                            )

                    self.next_topological_section_features.append(
                        next_boundary_section_feature
                    )

                    # Map the current to the next (boundary section feature).
                    self._map_curr_to_next_topological_section_feature[
                        curr_boundary_section_feature
                    ] = next_boundary_section_feature

                # Get the next boundary section feature associated with the current one.
                next_boundary_section_feature = (
                    self._map_curr_to_next_topological_section_feature[
                        curr_boundary_section_feature
                    ]
                )

                # The property name of the geometry referenced by the next boundary section.
                next_boundary_section_geometry_property_name = (
                    curr_boundary_sub_segment.get_topological_section()
                    .get_property()
                    .get_name()
                )
                # Whether to reverse the geometry referenced by the next boundary section.
                next_boundary_section_reverse_order = (
                    curr_boundary_sub_segment.was_geometry_reversed_in_topology()
                )

                # Create the *next* boundary section associated with the current one.
                next_boundary_section_property_value = pygplates.GpmlTopologicalSection.create(  # type: ignore
                    next_boundary_section_feature,
                    geometry_property_name=next_boundary_section_geometry_property_name,
                    reverse_order=next_boundary_section_reverse_order,
                    topological_geometry_type=pygplates.GpmlTopologicalPolygon,  # type: ignore
                )
                if not next_boundary_section_property_value:
                    raise RuntimeError(
                        "Unable to create topological section for a topological polygon."
                    )
                next_boundary_section_property_values.append(
                    next_boundary_section_property_value
                )

            # Create a topological boundary feature that enables the current topological boundary
            # to be resolved to the *next* time step.
            #
            # Note: If the current topology is a deforming network we still create a topological closed plate boundary
            #       (which is normally used for rigid plates) because we only need to detect if a seed point
            #       reconstructed to the next time step is contained within the network's *boundary*.
            #
            # Note: 'valid_time' argument is not specified which means valid for all time.
            #       This works for us since we only need it to be valid at the *next* time step.
            next_boundary_feature = pygplates.Feature.create_topological_feature(  # type: ignore
                pygplates.FeatureType.gpml_topological_closed_plate_boundary,  # type: ignore
                pygplates.GpmlTopologicalPolygon(  # type: ignore
                    next_boundary_section_property_values
                ),
            )
            self.next_topological_boundary_features.append(next_boundary_feature)

    class _NextTopologicalBoundary(object):
        """The next resolved topology boundary.

        If the next resolved topology is NOT consistent with the current resolved topology then we need to replace
        replace any inconsistent boundary sub-segments of the next resolved topology with the equivalent ones in the
        current resolved topology (but moved to the next time). For each inconsistent boundary sub-segment, this involves
        taking the boundary sub-segment of the *current* resolved topology and moving it to the *next* time, and using that
        (instead of the boundary sub-segment of the *next* resolved topology). This is not as ideal as actually resolving the
        all the *unclipped* boundary sections of the current topology at the *next* time, which would produce a more accurate
        plate boundary for the next resolved topology. But the inconsistency forces us to take an alternative approach.
        However this case shouldn't happen very often.

        The next resolved topology is not consistent with the current resolved topology if its boundary sections don't have
        the same number of intersections with neighbouring sections. For example, a boundary section of the current resolved topology
        intersects its neighbour once (as expected) but the same boundary section of the next resolved topology (which is really
        just the same topology resolved to the next time) intersects the same neighbour twice. In this case, the next resolved topology
        would likely have an unexpectedly different plate boundary shape than the current resolved topology, which might deactivate
        a bunch of points that it shouldn't. This happens because pyGPlates only considers the first intersection (if there are two or more).
        Usually the topological model is built such that each neighbouring section only intersects once. However, when you resolve
        that same topology at a different time (even if it's just 1 Myr different) then the intersections may not be what the
        model builder intended (especially if the current topology's time period does not include the *next* time).
        """

        # Distance from a polyline-polyline intersection such that if we clipped that much distance off the polylines
        # around that intersection then the clipped polylines are unlikely to intersect (unless they're almost parallel).
        #
        # We are a little conversative on this (ie, a little big) to ensure we don't get too close.
        INTERSECTION_THRESHOLD_RADIANS = 1e-4  # approx 600 metres

        # Maximum number of intersections to test between two intersecting polylines.
        #
        # If two polylines reach the max number of intersections then we've likely detected a partial overlap
        # between the two polylines (where a segment from each polyline are on top of each other and
        # hence essentially have an infinite number of intersections).
        MAX_NUM_INTERSECTIONS = 5

        def __init__(
            self, curr_resolved_topology, next_resolved_topology, rotation_model
        ):

            self.curr_resolved_topology = curr_resolved_topology
            self.next_resolved_topology = next_resolved_topology

            self.curr_boundary_sub_segments = (
                curr_resolved_topology.get_boundary_sub_segments()
            )
            self.next_boundary_sub_segments = (
                next_resolved_topology.get_boundary_sub_segments()
            )

            if len(self.curr_boundary_sub_segments) != len(
                self.next_boundary_sub_segments
            ):
                raise RuntimeError(
                    "Current and next topologies have a different number of boundary sub-segments."
                )
            self.num_boundary_sub_segments = len(self.curr_boundary_sub_segments)

            self.curr_time = curr_resolved_topology.get_reconstruction_time()
            self.next_time = next_resolved_topology.get_reconstruction_time()
            self.rotation_model = rotation_model

        def get_next_resolved_topology_boundary(self):

            # If the current and next resolved topologies are consistent then
            # just return the resolved boundary of the next resolved topology.
            inconsistent_boundary_sub_segment_indices = (
                self._get_inconsistent_boundary_sub_segments()
            )
            if not inconsistent_boundary_sub_segment_indices:
                return self.next_resolved_topology.get_resolved_boundary()

            # Instead of returning the boundary of the next resolved topology, replace any inconsistent boundary sub-segments
            # of the next resolved topology with the equivalent ones in the current resolved topology (but moved to the next time).
            return self._get_next_boundary_from_consistent_boundary_sub_segments(
                inconsistent_boundary_sub_segment_indices
            )

        def _get_next_boundary_from_consistent_boundary_sub_segments(
            self, inconsistent_boundary_sub_segment_indices
        ):

            #
            # Instead of the resolved boundary of the *next* resolved topology we'll replace the inconsistent boundary sub-segments
            # of the *next* resolved topology with the equivalent ones in the *current* resolved topology (but moved to the next time).
            #

            # Iterate over the next boundary sub-segments and join them to form a polygon boundary.
            next_boundary_polygon_points = []
            for boundary_sub_segment_index in range(self.num_boundary_sub_segments):

                # If the boundary sub-segment is consistent (ie, not inconsistent) then just add the
                # (potentially reversed) points of the boundary sub-segment of the *next* resolved topology.
                if (
                    boundary_sub_segment_index
                    not in inconsistent_boundary_sub_segment_indices
                ):
                    next_boundary_sub_segment = self.next_boundary_sub_segments[
                        boundary_sub_segment_index
                    ]
                    next_boundary_polygon_points.extend(
                        reversed(
                            next_boundary_sub_segment.get_resolved_geometry_points()
                        )
                        if next_boundary_sub_segment.was_geometry_reversed_in_topology()
                        else next_boundary_sub_segment.get_resolved_geometry_points()
                    )

                    continue

                #
                # The boundary sub-segment is inconsistent.
                #
                # So take the boundary sub-segment of the *current* resolved topology and move it to the *next* time,
                # and use that (instead of the boundary sub-segment of the *next* resolved topology).
                #

                curr_boundary_sub_segment = self.curr_boundary_sub_segments[
                    boundary_sub_segment_index
                ]
                curr_boundary_sub_segment_reversal = (
                    curr_boundary_sub_segment.was_geometry_reversed_in_topology()
                )

                # See if the boundary sub-segment is from a topological *line*.
                curr_boundary_sub_sub_segments = (
                    curr_boundary_sub_segment.get_sub_segments()
                )
                if not curr_boundary_sub_sub_segments:
                    #
                    # Boundary sub-segment is NOT from a topological *line*.
                    #
                    curr_boundary_sub_segment_feature = (
                        curr_boundary_sub_segment.get_resolved_feature()
                    )

                    # The sub-segment resolved feature is already reconstructed to 'curr_time'.
                    # So we need to reverse reconstruct it to present day (before we can reconstruct it to 'next_time').
                    #
                    # Note: We don't need to clone this feature before modifying it because it was already essentially cloned
                    #       when 'ResolvedTopologicalSubSegment.get_resolved_feature()' was called.
                    pygplates.reverse_reconstruct(  # type:ignore
                        curr_boundary_sub_segment_feature,
                        self.rotation_model,
                        self.curr_time,
                    )

                    # The resolved sub-segment feature is valid at the 'curr_time' but not necessarily at 'next_time'.
                    # If not valid at 'next_time' then set a valid time that includes 'next_time'.
                    if not curr_boundary_sub_segment_feature.is_valid_at_time(
                        self.next_time
                    ):
                        curr_boundary_sub_segment_feature.set_valid_time(
                            pygplates.GeoTimeInstant.create_distant_past(), 0.0  # type: ignore
                        )

                    # Reconstruct to 'next_time' (from present day).
                    next_boundary_sub_segment_reconstructed_geometries = (
                        pygplates.ReconstructSnapshot(  # type:ignore
                            curr_boundary_sub_segment_feature,
                            self.rotation_model,
                            self.next_time,
                        ).get_reconstructed_geometries()
                    )

                    # Should only be one reconstructed geometry.
                    if len(next_boundary_sub_segment_reconstructed_geometries) != 1:
                        raise RuntimeError("Expected a single reconstructed geometry.")
                    next_boundary_sub_segment_reconstructed_geometry = (
                        next_boundary_sub_segment_reconstructed_geometries[0]
                    )

                    # Add the (potentially reversed) points of the sub-segment to the next polygon boundary.
                    next_boundary_polygon_points.extend(
                        reversed(
                            next_boundary_sub_segment_reconstructed_geometry.get_reconstructed_geometry_points()
                        )
                        if curr_boundary_sub_segment_reversal
                        else next_boundary_sub_segment_reconstructed_geometry.get_reconstructed_geometry_points()
                    )

                    continue

                #
                # The boundary sub-segment IS from a topological *line*, so we need to iterate over its sub-segments.
                # This is because we need features that are *reconstructable*.
                # Topological lines are not reconstructable (because they're topologies linking reconstructable features).
                # Only the sub-segments of topological lines are reconstructable.
                #
                # Each sub-segment of the resolved topological *line* will contribute to the boundary sub-segment.
                #
                # Note: Only a sub-section of the full resolved topological *line* will contribute
                #       to the boundary of the current resolved topological boundary/network.
                #       Hence not all the sub-segments of the resolved topological *line* will contribute
                #       to the boundary of the current resolved topological boundary/network.
                #

                # Collect the sub-sub-segment features to use for the next boundary.
                #
                # Note: We need to reverse the order of sub-sub-segments if the resolved topological *line* is reversed
                #       in the current resolved boundary.
                curr_boundary_sub_sub_segment_features = [
                    sub_sub_segment.get_resolved_feature()
                    for sub_sub_segment in (
                        reversed(curr_boundary_sub_sub_segments)
                        if curr_boundary_sub_segment_reversal
                        else curr_boundary_sub_sub_segments
                    )
                ]
                # If a sub-sub-segment is reversed in the resolved topological *line* which is, in turn, reversed in the
                # current resolved boundary then the sub-sub-segment is NOT reversed in the current resolved boundary.
                curr_boundary_sub_sub_segment_reversals = [
                    curr_boundary_sub_segment_reversal
                    ^ sub_sub_segment.was_geometry_reversed_in_topology()
                    for sub_sub_segment in (
                        reversed(curr_boundary_sub_sub_segments)
                        if curr_boundary_sub_segment_reversal
                        else curr_boundary_sub_sub_segments
                    )
                ]

                # The sub-segment resolved features are already reconstructed to 'curr_time'.
                # So we need to reverse reconstruct them to present day (before we can reconstruct them to 'next_time').
                #
                # Note: We don't need to clone these features before modifying them because they were already essentially cloned
                #       when 'ResolvedTopologicalSubSegment.get_resolved_feature()' was called.
                pygplates.reverse_reconstruct(  # type:ignore
                    curr_boundary_sub_sub_segment_features,
                    self.rotation_model,
                    self.curr_time,
                )

                # The resolved sub-segment features are valid at the 'curr_time' but not necessarily at 'next_time'.
                # If not valid at 'next_time' then set a valid time that includes 'next_time'.
                for feature in curr_boundary_sub_sub_segment_features:
                    if not feature.is_valid_at_time(self.next_time):
                        feature.set_valid_time(
                            pygplates.GeoTimeInstant.create_distant_past(), 0.0  # type: ignore
                        )

                # Reconstruct to 'next_time' (from present day).
                next_boundary_sub_sub_segment_reconstructed_geometries = pygplates.ReconstructSnapshot(  # type:ignore
                    curr_boundary_sub_sub_segment_features,
                    self.rotation_model,
                    self.next_time,
                ).get_reconstructed_geometries(
                    # Make sure the order is retained...
                    same_order_as_reconstructable_features=True
                )

                # The order and length of reconstructed geometries should match that of the features.
                if len(next_boundary_sub_sub_segment_reconstructed_geometries) != len(
                    curr_boundary_sub_sub_segment_features
                ):
                    raise RuntimeError(
                        "Expected number of reconstructed geometries to match number of features."
                    )

                # Add the (potentially reversed) points of the sub-segments to the next polygon boundary.
                for (
                    feature_index,
                    sub_sub_segment_reconstructed_geometry,
                ) in enumerate(next_boundary_sub_sub_segment_reconstructed_geometries):
                    next_boundary_polygon_points.extend(
                        reversed(
                            sub_sub_segment_reconstructed_geometry.get_reconstructed_geometry_points()
                        )
                        if curr_boundary_sub_sub_segment_reversals[feature_index]
                        else sub_sub_segment_reconstructed_geometry.get_reconstructed_geometry_points()
                    )

            # Join the (potentially reversed) points of the boundary sub-segments and join them together into a polygon boundary.
            # Each sub-segment is either from the *current* or *next* resolved topology.
            return pygplates.PolygonOnSphere(  # type:ignore
                next_boundary_polygon_points
            )

        def _get_inconsistent_boundary_sub_segments(self):
            """Returns a set of indices of any boundary sub-segments that are inconsistent, otherwise None."""

            # If there's only 1 sub-segment then it cannot intersect with itself (but it can still be converted to a polygon).
            #
            # Note: If there's 2 sub-segments then they can intersect each other twice (to form a closed polygon).
            #       We will still process them to make sure they stay consistent (same number of intersections at current and next times).
            if self.num_boundary_sub_segments <= 1:
                return None

            inconsistent_boundary_sub_segment_indices = None

            # Compare number of neighbour intersections for the current and next boundary sub-segments.
            for boundary_sub_segment_index in range(self.num_boundary_sub_segments):
                curr_boundary_sub_segment = self.curr_boundary_sub_segments[
                    boundary_sub_segment_index
                ]
                next_boundary_sub_segment = self.next_boundary_sub_segments[
                    boundary_sub_segment_index
                ]

                prev_boundary_sub_segment_index = boundary_sub_segment_index - 1
                if prev_boundary_sub_segment_index < 0:
                    prev_boundary_sub_segment_index += self.num_boundary_sub_segments

                # Previous neighbour sub-segment (for the current and next resolved topologies).
                prev_curr_boundary_sub_segment = self.curr_boundary_sub_segments[
                    prev_boundary_sub_segment_index
                ]
                prev_next_boundary_sub_segment = self.next_boundary_sub_segments[
                    prev_boundary_sub_segment_index
                ]

                # Number of neighbour intersections of sub-segment (for the current and next resolved topologies).
                curr_num_intersections = self._get_num_intersections(
                    prev_curr_boundary_sub_segment.get_topological_section_geometry(),
                    curr_boundary_sub_segment.get_topological_section_geometry(),
                )
                next_num_intersections = self._get_num_intersections(
                    prev_next_boundary_sub_segment.get_topological_section_geometry(),
                    next_boundary_sub_segment.get_topological_section_geometry(),
                )

                # If the number of intersections differ then the current and next resolved topologies are inconsistent.
                if curr_num_intersections != next_num_intersections:
                    # We only create a set when we actually need one
                    # (most often we won't since most resolved topologies are consistent).
                    if inconsistent_boundary_sub_segment_indices is None:
                        inconsistent_boundary_sub_segment_indices = set()
                    # Add the previous and current sub-segments affected by the intersection(s).
                    inconsistent_boundary_sub_segment_indices.add(
                        prev_boundary_sub_segment_index
                    )
                    inconsistent_boundary_sub_segment_indices.add(
                        boundary_sub_segment_index
                    )

            return inconsistent_boundary_sub_segment_indices

        def _get_num_intersections(self, polyline1, polyline2, num_intersections=0):

            distance_tuple = pygplates.GeometryOnSphere.distance(  # type: ignore
                polyline1,
                polyline2,
                # We're only detecting zero distance (an intersection), so as an optimisation we
                # don't need to calculate an accurate non-zero distance (instead just return None)...
                distance_threshold_radians=self.INTERSECTION_THRESHOLD_RADIANS,
                return_closest_positions=True,
                return_closest_indices=True,
            )

            if distance_tuple is not None:
                (
                    distance,
                    polyline1_intersection,
                    polyline2_intersection,
                    polyline1_segment_index,
                    polyline2_segment_index,
                ) = distance_tuple

                if distance == 0.0:
                    # We found an intersection (between 'polyline1' and 'polyline2').
                    num_intersections += 1

                    # If we reached the max number of intersections then we've likely detected a partial overlap
                    # between the two polylines (where a segment from each polyline are on top of each other and
                    # hence essentially have an infinite number of intersections).
                    #
                    # In this case we'll just return the maximum number of intersections. This means if both the
                    # polylines overlap at both the current and next times then they'll return the same number
                    # of intersections (max) and compare equal and hence be considered consistent.
                    if num_intersections == self.MAX_NUM_INTERSECTIONS:
                        return num_intersections

                    # Split 'polyline1' into two polylines at its intersection point.
                    polyline1a = self._pre_split_polyline(
                        polyline1, polyline1_segment_index, polyline1_intersection
                    )
                    polyline1b = self._post_split_polyline(
                        polyline1, polyline1_segment_index, polyline1_intersection
                    )

                    # Split 'polyline2' into two polylines at its intersection point.
                    polyline2a = self._pre_split_polyline(
                        polyline2, polyline2_segment_index, polyline2_intersection
                    )
                    polyline2b = self._post_split_polyline(
                        polyline2, polyline2_segment_index, polyline2_intersection
                    )

                    # Test for intersections among the four combinations of split polylines.
                    #
                    # Note: It's possible there may be less than four split polylines if
                    #       either (or both) polylines were intersected at one of their end points.
                    if polyline1a:
                        if polyline2a:
                            num_intersections = self._get_num_intersections(
                                polyline1a, polyline2a, num_intersections
                            )
                            if num_intersections == self.MAX_NUM_INTERSECTIONS:
                                return num_intersections
                        if polyline2b:
                            num_intersections = self._get_num_intersections(
                                polyline1a, polyline2b, num_intersections
                            )
                            if num_intersections == self.MAX_NUM_INTERSECTIONS:
                                return num_intersections
                    if polyline1b:
                        if polyline2a:
                            num_intersections = self._get_num_intersections(
                                polyline1b, polyline2a, num_intersections
                            )
                            if num_intersections == self.MAX_NUM_INTERSECTIONS:
                                return num_intersections
                        if polyline2b:
                            num_intersections = self._get_num_intersections(
                                polyline1b, polyline2b, num_intersections
                            )
                            if num_intersections == self.MAX_NUM_INTERSECTIONS:
                                return num_intersections

            return num_intersections

        def _pre_split_polyline(
            self, polyline, polyline_segment_index, polyline_intersection
        ):
            # Get the points *before* the intersection (including the intersection).
            # This will have at least two points (needed for creating a polyline).
            pre_split_points = list(polyline[: polyline_segment_index + 1])
            pre_split_points.append(polyline_intersection)

            # Remove any zero length segments next to the intersected end of the polyline.
            # Because we've already detected the intersection and don't want it detected in subsequent intersection tests.
            while True:
                pre_split_last_segment = pygplates.GreatCircleArc(  # type: ignore
                    pre_split_points[-2],
                    pre_split_points[-1],
                )
                if not pre_split_last_segment.is_zero_length():
                    break

                # Remove the last segment (by removing last point).
                del pre_split_points[-1]
                if len(pre_split_points) < 2:
                    # All segments before the intersection are zero length.
                    # So there's no pre-split polyline.
                    return None

            # We've encountered a non-zero length segment. Let's shorten it such that it will no longer intersect the intersection point.
            # Because we've already detected the intersection and don't want it detected in subsequent intersection tests.
            if (
                self.INTERSECTION_THRESHOLD_RADIANS
                < pre_split_last_segment.get_arc_length()
            ):
                # Last segment is *longer* than the intersection threshold, so shorten it by that amount.
                # Replace last point with rotated point.
                pre_split_points[-1] = (
                    pygplates.FiniteRotation(  # type: ignore
                        pre_split_last_segment.get_rotation_axis(),
                        -self.INTERSECTION_THRESHOLD_RADIANS,  # negative rotates away from segment *end* point
                    )
                    * pre_split_points[-1]
                )
            else:
                # Last segment is *shorter* than the intersection threshold. But we know it's not a zero length segment,
                # so it should be far enough away from the intersection point. So just remove the segment (by removing last point).
                del pre_split_points[-1]
                if len(pre_split_points) < 2:
                    # There are no segments before the intersection.
                    # So there's no pre-split polyline.
                    return None

            return pygplates.PolylineOnSphere(pre_split_points)  # type: ignore

        def _post_split_polyline(
            self, polyline, polyline_segment_index, polyline_intersection
        ):
            # Get the points *after* the intersection (including the intersection).
            # This will have at least two points (needed for creating a polyline).
            post_split_points = [polyline_intersection]
            post_split_points.extend(polyline[polyline_segment_index + 1 :])

            # Remove any zero length segments next to the intersected start of the polyline.
            # Because we've already detected the intersection and don't want it detected in subsequent intersection tests.
            while True:
                post_split_first_segment = pygplates.GreatCircleArc(  # type: ignore
                    post_split_points[0],
                    post_split_points[1],
                )
                if not post_split_first_segment.is_zero_length():
                    break

                # Remove the first segment (by removing first point).
                del post_split_points[0]
                if len(post_split_points) < 2:
                    # All segments after the intersection are zero length.
                    # So there's no post-split polyline.
                    return None

            # We've encountered a non-zero length segment. Let's shorten it such that it will no longer intersect the intersection point.
            # Because we've already detected the intersection and don't want it detected in subsequent intersection tests.
            if (
                self.INTERSECTION_THRESHOLD_RADIANS
                < post_split_first_segment.get_arc_length()
            ):
                # First segment is *longer* than the intersection threshold, so shorten it by that amount.
                # Replace first point with rotated point.
                post_split_points[0] = (
                    pygplates.FiniteRotation(  # type: ignore
                        post_split_first_segment.get_rotation_axis(),
                        self.INTERSECTION_THRESHOLD_RADIANS,
                    )
                    * post_split_points[0]
                )
            else:
                # First segment is *shorter* than the intersection threshold. But we know it's not a zero length segment,
                # so it should be far enough away from the intersection point. So just remove the segment (by removing first point).
                del post_split_points[0]
                if len(post_split_points) < 2:
                    # There are no segments after the intersection.
                    # So there's no post-split polyline.
                    return None

            return pygplates.PolylineOnSphere(post_split_points)  # type: ignore


class _ReconstructByTopologicalModelImpl(_ReconstructByTopologies):
    """Similar to ``_ReconstructByTopologiesImpl`` except uses the ``pygplates.TopologicalModel`` class to reconstruct seed points.

    This is currently just a transition towards using ``pygplates.TopologicalModel``.
    But note that this still uses Python code, like in class ``_DefaultCollision``, to do the collision detection
    (as opposed to collision detection built into pyGPlates, which needs to be updated to better support GPlately).
    Also it's not as fast as ``_ReconstructByTopologiesImpl`` which has faster point-in-polygon testing (since it uses a spatial tree).

    So, for the time being it's probably still better to use class ``_ReconstructByTopologiesImpl`` instead.
    """

    class CollisionDelegator(pygplates.ReconstructedGeometryTimeSpan.DeactivatePoints):
        """Delegates collision detection from pyGPlates (pygplates.TopologicalModel) to the collision classes accessed by ``_ReconstructByTopologies``
        (like class ``_DefaultCollision``).

        Ultimately there should be less (or no) need to delegate once the default pyGPlates collision handler (``DefaultDeactivatePoints``) supports the
        funcionality of classes like ``_DefaultCollision``. Although ``_ContinentCollision`` might still need to be handled outside pyGPlates (ie, delegated to).
        """

        def __init__(
            self, detect_collisions, rotation_model, reconstruction_time_interval
        ):
            super().__init__()
            self.detect_collisions = detect_collisions
            self.rotation_model = rotation_model
            self.reconstruction_time_interval = reconstruction_time_interval

        def deactivate(
            self,
            prev_point,
            prev_location,
            prev_time,
            curr_point,
            curr_location,
            curr_time,
        ):
            # Get plate IDs of resolved topology containing previous and current point.
            #
            # Note that could be None, so the collision detection needs to handle that.
            prev_plate_id = None
            curr_plate_id = None

            # See if previous point is located in a resolved rigid plate or deforming network (or neither).
            prev_resolved_topology = prev_location.located_in_resolved_boundary()
            if not prev_resolved_topology:
                prev_resolved_topology = prev_location.located_in_resolved_network()
            if prev_resolved_topology:
                prev_plate_id = (
                    prev_resolved_topology.get_feature().get_reconstruction_plate_id()
                )

            # See if current point is located in a resolved rigid plate or deforming network (or neither).
            curr_resolved_topology = curr_location.located_in_resolved_boundary()
            if not curr_resolved_topology:
                curr_resolved_topology = curr_location.located_in_resolved_network()
            if curr_resolved_topology:
                curr_plate_id = (
                    curr_resolved_topology.get_feature().get_reconstruction_plate_id()
                )

            # Delegate to the Python implementation of collision detection (eg, classes '_DefaultCollision' and '_ContinentCollision').
            return self.detect_collisions(
                self.rotation_model,
                curr_time,
                self.reconstruction_time_interval,
                prev_point,
                curr_point,
                prev_plate_id,
                prev_resolved_topology,
                curr_plate_id,
                curr_resolved_topology,
            )

    def __init__(
        self,
        rotation_features_or_model,
        topology_features,
        reconstruction_begin_time,
        reconstruction_end_time,
        reconstruction_time_interval,
        points,
        *,
        point_begin_times: Union[np.ndarray, list, None] = None,
        point_end_times: Union[np.ndarray, list, None] = None,
        detect_collisions=_DEFAULT_COLLISION,
    ):
        """
        rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).

        topology_features: Topology feature collection(s), or list of features, or filename(s) or any combination of those.

        detect_collisions: Collision detection function, or None. Defaults to _DEFAULT_COLLISION.
        """

        super().__init__(
            reconstruction_begin_time,
            reconstruction_end_time,
            reconstruction_time_interval,
            points,
            point_begin_times,
            point_end_times,
        )

        self.topological_model = pygplates.TopologicalModel(
            topology_features,
            rotation_features_or_model,
            # Only really need to cache 2 topological snapshots since we progressively move through time (and hence never return to previous times).
            # Might only need to cache 1 topological snapshot (not sure). Best to be safe with 2...
            topological_snapshot_cache_size=2,
        )

        self.detect_collisions = self.CollisionDelegator(
            detect_collisions,
            self.topological_model.get_rotation_model(),
            self.reconstruction_time_interval,
        )

    def _begin_reconstruction_impl(self):
        # Activate any original points whose begin time is equal to (or older than) the newly updated current time.
        # Deactivate any activated points whose end time is equal to (or younger than) the newly updated current time.
        self._activate_deactivate_points_using_their_begin_end_times()

    def _reconstruct_to_next_time_impl(self):

        current_time = self.get_current_time()
        next_time = current_time + self.reconstruction_time_step

        # Find the active points.
        # These will be reconstructed to the next time step.
        current_active_points = []
        current_active_point_indices = []
        for point_index in range(self.num_points):
            current_point = self.curr_points[point_index]
            if not current_point:
                # Current point is not currently active.
                # So we cannot reconstruct to next time.
                self.next_points[point_index] = None
                continue

            # Keep track of the currently active points.
            # And their original indices into ALL points (active and inactive).
            current_active_points.append(current_point)
            current_active_point_indices.append(point_index)

        # Reconstruct active points to the next time step.
        reconstructed_time_span = self.topological_model.reconstruct_geometry(
            current_active_points,
            initial_time=current_time,
            oldest_time=max(current_time, next_time),
            youngest_time=min(current_time, next_time),
            time_increment=self.reconstruction_time_interval,  # must be positive
            deactivate_points=self.detect_collisions,
        )

        # Store the next points back to their original locations in ALL points (active and inactive).
        next_points = reconstructed_time_span.get_geometry_points(
            next_time, return_inactive_points=True
        )
        for next_point_index, next_point in enumerate(next_points):
            # Get index into ALL points (active and inactive).
            point_index = current_active_point_indices[next_point_index]
            # Update next point (note that this can be None if a collision was detected).
            self.next_points[point_index] = next_point

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
        #
        # Move the current time to the next time.
        self.current_time_index += 1

        # Activate any original points whose begin time is equal to (or older than) the newly updated current time.
        # Deactivate any activated points whose end time is equal to (or younger than) the newly updated current time.
        self._activate_deactivate_points_using_their_begin_end_times()


def reconstruct_points_by_topologies(
    rotation_features_or_model,
    topology_features,
    reconstruction_begin_time,
    reconstruction_end_time,
    reconstruction_time_interval,
    points,
    point_begin_times=None,
    point_end_times=None,
    point_plate_ids=None,
    detect_collisions=_DEFAULT_COLLISION,
):
    """Reconstruct points using the topological polygons."""

    topology_reconstruction = _ReconstructByTopologiesImpl(
        rotation_features_or_model,
        topology_features,
        reconstruction_begin_time,
        reconstruction_end_time,
        reconstruction_time_interval,
        points,
        point_begin_times=point_begin_times,
        point_end_times=point_end_times,
        point_plate_ids=point_plate_ids,
        detect_collisions=detect_collisions,
    )

    return topology_reconstruction.reconstruct()
