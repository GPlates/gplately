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
import logging
import math
import os
from typing import Union

import numpy as np
import pygplates

from .. import ptt as _ptt
from ..grids import read_netcdf_grid

logger = logging.getLogger("gplately")


class ReconstructByTopologies(object):
    """Reconstruct points using topologies and deactivates those that subduct at converging plate boundaries.

    This uses a new approach to collision detection of seed points
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
        plate_reconstruction,
        reconstruction_begin_time,
        reconstruction_end_time,
        reconstruction_time_interval,
        points,
        *,
        point_begin_times: Union[np.ndarray, list, None] = None,
        point_end_times: Union[np.ndarray, list, None] = None,
        continent_mask_filepath_format=None,
    ):
        """
        continent_mask_filepath_format: str
            The format of the file path of the continental mask grids that is converted to a
            file path using ``continent_mask_filepath_format.format(time)``.
        """

        self.plate_reconstruction = plate_reconstruction

        # Set up an array of reconstruction times covering the reconstruction time span.
        self.reconstruction_begin_time = reconstruction_begin_time
        self.reconstruction_end_time = reconstruction_end_time
        if reconstruction_time_interval <= 0.0:
            raise ValueError("'reconstruction_time_interval' must be positive.")
        # Reconstruction can go forward or backward in time.
        if self.reconstruction_begin_time > self.reconstruction_end_time:
            self._reconstruction_time_step = -reconstruction_time_interval
        else:
            self._reconstruction_time_step = reconstruction_time_interval
        # Get number of times including end time if time span is a multiple of time step.
        # The '1' is because, for example, 2 time intervals is 3 times.
        # The '1e-6' deals with limited floating-point precision, eg, we want (3.0 - 0.0) / 1.0 to be 3.0 and not 2.999999 (which gets truncated to 2).
        self._num_times = 1 + int(
            math.floor(
                1e-6
                + float(self.reconstruction_end_time - self.reconstruction_begin_time)
                / self._reconstruction_time_step
            )
        )
        # It's possible the time step is larger than the time span, in which case we change it to equal the time span
        # unless the reconstruction begin and end times are equal (in which case there'll be only one reconstruction snapshot).
        # This guarantees there'll be at least one time step (which has two times; one at either end of interval).
        if (
            self._num_times == 1
            and self.reconstruction_end_time != self.reconstruction_begin_time
        ):
            self._num_times = 2
            self._reconstruction_time_step = (
                self.reconstruction_end_time - self.reconstruction_begin_time
            )
        self._reconstruction_time_interval = math.fabs(self._reconstruction_time_step)

        self._last_time_index = self._num_times - 1

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

        self._continent_mask_filepath_format = continent_mask_filepath_format

    def reconstruct(self):
        # Initialise the reconstruction.
        self.begin_reconstruction()

        # Loop over the reconstruction times until reached end of the reconstruction time span, or
        # all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        while self.reconstruct_to_next_time():
            pass

        return self.get_current_active_points()

    def begin_reconstruction(self):
        self._current_time_index = 0

        # Store active and inactive points here.
        #
        # NOTE: Deactivated points will NOT get reset to None - we'll just ignore them.
        self._all_current_points = [None] * self.num_points

        # A boolean array indicating which points (in 'self._all_current_points') are active.
        self._point_is_currently_active = np.zeros(self.num_points, dtype=bool)

        # Each point can only get activated once (after deactivation it cannot be reactivated).
        self._point_has_been_activated = np.zeros(self.num_points, dtype=bool)

        # Activate any original points whose begin time is equal to (or older than) the newly updated current time.
        # Deactivate any activated points whose end time is equal to (or younger than) the newly updated current time.
        self._activate_deactivate_points_using_their_begin_end_times()

        # Deactivate any newly actived points that are masked by the continents.
        self._deactivate_continent_masked_points()

    def get_current_time_index(self):
        return self._current_time_index

    def get_current_time(self):
        return (
            self.reconstruction_begin_time
            + self._current_time_index * self._reconstruction_time_step
        )

    def get_current_active_points(self):
        """Return those points that are currently active."""
        return [
            self._all_current_points[point_index]
            for point_index in self.get_current_active_point_indices()
        ]

    def get_current_active_point_indices(self):
        """Return the indices of the currently active points (into the original points)."""
        return np.where(self._point_is_currently_active)[0]

    def reconstruct_to_next_time(self):
        # If we're at the last time then there is no next time to reconstruct to.
        if self._current_time_index == self._last_time_index:
            return False

        # If all points have been previously activated, but none are currently active then we're finished.
        # This means all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        if np.all(self._point_has_been_activated) and not np.any(
            self._point_is_currently_active
        ):
            return False

        # Call the main implementation.
        self._reconstruct_to_next_time_impl()

        # Move the current time to the next time.
        self._current_time_index += 1

        # Activate any original points whose begin time is equal to (or older than) the newly updated current time.
        # Deactivate any activated points whose end time is equal to (or younger than) the newly updated current time.
        self._activate_deactivate_points_using_their_begin_end_times()

        # Deactivate any active points that are masked by the continents.
        #
        # This includes deactivating any newly activated points.
        #
        # Note: This is done at the newly updated current time.
        self._deactivate_continent_masked_points()

        # We successfully reconstructed to the next time.
        return True

    def _reconstruct_to_next_time_impl(self):

        current_time = self.get_current_time()
        # Positive/negative time step means reconstructing backward/forward in time.
        next_time = current_time + self._reconstruction_time_step

        # Resolved the topologies at the current time.
        curr_topological_snapshot = self.plate_reconstruction.topological_snapshot(
            current_time, include_topological_slab_boundaries=False
        )

        # Get the resolved networks and rigid plates.
        curr_resolved_networks = curr_topological_snapshot.get_resolved_topologies(
            pygplates.ResolveTopologyType.network  # type:ignore
        )
        curr_resolved_boundaries = curr_topological_snapshot.get_resolved_topologies(
            pygplates.ResolveTopologyType.boundary  # type:ignore
        )

        # Get the currently active points and their indices (into the original points).
        curr_active_point_indices = self.get_current_active_point_indices()
        curr_active_points = [
            self._all_current_points[point_index]
            for point_index in curr_active_point_indices
        ]

        # For each active current point find the resolved topology network or rigid plate containing it.
        curr_active_points_spatial_tree = (
            _ptt.utils.points_spatial_tree.PointsSpatialTree(curr_active_points)
        )
        curr_active_point_resolved_networks = (
            _ptt.utils.points_in_polygons.find_polygons_using_points_spatial_tree(
                curr_active_points,
                curr_active_points_spatial_tree,
                [
                    resolved_topology.get_resolved_boundary()
                    for resolved_topology in curr_resolved_networks
                ],
                curr_resolved_networks,
            )
        )
        curr_active_point_resolved_boundaries = (
            _ptt.utils.points_in_polygons.find_polygons_using_points_spatial_tree(
                curr_active_points,
                curr_active_points_spatial_tree,
                [
                    resolved_topology.get_resolved_boundary()
                    for resolved_topology in curr_resolved_boundaries
                ],
                curr_resolved_boundaries,
            )
        )

        if len(curr_active_point_resolved_networks) != len(
            curr_active_point_resolved_boundaries
        ):
            raise RuntimeError(
                "Unexpected length mismatch of arrays of resolved topologies containing points."
            )

        # Map of current resolved topologies to the active points contained within them (their point indices).
        map_curr_resolved_topology_to_active_point_indices = {}

        # Iterate over the resolved topologies containing all currently active points and get a list of the unique resolved topologies.
        #
        # Note: A resolved topology for a currently active point can be None if it fell outside all resolved topologies.
        #       In this case we'll just deactivate the point. Previously we would keep it active but just not reconstruct it
        #       (ie, the current and next positions would be the same). However the point might end up in weird locations.
        #       So it's best to just remove it. And this shouldn't happen very often at all (for topologies with global coverage).
        for index in range(len(curr_active_points)):

            # Index into the active and inactive points 'self._all_current_points'.
            curr_point_index = curr_active_point_indices[index]

            # If point is inside a resolved *network* then prefer that, otherwise see if it's inside a rigid plate.
            # This is because deforming networks can overlay rigid plates.
            curr_active_point_resolved_topology = curr_active_point_resolved_networks[
                index
            ]
            if curr_active_point_resolved_topology is None:
                curr_active_point_resolved_topology = (
                    curr_active_point_resolved_boundaries[index]
                )

            # See if current active point fell outside all current resolved topologies.
            if curr_active_point_resolved_topology is None:
                # Current point is currently active but it fell outside all resolved topologies.
                # So we deactivate it.
                self._point_is_currently_active[curr_point_index] = False

                # Continue to next active point.
                continue

            # If resolved topology has not been encountered yet then give it an empty list of point indices.
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
            ].append(curr_point_index)

        # Prevent usage since might no longer be representative of the active points because
        # we might've just deactivated some points that fell outside all resolved topologies.
        del curr_active_points
        del curr_active_point_indices

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
            self.plate_reconstruction.rotation_model,
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
                self._point_is_currently_active[
                    curr_resolved_topology_active_point_indices
                ] = False
                continue

            # Get the boundary polygon of the next resolved topology.
            #
            # If the next resolved topology is consistent with the current resolved topology then we're fine,
            # otherwise we need to replace the boundary of the next resolved topology with something more consistent.
            next_resolved_topology_boundary = self._NextTopologicalBoundary(
                curr_resolved_topology,
                next_resolved_topology,
                self.plate_reconstruction.rotation_model,
            ).get_next_resolved_topology_boundary()

            # Iterate over the currently active points contained by the current resolved topology.
            for point_index in curr_resolved_topology_active_point_indices:

                # Reconstruct the currently active point from its position at current time to its position at the next time step.
                curr_active_point = self._all_current_points[point_index]
                next_active_point = curr_resolved_topology.reconstruct_point(
                    curr_active_point, next_time
                )

                # See if the location (of the current active point) reconstructed to the *next* time step
                # is inside the current topology resolved to the *next* time step.
                #
                # If it is outside then it is subducting going forward in time (or consumed by a mid-ocean ridge going backward in time).
                # It doesn't necessarily have to happen at a subduction zone (or mid-ocean ridge). It can be any part of a plate boundary
                # that is convergent forward in time (or divergent forward in time; which is convergent backward in time).
                if not next_resolved_topology_boundary.is_point_in_polygon(
                    next_active_point
                ):
                    self._point_is_currently_active[point_index] = False
                    continue

                # The current active point was not consumed by a plate boundary.
                # So update its position.
                self._all_current_points[point_index] = next_active_point

    def _activate_deactivate_points_using_their_begin_end_times(self):
        current_time = self.get_current_time()

        # Get indices of points that are not active and have never been activated.
        curr_inactive_and_not_yet_activated_point_indices = np.where(
            ~self._point_is_currently_active & ~self._point_has_been_activated
        )[0]
        # See which of those points we can activate.
        activate_point_indices = curr_inactive_and_not_yet_activated_point_indices[
            (
                current_time
                <= self.point_begin_times[
                    curr_inactive_and_not_yet_activated_point_indices
                ]
            )
            & (
                current_time
                >= self.point_end_times[
                    curr_inactive_and_not_yet_activated_point_indices
                ]
            )
        ]
        # Activate those points.
        self._point_is_currently_active[activate_point_indices] = True
        # Copy original point into currently active points array.
        for point_index in activate_point_indices:
            self._all_current_points[point_index] = self.points[point_index]
        # Mark those points as having been activated.
        self._point_has_been_activated[activate_point_indices] = True

        # Get indices of points that are active (and hence must have been previously activated).
        curr_active_point_indices = np.where(self._point_is_currently_active)[0]
        # See which of those points we can deactivate.
        deactivate_point_indices = curr_active_point_indices[
            ~(
                (current_time <= self.point_begin_times[curr_active_point_indices])
                & (current_time >= self.point_end_times[curr_active_point_indices])
            )
        ]
        # Deactivate those points.
        self._point_is_currently_active[deactivate_point_indices] = False
        # Leave deactivated points in the currently active points array (we'll just ignore them).

    def _deactivate_continent_masked_points(self):
        # If no continental collision detection requested then just return.
        if self._continent_mask_filepath_format is None:
            return

        # Get the currently active points and their indices (into the original points).
        curr_active_point_indices = self.get_current_active_point_indices()
        curr_active_points = [
            self._all_current_points[point_index]
            for point_index in curr_active_point_indices
        ]

        # Convert points to lat/lon.
        curr_points_lat = np.empty(len(curr_active_points))
        curr_points_lon = np.empty(len(curr_active_points))
        for index, point in enumerate(curr_active_points):
            point_lat, point_lon = point.to_lat_lon()
            curr_points_lat[index] = point_lat
            curr_points_lon[index] = point_lon

        # Read the continent mask at the current time.
        continent_mask_filepath = self._continent_mask_filepath_format.format(
            self.get_current_time()
        )
        gridZ, gridX, gridY = read_netcdf_grid(
            continent_mask_filepath, return_grids=True
        )
        ni, nj = gridZ.shape
        xmin = np.nanmin(gridX)
        xmax = np.nanmax(gridX)
        ymin = np.nanmin(gridY)
        ymax = np.nanmax(gridY)

        # Sample continent mask grid, which is one over continents and zero over oceans.
        points_i = (ni - 1) * ((curr_points_lat - ymin) / (ymax - ymin))
        points_j = (nj - 1) * ((curr_points_lon - xmin) / (xmax - xmin))
        points_i_uint = np.rint(points_i).astype(np.uint)
        points_j_uint = np.rint(points_j).astype(np.uint)
        try:
            mask_values = gridZ[points_i_uint, points_j_uint]
        except IndexError:
            points_i = np.clip(np.rint(points_i), 0, ni - 1).astype(np.int_)
            points_j = np.clip(np.rint(points_j), 0, nj - 1).astype(np.int_)
            mask_values = gridZ[points_i, points_j]

        # Deactivate any points that sampled inside the continent mask.
        self._point_is_currently_active[
            curr_active_point_indices[np.where(mask_values >= 0.5)[0]]
        ] = False

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


def reconstruct_points_by_topologies(
    plate_reconstruction,
    reconstruction_begin_time,
    reconstruction_end_time,
    reconstruction_time_interval,
    points,
    *,
    point_begin_times=None,
    point_end_times=None,
    continent_mask_filepath_format=None,
):
    """Reconstruct points using topologies."""

    topology_reconstruction = ReconstructByTopologies(
        plate_reconstruction,
        reconstruction_begin_time,
        reconstruction_end_time,
        reconstruction_time_interval,
        points,
        point_begin_times=point_begin_times,
        point_end_times=point_end_times,
        continent_mask_filepath_format=continent_mask_filepath_format,
    )

    return topology_reconstruction.reconstruct()
