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

import itertools
import math
import numpy as np
import pygplates
from functools import partial
from multiprocessing import Pool

from ..gpml import _load_FeatureCollection
from ..parallel import get_num_cpus
from ..ptt.utils import points_in_polygons, points_spatial_tree


class ReconstructContinents(object):

    def __init__(
        self,
        continent_features,
        max_time,
        min_time,
        time_step,
    ):
        """Create a ``ReconstructContinents`` object.

        Parameters
        ----------
        continent_features : str/`os.PathLike`, or a sequence (eg, `list` or `tuple`) of instances of `pygplates.Feature`_, or a single instance of `pygplates.Feature`_, or an instance of `pygplates.FeatureCollection`_, or a sequence of any combination of those four types
            These are any features on continental crust that are polygons (eg, continental polygons).
            Can be provided as a filename, or a sequence of features, or a single feature, or a feature collection, or a sequence (eg, a list or tuple) of any combination of those four types.
        max_time : float
            The maximum time for reconstructing continent features.
        min_time : float
            The minimum time for reconstructing continent features.
        time_step : float
            The delta time for reconstructing continent features.



        .. _pygplates.Feature: https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature
        .. _pygplates.FeatureCollection: https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection
        """
        self._continent_features_pickle = continent_features  # Pickle the __init__ argument (see __getstate__/__setstate__ for details).
        self.continent_features = _load_FeatureCollection(continent_features)

        self.max_time = max_time
        self.min_time = min_time
        self.time_step = time_step

        # The times from min to max time.
        self._times = np.arange(min_time, max_time + 1e-6, time_step)
        # A dict mapping times to their indices (into 'self._times').
        self._time_indices = {
            time: time_index for time_index, time in enumerate(self._times)
        }

        # Make sure the continent features are polygons and that they reconstruct by plate ID only.
        if not self._all_continent_features_are_polygons_that_reconstruct_by_plate_id():
            raise ValueError(
                "Each continent feature must contain a polygon and must reconstruct by plate ID only (excludes by half-stage, motion paths, flowlines, etc)"
            )

    def __getstate__(self):
        # Save the instance data variables.
        state = self.__dict__.copy()

        # Set 'continent_features' to None to avoid pickling it.
        #
        # Instead we're pickling '_continent_features_pickle'.
        # If they are just filenames then it's faster to rebuild FeatureCollection's (from files when unpickling)
        # than it is to pickle/unpickle FeatureCollection's.
        state["continent_features"] = None

        # Remove the unpicklable entries.
        #
        # This includes pygplates reconstructed feature geometries and resolved topological geometries.
        # Note: PyGPlates features and features collections (and rotation models) can be pickled though.
        #

        return state

    def __setstate__(self, state):
        # Restore the instance data variables.
        self.__dict__.update(state)

        # Load the continent features.
        #
        # The pickled attributes could be features and/or filesnames.
        # If they are just filenames then it's faster to reload FeatureCollection's
        # from files than it is to pickle/unpickle FeatureCollection's.
        self.continent_features = _load_FeatureCollection(
            self._continent_features_pickle
        )

        # Restore the unpicklable entries.
        #
        # This includes pygplates reconstructed feature geometries and resolved topological geometries.
        # Note: PyGPlates features and features collections (and rotation models) can be pickled though.
        #

    def reconstruct_continent_features(self, plate_reconstruction, nprocs=-2):
        """Reconstruct the continent features.

        The continent features are also deformed if ``plate_reconstruction`` contains a *deforming* topological model.

        Parameters
        ----------
        plate_reconstruction : PlateReconstruction
            A :class:`PlateReconstruction` object for reconstruction.
        nprocs : int, default=-2
            The number of CPUs to use for parts of the code that are parallelized.
            Must be an integer or convertible to an integer (eg, float is rounded towards zero).
            If positive then uses that many CPUs.
            If ``1`` then executes in serial (ie, is not parallelized).
            If ``0`` then a ``ValueError`` is raised.
            If ``-1`` then all available CPUs are used.
            If ``-2`` then all available CPUs except one are used, etc.
            Defaults to ``-2`` (ie, uses all available CPUs except one to keep system responsive).
        """
        self._reconstructed_continents = None

        if not self.continent_features:
            return

        # If the topological model has deforming networks then reconstruct and deform the continent features from min to max time.
        # Otherwise don't do anything - we'll just rigidly reconstruct them when the user asks for them.
        if self._has_deforming_network_features(plate_reconstruction):
            # Determine number of CPUs to use.
            num_cpus = get_num_cpus(nprocs)

            self._reconstructed_continents = (
                self._reconstruct_and_deform_continent_features(
                    plate_reconstruction,
                    num_cpus,
                )
            )

            all_reconstructed_continents = []
            for time in self._times:
                reconstructed_continents = self.get_reconstructed_continents(
                    time, plate_reconstruction
                )
                reconstructed_continent_features = [
                    reconstructed_continent.get_feature()
                    for reconstructed_continent in reconstructed_continents
                ]
                for feature in reconstructed_continent_features:
                    feature.set_valid_time(
                        time + 0.5 * self.time_step, time - 0.5 * self.time_step
                    )
                all_reconstructed_continents.extend(reconstructed_continent_features)
            pygplates.FeatureCollection(all_reconstructed_continents).write(
                "deformed_continents.gpmlz"
            )
        else:
            self._reconstructed_continents = None

    def get_reconstructed_continents(self, time, plate_reconstruction):
        """Retrieve the reconstructed continent features.

        Parameters
        ----------
        time : float
            The reconstruction time.
        plate_reconstruction : PlateReconstruction
            A :class:`PlateReconstruction` object for reconstruction.
        """
        if not self.continent_features:
            return []

        # Ensure 'time' is one of the valid reconstruction times.
        time_index = self._time_indices.get(time, None)
        if time_index is None:
            # The requested 'time' was not one of the valid times.
            raise ValueError(
                f"{time} is not a valid time in the range [{self.max_time}, {self.min_time}] with interval {self.time_step}"
            )

        if self._reconstructed_continents:
            continent_features = []

            # Get the already-reconstructed-and-deformed continental polygons.
            for feature_index, reconstructed_continent in enumerate(
                self._reconstructed_continents
            ):
                continent_feature = self.continent_features[feature_index].clone()

                for (
                    property_name,
                    reconstructed_polygon_time_spans,
                ) in reconstructed_continent.items():
                    reconstructed_geometries = [
                        pygplates.PolygonOnSphere(time_span[time_index])  # type:ignore
                        for time_span in reconstructed_polygon_time_spans
                    ]

                    continent_feature.set_geometry(
                        reconstructed_geometries,
                        property_name,
                        verify_information_model=pygplates.VerifyInformationModel.no,  # type:ignore
                    )

                continent_features.append(continent_feature)

            pygplates.reverse_reconstruct(  # type:ignore
                continent_features, plate_reconstruction.rotation_model, time
            )

        else:
            continent_features = self.continent_features

        # If the time is one of the valid times then rigidly reconstruct the present-day continent features to 'time'.
        return plate_reconstruction.reconstruct(continent_features, time)

    def _all_continent_features_are_polygons_that_reconstruct_by_plate_id(self):
        if not self.continent_features:
            return True

        # Make sure all the continent features reconstruct by plate ID only and contain a polygon.
        for feature in self.continent_features:
            # Exclude features that reconstruct by half-stage rotation.
            if feature.get_reconstruction_method().startswith("HalfStageRotation"):
                return False

            # Exclude features with feature type motion path, flowline, VGP.
            feature_type = feature.get_feature_type()
            if (
                feature_type == pygplates.FeatureType.gpml_motion_path  # type:ignore
                or feature_type == pygplates.FeatureType.gpml_flowline  # type:ignore
                or feature_type
                == pygplates.FeatureType.gpml_virtual_geomagnetic_pole  # type:ignore
            ):
                return False

            # Exclude features that don't have any polygon geometries.
            if not any(
                isinstance(geometry, pygplates.PolygonOnSphere)  # type:ignore
                for geometry in feature.get_all_geometries()
            ):
                return False

        return True

    def _has_deforming_network_features(self, plate_reconstruction):
        # See if topological model has any deforming network features.
        return any(
            feature.get_feature_type()
            == pygplates.FeatureType.gpml_topological_network  # type:ignore
            for feature in plate_reconstruction.topology_features
        )

    def _reconstruct_and_deform_continent_features(
        self, plate_reconstruction, num_cpus
    ):
        if not self.continent_features:
            return None

        if num_cpus > 1:

            continent_features = list(self.continent_features)
            continent_features_lists = []
            num_continent_features_per_task = math.ceil(
                len(continent_features) / (2 * num_cpus)
            )
            task_start_feature_index = 0
            while task_start_feature_index < len(continent_features):
                continent_features_lists.append(
                    continent_features[
                        task_start_feature_index : task_start_feature_index
                        + num_continent_features_per_task
                    ]
                )
                task_start_feature_index += num_continent_features_per_task

            with Pool(num_cpus) as pool:
                reconstructed_continents_list = pool.map(
                    partial(
                        _reconstruct_and_deform_continent_features_impl,
                        plate_reconstruction=plate_reconstruction,
                        times=self._times,
                        min_time=self.min_time,
                        max_time=self.max_time,
                        time_step=self.time_step,
                    ),
                    continent_features_lists,
                )

            # Merge output lists back into one list.
            return list(itertools.chain.from_iterable(reconstructed_continents_list))

        else:
            return _reconstruct_and_deform_continent_features_impl(
                self.continent_features,
                plate_reconstruction,
                self._times,
                self.min_time,
                self.max_time,
                self.time_step,
            )


class _ContinentalPolygon(object):

    def __init__(self, feature_index, property_name, present_day_polygon, times):
        self.feature_index = feature_index
        self.property_name = property_name
        self.reconstructed_polygon = present_day_polygon

        # Get the polygon boundary points (exterior ring).
        self.reconstructed_points = list(
            self.reconstructed_polygon.get_exterior_ring_points()
        )

        # Storage for the reconstructed geometry points (as lat/lon) over the time range.
        self.reconstructed_polygon_time_span = np.empty(
            (len(times), len(self.reconstructed_points), 2),
            dtype=float,
        )

    def update_polygon(self):
        self.reconstructed_polygon = pygplates.PolygonOnSphere(  # type:ignore
            self.reconstructed_points
        )


class _ContinentAggregate(object):
    def __init__(self, stage_rotation, continental_polygon):
        self.stage_rotation = stage_rotation
        self.continental_polygons = [continental_polygon]
        self.attached_resolved_networks = []

    def add(self, new_continental_polygon, new_stage_rotation):
        # new_stage_rotation_pole, new_stage_rotation_angle_radians = (
        #    new_stage_rotation.get_euler_pole_and_angle()
        # )
        # stage_rotation_pole, stage_rotation_angle_radians = (
        #    self.stage_rotation.get_euler_pole_and_angle()
        # )
        # if (
        #    pygplates.Vector3D.angle_between(  # type:ignore
        #        new_stage_rotation_pole.to_xyz(), stage_rotation_pole.to_xyz()
        #    )
        #    < math.radians(0.01)
        #    and new_stage_rotation == self.stage_rotation
        # ):
        if pygplates.FiniteRotation.are_equal(  # type:ignore
            new_stage_rotation, self.stage_rotation, threshold_degrees=0.01
        ):
            for continental_polygon in self.continental_polygons:
                if (
                    pygplates.GeometryOnSphere.distance(  # type:ignore
                        new_continental_polygon.reconstructed_polygon,
                        continental_polygon.reconstructed_polygon,
                        distance_threshold_radians=1e-4,
                        geometry1_is_solid=True,
                        geometry2_is_solid=True,
                    )
                    is None
                ):
                    self.continental_polygons.append(new_continental_polygon)
                    return True

        return False

    def find_attached_networks(
        self,
        resolved_networks,
        all_resolved_network_points,
        all_resolved_network_points_spatial_tree,
        all_resolved_network_point_velocities,
        resolved_network_point_ranges,
    ):
        #
        # Find all networks that are attached to this continent aggregrate:
        # - For each network find any boundary points that are inside any continent polygons in this aggregrate.
        # - If any two consecutive boundary points have the same velocity (at same location) as rigid continent velocity:
        #   + then that network is attached to this aggregrate.
        #

        resolved_network_boundary_points_inside_aggregate = (
            points_in_polygons.find_polygons_using_points_spatial_tree(
                all_resolved_network_points,
                all_resolved_network_points_spatial_tree,
                [
                    continental_polygon.reconstructed_polygon
                    for continental_polygon in self.continental_polygons
                ],
            )
        )

        resolved_network_index = 0
        point_index = 0
        while point_index < len(all_resolved_network_points):
            _, resolved_network_point_end_index = resolved_network_point_ranges[
                resolved_network_index
            ]

            resolved_network_is_attached = False
            if (
                resolved_network_boundary_points_inside_aggregate[point_index]
                is not None
            ):
                resolved_network_point = all_resolved_network_points[point_index]
                resolved_network_point_velocity = all_resolved_network_point_velocities[
                    point_index
                ]
                aggregrate_velocity_at_point = pygplates.calculate_velocities(  # type:ignore
                    (resolved_network_point,),
                    # Stage rotation goes backward in time but velocity needs to go forward...
                    self.stage_rotation.get_inverse(),
                    1.0,  # time interval
                )[
                    0
                ]
                resolved_network_is_attached = (
                    aggregrate_velocity_at_point - resolved_network_point_velocity
                ).get_magnitude() < 1.0

            if resolved_network_is_attached:
                # Current resolved network boundary point is *inside* this aggregate.
                # So add current network as an attached network.
                self.attached_resolved_networks.append(
                    resolved_networks[resolved_network_index]
                )
                # Move to the next network.
                point_index = (
                    resolved_network_point_end_index  # same as start of next network
                )
                resolved_network_index += 1
            else:
                # Current resolved network boundary point is *outside* this aggregate.
                # So just move to the next point.
                point_index += 1
                if point_index == resolved_network_point_end_index:
                    resolved_network_index += 1

    def reconstruct_continental_polygons(self, time):
        #
        # For each continental polygon, reconstruct using attached networks (or rigid stage rotations):
        # - For each point in a continental polygon determine if inside an attached network.
        # - If so then reconstruct using it, otherwise reconstruct using this aggregate's rigid rotation.
        #

        continent_aggregrate_points = []
        continental_polygon_point_ranges = []
        for continental_polygon in self.continental_polygons:
            continental_polygon_start_point_index = len(continent_aggregrate_points)
            continent_aggregrate_points.extend(continental_polygon.reconstructed_points)
            continental_polygon_end_point_index = len(continent_aggregrate_points)

            continental_polygon_point_ranges.append(
                (
                    continental_polygon_start_point_index,
                    continental_polygon_end_point_index,
                )
            )

        continent_aggregrate_points_inside_attached_resolved_networks = (
            points_in_polygons.find_polygons(
                continent_aggregrate_points,
                [
                    resolved_network.get_resolved_boundary()
                    for resolved_network in self.attached_resolved_networks
                ],
                self.attached_resolved_networks,
            )
        )

        continental_polygon_index = 0
        point_index = 0
        while point_index < len(continent_aggregrate_points):
            continental_polygon = self.continental_polygons[continental_polygon_index]

            (
                continental_polygon_point_start_index,
                continental_polygon_point_end_index,
            ) = continental_polygon_point_ranges[continental_polygon_index]

            continental_polygon_point_index = (
                point_index - continental_polygon_point_start_index
            )

            resolved_network_containing_point = (
                continent_aggregrate_points_inside_attached_resolved_networks[
                    point_index
                ]
            )

            reconstructed_point = continental_polygon.reconstructed_points[
                continental_polygon_point_index
            ]

            if resolved_network_containing_point is not None:
                # Current continental polygon point is *inside* an attached resolved network.
                # So reconstruct it using the attached resolved network.
                reconstructed_point = (
                    resolved_network_containing_point.reconstruct_point(
                        reconstructed_point,
                        time,
                    )
                )
            else:
                reconstructed_point = self.stage_rotation * reconstructed_point

            continental_polygon.reconstructed_points[
                continental_polygon_point_index
            ] = reconstructed_point

            # Move to the next point.
            point_index += 1
            if point_index == continental_polygon_point_end_index:
                continental_polygon_index += 1

        # Now that we've reconstructed all the continental polygon points
        # we can create reconstructed polygons from those reconstructed points.
        for continental_polygon in self.continental_polygons:
            continental_polygon.update_polygon()

    def record_reconstructed_points(self, time_index):
        for continental_polygon in self.continental_polygons:
            continental_polygon.reconstructed_polygon_time_span[time_index] = (
                np.fromiter(
                    (
                        point.to_lat_lon()
                        for point in continental_polygon.reconstructed_points
                    ),
                    dtype=np.dtype((float, 2)),
                    count=len(continental_polygon.reconstructed_points),
                )
            )


def _reconstruct_and_deform_continent_features_impl(
    continent_features, plate_reconstruction, times, min_time, max_time, time_step
):

    continental_polygons = []
    for feature_index, feature in enumerate(continent_features):
        continent_plate_id = continent_features[
            feature_index
        ].get_reconstruction_plate_id()
        # Adjustment for present day geometries in case present day rotation is non-zero
        # (generally it shouldn't be though).
        present_day_rotation = plate_reconstruction.rotation_model.get_rotation(
            0.0, continent_plate_id
        )
        for property in feature:
            property_value = property.get_value()
            if property_value:
                geometry = property_value.get_geometry()
                if geometry and isinstance(
                    geometry, pygplates.PolygonOnSphere  # type:ignore
                ):
                    continental_polygon = _ContinentalPolygon(
                        feature_index,
                        property.get_name(),
                        present_day_rotation * geometry,
                        times,
                    )
                    continental_polygons.append(continental_polygon)

    initial_times = np.arange(time_step, min_time - 1e-6, time_step)
    reconstruction_times = np.concatenate((initial_times, times))

    previous_time = 0.0
    time_index = 0
    for time in reconstruction_times:

        topological_snapshot = pygplates.TopologicalSnapshot(  # type:ignore
            plate_reconstruction.topology_features,
            plate_reconstruction.rotation_model,
            previous_time,
        )
        resolved_networks = topological_snapshot.get_resolved_topologies(
            pygplates.ResolveTopologyType.network  # type:ignore
        )

        #
        # 1. Find rigid aggregates:
        #    - Find groups of continental polygons with stage rotations that are equivalent.
        #    - Within each group further divide into those that are touching/overlapping each other.
        # 2. For each continent group find all networks that are attached to each group:
        #    - For each network find any boundary points that are inside any continent polygons in a group.
        #    - If any two consecutive boundary points have the same velocity (at same location) as rigid continent velocity:
        #      + then that network is attached to the group.
        # 3. For each continental polygon in a group, reconstruct using attached networks (or rigid stage rotations):
        #    - For each point in a continental polygon determine if inside an attached network.
        #    - If so then reconstruct using it, otherwise reconstruct using this aggregate's rigid rotation.
        #

        continent_aggregates = _find_continent_aggregates(
            previous_time,
            time,
            continental_polygons,
            continent_features,
            plate_reconstruction.rotation_model,
        )

        all_resolved_network_points = []
        all_resolved_network_point_velocities = []
        resolved_network_point_ranges = []
        for resolved_network in resolved_networks:

            resolved_network_start_point_index = len(all_resolved_network_points)

            all_resolved_network_points.extend(
                resolved_network.get_resolved_geometry_points()
            )
            all_resolved_network_point_velocities.extend(
                resolved_network.get_resolved_geometry_point_velocities()
            )

            resolved_network_end_point_index = len(all_resolved_network_points)
            resolved_network_point_ranges.append(
                (resolved_network_start_point_index, resolved_network_end_point_index)
            )
        all_resolved_network_points_spatial_tree = (
            points_spatial_tree.PointsSpatialTree(all_resolved_network_points)
        )

        for continent_aggregate in continent_aggregates:
            continent_aggregate.find_attached_networks(
                resolved_networks,
                all_resolved_network_points,
                all_resolved_network_points_spatial_tree,
                all_resolved_network_point_velocities,
                resolved_network_point_ranges,
            )

            continent_aggregate.reconstruct_continental_polygons(time)

            # Record the reconstructed points if we're in the [min_time, max_time] time range.
            if time >= min_time:
                continent_aggregate.record_reconstructed_points(time_index)

        if time >= min_time:
            time_index += 1

        previous_time = time

    reconstructed_continents = [{} for _ in range(len(continent_features))]

    for continental_polygon in continental_polygons:
        reconstructed_continent = reconstructed_continents[
            continental_polygon.feature_index
        ]
        if continental_polygon.property_name not in reconstructed_continent:
            reconstructed_continent[continental_polygon.property_name] = []
        reconstructed_continent[continental_polygon.property_name].append(
            continental_polygon.reconstructed_polygon_time_span
        )

    return reconstructed_continents


def _find_continent_aggregates(
    previous_time, time, continental_polygons, continent_features, rotation_model
):
    #
    # Find rigid aggregates:
    # - Find groups of continental polygons with stage rotations that are equivalent.
    # - Within each group further divide into those that are touching/overlapping each other.
    #

    continent_aggregates = []

    for continental_polygon in continental_polygons:
        continent_plate_id = continent_features[
            continental_polygon.feature_index
        ].get_reconstruction_plate_id()

        stage_rotation = rotation_model.get_rotation(
            time, continent_plate_id, previous_time
        )

        # Add continental polygon to an existing continent aggegrate, if possible.
        was_added_to_aggregrate = False
        for continent_aggregate in continent_aggregates:
            if continent_aggregate.add(continental_polygon, stage_rotation):
                was_added_to_aggregrate = True
                break

        if not was_added_to_aggregrate:
            # Create a new continent aggregrate containing the current continental polygon.
            continent_aggregate = _ContinentAggregate(
                stage_rotation, continental_polygon
            )
            continent_aggregates.append(continent_aggregate)

    return continent_aggregates
