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
from ..ptt.utils import points_in_polygons


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
            These are any features on continental crust (eg, continental polygons).
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

        # Make sure the continent features all reconstruct by plate ID only.
        if not self._do_all_continent_features_reconstruct_by_plate_id():
            raise ValueError(
                "All continent features must reconstruct by plate ID only (excludes by half-stage, motion paths, flowlines, etc)"
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
                    reconstructed_geometry_time_spans,
                ) in reconstructed_continent.items():
                    reconstructed_geometries = [
                        pygplates.PolygonOnSphere(time_span[time_index])
                        for time_span in reconstructed_geometry_time_spans
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

    def _do_all_continent_features_reconstruct_by_plate_id(self):
        if not self.continent_features:
            return True

        # Make sure all the continent features reconstruct by plate ID only.
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


def _reconstruct_and_deform_continent_features_impl(
    continent_features, plate_reconstruction, times, min_time, max_time, time_step
):

    reconstructed_continents = [{} for _ in range(len(continent_features))]
    for feature_index, feature in enumerate(continent_features):
        reconstructed_continent = reconstructed_continents[feature_index]
        for property in feature:
            property_value = property.get_value()
            if property_value:
                geometry = property_value.get_geometry()
                if geometry:
                    property_name = property.get_name()
                    if property_name not in reconstructed_continent:
                        reconstructed_continent[property_name] = []
                    reconstructed_continent[property_name].append(geometry)

    topological_model = pygplates.TopologicalModel(  # type:ignore
        plate_reconstruction.topology_features, plate_reconstruction.rotation_model
    )

    initial_times = np.arange(time_step, min_time - 1e-6, time_step)
    reconstruction_times = np.concatenate((initial_times, times))

    for feature_index, reconstructed_continent in enumerate(reconstructed_continents):
        continent_plate_id = continent_features[
            feature_index
        ].get_reconstruction_plate_id()

        for property_name, geometries in reconstructed_continent.items():

            reconstructed_geometry_time_spans = []

            for geometry in geometries:
                # Adjust present day geometry if present day rotation is non-zero
                # (generally it shouldn't be though).
                present_day_rotation = plate_reconstruction.rotation_model.get_rotation(
                    0.0, continent_plate_id
                )
                geometry = present_day_rotation * geometry

                if isinstance(geometry, pygplates.PolygonOnSphere):  # type:ignore
                    reconstructed_geometry_points = list(
                        geometry.get_exterior_ring_points()
                    )
                else:
                    reconstructed_geometry_points = list(geometry.get_points())

                reconstructed_geometry_time_span = np.empty(
                    (len(times), len(reconstructed_geometry_points), 2),
                    dtype=float,
                )

                previous_time = 0.0
                time_index = 0
                for time in reconstruction_times:

                    topological_snapshot = topological_model.topological_snapshot(
                        previous_time
                    )
                    resolved_networks = topological_snapshot.get_resolved_topologies(
                        pygplates.ResolveTopologyType.network  # type:ignore
                    )

                    reconstructed_geometry_point_resolved_networks = (
                        points_in_polygons.find_polygons(
                            reconstructed_geometry_points,
                            [
                                resolved_topology.get_resolved_boundary()
                                for resolved_topology in resolved_networks
                            ],
                            resolved_networks,
                        )
                    )

                    stage_rotation = plate_reconstruction.rotation_model.get_rotation(
                        time, continent_plate_id, previous_time
                    )

                    for point_index, resolved_network in enumerate(
                        reconstructed_geometry_point_resolved_networks
                    ):
                        if resolved_network:
                            reconstructed_geometry_points[point_index] = (
                                resolved_network.reconstruct_point(
                                    reconstructed_geometry_points[point_index], time
                                )
                            )
                        else:
                            reconstructed_geometry_points[point_index] = (
                                stage_rotation
                                * reconstructed_geometry_points[point_index]
                            )

                    # Record the reconstructed geometry points if we're in the [min_time, max_time] time range.
                    if time >= min_time:
                        reconstructed_geometry_time_span[time_index] = np.fromiter(
                            (
                                point.to_lat_lon()
                                for point in reconstructed_geometry_points
                            ),
                            dtype=np.dtype((float, 2)),
                            count=len(reconstructed_geometry_points),
                        )
                        time_index += 1

                    previous_time = time

                reconstructed_geometry_time_spans.append(
                    reconstructed_geometry_time_span
                )

            reconstructed_continent[property_name] = reconstructed_geometry_time_spans

    return reconstructed_continents
