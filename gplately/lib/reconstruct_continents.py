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

from ..gpml import _load_FeatureCollection
from ..parallel import get_num_cpus

import numpy as np
import pygplates


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
        # See if model has any topological network features.
        self._topological_model_is_deforming = self._is_topological_model_deforming(
            plate_reconstruction
        )

        # If the topological model has deforming networks then reconstruct and deform the continent features from min to max time.
        # Otherwise don't do anything - we'll just rigidly reconstruct them when the user asks for them.
        if self._topological_model_is_deforming:
            # Determine number of CPUs to use.
            num_cpus = get_num_cpus(nprocs)

            self._valid_reconstructions = (
                self._reconstruct_and_deform_continent_features(
                    plate_reconstruction, num_cpus
                )
            )
        else:
            self._valid_reconstruction_times = set(
                np.arange(self.max_time, self.min_time - 1e-4, -self.time_step)
            )

    def get_reconstructed_continents(self, time, plate_reconstruction):
        """Retrieve the reconstructed continent features.

        Parameters
        ----------
        time : float
            The reconstruction time.
        plate_reconstruction : PlateReconstruction
            A :class:`PlateReconstruction` object for reconstruction.
        """
        if self._topological_model_is_deforming:
            # Get the already-reconstructed-and-deformed continents (indexed by time).
            reconstructed_continents = self._valid_reconstructions.get(time, None)
            if reconstructed_continents is not None:
                return reconstructed_continents
        else:
            # If the time is one of the valid times then rigidly reconstruct the present-day continent features to 'time'.
            if time in self._valid_reconstruction_times:
                if not self.continent_features:
                    return []
                return plate_reconstruction.reconstruct(self.continent_features, time)

        # The requested 'time' was not one of the valid times.
        raise ValueError(
            f"{time} is not a valid time in the range [{self.max_time}, {self.min_time}] with interval {self.time_step}"
        )

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

    def _is_topological_model_deforming(self, plate_reconstruction):
        # See if model has any topological network features.
        return any(
            feature.get_feature_type()
            == pygplates.FeatureType.gpml_topological_network  # type:ignore
            for feature in plate_reconstruction.topology_features
        )

    def _reconstruct_and_deform_continent_features(
        self, plate_reconstruction, num_cpus
    ):
        # If there are no continent features then just return an empty list for each reconstruction time.
        if not self.continent_features:
            return {
                time: []
                for time in np.arange(
                    self.max_time, self.min_time - 1e-4, -self.time_step
                )
            }

        # TODO: Implement reconstruction and deformation of continental polygons.

        return {}
