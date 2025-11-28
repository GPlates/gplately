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


class ReconstructContinents(object):

    def __init__(self, plate_reconstruction, continent_features):
        self.plate_reconstruction = plate_reconstruction

        self._continent_features_pickle = continent_features  # Pickle the __init__ argument (see __getstate__/__setstate__ for details).
        self.continent_features = _load_FeatureCollection(continent_features)

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

    def get_reconstructed_continents(self, time):
        return self.plate_reconstruction.reconstruct(self.continent_features, time)
