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

"""
This sub-module contains tools that wrap up pyGPlates and Plate Tectonic Tools functionalities for reconstructing features,
working with point data, and calculating plate velocities at specific geological times.
"""

# pyright: reportMissingTypeStubs=true

import logging
import warnings
from typing import Union

import numpy as np
import pygplates
from plate_model_manager import PlateModel  # type: ignore

from . import tools as _tools
from .gpml import _load_FeatureCollection
from .ptt import separate_ridge_transform_segments

logger = logging.getLogger("gplately")


class PlateReconstruction(object):
    """Reconstruct topology features to specific geological times given a :py:attr:`~rotation_model`,
    a set of :py:attr:`~topology_features` and a set of :py:attr:`~static_polygons`.
    Topological plate velocity data at specific geological times can also be
    calculated from these reconstructed features.
    """

    def __init__(
        self,
        rotation_model,
        topology_features=None,
        static_polygons=None,
        anchor_plate_id: Union[int, None] = None,
        plate_model: Union[PlateModel, None] = None,
    ):
        """
        Parameters
        ----------
        rotation_model : str/`os.PathLike`, or instance of `pygplates.FeatureCollection`_, or `pygplates.Feature`_, or sequence of `pygplates.Feature`_, or instance of `pygplates.RotationModel`_
            A rotation model to query equivalent and/or relative topological plate rotations
            from a time in the past relative to another time in the past or to present day. Can be
            provided as a rotation filename, or rotation feature collection, or rotation feature, or
            sequence of rotation features, or a sequence (eg, a list or tuple) of any combination of
            those four types.
        topology_features : str/`os.PathLike`, or a sequence (eg, `list` or `tuple`) of instances of `pygplates.Feature`_, or a single instance of `pygplates.Feature`_, or an instance of `pygplates.FeatureCollection`_, default None
            Reconstructable topological features like trenches, ridges and transforms. Can be provided
            as an optional topology-feature filename, or sequence of features, or a single feature.
        static_polygons : str/`os.PathLike`, or instance of `pygplates.Feature`_, or sequence of `pygplates.Feature`_, or an instance of `pygplates.FeatureCollection`_, default None
            Present-day polygons whose shapes do not change through geological time. They are
            used to cookie-cut dynamic polygons into identifiable topological plates (assigned
            an ID) according to their present-day locations. Can be provided as a static polygon feature
            collection, or optional filename, or a single feature, or a sequence of
            features.
        anchor_plate_id : int, optional
            Default anchor plate ID for reconstruction.
            If not specified then uses the default anchor plate of :py:attr:`~rotation_model`.
        plate_model : PlateModel, optional
            Only if users would like the :py:class:`PlateReconstruction` object tracks the :class:`PlateModel`.


        .. _pygplates.RotationModel: https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel
        .. _pygplates.Feature: https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature
        .. _pygplates.FeatureCollection: https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection
        """
        #: A `pygplates.RotationModel <https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel>`__ object
        #: to query equivalent and/or relative topological plate rotations
        #: from a time in the past relative to another time in the past or to present day.
        self.rotation_model: Union[pygplates.RotationModel, None] = None

        # Add a warning if the rotation_model is empty
        if not rotation_model:
            logger.warning(
                "No rotation features were passed to the constructor of PlateReconstruction. The reconstruction will not work. Check your rotation file(s)."
            )

        if hasattr(rotation_model, "reconstruction_identifier"):
            self.name = rotation_model.reconstruction_identifier
        else:
            self.name = None

        if anchor_plate_id is None:
            if isinstance(rotation_model, pygplates.RotationModel):  # type: ignore
                # Use the default anchor plate of 'rotation_model'.
                self.rotation_model = rotation_model
            else:
                # Using rotation features/files, so default anchor plate is 0.
                self.rotation_model = pygplates.RotationModel(rotation_model)
        else:
            # User has explicitly specified an anchor plate ID, so let's check it.
            anchor_plate_id = self._check_anchor_plate_id(anchor_plate_id)
            # This works when 'rotation_model' is a RotationModel or rotation features/files.
            self.rotation_model = pygplates.RotationModel(  # type: ignore
                rotation_model, default_anchor_plate_id=anchor_plate_id
            )

        #: A `pygplates.FeatureCollection <https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection>`__
        #: object containing topological features like trenches, ridges and transforms.
        self.topology_features: Union[pygplates.FeatureCollection, None] = (
            _load_FeatureCollection(topology_features)
        )

        #: A `pygplates.FeatureCollection <https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection>`__
        #: object containing the present-day static polygons whose shapes do not change through geological time when reconstructed.
        self.static_polygons: Union[pygplates.FeatureCollection, None] = (
            _load_FeatureCollection(static_polygons)
        )

        #: Optional plate model name
        self.plate_model = plate_model

        # Keep a snapshot of the resolved topologies at its last requested snapshot time (and anchor plate).
        # Also keep a snapshot of the reconstructed static polygons at its the last requested snapshot time (and anchor plate)
        # which, by the way, could be a different snapshot time and anchor plate than the topological snapshot.
        #
        # This avoids having to do unnessary work if the same snapshot time (and anchor plate) is requested again.
        # But if the requested time (or anchor plate) changes then we'll create a new snapshot.
        #
        # Note: Both pygplates.TopologicalSnapshot and pygplates.ReconstructSnapshot can be pickled.
        self._topological_snapshot = None
        self._static_polygons_snapshot = None

    def __getstate__(self):
        state = self.__dict__.copy()

        # Remove the unpicklable entries.
        #
        # This includes pygplates reconstructed feature geometries and resolved topological geometries.
        # Note: PyGPlates features and features collections (and rotation models) can be pickled though.
        #

        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

        # Restore the unpicklable entries.
        #
        # This includes pygplates reconstructed feature geometries and resolved topological geometries.
        # Note: PyGPlates features and features collections (and rotation models) can be pickled though.
        #

    @property
    def anchor_plate_id(self):
        """Default anchor plate ID for reconstruction. Must be an integer >= 0."""
        # The default anchor plate comes from the RotationModel.
        if self.rotation_model:
            return self.rotation_model.get_default_anchor_plate_id()

    @anchor_plate_id.setter
    def anchor_plate_id(self, anchor_plate):
        # Note: Caller cannot specify None when setting the anchor plate.
        anchor_plate = self._check_anchor_plate_id(anchor_plate)
        # Only need to update if the anchor plate changed.
        if anchor_plate != self.anchor_plate_id:
            # Update the RotationModel (which is where the anchor plate is stored).
            # This keeps the same rotation model but just changes the anchor plate.
            self.rotation_model = pygplates.RotationModel(  # type: ignore
                self.rotation_model, default_anchor_plate_id=anchor_plate
            )

    @staticmethod
    def _check_anchor_plate_id(id):
        id = int(id)
        if id < 0:
            raise ValueError("Invalid anchor plate ID: {}".format(id))
        return id

    def _check_topology_features(self, *, include_topological_slab_boundaries=True):
        if self.topology_features is None:
            raise ValueError(
                "Topology features have not been set in this PlateReconstruction."
            )

        # If not including topological slab boundaries then remove them.
        if not include_topological_slab_boundaries:
            return [
                feature
                for feature in self.topology_features
                if feature.get_feature_type()
                != pygplates.FeatureType.gpml_topological_slab_boundary  # type: ignore
            ]

        return self.topology_features

    def topological_snapshot(
        self, time, *, anchor_plate_id=None, include_topological_slab_boundaries=True
    ):
        """Create a snapshot of resolved topologies at the specified reconstruction time.

        This returns a `pygplates.TopologicalSnapshot <https://www.gplates.org/docs/pygplates/generated/pygplates.TopologicalSnapshot>`__
        from which you can extract resolved topologies, calculate velocities at point locations, calculate plate boundary statistics, etc.

        Parameters
        ----------
        time : float, int or pygplates.GeoTimeInstant
            The geological time at which to create the topological snapshot.
        anchor_plate_id : int, optional
            The anchored plate id to use when resolving topologies.
            If not specified then uses the current anchor plate (:py:attr:`PlateReconstruction.anchor_plate_id` attribute).
        include_topological_slab_boundaries : bool, default=True
            Include topological boundary features of type ``gpml:TopologicalSlabBoundary``.
            By default all features passed into constructor ``__init__()`` are included in the snapshot.
            However setting this to False is useful when you're only interested in plate boundaries.

        Returns
        -------
        topological_snapshot : pygplates.TopologicalSnapshot
            The `topological snapshot <https://www.gplates.org/docs/pygplates/generated/pygplates.TopologicalSnapshot>`__
            at the specified ``time`` (and anchor plate).

        Raises
        ------
        ValueError
            If topology features have not been set in this :py:class:`PlateReconstruction` object.
        """
        if anchor_plate_id is None:
            anchor_plate_id = self.anchor_plate_id

        # Only need to create a new snapshot if we don't have one, or if any of the following have changed since the last snapshot:
        # - the reconstruction time,
        # - the anchor plate,
        # - whether to include topological slab boundaries or not.
        if (
            self._topological_snapshot is None
            # last snapshot time...
            or self._topological_snapshot.get_reconstruction_time()
            # use pygplates.GeoTimeInstant to get a numerical tolerance in floating-point time comparison...
            != pygplates.GeoTimeInstant(time)  # type: ignore
            # last snapshot anchor plate...
            or self._topological_snapshot.get_rotation_model().get_default_anchor_plate_id()
            != anchor_plate_id
            # whether last snapshot included topological slab boundaries...
            or self._topological_snapshot_includes_topological_slab_boundaries
            != include_topological_slab_boundaries
        ):
            # Create snapshot for current parameters.
            self._topological_snapshot = pygplates.TopologicalSnapshot(  # type: ignore
                self._check_topology_features(
                    include_topological_slab_boundaries=include_topological_slab_boundaries
                ),
                self.rotation_model,
                time,
                anchor_plate_id=anchor_plate_id,
            )

            # Parameters used for the last snapshot.
            #
            # The snapshot time and anchor plate are stored in the snapshot itself (so not added here).
            #
            # Note: These don't need to be initialised in '__init__()' as long as it sets "self._topological_snapshot = None".
            #
            # Note: If we add more parameters then perhaps create a single nested private (leading '_') class for them.
            self._topological_snapshot_includes_topological_slab_boundaries = (
                include_topological_slab_boundaries
            )

        return self._topological_snapshot

    def _check_static_polygons(self):
        # Check we have static polygons.
        #
        # Currently all available models have them, but it's possible for a user to create a PlateReconstruction without them.
        if self.static_polygons is None:
            raise ValueError(
                "Static polygons have not been set in this PlateReconstruction."
            )

        return self.static_polygons

    def static_polygons_snapshot(self, time, *, anchor_plate_id=None):
        """Create a reconstructed snapshot of the static polygons at the specified reconstruction time.

        Return a `pygplates.ReconstructSnapshot <https://www.gplates.org/docs/pygplates/generated/pygplates.ReconstructSnapshot>`__
        from which you can extract reconstructed static polygons, find reconstructed polygons containing points and calculate velocities at point locations, etc.

        Parameters
        ----------
        time : float, int or pygplates.GeoTimeInstant
            The geological time at which to create the reconstructed static polygons snapshot.
        anchor_plate_id : int, optional
            The anchored plate id to use when reconstructing the static polygons.
            If not specified then uses the current anchor plate (:py:attr:`anchor_plate_id` attribute).

        Returns
        -------
        static_polygons_snapshot : pygplates.ReconstructSnapshot
            The reconstructed static polygons `snapshot <https://www.gplates.org/docs/pygplates/generated/pygplates.ReconstructSnapshot>`__
            at the specified ``time`` (and anchor plate).

        Raises
        ------
        ValueError
            If static polygons have not been set in this :py:class:`PlateReconstruction` object.
        """
        if anchor_plate_id is None:
            anchor_plate_id = self.anchor_plate_id

        # Only need to create a new snapshot if we don't have one, or if any of the following have changed since the last snapshot:
        # - the reconstruction time,
        # - the anchor plate.
        if (
            self._static_polygons_snapshot is None
            # last snapshot time...
            or self._static_polygons_snapshot.get_reconstruction_time()
            # use pygplates.GeoTimeInstant to get a numerical tolerance in floating-point time comparison...
            != pygplates.GeoTimeInstant(time)  # type: ignore
            # last snapshot anchor plate...
            or self._static_polygons_snapshot.get_rotation_model().get_default_anchor_plate_id()
            != anchor_plate_id
        ):
            # Create snapshot for current parameters.
            self._static_polygons_snapshot = pygplates.ReconstructSnapshot(  # type: ignore
                self._check_static_polygons(),
                self.rotation_model,
                time,
                anchor_plate_id=anchor_plate_id,
            )

        return self._static_polygons_snapshot

    def divergent_convergent_plate_boundaries(
        self,
        time,
        uniform_point_spacing_radians=0.001,
        divergence_velocity_threshold=0.0,
        convergence_velocity_threshold=0.0,
        *,
        first_uniform_point_spacing_radians=None,
        anchor_plate_id=None,
        velocity_delta_time=1.0,
        velocity_delta_time_type=pygplates.VelocityDeltaTimeType.t_plus_delta_t_to_t,  # type: ignore
        velocity_units=pygplates.VelocityUnits.cms_per_yr,  # type: ignore
        earth_radius_in_kms=pygplates.Earth.mean_radius_in_kms,  # type: ignore
        include_network_boundaries=False,
        include_topological_slab_boundaries=False,
        boundary_section_filter=None,
    ):
        """Samples points uniformly along plate boundaries and calculates statistics at diverging/converging locations at a particular geological time.

        Resolves topologies at ``time``, uniformly samples all plate boundaries into points and returns two lists of
        `pygplates.PlateBoundaryStatistic <https://www.gplates.org/docs/pygplates/generated/pygplates.PlateBoundaryStatistic>`__.
        The first list represents sample points where the plates are diverging, and the second where plates are converging.

        Parameters
        ----------
        time : float
            The reconstruction time (Ma) at which to query divergent/convergent plate boundaries.
        uniform_point_spacing_radians : float, default=0.001
            The spacing between uniform points along plate boundaries (in radians).
        divergence_velocity_threshold : float, default=0.0
            Orthogonal (ie, in the direction of boundary normal) velocity threshold for ``diverging`` sample points.
            Points with an orthogonal ``diverging`` velocity above this value will be returned in ``diverging_data``.
            The default is 0.0 which removes all converging sample points (leaving only diverging points).
            This value can be negative which means a small amount of convergence is allowed for the diverging points.
            The units should match the units of ``velocity_units`` (eg, if that's cm/yr then this threshold should also be in cm/yr).
        convergence_velocity_threshold : float, default=0.0
            Orthogonal (ie, in the direction of boundary normal) velocity threshold for ``converging`` sample points.
            Points with an orthogonal ``converging`` velocity above this value will be returned in ``converging_data``.
            The default is 0.0 which removes all diverging sample points (leaving only converging points).
            This value can be negative which means a small amount of divergence is allowed for the converging points.
            The units should match the units of ``velocity_units`` (eg, if that's cm/yr then this threshold should also be in cm/yr).
        first_uniform_point_spacing_radians : float, optional
            Spacing of first uniform point in each resolved topological section (in radians) - see
            `pygplates.TopologicalSnapshot.calculate_plate_boundary_statistics() <https://www.gplates.org/docs/pygplates/generated/pygplates.topologicalsnapshot#pygplates.TopologicalSnapshot.calculate_plate_boundary_statistics>`__
            for more details. Defaults to half of ``uniform_point_spacing_radians``.
        anchor_plate_id : int, optional
            Anchor plate ID. Defaults to the current anchor plate ID (:py:attr:`anchor_plate_id` attribute).
        velocity_delta_time : float, default=1.0
            The time delta used to calculate velocities (defaults to 1 Myr).
        velocity_delta_time_type : pygplates.VelocityDeltaTimeType, default=pygplates.VelocityDeltaTimeType.t_plus_delta_t_to_t
            How the two velocity times are calculated relative to ``time`` (defaults to ``[time + velocity_delta_time, time]``).
        velocity_units : pygplates.VelocityUnits, default=pygplates.VelocityUnits.cms_per_yr
            Whether to return velocities in centimetres per year or kilometres per million years (defaults to centimetres per year).
        earth_radius_in_kms : float, default=pygplates.Earth.mean_radius_in_kms
            Radius of the Earth in kilometres.
            This is only used to calculate velocities (strain rates always use ``pygplates.Earth.equatorial_radius_in_kms``).
        include_network_boundaries : bool, default=False
            Whether to sample along network boundaries that are not also plate boundaries (defaults to False).
            If a deforming network shares a boundary with a plate then it'll get included regardless of this option.
        include_topological_slab_boundaries : bool, default=False
            Whether to sample along slab boundaries (features of type ``gpml:TopologicalSlabBoundary``).
            By default they are not sampled since they are not plate boundaries.
        boundary_section_filter
            Same as the ``boundary_section_filter`` argument in
            `pygplates.TopologicalSnapshot.calculate_plate_boundary_statistics() <https://www.gplates.org/docs/pygplates/generated/pygplates.topologicalsnapshot#pygplates.TopologicalSnapshot.calculate_plate_boundary_statistics>`__.
            Defaults to None (meaning all plate boundaries are included by default).

        Returns
        -------
        diverging_data : list of pygplates.PlateBoundaryStatistic
            The results for all uniformly sampled points along plate boundaries that are ``diverging`` relative to ``divergence_threshold``.
            The size of the returned list is equal to the number of sampled points that are ``diverging``.
            Each `pygplates.PlateBoundaryStatistic <https://www.gplates.org/docs/pygplates/generated/pygplates.PlateBoundaryStatistic>`__ is guaranteed to have a valid (ie, not None)
            `convergence velocity <https://www.gplates.org/docs/pygplates/generated/pygplates.PlateBoundaryStatistic.html#pygplates.PlateBoundaryStatistic.convergence_velocity>`__.
        converging_data : list of pygplates.PlateBoundaryStatistic
            The results for all uniformly sampled points along plate boundaries that are ``converging`` relative to ``convergence_threshold``.
            The size of the returned list is equal to the number of sampled points that are ``converging``.
            Each `pygplates.PlateBoundaryStatistic <https://www.gplates.org/docs/pygplates/generated/pygplates.PlateBoundaryStatistic>`__ is guaranteed to have a valid (ie, not None)
            `convergence velocity <https://www.gplates.org/docs/pygplates/generated/pygplates.PlateBoundaryStatistic.html#pygplates.PlateBoundaryStatistic.convergence_velocity>`__.

        Raises
        ------
        ValueError
            If topology features have not been set in this :py:class:`PlateReconstruction` object.

        Examples
        --------
        To sample diverging/converging points along plate boundaries at 50Ma:

        .. code-block:: python
            :linenos:

            diverging_data, converging_data = (
                plate_reconstruction.divergent_convergent_plate_boundaries(50)
            )

        To do the same, but restrict converging data to points where orthogonal converging velocities are greater than 0.2 cm/yr
        (with diverging data remaining unchanged with the default 0.0 threshold):

        .. code-block:: python
            :linenos:

            diverging_data, converging_data = (
                plate_reconstruction.divergent_convergent_plate_boundaries(
                    50, convergence_velocity_threshold=0.2
                )
            )

        Notes
        -----
        If you want to access all sampled points regardless of their convergence/divergence you can call :meth:`topological_snapshot()`
        and then use it to directly call `pygplates.TopologicalSnapshot.calculate_plate_boundary_statistics() <https://www.gplates.org/docs/pygplates/generated/pygplates.topologicalsnapshot#pygplates.TopologicalSnapshot.calculate_plate_boundary_statistics>`__.
        Then you can do your own analysis on the returned data:

        .. code-block:: python
            :linenos:

            plate_boundary_statistics = plate_reconstruction.topological_snapshot(
                time, include_topological_slab_boundaries=False
            ).calculate_plate_boundary_statistics(uniform_point_spacing_radians=0.001)

            for stat in plate_boundary_statistics:
                if np.isnan(stat.convergence_velocity_orthogonal):
                    continue  # missing left or right plate
                latitude, longitude = stat.boundary_point.to_lat_lon()

        """

        # Generate statistics at uniformly spaced points along plate boundaries.
        plate_boundary_statistics = self.topological_snapshot(
            time,
            anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id'
            include_topological_slab_boundaries=include_topological_slab_boundaries,
        ).calculate_plate_boundary_statistics(
            uniform_point_spacing_radians,
            first_uniform_point_spacing_radians=first_uniform_point_spacing_radians,
            velocity_delta_time=velocity_delta_time,
            velocity_delta_time_type=velocity_delta_time_type,
            velocity_units=velocity_units,
            earth_radius_in_kms=earth_radius_in_kms,
            include_network_boundaries=include_network_boundaries,
            boundary_section_filter=boundary_section_filter,
        )

        diverging_point_stats = []
        converging_point_stats = []

        for stat in plate_boundary_statistics:

            # Convergence velocity.
            #
            # Note: We use the 'orthogonal' component of velocity vector.
            convergence_velocity_orthogonal = stat.convergence_velocity_orthogonal
            # Skip current point if missing left or right plate (cannot calculate convergence).
            if np.isnan(convergence_velocity_orthogonal):
                continue

            # Add to diverging points if within the specified divergence velocity threshold.
            if -convergence_velocity_orthogonal >= divergence_velocity_threshold:
                diverging_point_stats.append(stat)

            # Add to converging points if within the specified convergence velocity threshold.
            if convergence_velocity_orthogonal >= convergence_velocity_threshold:
                converging_point_stats.append(stat)

        return diverging_point_stats, converging_point_stats

    def crustal_production_destruction_rate(
        self,
        time,
        uniform_point_spacing_radians=0.001,
        divergence_velocity_threshold_in_cms_per_yr=0.0,
        convergence_velocity_threshold_in_cms_per_yr=0.0,
        *,
        first_uniform_point_spacing_radians=None,
        velocity_delta_time=1.0,
        velocity_delta_time_type=pygplates.VelocityDeltaTimeType.t_plus_delta_t_to_t,  # type: ignore
        include_network_boundaries=False,
        include_topological_slab_boundaries=False,
        boundary_section_filter=None,
    ):
        """Calculates the total crustal production and destruction rates (in km\\ :sup:`2`/yr) of divergent and convergent
        plate boundaries at the specified geological time (Ma).

        Resolves topologies at ``time`` and uniformly samples all plate boundaries into divergent and convergent boundary points.

        Total crustal production (and destruction) rate is then calculated by accumulating divergent (and convergent) orthogonal
        velocities multiplied by their local boundary lengths.
        Velocities and lengths are scaled using the geocentric radius (at each divergent and convergent sampled point).

        Parameters
        ----------
        time : float
            The reconstruction time (Ma) at which to query divergent/convergent plate boundaries.
        uniform_point_spacing_radians : float, default=0.001
            The spacing between uniform points along plate boundaries (in radians).
        divergence_velocity_threshold_in_cms_per_yr : float, default=0.0
            Orthogonal (ie, in the direction of boundary normal) velocity threshold for *diverging* sample points.
            Points with an orthogonal *diverging* velocity above this value will accumulate crustal **production**.
            The default is 0.0 which removes all converging sample points (leaving only diverging points).
            This value can be negative which means a small amount of convergence is allowed for the diverging points.
            The units should be in cm/yr.
        convergence_velocity_threshold_in_cms_per_yr : float, default=0.0
            Orthogonal (ie, in the direction of boundary normal) velocity threshold for **converging** sample points.
            Points with an orthogonal *converging* velocity above this value will accumulate crustal *destruction*.
            The default is 0.0 which removes all diverging sample points (leaving only converging points).
            This value can be negative which means a small amount of divergence is allowed for the converging points.
            The units should be in cm/yr.
        first_uniform_point_spacing_radians : float, optional
            Spacing of first uniform point in each resolved topological section (in radians) - see
            :meth:`divergent_convergent_plate_boundaries()` for more details. Defaults to half of ``uniform_point_spacing_radians``.
        velocity_delta_time : float, default=1.0
            The time delta used to calculate velocities (defaults to 1 Myr).
        velocity_delta_time_type : pygplates.VelocityDeltaTimeType, default=pygplates.VelocityDeltaTimeType.t_plus_delta_t_to_t
            How the two velocity times are calculated relative to ``time`` (defaults to ``[time + velocity_delta_time, time]``).
        include_network_boundaries : bool, default=False
            Whether to sample along network boundaries that are not also plate boundaries (defaults to False).
            If a deforming network shares a boundary with a plate then it'll get included regardless of this option.
        include_topological_slab_boundaries : bool, default=False
            Whether to sample along slab boundaries (features of type ``gpml:TopologicalSlabBoundary``).
            By default they are **not** sampled since they are *not* plate boundaries.
        boundary_section_filter
            Same as the ``boundary_section_filter`` argument in :meth:`divergent_convergent_plate_boundaries()`.
            Defaults to ``None`` (meaning all plate boundaries are included by default).

        Returns
        -------
        total_crustal_production_rate_in_km_2_per_yr : float
            The total rate of crustal *production* at divergent plate boundaries (in km\\ :sup:`2`/yr) at the specified ``time``.
        total_crustal_destruction_rate_in_km_2_per_yr : float
            The total rate of crustal *destruction* at convergent plate boundaries (in km\\ :sup:`2`/yr) at the specified ``time``.

        Raises
        ------
        ValueError
            If topology features have not been set in this :py:class:`PlateReconstruction` object.

        Examples
        --------
        To calculate total crustal production/destruction along plate boundaries at 50Ma:

        .. code-block:: python
            :linenos:

            (
                total_crustal_production_rate_in_km_2_per_yr,
                total_crustal_destruction_rate_in_km_2_per_yr,
            ) = plate_reconstruction.crustal_production_destruction_rate(50)

        To do the same, but restrict convergence to points where orthogonal converging velocities are greater than 0.2 cm/yr
        (with divergence remaining unchanged with the default 0.0 threshold):

        .. code-block:: python
            :linenos:

            (
                total_crustal_production_rate_in_km_2_per_yr,
                total_crustal_destruction_rate_in_km_2_per_yr,
            ) = plate_reconstruction.crustal_production_destruction_rate(
                50, convergence_velocity_threshold_in_cms_per_yr=0.2
            )
        """

        # Generate statistics at uniformly spaced points along plate boundaries.
        diverging_data, converging_data = self.divergent_convergent_plate_boundaries(
            time,
            uniform_point_spacing_radians=uniform_point_spacing_radians,
            divergence_velocity_threshold=divergence_velocity_threshold_in_cms_per_yr,
            convergence_velocity_threshold=convergence_velocity_threshold_in_cms_per_yr,
            first_uniform_point_spacing_radians=first_uniform_point_spacing_radians,
            velocity_delta_time=velocity_delta_time,
            velocity_delta_time_type=velocity_delta_time_type,
            velocity_units=pygplates.VelocityUnits.cms_per_yr,
            earth_radius_in_kms=pygplates.Earth.mean_radius_in_kms,
            include_network_boundaries=include_network_boundaries,
            include_topological_slab_boundaries=include_topological_slab_boundaries,
            boundary_section_filter=boundary_section_filter,
        )

        # Total crustal production rate at divergent plate boundaries.
        total_crustal_production_rate = 0.0
        for stat in diverging_data:
            # Get actual Earth radius at current latitude.
            boundary_lat, _ = stat.boundary_point.to_lat_lon()
            earth_radius_kms = _tools.geocentric_radius(boundary_lat) / 1e3

            # Convergence velocity was calculated using pygplates.Earth.mean_radius_in_kms,
            # so adjust for actual Earth radius 'earth_radius_kms' at current latitude.
            convergence_velocity_orthogonal = stat.convergence_velocity_orthogonal * (
                earth_radius_kms / pygplates.Earth.mean_radius_in_kms
            )

            # Calculate crustal production rate at current location (in km^2/yr).
            #
            # Note: Orthogonal convergence velocity is guaranteed to be non-NaN.
            crustal_production_rate = (
                -convergence_velocity_orthogonal  # negate for divergence
                * 1e-5  # convert cm/yr to km/yr
                * stat.boundary_length  # radians
                * earth_radius_kms  # km
            )

            total_crustal_production_rate += crustal_production_rate

        # Total crustal destruction rate at convergent plate boundaries.
        total_crustal_destruction_rate = 0.0
        for stat in converging_data:
            # Get actual Earth radius at current latitude.
            boundary_lat, _ = stat.boundary_point.to_lat_lon()
            earth_radius_kms = _tools.geocentric_radius(boundary_lat) / 1e3

            # Convergence velocity was calculated using pygplates.Earth.mean_radius_in_kms,
            # so adjust for actual Earth radius 'earth_radius_kms' at current latitude.
            convergence_velocity_orthogonal = stat.convergence_velocity_orthogonal * (
                earth_radius_kms / pygplates.Earth.mean_radius_in_kms
            )

            # Calculate crustal destruction rate at current location (in km^2/yr).
            #
            # Note: Orthogonal convergence velocity is guaranteed to be non-NaN.
            crustal_destruction_rate = (
                convergence_velocity_orthogonal
                * 1e-5  # convert cm/yr to km/yr
                * stat.boundary_length  # radians
                * earth_radius_kms  # km
            )

            total_crustal_destruction_rate += crustal_destruction_rate

        return total_crustal_production_rate, total_crustal_destruction_rate

    def _subduction_convergence(
        self,
        time,
        uniform_point_spacing_radians,
        velocity_delta_time,
        anchor_plate_id,
        include_network_boundaries,
        convergence_threshold_in_cm_per_yr,
        output_distance_to_nearest_edge_of_trench=False,
        output_distance_to_start_edge_of_trench=False,
        output_convergence_velocity_components=False,
        output_trench_absolute_velocity_components=False,
        output_subducting_absolute_velocity=False,
        output_subducting_absolute_velocity_components=False,
        output_trench_normal=False,
    ):
        #
        # This is essentially a replacement for 'ptt.subduction_convergence.subduction_convergence()'.
        #
        # Instead of calculating convergence along subduction zones using subducting and overriding plate IDs,
        # it uses pyGPlates 1.0 functionality that calculates statistics along plate boundaries
        # (such as plate velocities, from which convergence velocity can be obtained).
        #
        # Note that this function has an advantage over 'ptt.subduction_convergence.subduction_convergence()':
        #   It does not reject subducting boundaries that have more than one (or even zero) subducting plates (or subducting networks),
        #   which can happen if the topological model was built incorrectly (eg, mislabelled plate boundaries).
        #   As long as there's at least one plate (or network) on the subducting side then it can find it
        #   (even if the plate is not directly attached to the subduction zone, ie, doesn't specify it as part of its boundary).
        # However, like 'ptt.subduction_convergence.subduction_convergence()', it only samples plate boundaries that have a
        # subduction polarity (eg, subduction zones) since we still need to know which plates are subducting and overriding,
        # and hence cannot calculate convergence over all plate boundaries.

        # Restrict plate boundaries to those that have a subduction polarity.
        # This is just an optimisation to avoid unnecessarily sampling all plate boundaries.
        def _boundary_section_filter_function(resolved_topological_section):
            return (
                resolved_topological_section.get_feature().get_enumeration(
                    pygplates.PropertyName.gpml_subduction_polarity
                )
                is not None
            )

        # Generate statistics at uniformly spaced points along plate boundaries.
        plate_boundary_statistics_dict = self.topological_snapshot(
            time,
            anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
            # Ignore topological slab boundaries since they are not *plate* boundaries
            # (a slab edge could have a subduction polarity, and would otherwise get included)...
            include_topological_slab_boundaries=False,
        ).calculate_plate_boundary_statistics(
            uniform_point_spacing_radians,
            first_uniform_point_spacing_radians=0,
            velocity_delta_time=velocity_delta_time,
            velocity_units=pygplates.VelocityUnits.cms_per_yr,
            include_network_boundaries=include_network_boundaries,
            boundary_section_filter=_boundary_section_filter_function,
            return_shared_sub_segment_dict=True,
        )

        subduction_data = []

        # Iterate over the shared boundary sub-segments (each one will have a list of uniform points).
        for (
            shared_sub_segment,
            shared_sub_segment_stats,
        ) in plate_boundary_statistics_dict.items():

            # Find the subduction plate of the current shared boundary sub-segment.
            subducting_plate_and_polarity = shared_sub_segment.get_subducting_plate(
                return_subduction_polarity=True,
                enforce_single_plate=False,
            )
            # Skip current shared boundary sub-segment if it doesn't have a valid subduction polarity.
            #
            # Note: There might not even be a subducting plate directly attached, but that's fine because
            #       we're only interested in the subduction polarity. Later we'll get the subducting plate
            #       from the plate boundary statistics instead (since that's more reliable).
            if not subducting_plate_and_polarity:
                continue
            _, subduction_polarity = subducting_plate_and_polarity

            if subduction_polarity == "Left":
                overriding_plate_is_on_left = True
            else:
                overriding_plate_is_on_left = False

            # Iterate over the uniform points of the current shared boundary sub-segment.
            for stat in shared_sub_segment_stats:
                # Find subducting plate velocity (opposite to overriding plate).
                if overriding_plate_is_on_left:
                    subducting_plate_velocity = stat.right_plate_velocity
                else:
                    subducting_plate_velocity = stat.left_plate_velocity
                # Reject point if there's no subducting plate (or network).
                if subducting_plate_velocity is None:
                    continue

                # The convergence velocity is actually that of the subducting plate relative to the trench line.
                # It's not the right plate relative to the left (or vice versa).
                convergence_velocity = (
                    subducting_plate_velocity - stat.boundary_velocity
                )

                # Get the trench normal (and azimuth).
                trench_normal = stat.boundary_normal
                trench_normal_azimuth = stat.boundary_normal_azimuth
                # If the trench normal (in direction of overriding plate) is opposite the boundary line normal
                # (which is to the left) then flip it.
                if not overriding_plate_is_on_left:
                    trench_normal = -trench_normal
                    trench_normal_azimuth -= np.pi
                    # Keep in the range [0, 2*pi].
                    if trench_normal_azimuth < 0:
                        trench_normal_azimuth += 2 * np.pi

                # If requested, reject point if it's not converging within specified threshold.
                if convergence_threshold_in_cm_per_yr is not None:
                    # Note that we use the 'orthogonal' component of velocity vector.
                    if (
                        pygplates.Vector3D.dot(convergence_velocity, trench_normal)
                        < convergence_threshold_in_cm_per_yr
                    ):
                        continue

                # Convergence velocity magnitude and obliquity.
                if convergence_velocity.is_zero_magnitude():
                    convergence_velocity_magnitude = 0
                    convergence_obliquity = 0
                else:
                    convergence_velocity_magnitude = (
                        convergence_velocity.get_magnitude()
                    )
                    convergence_obliquity = pygplates.Vector3D.angle_between(
                        convergence_velocity, trench_normal
                    )

                    # The direction towards which we rotate from the trench normal in a clockwise fashion.
                    clockwise_direction = pygplates.Vector3D.cross(
                        trench_normal, stat.boundary_point.to_xyz()
                    )
                    # Anti-clockwise direction has range (0, -pi) instead of (0, pi).
                    if (
                        pygplates.Vector3D.dot(
                            convergence_velocity, clockwise_direction
                        )
                        < 0
                    ):
                        convergence_obliquity = -convergence_obliquity

                    # See if plates are diverging (moving away from each other).
                    # If plates are diverging (moving away from each other) then make the
                    # velocity magnitude negative to indicate this. This could be inferred from
                    # the obliquity but it seems this is the standard way to output convergence rate.
                    #
                    # Note: This is the same as done in 'ptt.subduction_convergence.subduction_convergence()'.
                    if pygplates.Vector3D.dot(convergence_velocity, trench_normal) < 0:
                        convergence_velocity_magnitude = -convergence_velocity_magnitude

                # Trench absolute velocity magnitude and obliquity.
                trench_absolute_velocity_magnitude = stat.boundary_velocity_magnitude
                trench_absolute_velocity_obliquity = stat.boundary_velocity_obliquity

                # If the trench normal (in direction of overriding plate) is opposite the boundary line normal (which is to the left)
                # then we need to flip the obliquity of the trench absolute velocity vector. This is because it's currently relative
                # to the boundary line normal but needs to be relative to the trench normal.
                if not overriding_plate_is_on_left:
                    trench_absolute_velocity_obliquity -= np.pi
                    # Keep obliquity in the range [-pi, pi].
                    if trench_absolute_velocity_obliquity < -np.pi:
                        trench_absolute_velocity_obliquity += 2 * np.pi

                # See if the trench absolute motion is heading in the direction of the
                # overriding plate. If it is then make the velocity magnitude negative to
                # indicate this. This could be inferred from the obliquity but it seems this
                # is the standard way to output trench velocity magnitude.
                #
                # Note that we are not calculating the motion of the trench
                # relative to the overriding plate - they are usually attached to each other
                # and hence wouldn't move relative to each other.
                #
                # Note: This is the same as done in 'ptt.subduction_convergence.subduction_convergence()'.
                if np.abs(trench_absolute_velocity_obliquity) < 0.5 * np.pi:
                    trench_absolute_velocity_magnitude = (
                        -trench_absolute_velocity_magnitude
                    )

                lat, lon = stat.boundary_point.to_lat_lon()

                # The plate ID along the trench line.
                #
                # Note: The plate IDs along the trench line and overriding plate ID can differ even in a non-deforming model
                #       due to smaller plates, not modelled by topologies, moving differently than the larger topological
                #       plate being modelled - and the trench line having plate IDs of the smaller plates near them.
                #       For that reason we use the plate IDs of the trench line (rather than the overriding plate ID).
                #
                # Note: Using 'pygplates.PlateBoundaryStatistic.boundary_feature' means that if the current shared sub-segment
                #       is part of a topological line then this will obtain a plate ID from whichever sub-segment, of the current
                #       shared sub-segment, the current boundary point lies on (rather than obtaining a plate ID from the
                #       shared sub-segment itself). This is because a trench line that is a topological line might actually be
                #       deforming (or intended to be deforming) and hence its plate ID is not meaningful, or at least we can't
                #       be sure whether it will be zero or the overriding plate (or something else).
                #
                trench_plate_id = stat.boundary_feature.get_reconstruction_plate_id()

                if overriding_plate_is_on_left:
                    subducting_plate = stat.right_plate
                else:
                    subducting_plate = stat.left_plate

                # Get the subducting plate ID from resolved topological boundary (or network).
                if subducting_plate.located_in_resolved_boundary():
                    subducting_plate_id = (
                        subducting_plate.located_in_resolved_boundary()
                        .get_feature()
                        .get_reconstruction_plate_id()
                    )
                else:
                    subducting_plate_id = (
                        subducting_plate.located_in_resolved_network()
                        .get_feature()
                        .get_reconstruction_plate_id()
                    )

                output_tuple = (
                    lon,
                    lat,
                    convergence_velocity_magnitude,
                    np.degrees(convergence_obliquity),
                    trench_absolute_velocity_magnitude,
                    np.degrees(trench_absolute_velocity_obliquity),
                    np.degrees(stat.boundary_length),
                    np.degrees(trench_normal_azimuth),
                    subducting_plate_id,
                    trench_plate_id,
                )

                if output_distance_to_nearest_edge_of_trench:
                    distance_to_nearest_edge_of_trench = min(
                        stat.distance_from_start_of_topological_section,
                        stat.distance_to_end_of_topological_section,
                    )
                    output_tuple += (np.degrees(distance_to_nearest_edge_of_trench),)

                if output_distance_to_start_edge_of_trench:
                    # We want the distance to be along the clockwise direction around the overriding plate.
                    if overriding_plate_is_on_left:
                        # The overriding plate is on the left of the trench.
                        # So the clockwise direction starts at the end of the trench.
                        distance_to_start_edge_of_trench = (
                            stat.distance_to_end_of_topological_section
                        )
                    else:
                        # The overriding plate is on the right of the trench.
                        # So the clockwise direction starts at the beginning of the trench.
                        distance_to_start_edge_of_trench = (
                            stat.distance_from_start_of_topological_section
                        )
                    output_tuple += (np.degrees(distance_to_start_edge_of_trench),)

                if output_convergence_velocity_components:
                    # The orthogonal and parallel components are just magnitude multiplied by cosine and sine.
                    convergence_velocity_orthogonal = np.cos(
                        convergence_obliquity
                    ) * np.abs(convergence_velocity_magnitude)
                    convergence_velocity_parallel = np.sin(
                        convergence_obliquity
                    ) * np.abs(convergence_velocity_magnitude)
                    output_tuple += (
                        convergence_velocity_orthogonal,
                        convergence_velocity_parallel,
                    )

                if output_trench_absolute_velocity_components:
                    # The orthogonal and parallel components are just magnitude multiplied by cosine and sine.
                    trench_absolute_velocity_orthogonal = np.cos(
                        trench_absolute_velocity_obliquity
                    ) * np.abs(trench_absolute_velocity_magnitude)
                    trench_absolute_velocity_parallel = np.sin(
                        trench_absolute_velocity_obliquity
                    ) * np.abs(trench_absolute_velocity_magnitude)
                    output_tuple += (
                        trench_absolute_velocity_orthogonal,
                        trench_absolute_velocity_parallel,
                    )

                if (
                    output_subducting_absolute_velocity
                    or output_subducting_absolute_velocity_components
                ):
                    # Subducting absolute velocity magnitude and obliquity.
                    #
                    # Note: Subducting plate is opposite the overriding plate.
                    if overriding_plate_is_on_left:
                        subducting_absolute_velocity_magnitude = (
                            stat.right_plate_velocity_magnitude
                        )
                        subducting_absolute_velocity_obliquity = (
                            stat.right_plate_velocity_obliquity
                        )
                    else:
                        subducting_absolute_velocity_magnitude = (
                            stat.left_plate_velocity_magnitude
                        )
                        subducting_absolute_velocity_obliquity = (
                            stat.left_plate_velocity_obliquity
                        )
                        # Flip obliquity since trench normal (towards overidding plate on right)
                        # is opposite the boundary line normal (towards left).
                        subducting_absolute_velocity_obliquity -= np.pi
                        # Keep obliquity in the range [-pi, pi].
                        if subducting_absolute_velocity_obliquity < -np.pi:
                            subducting_absolute_velocity_obliquity += 2 * np.pi

                    # Similar to the trench absolute motion, if subducting absolute motion is heading
                    # in the direction of the overriding plate then make the velocity magnitude negative.
                    if np.abs(subducting_absolute_velocity_obliquity) < 0.5 * np.pi:
                        subducting_absolute_velocity_magnitude = (
                            -subducting_absolute_velocity_magnitude
                        )

                    if output_subducting_absolute_velocity:
                        output_tuple += (
                            subducting_absolute_velocity_magnitude,
                            np.degrees(subducting_absolute_velocity_obliquity),
                        )
                    if output_subducting_absolute_velocity_components:
                        # The orthogonal and parallel components are just magnitude multiplied by cosine and sine.
                        subducting_absolute_velocity_orthogonal = np.cos(
                            subducting_absolute_velocity_obliquity
                        ) * np.abs(subducting_absolute_velocity_magnitude)
                        subducting_absolute_velocity_parallel = np.sin(
                            subducting_absolute_velocity_obliquity
                        ) * np.abs(subducting_absolute_velocity_magnitude)
                        output_tuple += (
                            subducting_absolute_velocity_orthogonal,
                            subducting_absolute_velocity_parallel,
                        )

                if output_trench_normal:
                    output_tuple += trench_normal.to_xyz()

                subduction_data.append(output_tuple)

        return subduction_data

    def tessellate_subduction_zones(
        self,
        time,
        tessellation_threshold_radians=0.001,
        ignore_warnings=False,
        return_geodataframe=False,
        *,
        use_ptt=False,
        include_network_boundaries=False,
        convergence_threshold_in_cm_per_yr=None,
        anchor_plate_id=None,
        velocity_delta_time=1.0,
        output_distance_to_nearest_edge_of_trench=False,
        output_distance_to_start_edge_of_trench=False,
        output_convergence_velocity_components=False,
        output_trench_absolute_velocity_components=False,
        output_subducting_absolute_velocity=False,
        output_subducting_absolute_velocity_components=False,
        output_trench_normal=False,
    ):
        """Samples points along subduction zone trenches and obtains subduction data at a particular geological time.

        Resolves topologies at ``time`` and tessellates all resolved subducting features into points.

        Returns a 10-column vertically-stacked tuple with the following data per sampled trench point:

        * Col. 0 - longitude of sampled trench point
        * Col. 1 - latitude of sampled trench point
        * Col. 2 - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
        * Col. 3 - subducting convergence velocity obliquity angle in degrees (angle between trench normal vector and convergence velocity vector)
        * Col. 4 - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
        * Col. 5 - trench absolute velocity obliquity angle in degrees (angle between trench normal vector and trench absolute velocity vector)
        * Col. 6 - length of arc segment (in degrees) that current point is on
        * Col. 7 - trench normal (in subduction direction, ie, towards overriding plate) azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
        * Col. 8 - subducting plate ID
        * Col. 9 - trench plate ID

        The optional ``output_*`` parameters can be used to append extra data to the output tuple of each sampled trench point.
        The order of any extra data is the same order in which the parameters are listed below.

        Parameters
        ----------
        time : float
            The reconstruction time (Ma) at which to query subduction convergence.
        tessellation_threshold_radians : float, default=0.001
            The threshold sampling distance along the plate boundaries (in radians).
        ignore_warnings : bool, default=False
            Choose to ignore warnings from :py:func:`gplately.subduction_convergence` (if ``use_ptt`` is ``True``).
        return_geodataframe : bool, default=False
            Choose to return data in a ``geopandas.GeoDataFrame``.
        use_ptt : bool, default=False
            If set to ``True`` then uses :py:func:`gplately.subduction_convergence` to calculate subduction convergence
            (which uses the subducting stage rotation of the subduction/trench plate IDs calculate subducting velocities).
            If set to ``False`` then uses plate convergence to calculate subduction convergence
            (which samples velocities of the two adjacent boundary plates at each sampled point to calculate subducting velocities).
            Both methods ignore plate boundaries that do not have a subduction polarity (feature property), which essentially means
            they only sample subduction zones.
        include_network_boundaries : bool, default=False
            Whether to calculate subduction convergence along network boundaries that are not also plate boundaries (defaults to False).
            If a deforming network shares a boundary with a plate then it'll get included regardless of this option.
            Since subduction zones occur along *plate* boundaries this would only be an issue if an intra-plate network boundary was incorrectly labelled as subducting.
        convergence_threshold_in_cm_per_yr : float, optional
            Only return sample points with an orthogonal (ie, in the subducting geometry's normal direction) converging velocity above this value (in cm/yr).
            For example, setting this to ``0.0`` would remove all diverging sample points (leaving only converging points).
            This value can be negative which means a small amount of divergence is allowed.
            If ``None`` then all (converging and diverging) sample points are returned. This is the default.
            Note that this parameter can only be specified if ``use_ptt`` is ``False``.
        anchor_plate_id : int, optional
            Anchor plate ID. Defaults to the current anchor plate ID (:py:attr:`PlateReconstruction.anchor_plate_id` attribute).
        velocity_delta_time : float, default=1.0
            Velocity delta time used in convergence velocity calculations (defaults to 1 Myr).
        output_distance_to_nearest_edge_of_trench : bool, default=False
            Append the distance (in degrees) along the trench line to the nearest trench edge to each returned sample point.
            A trench edge is the farthermost location on the current trench feature that contributes to a plate boundary.
        output_distance_to_start_edge_of_trench : bool, default=False
            Append the distance (in degrees) along the trench line from the start edge of the trench to each returned sample point.
            The start of the trench is along the clockwise direction around the overriding plate.
        output_convergence_velocity_components : bool, default=False
            Append the convergence velocity orthogonal and parallel components (in cm/yr) to each returned sample point.
            Orthogonal is normal to trench (in direction of overriding plate when positive).
            Parallel is along trench (90 degrees clockwise from trench normal when positive).
        output_trench_absolute_velocity_components : bool, default=False
            Append the trench absolute velocity orthogonal and parallel components (in cm/yr) to each returned sample point.
            Orthogonal is normal to trench (in direction of overriding plate when positive).
            Parallel is along trench (90 degrees clockwise from trench normal when positive).
        output_subducting_absolute_velocity : bool, default=False
            Append the subducting plate absolute velocity magnitude (in cm/yr) and obliquity angle (in degrees) to each returned sample point.
        output_subducting_absolute_velocity_components : bool, default=False
            Append the subducting plate absolute velocity orthogonal and parallel components (in cm/yr) to each returned sample point.
            Orthogonal is normal to trench (in direction of overriding plate when positive).
            Parallel is along trench (90 degrees clockwise from trench normal when positive).
        output_trench_normal : bool, default=False
            Append the x, y and z components of the trench normal unit-length 3D vectors.
            These vectors are normal to the trench in the direction of subduction (towards overriding plate).
            These are global 3D vectors which differ from trench normal azimuth angles (ie, angles relative to North).

        Returns
        -------
        subduction_data : a list of vertically-stacked tuples
            The results for all tessellated points sampled along the trench.
            The size of the returned list is equal to the number of tessellated points.
            Each tuple in the list corresponds to a tessellated point and has the following tuple items:

            * Col. 0 - longitude of sampled trench point
            * Col. 1 - latitude of sampled trench point
            * Col. 2 - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
            * Col. 3 - subducting convergence velocity obliquity angle in degrees (angle between trench normal vector and convergence velocity vector)
            * Col. 4 - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
            * Col. 5 - trench absolute velocity obliquity angle in degrees (angle between trench normal vector and trench absolute velocity vector)
            * Col. 6 - length of arc segment (in degrees) that current point is on
            * Col. 7 - trench normal (in subduction direction, ie, towards overriding plate) azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
            * Col. 8 - subducting plate ID
            * Col. 9 - trench plate ID

            The optional ``output_*`` parameters can be used to append extra data to the tuple of each sampled trench point.
            The order of any extra data is the same order in which the parameters are listed in this function.

        Raises
        ------
        ValueError
            If topology features have not been set in this :py:class:`PlateReconstruction` object.
        ValueError
            If ``use_ptt`` is ``True`` and ``convergence_threshold_in_cm_per_yr`` is not ``None``.

        Notes
        -----
        If ``use_ptt`` is False then each trench is sampled at **exactly** uniform intervals along its length such that the sampled points
        have a uniform spacing (along each trench polyline) that is **equal** to ``tessellation_threshold_radians``.
        If ``use_ptt`` is True then each trench is sampled at **approximately** uniform intervals along its length such that the sampled points
        have a uniform spacing (along each trench polyline) that is **less than or equal to** ``tessellation_threshold_radians``.

        The trench normal (at each sampled trench point) always points **towards** the overriding plate.
        The obliquity angles are in the range (-180, 180). The range (0, 180) goes clockwise (when viewed from above the Earth)
        from the trench normal direction to the velocity vector. The range (0, -180) goes counter-clockwise.
        You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.
        The trench normal is perpendicular to the trench and pointing toward the overriding plate.

        Note that the convergence velocity magnitude is negative if the plates are diverging (if convergence obliquity angle
        is greater than 90 or less than -90). And note that the trench absolute velocity magnitude is negative if the trench
        (subduction zone) is moving towards the overriding plate (if trench absolute obliquity angle is less than 90 and greater
        than -90) - note that this ignores the kinematics of the subducting plate. Similiarly for the subducting plate absolute
        velocity magnitude (if keyword argument ``output_subducting_absolute_velocity`` is True).

        The trench plate ID at each sample point can differ from the overriding plate ID.
        This is because, even in a non-deforming model, the smaller plates (not modelled by topologies) can move differently
        than the larger topological plate. So the trench line has the plate IDs of the smaller plates.

        Examples
        --------
        To sample points along subduction zones at 50Ma:

        .. code-block:: python
            :linenos:

            subduction_data = plate_reconstruction.tessellate_subduction_zones(50)

        To sample points along subduction zones at 50Ma, but only where there's convergence:

        .. code-block:: python
            :linenos:

            subduction_data = plate_reconstruction.tessellate_subduction_zones(
                50, convergence_threshold_in_cm_per_yr=0.0
            )
        """

        if use_ptt:
            from . import ptt as _ptt

            if convergence_threshold_in_cm_per_yr is not None:
                raise ValueError(
                    "Can only specify 'convergence_threshold_in_cm_per_yr' if 'use_ptt' is False."
                )

            with warnings.catch_warnings():
                if ignore_warnings:
                    warnings.simplefilter("ignore")

                subduction_data = _ptt.subduction_convergence.subduction_convergence(
                    self.rotation_model,
                    self._check_topology_features(
                        # Ignore topological slab boundaries since they are not *plate* boundaries
                        # (actually they get ignored by default in 'ptt.subduction_convergence' anyway)...
                        include_topological_slab_boundaries=False
                    ),
                    tessellation_threshold_radians,
                    time,
                    velocity_delta_time=velocity_delta_time,
                    anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
                    include_network_boundaries=include_network_boundaries,
                    output_distance_to_nearest_edge_of_trench=output_distance_to_nearest_edge_of_trench,
                    output_distance_to_start_edge_of_trench=output_distance_to_start_edge_of_trench,
                    output_convergence_velocity_components=output_convergence_velocity_components,
                    output_trench_absolute_velocity_components=output_trench_absolute_velocity_components,
                    output_subducting_absolute_velocity=output_subducting_absolute_velocity,
                    output_subducting_absolute_velocity_components=output_subducting_absolute_velocity_components,
                    output_trench_normal=output_trench_normal,
                )

        else:
            subduction_data = self._subduction_convergence(
                time,
                uniform_point_spacing_radians=tessellation_threshold_radians,
                velocity_delta_time=velocity_delta_time,
                anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
                include_network_boundaries=include_network_boundaries,
                convergence_threshold_in_cm_per_yr=convergence_threshold_in_cm_per_yr,
                output_distance_to_nearest_edge_of_trench=output_distance_to_nearest_edge_of_trench,
                output_distance_to_start_edge_of_trench=output_distance_to_start_edge_of_trench,
                output_convergence_velocity_components=output_convergence_velocity_components,
                output_trench_absolute_velocity_components=output_trench_absolute_velocity_components,
                output_subducting_absolute_velocity=output_subducting_absolute_velocity,
                output_subducting_absolute_velocity_components=output_subducting_absolute_velocity_components,
                output_trench_normal=output_trench_normal,
            )

        if subduction_data:
            subduction_data = np.vstack(subduction_data)
        else:
            # No subduction data.
            num_columns = 10
            if output_distance_to_nearest_edge_of_trench:
                num_columns += 1
            if output_distance_to_start_edge_of_trench:
                num_columns += 1
            if output_convergence_velocity_components:
                num_columns += 2
            if output_trench_absolute_velocity_components:
                num_columns += 2
            if output_subducting_absolute_velocity:
                num_columns += 2
            if output_subducting_absolute_velocity_components:
                num_columns += 2
            if output_trench_normal:
                num_columns += 3
            subduction_data = np.empty((0, num_columns))

        if return_geodataframe:
            import geopandas as gpd
            from shapely import geometry

            points = [
                geometry.Point(lon, lat)
                for lon, lat in zip(subduction_data[:, 0], subduction_data[:, 1])
            ]
            # Required data.
            gdf_data = {
                "geometry": points,
                "convergence velocity (cm/yr)": subduction_data[:, 2],
                "convergence obliquity angle (degrees)": subduction_data[:, 3],
                "trench velocity (cm/yr)": subduction_data[:, 4],
                "trench obliquity angle (degrees)": subduction_data[:, 5],
                "length (degrees)": subduction_data[:, 6],
                "trench normal angle (degrees)": subduction_data[:, 7],
                "subducting plate ID": subduction_data[:, 8],
                "trench plate ID": subduction_data[:, 9],
            }

            # Optional data.
            #
            # Note: The order must match the output order.
            optional_gdf_data_index = 10
            if output_distance_to_nearest_edge_of_trench:
                gdf_data["distance to nearest trench edge (degrees)"] = subduction_data[
                    :, optional_gdf_data_index
                ]
                optional_gdf_data_index += 1
            if output_distance_to_start_edge_of_trench:
                gdf_data["distance to start of trench edge (degrees)"] = (
                    subduction_data[:, optional_gdf_data_index]
                )
                optional_gdf_data_index += 1
            if output_convergence_velocity_components:
                gdf_data["convergence velocity orthogonal component (cm/yr)"] = (
                    subduction_data[:, optional_gdf_data_index]
                )
                gdf_data["convergence velocity parallel component (cm/yr)"] = (
                    subduction_data[:, optional_gdf_data_index + 1]
                )
                optional_gdf_data_index += 2
            if output_trench_absolute_velocity_components:
                gdf_data["trench absolute velocity orthogonal component (cm/yr)"] = (
                    subduction_data[:, optional_gdf_data_index]
                )
                gdf_data["trench absolute velocity parallel component (cm/yr)"] = (
                    subduction_data[:, optional_gdf_data_index + 1]
                )
                optional_gdf_data_index += 2
            if output_subducting_absolute_velocity:
                gdf_data["subducting absolute velocity (cm/yr)"] = subduction_data[
                    :, optional_gdf_data_index
                ]
                gdf_data["subducting absolute obliquity angle (degrees)"] = (
                    subduction_data[:, optional_gdf_data_index + 1]
                )
                optional_gdf_data_index += 2
            if output_subducting_absolute_velocity_components:
                gdf_data[
                    "subducting absolute velocity orthogonal component (cm/yr)"
                ] = subduction_data[:, optional_gdf_data_index]
                gdf_data["subducting absolute velocity parallel component (cm/yr)"] = (
                    subduction_data[:, optional_gdf_data_index + 1]
                )
                optional_gdf_data_index += 2
            if output_trench_normal:
                gdf_data["trench normal (unit-length 3D vector) x component"] = (
                    subduction_data[:, optional_gdf_data_index]
                )
                gdf_data["trench normal (unit-length 3D vector) y component"] = (
                    subduction_data[:, optional_gdf_data_index + 1]
                )
                gdf_data["trench normal (unit-length 3D vector) z component"] = (
                    subduction_data[:, optional_gdf_data_index + 2]
                )
                optional_gdf_data_index += 3

            gdf = gpd.GeoDataFrame(gdf_data, geometry="geometry")
            return gdf

        else:
            return subduction_data

    def total_subduction_zone_length(
        self,
        time,
        use_ptt=False,
        ignore_warnings=False,
        *,
        include_network_boundaries=False,
        convergence_threshold_in_cm_per_yr=None,
    ):
        """Calculates the total length of all subduction zones (km) at the specified geological time (Ma).

        Resolves topologies at ``time`` and tessellates all resolved subducting features into points (see :py:meth:`tessellate_subduction_zones`).

        Total length is calculated by sampling points along the resolved subducting features (e.g. subduction zones) and accumulating their lengths
        (see :py:meth:`tessellate_subduction_zones`). Scales lengths to kilometres using the geocentric radius (at each sampled point).

        Parameters
        ----------
        time : int
            The geological time at which to calculate total subduction zone lengths.
        use_ptt : bool, default=False
            If set to ``True`` then uses :py:func:`gplately.subduction_convergence` to calculate total subduction zone length.
            If set to ``False`` then uses plate convergence instead.
            Plate convergence is the more general approach that works along all plate boundaries (not just subduction zones).
        ignore_warnings : bool, default=False
            Choose to ignore warnings from :py:func:`gplately.subduction_convergence` (if ``use_ptt`` is ``True``).
        include_network_boundaries : bool, default=False
            Whether to count lengths along network boundaries that are not also plate boundaries (defaults to False).
            If a deforming network shares a boundary with a plate then it'll get included regardless of this option.
            Since subduction zones occur along plate boundaries this would only be an issue if an intra-plate network boundary was incorrectly labelled as subducting.
        convergence_threshold_in_cm_per_yr : float, optional
            Only count lengths associated with sample points that have an orthogonal (ie, in the subducting geometry's normal direction) converging velocity above this value (in cm/yr).
            For example, setting this to ``0.0`` would remove all diverging sample points (leaving only converging points).
            This value can be negative which means a small amount of divergence is allowed.
            If ``None`` then all (converging and diverging) sample points are counted. This is the default.
            Note that this parameter can only be specified if ``use_ptt`` is ``False``.

        Returns
        -------
        total_subduction_zone_length_kms : float
            The total subduction zone length (in km) at the specified ``time``.

        Raises
        ------
        ValueError
            If topology features have not been set in this :py:class:`PlateReconstruction` object.
        ValueError
            If ``use_ptt`` is ``True`` and ``convergence_threshold_in_cm_per_yr`` is not ``None``.

        Examples
        --------
        To calculate the total length of subduction zones at 50Ma:

        .. code-block:: python
            :linenos:

            total_subduction_zone_length_kms = plate_reconstruction.total_subduction_zone_length(50)

        To calculate the total length of subduction zones at 50Ma, but only where there's actual convergence:

        .. code-block:: python
            :linenos:

            total_subduction_zone_length_kms = plate_reconstruction.total_subduction_zone_length(50,
                    convergence_threshold_in_cm_per_yr=0.0)
        """
        subduction_data = self.tessellate_subduction_zones(
            time,
            ignore_warnings=ignore_warnings,
            use_ptt=use_ptt,
            include_network_boundaries=include_network_boundaries,
            convergence_threshold_in_cm_per_yr=convergence_threshold_in_cm_per_yr,
        )

        trench_arcseg = subduction_data[:, 6]
        trench_pt_lat = subduction_data[:, 1]

        total_subduction_zone_length_kms = 0
        for i, segment in enumerate(trench_arcseg):
            earth_radius = _tools.geocentric_radius(trench_pt_lat[i]) / 1e3
            total_subduction_zone_length_kms += np.deg2rad(segment) * earth_radius

        return total_subduction_zone_length_kms

    def total_continental_arc_length(
        self,
        time,
        continental_grid,
        trench_arc_distance,
        ignore_warnings=True,
        *,
        use_ptt=False,
        include_network_boundaries=False,
        convergence_threshold_in_cm_per_yr=None,
    ):
        """Calculates the total length of all global continental arcs (km) at the specified geological time (Ma).

        Resolves topologies at ``time`` and tessellates all resolved subducting features into points (see :py:meth:`tessellate_subduction_zones`).
        The resolved points then are projected out by the ``trench_arc_distance`` (towards overriding plate) and their new locations are
        linearly interpolated onto the supplied ``continental_grid``. If the projected trench points lie in the grid, they are considered
        continental arc points, and their arc segment lengths are appended to the total continental arc length for the specified ``time``.
        The total length is scaled to kilometres using the geocentric radius (at each sampled point).

        Parameters
        ----------
        time : int
            The geological time at which to calculate total continental arc lengths.
        continental_grid: Raster, array_like, or str
            The continental grid used to identify continental arc points. Must
            be convertible to :class:`Raster`. For an array, a global extent is
            assumed [-180,180,-90,90]. For a filename, the extent is obtained
            from the file.
        trench_arc_distance : float
            The trench-to-arc distance (in kilometres) to project sampled trench points out by in the direction of the overriding plate.
        ignore_warnings : bool, default=True
            Choose whether to ignore warning messages from :py:func:`gplately.subduction_convergence` (if ``use_ptt`` is ``True``)
            that alerts the user of subduction sub-segments that are ignored due to unidentified polarities and/or subducting plates.
        use_ptt : bool, default=False
            If set to ``True`` then uses :py:func:`gplately.subduction_convergence` to sample subducting features and their subduction polarities.
            If set to ``False`` then uses plate convergence instead.
            Plate convergence is the more general approach that works along all plate boundaries (not just subduction zones).
        include_network_boundaries : bool, default=False
            Whether to sample subducting features along network boundaries that are not also plate boundaries (defaults to False).
            If a deforming network shares a boundary with a plate then it'll get included regardless of this option.
            Since subduction zones occur along plate boundaries this would only be an issue if an intra-plate network boundary was incorrectly labelled as subducting.
        convergence_threshold_in_cm_per_yr : float, optional
            Only sample points with an orthogonal (ie, in the subducting geometry's normal direction) converging velocity above this value (in cm/yr).
            For example, setting this to ``0.0`` would remove all diverging sample points (leaving only converging points).
            This value can be negative which means a small amount of divergence is allowed.
            If ``None`` then all (converging and diverging) points are sampled. This is the default.
            Note that this parameter can only be specified if ``use_ptt`` is ``False``.

        Returns
        -------
        total_continental_arc_length_kms : float
            The continental arc length (in km) at the specified time.

        Raises
        ------
        ValueError
            If topology features have not been set in this :py:class:`gplately.PlateReconstruction` object.
        ValueError
            If ``use_ptt`` is ``True`` and ``convergence_threshold_in_cm_per_yr`` is not ``None``.

        Examples
        --------
        To calculate the total length of continental arcs at 50Ma:

        .. code-block:: python
            :linenos:

            total_continental_arc_length_kms = plate_reconstruction.total_continental_arc_length(50)

        To calculate the total length of subduction zones adjacent to continents at 50Ma, but only where there's actual convergence:

        .. code-block:: python
            :linenos:

            total_continental_arc_length_kms = (
                plate_reconstruction.total_continental_arc_length(
                    50, convergence_threshold_in_cm_per_yr=0.0
                )
            )
        """
        from . import grids as _grids

        if isinstance(continental_grid, _grids.Raster):
            graster = continental_grid
        elif isinstance(continental_grid, str):
            # Process the continental grid directory
            graster = _grids.Raster(
                data=continental_grid,
                realign=True,
                time=float(time),
            )
        else:
            # Process the masked continental grid
            try:
                continental_grid = np.array(continental_grid)
                graster = _grids.Raster(
                    data=continental_grid,
                    extent=(-180, 180, -90, 90),
                    time=float(time),
                )
            except Exception as e:
                raise TypeError(
                    "Invalid type for `continental_grid` (must be Raster,"
                    + " str, or array_like)"
                ) from e
        if (time != graster.time) and (not ignore_warnings):
            raise RuntimeWarning(
                "`continental_grid.time` ({}) ".format(graster.time)
                + "does not match `time` ({})".format(time)
            )

        # Obtain trench data.
        trench_data = self.tessellate_subduction_zones(
            time,
            ignore_warnings=ignore_warnings,
            use_ptt=use_ptt,
            include_network_boundaries=include_network_boundaries,
            convergence_threshold_in_cm_per_yr=convergence_threshold_in_cm_per_yr,
        )

        # Extract trench data
        trench_normal_azimuthal_angle = trench_data[:, 7]
        trench_arcseg = trench_data[:, 6]
        trench_pt_lon = trench_data[:, 0]
        trench_pt_lat = trench_data[:, 1]

        # Modify the trench-arc distance using the geocentric radius
        arc_distance = trench_arc_distance / (
            _tools.geocentric_radius(trench_pt_lat) / 1000
        )

        # Project trench points out along trench-arc distance, and obtain their new lat-lon coordinates
        dlon = arc_distance * np.sin(np.radians(trench_normal_azimuthal_angle))
        dlat = arc_distance * np.cos(np.radians(trench_normal_azimuthal_angle))
        ilon = trench_pt_lon + np.degrees(dlon)
        ilat = trench_pt_lat + np.degrees(dlat)

        # Linearly interpolate projected points onto continental grids, and collect the indices of points that lie
        # within the grids.
        sampled_points = graster.interpolate(
            ilon,
            ilat,
            method="linear",
            return_indices=False,
        )
        # the code below will not work if graster.interpolate() returned a tuple
        assert isinstance(sampled_points, np.ndarray)
        continental_indices = np.where(sampled_points > 0)
        point_lats = ilat[continental_indices]
        point_radii = _tools.geocentric_radius(point_lats) * 1.0e-3  # km
        # the code below will not work if trench_arcseg is GeoDataFrame
        assert isinstance(trench_arcseg, np.ndarray)
        segment_arclens = np.deg2rad(trench_arcseg[continental_indices])
        segment_lengths = point_radii * segment_arclens
        return np.sum(segment_lengths)

    def _ridge_spreading_rates(
        self,
        time,
        uniform_point_spacing_radians,
        velocity_delta_time,
        anchor_plate_id,
        spreading_feature_types,
        transform_segment_deviation_in_radians,
        include_network_boundaries,
        divergence_threshold_in_cm_per_yr,
        output_obliquity_and_normal_and_left_right_plates,
    ):
        #
        # This is essentially a replacement for 'ptt.ridge_spreading_rate.spreading_rates()'.
        #
        # Instead of calculating spreading rates along mid-ocean ridges using left/right plate IDs,
        # it uses pyGPlates 1.0 functionality that calculates statistics along plate boundaries
        # (such as plate velocities, from which divergence spreading velocity can be obtained).
        #
        # Note that this function has an advantage over 'ptt.ridge_spreading_rate.spreading_rates()'.
        # It can work on all plate boundaries, not just those that are spreading (eg, have left/right plate IDs).
        # This is because it uses plate velocities to calculate divergence (and hence spreading rates).
        #

        # Generate statistics at uniformly spaced points along plate boundaries.
        plate_boundary_statistics = self.topological_snapshot(
            time,
            anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
            # Ignore topological slab boundaries since they are not *plate* boundaries
            # (useful when 'spreading_feature_types' is None, and hence all plate boundaries are considered)...
            include_topological_slab_boundaries=False,
        ).calculate_plate_boundary_statistics(
            uniform_point_spacing_radians,
            first_uniform_point_spacing_radians=0,
            velocity_delta_time=velocity_delta_time,
            velocity_units=pygplates.VelocityUnits.cms_per_yr,
            include_network_boundaries=include_network_boundaries,
            boundary_section_filter=spreading_feature_types,
        )

        ridge_data = []

        for stat in plate_boundary_statistics:
            # Reject point if there's not a plate (or network) on both the left and right sides.
            if not stat.convergence_velocity:
                continue

            spreading_obliquity = stat.convergence_velocity_obliquity

            # If requested, reject point if it's not diverging within specified threshold.
            if divergence_threshold_in_cm_per_yr is not None:
                # Note that we use the 'orthogonal' component of velocity vector.
                if (
                    -stat.convergence_velocity_orthogonal
                    < divergence_threshold_in_cm_per_yr
                ):
                    continue

            if (
                output_obliquity_and_normal_and_left_right_plates
                or transform_segment_deviation_in_radians is not None
            ):
                # Convert obliquity from the range [-pi, pi] to [0, pi/2].
                # We're only interested in the deviation angle from the normal line (positive or negative normal direction).
                # not interested in clockwise vs anti-clockwise
                spreading_obliquity = np.abs(stat.convergence_velocity_obliquity)
                # angle relative to negative normal direction
                if spreading_obliquity > 0.5 * np.pi:
                    spreading_obliquity = np.pi - spreading_obliquity

                # If a transform segment deviation was specified then we need to reject transform segments.
                if transform_segment_deviation_in_radians is not None:
                    # Reject if spreading direction is too oblique compared to the plate boundary normal.
                    #
                    # Note: If there is zero spreading then we don't actually have an obliquity.
                    #       In which case we reject the current point to match the behaviour of
                    #       'ptt.ridge_spreading_rate.spreading_rates()' which rejects zero spreading stage rotations.
                    if (
                        stat.convergence_velocity.is_zero_magnitude()
                        or spreading_obliquity > transform_segment_deviation_in_radians
                    ):
                        continue

            lat, lon = stat.boundary_point.to_lat_lon()
            spreading_velocity = stat.convergence_velocity_magnitude

            if output_obliquity_and_normal_and_left_right_plates:
                # Get the left plate ID from resolved topological boundary (or network).
                if stat.left_plate.located_in_resolved_boundary():
                    left_plate_id = (
                        stat.left_plate.located_in_resolved_boundary()
                        .get_feature()
                        .get_reconstruction_plate_id()
                    )
                else:
                    left_plate_id = (
                        stat.left_plate.located_in_resolved_network()
                        .get_feature()
                        .get_reconstruction_plate_id()
                    )
                # Get the right plate ID from resolved topological boundary (or network).
                if stat.right_plate.located_in_resolved_boundary():
                    right_plate_id = (
                        stat.right_plate.located_in_resolved_boundary()
                        .get_feature()
                        .get_reconstruction_plate_id()
                    )
                else:
                    right_plate_id = (
                        stat.right_plate.located_in_resolved_network()
                        .get_feature()
                        .get_reconstruction_plate_id()
                    )

                ridge_data.append(
                    (
                        lon,
                        lat,
                        spreading_velocity,
                        np.degrees(spreading_obliquity),
                        np.degrees(stat.boundary_length),
                        np.degrees(stat.boundary_normal_azimuth),
                        left_plate_id,
                        right_plate_id,
                    )
                )
            else:
                ridge_data.append(
                    (
                        lon,
                        lat,
                        spreading_velocity,
                        np.degrees(stat.boundary_length),
                    )
                )

        return ridge_data

    def tessellate_mid_ocean_ridges(
        self,
        time,
        tessellation_threshold_radians=0.001,
        ignore_warnings=False,
        return_geodataframe=False,
        *,
        use_ptt=False,
        spreading_feature_types=[pygplates.FeatureType.gpml_mid_ocean_ridge],
        transform_segment_deviation_in_radians=separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS,
        include_network_boundaries=False,
        divergence_threshold_in_cm_per_yr=None,
        output_obliquity_and_normal_and_left_right_plates=False,
        anchor_plate_id=None,
        velocity_delta_time=1.0,
    ):
        """Samples points along resolved spreading features (e.g. mid-ocean ridges) and calculates spreading rates and
        lengths of ridge segments at a particular geological time.

        Resolves topologies at time and tessellates all resolved spreading features into points.

        The transform segments of spreading features are ignored (unless ``transform_segment_deviation_in_radians`` is None).

        Returns a 4-column vertically stacked tuple with the following data per sampled ridge point
        (depending on ``output_obliquity_and_normal_and_left_right_plates``):

        If ``output_obliquity_and_normal_and_left_right_plates`` is False (the default):

        * Col. 0 - longitude of sampled ridge point
        * Col. 1 - latitude of sampled ridge point
        * Col. 2 - spreading velocity magnitude (in cm/yr)
        * Col. 3 - length of arc segment (in degrees) that current point is on

        If ``output_obliquity_and_normal_and_left_right_plates`` is True:

        * Col. 0 - longitude of sampled ridge point
        * Col. 1 - latitude of sampled ridge point
        * Col. 2 - spreading velocity magnitude (in cm/yr)
        * Col. 3 - spreading obliquity in degrees (deviation from normal line in range 0 to 90 degrees)
        * Col. 4 - length of arc segment (in degrees) that current point is on
        * Col. 5 - azimuth of vector normal to the arc segment in degrees (clockwise starting at North, ie, 0 to 360 degrees)
        * Col. 6 - left plate ID
        * Col. 7 - right plate ID

        Parameters
        ----------
        time : float
            The reconstruction time (Ma) at which to query spreading rates.
        tessellation_threshold_radians : float, default=0.001
            The threshold sampling distance along the plate boundaries (in radians).
        ignore_warnings : bool, default=False
            Choose to ignore warnings from Plate Tectonic Tools' :py:func:`gplately.ridge_spreading_rate` workflow (if ``use_ptt`` is True).
        return_geodataframe : bool, default=False
            Choose to return data in a ``geopandas.GeoDataFrame``.
        use_ptt : bool, default=False
            If set to True then uses Plate Tectonic Tools' :py:func:`gplately.ridge_spreading_rate` workflow to calculate ridge spreading rates
            (which uses the spreading stage rotation of the left/right plate IDs calculate spreading velocities).
            If set to False then uses plate divergence to calculate ridge spreading rates
            (which samples velocities of the two adjacent boundary plates at each sampled point to calculate spreading velocities).
            Plate divergence is the more general approach that works along all plate boundaries (not just mid-ocean ridges).
        spreading_feature_types : `pygplates.FeatureType`_ or sequence of `pygplates.FeatureType`_, default=pygplates.FeatureType.gpml_mid_ocean_ridge
            Only sample points along plate boundaries of the specified feature types.
            Default is to only sample mid-ocean ridges.
            You can explicitly specify None to sample all plate boundaries, but note that if ``use_ptt`` is True
            then only plate boundaries that are spreading feature types are sampled
            (since Plate Tectonic Tools only works on spreading plate boundaries, eg, mid-ocean ridges).
        transform_segment_deviation_in_radians : float, default=<implementation-defined>
            How much a spreading direction can deviate from the segment normal before it's considered a transform segment (in radians).
            The default value has been empirically determined to give the best results for typical models.
            If ``None`` then the full feature geometry is used (ie, it is not split into ridge and transform segments with the transform segments getting ignored).
        include_network_boundaries : bool, default=False
            Whether to calculate spreading rate along network boundaries that are not also plate boundaries (defaults to False).
            If a deforming network shares a boundary with a plate then it'll get included regardless of this option.
            Since spreading features occur along *plate* boundaries this would only be an issue if an intra-plate network boundary was incorrectly labelled as spreading.
        divergence_threshold_in_cm_per_yr : float, optional
            Only return sample points with an orthogonal (ie, in the spreading geometry's normal direction) diverging velocity above this value (in cm/yr).
            For example, setting this to ``0.0`` would remove all converging sample points (leaving only diverging points).
            This value can be negative which means a small amount of convergence is allowed.
            If ``None`` then all (diverging and converging) sample points are returned.
            This is the default since ``spreading_feature_types`` is instead used (by default) to include only plate boundaries that are typically diverging (eg, mid-ocean ridges).
            However, setting ``spreading_feature_types`` to None (and ``transform_segment_deviation_in_radians`` to None) and explicitly specifying this parameter (eg, to 0.0)
            can be used to find points along all plate boundaries that are diverging.
            However, this parameter can only be specified if ``use_ptt`` is False.
        output_obliquity_and_normal_and_left_right_plates : bool, default=False
            Whether to also return spreading obliquity, normal azimuth and left/right plates.
        anchor_plate_id : int, optional
            Anchor plate ID. Defaults to the current anchor plate ID (:py:attr:`anchor_plate_id` attribute)..
        velocity_delta_time : float, default=1.0
            Velocity delta time used in spreading velocity calculations (defaults to 1 Myr).

        Returns
        -------
        ridge_data : a list of vertically-stacked tuples
            The results for all tessellated points sampled along the mid-ocean ridges.
            The size of the returned list is equal to the number of tessellated points.
            Each tuple in the list corresponds to a tessellated point and has the following tuple items
            (depending on ``output_obliquity_and_normal_and_left_right_plates``):

            If ``output_obliquity_and_normal_and_left_right_plates`` is False (the default):

            * longitude of sampled point
            * latitude of sampled point
            * spreading velocity magnitude (in cm/yr)
            * length of arc segment (in degrees) that sampled point is on

            If ``output_obliquity_and_normal_and_left_right_plates`` is True:

            * longitude of sampled point
            * latitude of sampled point
            * spreading velocity magnitude (in cm/yr)
            * spreading obliquity in degrees (deviation from normal line in range 0 to 90 degrees)
            * length of arc segment (in degrees) that sampled point is on
            * azimuth of vector normal to the arc segment in degrees (clockwise starting at North, ie, 0 to 360 degrees)
            * left plate ID
            * right plate ID

        Raises
        ------
        ValueError
            If topology features have not been set in this :py:class:`PlateReconstruction` object.
        ValueError
            If ``use_ptt`` is True and ``divergence_threshold_in_cm_per_yr`` is not None.

        Notes
        -----
        If ``use_ptt`` is False then each ridge segment is sampled at exactly uniform intervals along its length such that the sampled points
        have a uniform spacing (along each ridge segment polyline) that is equal to ``tessellation_threshold_radians``.
        If ``use_ptt`` is True then each ridge segment is sampled at approximately uniform intervals along its length such that the sampled points
        have a uniform spacing (along each ridge segment polyline) that is less than or equal to ``tessellation_threshold_radians``.

        Examples
        --------
        To sample points along mid-ocean ridges at 50Ma, but ignoring the transform segments (of the ridges):

        .. code-block:: python
            :linenos:

            ridge_data = plate_reconstruction.tessellate_mid_ocean_ridges(50)

        To do the same, but instead of ignoring transform segments include both ridge and transform segments,
        but only where orthogonal diverging velocities are greater than 0.2 cm/yr:

        .. code-block:: python
            :linenos:

            ridge_data = plate_reconstruction.tessellate_mid_ocean_ridges(
                50,
                transform_segment_deviation_in_radians=None,
                divergence_threshold_in_cm_per_yr=0.2,
            )

        .. _pygplates.FeatureType: https://www.gplates.org/docs/pygplates/generated/pygplates.featuretype
        """

        if use_ptt:
            from . import ptt as _ptt

            if divergence_threshold_in_cm_per_yr is not None:
                raise ValueError(
                    "Can only specify 'divergence_threshold_in_cm_per_yr' if 'use_ptt' is False."
                )

            with warnings.catch_warnings():
                if ignore_warnings:
                    warnings.simplefilter("ignore")

                ridge_data = _ptt.ridge_spreading_rate.spreading_rates(
                    self.rotation_model,
                    self._check_topology_features(
                        # Ignore topological slab boundaries since they are not *plate* boundaries
                        # (not really needed since only *spreading* feature types are considered, and
                        # they typically wouldn't get used for a slab's boundary)...
                        include_topological_slab_boundaries=False
                    ),
                    time,
                    tessellation_threshold_radians,
                    spreading_feature_types=spreading_feature_types,
                    transform_segment_deviation_in_radians=transform_segment_deviation_in_radians,
                    velocity_delta_time=velocity_delta_time,
                    anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
                    include_network_boundaries=include_network_boundaries,
                    output_obliquity_and_normal_and_left_right_plates=output_obliquity_and_normal_and_left_right_plates,
                )

        else:
            ridge_data = self._ridge_spreading_rates(
                time,
                uniform_point_spacing_radians=tessellation_threshold_radians,
                velocity_delta_time=velocity_delta_time,
                anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
                spreading_feature_types=spreading_feature_types,
                transform_segment_deviation_in_radians=transform_segment_deviation_in_radians,
                include_network_boundaries=include_network_boundaries,
                divergence_threshold_in_cm_per_yr=divergence_threshold_in_cm_per_yr,
                output_obliquity_and_normal_and_left_right_plates=output_obliquity_and_normal_and_left_right_plates,
            )

        if ridge_data:
            ridge_data = np.vstack(ridge_data)
        else:
            # No ridge data.
            if output_obliquity_and_normal_and_left_right_plates:
                ridge_data = np.empty((0, 8))
            else:
                ridge_data = np.empty((0, 4))

        if return_geodataframe:
            import geopandas as gpd
            from shapely import geometry

            points = [
                geometry.Point(lon, lat)
                for lon, lat in zip(ridge_data[:, 0], ridge_data[:, 1])
            ]
            gdf_data = {
                "geometry": points,
                "velocity (cm/yr)": ridge_data[:, 2],
            }
            if output_obliquity_and_normal_and_left_right_plates:
                gdf_data["obliquity (degrees)"] = ridge_data[:, 3]
                gdf_data["length (degrees)"] = ridge_data[:, 4]
                gdf_data["normal azimuth (degrees)"] = ridge_data[:, 5]
                gdf_data["left plate ID"] = ridge_data[:, 6]
                gdf_data["right plate ID"] = ridge_data[:, 7]
            else:
                gdf_data["length (degrees)"] = ridge_data[:, 3]
            return gpd.GeoDataFrame(gdf_data, geometry="geometry")

        else:
            return ridge_data

    def total_ridge_length(
        self,
        time,
        use_ptt=False,
        ignore_warnings=False,
        *,
        spreading_feature_types=[pygplates.FeatureType.gpml_mid_ocean_ridge],
        transform_segment_deviation_in_radians=separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS,
        include_network_boundaries=False,
        divergence_threshold_in_cm_per_yr=None,
    ):
        """Calculates the total length of all resolved spreading features (e.g. mid-ocean ridges) at the specified geological time (Ma).

        Resolves topologies at ``time`` and tessellates all resolved spreading features into points (see :py:meth:`tessellate_mid_ocean_ridges`).

        The transform segments of spreading features are ignored (unless ``transform_segment_deviation_in_radians`` is ``None``).

        Total length is calculated by sampling points along the resolved spreading features (e.g. mid-ocean ridges) and accumulating their lengths
        (see :py:meth:`tessellate_mid_ocean_ridges`). Scales lengths to kilometres using the geocentric radius (at each sampled point).

        Parameters
        ----------
        time : int
            The geological time at which to calculate total mid-ocean ridge lengths.
        use_ptt : bool, default=False
            If set to ``True`` then uses :py:func:`gplately.ridge_spreading_rate()` to calculate total ridge length
            (which uses the spreading stage rotation of the left/right plate IDs to calculate spreading directions - see ``transform_segment_deviation_in_radians``).
            If set to ``False`` then uses plate divergence to calculate total ridge length (which samples velocities of the two adjacent
            boundary plates at each sampled point to calculate spreading directions - see ``transform_segment_deviation_in_radians``).
            Plate divergence is the more general approach that works along all plate boundaries (not just mid-ocean ridges).
        ignore_warnings : bool, default=False
            Choose to ignore warnings from :py:func:`gplately.ridge_spreading_rate()` (if ``use_ptt`` is ``True``).
        spreading_feature_types : `pygplates.FeatureType`_ or sequence of `pygplates.FeatureType`_, default=pygplates.FeatureType.gpml_mid_ocean_ridge
            Only count lengths along plate boundaries of the specified feature types.
            Default is to only sample mid-ocean ridges.
            You can explicitly specify ``None`` to sample all plate boundaries, but note that if ``use_ptt`` is ``True``
            then only plate boundaries that are spreading feature types are sampled
            (since Plate Tectonic Tools only works on **spreading** plate boundaries, eg, mid-ocean ridges).
        transform_segment_deviation_in_radians : float, default=<implementation-defined>
            How much a spreading direction can deviate from the segment normal before it's considered a transform segment (in radians).
            The default value has been empirically determined to give the best results for typical models.
            If ``None`` then the full feature geometry is used (ie, it is not split into ridge and transform segments with the transform segments getting ignored).
        include_network_boundaries : bool, default=False
            Whether to count lengths along network boundaries that are not also plate boundaries (defaults to False).
            If a deforming network shares a boundary with a plate then it'll get included regardless of this option.
            Since spreading features occur along plate boundaries this would only be an issue if an intra-plate network boundary was incorrectly labelled as spreading.
        divergence_threshold_in_cm_per_yr : float, optional
            Only count lengths associated with sample points that have an orthogonal (ie, in the spreading geometry's normal direction) diverging velocity above this value (in cm/yr).
            For example, setting this to ``0.0`` would remove all converging sample points (leaving only diverging points).
            This value can be negative which means a small amount of convergence is allowed.
            If ``None`` then all (diverging and converging) sample points are counted.
            This is the default since ``spreading_feature_types`` is instead used (by default) to include only plate boundaries that are typically diverging (eg, mid-ocean ridges).
            However, setting ``spreading_feature_types`` to ``None`` (and ``transform_segment_deviation_in_radians`` to ``None``) and explicitly specifying this parameter (eg, to ``0.0``)
            can be used to count points along all plate boundaries that are diverging.
            However, this parameter can only be specified if ``use_ptt`` is ``False``.

        Returns
        -------
        total_ridge_length_kms : float
            The total length of global mid-ocean ridges (in kilometres) at the specified time.

        Raises
        ------
        ValueError
            If topology features have not been set in this :py:class:`PlateReconstruction` object.
        ValueError
            If ``use_ptt`` is ``True`` and ``divergence_threshold_in_cm_per_yr`` is not ``None``.

        Examples
        --------
        To calculate the total length of mid-ocean ridges at 50Ma, but ignoring the transform segments (of the ridges):

        .. code-block:: python
            :linenos:

            total_ridge_length_kms = plate_reconstruction.total_ridge_length(50)

        To do the same, but instead of ignoring transform segments include both ridge and transform segments,
        but only where orthogonal diverging velocities are greater than 0.2 cm/yr:

        .. code-block:: python
            :linenos:

            total_ridge_length_kms = plate_reconstruction.total_ridge_length(50,
                    transform_segment_deviation_in_radians=None,
                    divergence_threshold_in_cm_per_yr=0.2)
        """
        ridge_data = self.tessellate_mid_ocean_ridges(
            time,
            ignore_warnings=ignore_warnings,
            use_ptt=use_ptt,
            spreading_feature_types=spreading_feature_types,
            transform_segment_deviation_in_radians=transform_segment_deviation_in_radians,
            include_network_boundaries=include_network_boundaries,
            divergence_threshold_in_cm_per_yr=divergence_threshold_in_cm_per_yr,
        )

        ridge_arcseg = ridge_data[:, 3]
        ridge_pt_lat = ridge_data[:, 1]

        total_ridge_length_kms = 0
        for i, segment in enumerate(ridge_arcseg):
            earth_radius = _tools.geocentric_radius(ridge_pt_lat[i]) / 1e3
            total_ridge_length_kms += np.deg2rad(segment) * earth_radius

        return total_ridge_length_kms

    def reconstruct_snapshot(
        self,
        reconstructable_features,
        time,
        *,
        anchor_plate_id=None,
        from_time=0,
    ):
        """Create a snapshot of reconstructed regular features (including motion paths and flowlines) at a specific geological time.

        Parameters
        ----------
        reconstructable_features : str/os.PathLike, or a sequence (eg, list or tuple) of instances of `pygplates.Feature`_,
            or a single instance of `pygplates.Feature`_, or an instance of `pygplates.FeatureCollection`_
            Regular reconstructable features (including motion paths and flowlines). It can be provided as a feature collection, or
            filename, or feature, or sequence of features, or a sequence (eg, list or tuple) of any combination of those four types.

        time : float, or pygplates.GeoTimeInstant
            The specific geological time to reconstruct to.

        anchor_plate_id : int, optional
            Anchor plate ID. Defaults to the current anchor plate ID (:py:attr:`anchor_plate_id` attribute).

        from_time : float, default=0
            The specific geological time to reconstruct from. By default, this is set to present day.
            If not set to 0 Ma (present day) then the geometry in ``feature`` is assumed to be a reconstructed snapshot
            at ``from_time``, in which case it is reverse reconstructed to present day before reconstructing to ``to_time``.
            Usually features should contain present day geometry but might contain reconstructed geometry in some cases,
            such as those generated by the reconstruction export in GPlates.

        Returns
        -------
        reconstruct_snapshot : pygplates.ReconstructSnapshot
            A `pygplates.ReconstructSnapshot <https://www.gplates.org/docs/pygplates/generated/pygplates.ReconstructSnapshot>`__
            of the specified reconstructable features reconstructed using the internal rotation model to the specified reconstruction time.
        """

        # If the features represent a snapshot at a *past* geological time then we need to reverse reconstruct them
        # such that they contain present-day geometry (not reconstructed geometry).
        if from_time != 0:
            # Extract the reconstructed features and clone them so we don't modify the caller's features.
            reconstructable_features = [
                feature.clone()
                for feature in pygplates.FeaturesFunctionArgument(
                    reconstructable_features
                ).get_features()
            ]
            # Reverse reconstruct in-place (modifies each feature's geometry).

            pygplates.reverse_reconstruct(  # type: ignore
                reconstructable_features,
                self.rotation_model,
                from_time,
                anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
            )

        return pygplates.ReconstructSnapshot(
            reconstructable_features,
            self.rotation_model,
            time,
            anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
        )

    def reconstruct(
        self,
        feature,
        to_time,
        from_time=0,
        anchor_plate_id=None,
        *,
        reconstruct_type=pygplates.ReconstructType.feature_geometry,
        group_with_feature=False,
    ):
        """Reconstructs regular geological features, motion paths or flowlines to a specific geological time.

        Parameters
        ----------
        feature : str/os.PathLike, or `pygplates.FeatureCollection`_, or `pygplates.Feature`_, or sequence of `pygplates.Feature`_
            The geological features to reconstruct. It can be provided as a feature collection, or filename,
            or feature, or sequence of features, or a sequence (eg, a list or tuple) of any combination of those four types.

        to_time : float, or pygplates.GeoTimeInstant
            The specific geological time to reconstruct to.

        from_time : float, default=0
            The specific geological time to reconstruct from. By default, this is set to present day.
            If not set to 0 Ma (present day) then the geometry in ``feature`` is assumed to be a reconstructed snapshot
            at ``from_time``, in which case it is reverse reconstructed to present day before reconstructing to ``to_time``.
            Usually features should contain present day geometry but might contain reconstructed geometry in some cases,
            such as those generated by the reconstruction export in GPlates.

        anchor_plate_id : int, optional
            Anchor plate ID. Defaults to the current anchor plate ID (:py:attr:`anchor_plate_id` attribute).

        reconstruct_type : pygplates.ReconstructType, default=pygplates.ReconstructType.feature_geometry
            The specific reconstruction type to generate based on input feature geometry type. It can be provided as
            ``pygplates.ReconstructType.feature_geometry`` to only reconstruct regular feature geometries, or
            ``pygplates.ReconstructType.motion_path`` to only reconstruct motion path features, or
            ``pygplates.ReconstructType.flowline`` to only reconstruct flowline features.
            Generates ``pygplates.ReconstructedFeatureGeometry``, or ``pygplates.ReconstructedMotionPath``, or
            ``pygplates.ReconstructedFlowline`` respectively.

        group_with_feature : bool, default=False
            Used to group reconstructed geometries with their features. This can be useful when a feature has more than one
            geometry and hence more than one reconstructed geometry. The returned list then becomes a list of tuples where
            each tuple contains a pygplates.Feature and a list of reconstructed geometries.

        Returns
        -------
        reconstructed_features : list
            The reconstructed geological features.
            The reconstructed geometries are output in the same order as that of their respective input features (in the
            parameter ``features``). This includes the order across any input feature collections or files. If ``group_with_feature``
            is True then the list contains tuples that group each pygplates.Feature with a list of its reconstructed geometries.

        See Also
        --------
        reconstruct_snapshot
        """
        reconstruct_snapshot = self.reconstruct_snapshot(
            feature,
            to_time,
            anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
            from_time=from_time,
        )

        if group_with_feature:
            # These are always sorted in same order as the input features.
            return reconstruct_snapshot.get_reconstructed_features(reconstruct_type)
        else:
            return reconstruct_snapshot.get_reconstructed_geometries(
                reconstruct_type, same_order_as_reconstructable_features=True
            )

    def get_point_velocities(
        self,
        lons,
        lats,
        time,
        delta_time=1.0,
        *,
        velocity_delta_time_type=pygplates.VelocityDeltaTimeType.t_plus_delta_t_to_t,
        velocity_units=pygplates.VelocityUnits.kms_per_my,
        earth_radius_in_kms=pygplates.Earth.mean_radius_in_kms,
        include_networks=True,
        include_topological_slab_boundaries=False,
        anchor_plate_id=None,
        return_east_north_arrays=False,
    ):
        """Calculates the north and east components of the velocity vector (in kms/myr) for each specified point
        (from ``lons`` and ``lats``) at a particular geological ``time``.

        Parameters
        ----------
        lons : array
            A 1D array of point data's longitudes.

        lats : array
            A 1D array of point data's latitudes.

        time : float
            The specific geological time (Ma) at which to calculate plate velocities.

        delta_time : float, default=1.0
            The time interval used for velocity calculations. 1.0Ma by default.

        velocity_delta_time_type : pygplates.VelocityDeltaTimeType, default=pygplates.VelocityDeltaTimeType.t_plus_delta_t_to_t
            How the two velocity times are calculated relative to ``time`` (defaults to ``[time + velocity_delta_time, time]``).

        velocity_units : pygplates.VelocityUnits, default=pygplates.VelocityUnits.kms_per_my
            Whether to return velocities in centimetres per year or kilometres per million years (defaults to kilometres per million years).

        earth_radius_in_kms : float, default=pygplates.Earth.mean_radius_in_kms
            Radius of the Earth in kilometres.
            This is only used to calculate velocities (strain rates always use ``pygplates.Earth.equatorial_radius_in_kms``).

        include_networks : bool, default=True
            Whether to include deforming networks when calculating velocities.
            By default they are included (and also given precedence since they typically overlay a rigid plate).

        include_topological_slab_boundaries : bool, default=False
            Whether to include features of type ``gpml:TopologicalSlabBoundary`` when calculating velocities.
            By default they are not included (they tend to overlay a rigid plate which should instead be used to calculate plate velocity).

        anchor_plate_id : int, optional
            Anchor plate ID. Defaults to the current anchor plate ID (:py:attr:`anchor_plate_id` attribute).

        return_east_north_arrays : bool, default=False
            Return the velocities as arrays separately containing the east and north components of the velocities.
            Note that setting this to True matches the output of :py:meth:`Points.plate_velocity`.

        Returns
        -------
        north_east_velocities : 2D ndarray
            Only provided if ``return_east_north_arrays`` is False.
            Each array element contains the (north, east) velocity components of a single point.
        east_velocities, north_velocities : 1D ndarray
            Only provided if ``return_east_north_arrays`` is True.
            The east and north components of velocities as separate arrays.
            These are also ordered (east, north) instead of (north, east).

        Raises
        ------
        ValueError
            If topology features have not been set in this :py:class:`PlateReconstruction` object.


        .. note::

            The velocities are in ``kilometres per million years`` by default (not ``centimetres per year``,
            the default in :py:meth:`Points.plate_velocity`).
            This difference is maintained for backward compatibility.

            For each velocity, the ``north`` component is first followed by the ``east`` component.
            This is different to :py:meth:`Points.plate_velocity` where the ``east`` component is first.
            This difference is maintained for backward compatibility.
        """
        # Add points to a multipoint geometry

        points = [pygplates.PointOnSphere(lat, lon) for lat, lon in zip(lats, lons)]

        topological_snapshot = self.topological_snapshot(
            time,
            anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
            include_topological_slab_boundaries=include_topological_slab_boundaries,
        )

        # If requested, exclude resolved topological *networks*.
        resolve_topology_types = pygplates.ResolveTopologyType.boundary
        if include_networks:
            resolve_topology_types = (
                resolve_topology_types | pygplates.ResolveTopologyType.network
            )

        point_velocities = topological_snapshot.get_point_velocities(
            points,
            resolve_topology_types=resolve_topology_types,
            velocity_delta_time=delta_time,
            velocity_delta_time_type=velocity_delta_time_type,
            velocity_units=velocity_units,
            earth_radius_in_kms=earth_radius_in_kms,
        )

        # Replace any missing velocities with zero velocity.
        #
        # If a point does not intersect a topological plate (or network) then its velocity is None.
        for point_index in range(len(points)):
            if point_velocities[point_index] is None:
                point_velocities[point_index] = pygplates.Vector3D.zero

        # Convert global 3D velocity vectors to local (North, East, Down) vectors (one per point).
        point_velocities_north_east_down = (
            pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                points, point_velocities
            )
        )

        if return_east_north_arrays:
            # Extract the East and North velocity components into separate arrays.
            east_velocities = [ned.get_y() for ned in point_velocities_north_east_down]
            north_velocities = [ned.get_x() for ned in point_velocities_north_east_down]
            # Note: This is the opposite order (ie, (east,north) instead of (north,east)).
            return np.array(east_velocities), np.array(north_velocities)
        else:
            # Extract the North and East velocity components into a single array.
            north_east_velocities = [
                (ned.get_x(), ned.get_y()) for ned in point_velocities_north_east_down
            ]
            return np.array(north_east_velocities)

    def create_motion_path(
        self,
        lons,
        lats,
        time_array,
        plate_id=None,
        anchor_plate_id=None,
        return_rate_of_motion=False,
    ):
        """Create a path of points to mark the trajectory of a plate's motion through geological time.

        Parameters
        ----------
        lons : arr
            An array containing the longitudes of seed points on a plate in motion.
        lats : arr
            An array containing the latitudes of seed points on a plate in motion.
        time_array : arr
            An array of reconstruction times at which to determine the trajectory
            of a point on a plate. For example:

            .. code-block:: python
                :linenos:

                import numpy as np
                min_time = 30
                max_time = 100
                time_step = 2.5
                time_array = np.arange(min_time, max_time + time_step, time_step)

        plate_id : int, optional
            The ID of the moving plate. If this is not passed, the plate ID of the
            seed points are ascertained using pygplates' ``PlatePartitioner``.
        anchor_plate_id : int, optional
            The ID of the anchor plate. Defaults to the default anchor plate
            (specified in ``__init__`` or set with ``anchor_plate_id`` attribute).
        return_rate_of_motion : bool, default=False
            Choose whether to return the rate of plate motion through time for each step.

        Returns
        -------
        rlons : ndarray
            An n-dimensional array with columns containing the longitudes of
            the seed points at each timestep in ``time_array``. There are n columns for n seed points.
        rlats : ndarray
            An n-dimensional array with columns containing the latitudes of
            the seed points at each timestep in ``time_array``. There are n columns for n seed points.
        StepTimes : ndarray
            Time interval for each step.
        StepRates : ndarray
            The rate of plate motion for each step.

        Raises
        ------
        ValueError
            If ``plate_id`` is ``None`` and topology features have not been set in this :class:`PlateReconstruction`.

        Examples
        --------
        To access the latitudes and longitudes of each seed point's motion path:

        .. code-block:: python
            :linenos:

            for i in np.arange(0,len(seed_points)):
                current_lons = lon[:,i]
                current_lats = lat[:,i]
        """
        lons = np.atleast_1d(lons)
        lats = np.atleast_1d(lats)
        time_array = np.atleast_1d(time_array)

        # ndarrays to fill with reconstructed points and
        # rates of motion (if requested)
        rlons = np.empty((len(time_array), len(lons)))
        rlats = np.empty((len(time_array), len(lons)))

        if plate_id is None:
            query_plate_id = True
            plate_ids = []
        else:
            query_plate_id = False
            plate_ids = np.ones(len(lons), dtype=int) * plate_id

        seed_points = zip(lats, lons)
        if return_rate_of_motion is True:
            StepTimes = np.empty(((len(time_array) - 1) * 2, len(lons)))
            StepRates = np.empty(((len(time_array) - 1) * 2, len(lons)))
        else:
            StepTimes = np.array([])
            StepRates = np.array([])
        for i, lat_lon in enumerate(seed_points):
            seed_points_at_digitisation_time = pygplates.PointOnSphere(
                pygplates.LatLonPoint(float(lat_lon[0]), float(lat_lon[1]))
            )
            # Allocate the present-day plate ID to the PointOnSphere if
            # it was not given.
            if query_plate_id:
                plate_id = _tools.plate_partitioner_for_point(
                    lat_lon, self._check_topology_features(), self.rotation_model
                )
            else:
                plate_id = plate_ids[i]

            # Create the motion path feature. enforce float and int for C++ signature.
            motion_path_feature = pygplates.Feature.create_motion_path(
                seed_points_at_digitisation_time,
                time_array,
                valid_time=(time_array.max(), time_array.min()),
                relative_plate=(  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
                    anchor_plate_id
                    if anchor_plate_id is not None
                    else self.anchor_plate_id
                ),
                reconstruction_plate_id=int(plate_id),
            )

            reconstructed_motion_paths = self.reconstruct(
                motion_path_feature,
                to_time=0,
                reconstruct_type=pygplates.ReconstructType.motion_path,
                anchor_plate_id=anchor_plate_id,  # if None then uses 'self.anchor_plate_id' (default anchor plate of 'self.rotation_model')
            )
            # Turn motion paths in to lat-lon coordinates
            trail = None
            for reconstructed_motion_path in reconstructed_motion_paths:
                # not sure about this. always set the "trail" to the last one in reconstructed_motion_paths?
                # or there is only one path in reconstructed_motion_paths? -- Michael Chin
                trail = reconstructed_motion_path.get_motion_path().to_lat_lon_array()
            assert trail is not None
            lon, lat = np.flipud(trail[:, 1]), np.flipud(trail[:, 0])

            rlons[:, i] = lon
            rlats[:, i] = lat

            # Obtain step-plot coordinates for rate of motion
            if return_rate_of_motion is True:
                # Get timestep
                TimeStep = []
                for j in range(len(time_array) - 1):
                    diff = time_array[j + 1] - time_array[j]
                    TimeStep.append(diff)

                # Iterate over each segment in the reconstructed motion path, get the distance travelled by the moving
                # plate relative to the fixed plate in each time step
                Dist = []
                for reconstructed_motion_path in reconstructed_motion_paths:
                    for (
                        segment
                    ) in reconstructed_motion_path.get_motion_path().get_segments():
                        Dist.append(
                            segment.get_arc_length()
                            * _tools.geocentric_radius(
                                segment.get_start_point().to_lat_lon()[0]
                            )
                            / 1e3
                        )

                # Note that the motion path coordinates come out starting with the oldest time and working forwards
                # So, to match our 'times' array, we flip the order
                Dist = np.flipud(Dist)

                # Get rate of motion as distance per Myr
                Rate = np.asarray(Dist) / TimeStep

                # Manipulate arrays to get a step plot
                StepRate = np.zeros(len(Rate) * 2)
                StepRate[::2] = Rate
                StepRate[1::2] = Rate

                StepTime = np.zeros(len(Rate) * 2)
                StepTime[::2] = time_array[:-1]
                StepTime[1::2] = time_array[1:]

                # Append the nth point's step time and step rate coordinates to the ndarray
                StepTimes[:, i] = StepTime
                StepRates[:, i] = StepRate * 0.1  # cm/yr

                # Obseleted by Lauren's changes above (though it is more efficient)
                # multiply arc length of the motion path segment by a latitude-dependent Earth radius
                # use latitude of the segment start point
                # distance.append( segment.get_arc_length() * _tools.geocentric_radius(segment.get_start_point().to_lat_lon()[0]) / 1e3)
                # rate = np.asarray(distance)/np.diff(time_array)
                # rates[:,i] = np.flipud(rate)
                # rates *= 0.1 # cm/yr

        if return_rate_of_motion is True:
            return (
                np.squeeze(rlons),
                np.squeeze(rlats),
                np.squeeze(StepTimes),
                np.squeeze(StepRates),
            )
        else:
            return np.squeeze(rlons), np.squeeze(rlats)

    def create_flowline(
        self,
        lons,
        lats,
        time_array,
        left_plate_ID,
        right_plate_ID,
        return_rate_of_motion=False,
    ):
        """Create a path of points to track plate motion away from spreading ridges over time using half-stage rotations.

        Parameters
        ----------
        lons : arr
            An array of longitudes of points along spreading ridges.
        lats : arr
            An array of latitudes of points along spreading ridges.
        time_array : arr
            A list of times to obtain seed point locations at.
        left_plate_ID : int
            The plate ID of the polygon to the left of the spreading
            ridge.
        right_plate_ID : int
            The plate ID of the polygon to the right of the spreading
            ridge.
        return_rate_of_motion : bool, default False
            Choose whether to return a step time and step rate array
            for a step plot of motion.

        Returns
        -------
        left_lon : ndarray
            The longitudes of the left flowline for n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in time_array.
        left_lat : ndarray
            The latitudes of the left flowline of n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in time_array.
        right_lon : ndarray
            The longitudes of the right flowline of n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in time_array.
        right_lat : ndarray
            The latitudes of the right flowline of n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in time_array.

        Examples
        --------
        To access the i\\ :sup:`th` seed point's left and right latitudes and
        longitudes:

        .. code-block:: python
            :linenos:

            for i in np.arange(0,len(seed_points)):
                left_flowline_longitudes = left_lon[:,i]
                left_flowline_latitudes = left_lat[:,i]
                right_flowline_longitudes = right_lon[:,i]
                right_flowline_latitudes = right_lat[:,i]
        """
        lats = np.atleast_1d(lats)
        lons = np.atleast_1d(lons)
        time_array = np.atleast_1d(time_array)

        seed_points = list(zip(lats, lons))
        multi_point = pygplates.MultiPointOnSphere(seed_points)

        start = 0
        if time_array[0] != 0:
            start = 1
            time_array = np.hstack([[0], time_array])

        # Create the flowline feature
        flowline_feature = pygplates.Feature.create_flowline(
            multi_point,
            time_array.tolist(),
            valid_time=(time_array.max(), time_array.min()),
            left_plate=left_plate_ID,
            right_plate=right_plate_ID,
        )

        # reconstruct the flowline in present-day coordinates
        reconstructed_flowlines = self.reconstruct(
            flowline_feature,
            to_time=0,
            reconstruct_type=pygplates.ReconstructType.flowline,
        )

        # Wrap things to the dateline, to avoid plotting artefacts.
        date_line_wrapper = pygplates.DateLineWrapper()

        # Create lat-lon ndarrays to store the left and right lats and lons of flowlines
        left_lon = np.empty((len(time_array), len(lons)))
        left_lat = np.empty((len(time_array), len(lons)))
        right_lon = np.empty((len(time_array), len(lons)))
        right_lat = np.empty((len(time_array), len(lons)))
        StepTimes = np.empty(((len(time_array) - 1) * 2, len(lons)))
        StepRates = np.empty(((len(time_array) - 1) * 2, len(lons)))

        # Iterate over the reconstructed flowlines. Each seed point results in a 'left' and 'right' flowline
        for i, reconstructed_flowline in enumerate(reconstructed_flowlines):
            # Get the points for the left flowline only
            left_latlon = reconstructed_flowline.get_left_flowline().to_lat_lon_array()
            left_lon[:, i] = left_latlon[:, 1]
            left_lat[:, i] = left_latlon[:, 0]

            # Repeat for the right flowline points
            right_latlon = (
                reconstructed_flowline.get_right_flowline().to_lat_lon_array()
            )
            right_lon[:, i] = right_latlon[:, 1]
            right_lat[:, i] = right_latlon[:, 0]

        if return_rate_of_motion:
            for i, reconstructed_motion_path in enumerate(reconstructed_flowlines):
                distance = []
                for (
                    segment
                ) in reconstructed_motion_path.get_left_flowline().get_segments():
                    distance.append(
                        segment.get_arc_length()
                        * _tools.geocentric_radius(
                            segment.get_start_point().to_lat_lon()[0]
                        )
                        / 1e3
                    )

                # Get rate of motion as distance per Myr
                # Need to multiply rate by 2, since flowlines give us half-spreading rate
                time_step = time_array[1] - time_array[0]
                Rate = (
                    np.asarray(distance) / time_step
                ) * 2  # since we created the flowline at X increment

                # Manipulate arrays to get a step plot
                StepRate = np.zeros(len(Rate) * 2)
                StepRate[::2] = Rate
                StepRate[1::2] = Rate

                StepTime = np.zeros(len(Rate) * 2)
                StepTime[::2] = time_array[:-1]
                StepTime[1::2] = time_array[1:]

                # Append the nth point's step time and step rate coordinates to the ndarray
                StepTimes[:, i] = StepTime
                StepRates[:, i] = StepRate * 0.1  # cm/yr

            return (
                left_lon[start:],
                left_lat[start:],
                right_lon[start:],
                right_lat[start:],
                StepTimes,
                StepRates,
            )

        else:
            return (
                left_lon[start:],
                left_lat[start:],
                right_lon[start:],
                right_lat[start:],
            )
