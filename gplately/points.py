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
from typing import Union

import numpy as np
import pygplates

from . import tools as _tools

logger = logging.getLogger("gplately")


class Points(object):
    """Reconstruct and work with geological point data.

    The locations and plate velocities of point data can be calculated at a specific geological time.
    The :py:class:`Points` class depends on the :py:class:`PlateReconstruction` class,
    as it provides essential components for reconstructing plate motions.
    Specifically, it uses the :py:attr:`PlateReconstruction.rotation_model` to compute point rotations through time,
    and the :py:attr:`PlateReconstruction.static_polygons` to assign points to their respective tectonic plates.
    """

    def __init__(
        self,
        plate_reconstruction,
        lons,
        lats,
        time: float = 0,
        plate_id=None,
        age: Union[float, np.ndarray, None] = np.inf,
        *,
        anchor_plate_id=None,
        remove_unreconstructable_points=False,
    ):
        """
        Parameters
        ----------
        plate_reconstruction : PlateReconstruction
            Object to provide the following essential components for reconstructing points.

            * :py:attr:`PlateReconstruction.rotation_model`
            * :py:attr:`PlateReconstruction.topology_featues`
            * :py:attr:`PlateReconstruction.static_polygons`

        lons : float or 1D array
            Longitudes of the initial points at the initial ``time``.

        lats : float or 1D array
            Latitudes of the initial points at the initial ``time``.

        time : float, default=0
            The initial time (Ma) of the points.
            The ``lons`` and ``lats`` are the initial coordinates of the points at this time.
            By default, it is set to the present day (0 Ma).

        plate_id : int or 1D array or None, default=None
            Plate ID(s) of a particular tectonic plate on which point data lies, if known.
            If it is a single integer then all points will have the same plate ID. If it is a 1D array then length must match the number of points.
            If ``None`` then plate IDs are determined using the :py:attr:`PlateReconstruction.static_polygons`.
            By default, the plate IDs are determined using the static polygons.

        age : float or 1D array or None, default=numpy.inf
            Age(s) at which each point appears, if known.
            If it is a single float then all points will have the same age.
            If it is a 1D array then length must match the number of points.
            If ``None`` then ages are determined using the :py:attr:`PlateReconstruction.static_polygons`.
            For points on oceanic crust this is when they were created at a mid-ocean ridge.
            By default, all points exist for all time (ie, time of appearance is infinity). This default is for backward
            compatibility, but you'll typically only want this if all your points are on **continental** crust (not *oceanic*).

        anchor_plate_id : int, optional
            Anchor plate ID that the specified ``lons`` and ``lats`` are relative to.
            Defaults to the current anchor plate ID of ``plate_reconstruction`` (its ``anchor_plate_id`` attribute).

        remove_unreconstructable_points : bool or list, default=False
            Whether to remove points that cannot be reconstructed.
            By default, all unreconstructable points are retained.
            A point cannot be reconstructed if it cannot be assigned a plate ID, or cannot be assigned an age, because it did not
            intersect any reconstructed static polygons (note that this can only happen when ``plate_id`` and/or ``age`` is None).
            Also, a point cannot be reconstructed if point ages were **explicitly** provided (ie, ``age`` was **not** None) and
            a point's age was less than (younger than) ``time``, meaning it did not exist as far back as ``time``.
            Additionally, if this variable is a :py:class:`list` then the indices (into the supplied ``lons`` and ``lats`` arguments)
            of any removed points (ie, that are unreconstructable) are appended to that list.


        .. _points-note:
        .. note::

            If ``time`` is non-zero (ie, not present day) then ``lons`` and ``lats`` are assumed to be the **reconstructed** point
            locations at ``time``. And the reconstructed positions are assumed to be relative to the anchor plate
            (which is ``plate_reconstruction.anchor_plate_id`` if ``anchor_plate_id`` is None).

            If ``plate_id`` and/or ``age`` is None then the plate ID and/or age of each point is determined by reconstructing the static polygons
            of ``plate_reconstruction`` to ``time`` and reconstructing relative to the anchor plate (regardless of whether ``time`` is present day or not).
            And then, for each point, assigning the plate ID and/or time-of-appearance (begin time) of the static polygon containing the point.

            A point is considered unreconstructable if it does not exist at ``time``. This can happen if its age was explicitly provided (ie, ``age`` is **not** None)
            but is younger than ``time``. It can also happen if the point is automatically assigned a plate ID (ie, ``plate_id`` is None) or an age (ie, ``age`` is None)
            but does not intersect any reconstructed static polygons (at ``time``). In either of these cases it is marked as unreconstructable and will not be available
            for any method outputing a reconstruction, such as :meth:`reconstruct()`, or any method depending on a reconstruction, such as :meth:`plate_velocity`.
            However, all the initial locations and their associated plate IDs and ages will still be accessible as attributes, regardless of whether all the points
            are reconstructable or not. That is, unless ``remove_unreconstructable_points`` is True (or a :py:class:`list`),
            in which case only the reconstructable points are retained.
        """
        # If anchor plate is None then use default anchor plate of 'plate_reconstruction'.
        if anchor_plate_id is None:
            anchor_plate_id = plate_reconstruction.anchor_plate_id
        else:
            anchor_plate_id = self._check_anchor_plate_id(anchor_plate_id)

        point_ages = np.array([])
        point_plate_ids = np.array([])

        # The caller can specify a 'list' for the 'remove_unreconstructable_points' argument if they want us to
        # return the indices of any points that are NOT reconstructable.
        #
        # Otherwise 'remove_unreconstructable_points' must be true or false.
        if isinstance(remove_unreconstructable_points, list):
            unreconstructable_point_indices_list = remove_unreconstructable_points
            remove_unreconstructable_points = True
        else:
            unreconstructable_point_indices_list = None

        # Most common case first: both are sequences.
        if not np.isscalar(lons) and not np.isscalar(lats):
            # Make sure numpy arrays (if not already).
            lons = np.asarray(lons)
            lats = np.asarray(lats)
            if len(lons) != len(lats):
                raise ValueError(
                    "'lons' and 'lats' must be of equal length ({} != {})".format(
                        len(lons), len(lats)
                    )
                )
        elif np.isscalar(lons) and np.isscalar(lats):
            # Both are scalars. Convert to arrays with one element.
            lons = np.atleast_1d(lons)
            lats = np.atleast_1d(lats)
        else:
            raise ValueError(
                "Both 'lats' and 'lons' must both be a sequence or both a scalar"
            )

        num_points = len(lons)

        # If caller provided plate IDs.
        if plate_id is not None:
            # If plate ID is a scalar then all points have the same plate ID.
            if np.isscalar(plate_id):
                point_plate_ids = np.full(num_points, plate_id)
            else:
                point_plate_ids = np.asarray(plate_id)
                if len(point_plate_ids) != num_points:
                    raise ValueError(
                        "'plate_id' must be same length as 'lons' and 'lats' ({} != {})".format(
                            len(point_plate_ids), num_points
                        )
                    )

        # If caller provided begin ages.
        if age is not None:
            # If age is a scalar then all points have the same age.
            if np.isscalar(age):
                point_ages = np.full(num_points, age)
            else:
                point_ages = np.asarray(age)
                if len(point_ages) != num_points:
                    raise ValueError(
                        "'age' must be same length as 'lons' and 'lats' ({} != {})".format(
                            len(point_ages), num_points
                        )
                    )

        # Create pygplates points.
        points = [pygplates.PointOnSphere(lat, lon) for lon, lat in zip(lons, lats)]

        # If plate IDs and/or ages are automatically assigned using reconstructed static polygons then
        # some points might be outside all reconstructed static polygons, and hence not reconstructable.
        #
        # However, if the user provided both plate IDs and ages then all points will be reconstructable.
        points_are_reconstructable = np.full(num_points, True)

        # If caller did not provide plate IDs or begin ages then
        # we need to determine them using the static polygons.
        if plate_id is None or age is None:
            if plate_id is None:
                point_plate_ids = np.empty(num_points, dtype=int)
            if age is None:
                point_ages = np.empty(num_points)

            # Assign a plate ID to each point based on which reconstructed static polygon it's inside.
            static_polygons_snapshot = plate_reconstruction.static_polygons_snapshot(
                time,
                anchor_plate_id=anchor_plate_id,
            )
            reconstructed_static_polygons_containing_points = (
                static_polygons_snapshot.get_point_locations(points)
            )
            for point_index in range(num_points):
                reconstructed_static_polygon = (
                    reconstructed_static_polygons_containing_points[point_index]
                )

                # If current point is inside a reconstructed static polygon then assign its plate ID to the point,
                # otherwise assign the anchor plate to the point.
                if reconstructed_static_polygon is not None:
                    reconstructed_static_polygon_feature = (
                        reconstructed_static_polygon.get_feature()
                    )

                    if plate_id is None:
                        point_plate_ids[point_index] = (
                            reconstructed_static_polygon_feature.get_reconstruction_plate_id()
                        )
                    if age is None:
                        point_ages[point_index], _ = (
                            reconstructed_static_polygon_feature.get_valid_time()
                        )

                else:  # current point did NOT intersect a reconstructed static polygon ...

                    # We're trying to assign a plate ID or assign an age (or both), neither of which we can assign.
                    # That essentially makes the current point unreconstructable.
                    #
                    # Mark the current point as unreconstructable.
                    points_are_reconstructable[point_index] = False

                    if plate_id is None:
                        # Assign the anchor plate ID to indicate we could NOT assign a proper plate ID.
                        point_plate_ids[point_index] = anchor_plate_id
                    if age is None:
                        # Assign the distant future (not distant past) to indicate we could NOT assign a proper age.
                        point_ages[point_index] = -np.inf  # distant future

        # If point ages were explicitly provided by the caller then we need to check if points existed at 'time'.
        if age is not None:
            # Any point with an age younger than 'time' did not exist at 'time' and hence is not reconstructable.
            points_are_reconstructable[point_ages < time] = False

        # If requested, remove any unreconstructable points.
        if remove_unreconstructable_points and not points_are_reconstructable.all():
            if unreconstructable_point_indices_list is not None:
                # Caller requested the indices of points that are NOT reconstructable.
                unreconstructable_point_indices_list.extend(
                    np.where(points_are_reconstructable == False)[0]
                )
            lons = lons[points_are_reconstructable]
            lats = lats[points_are_reconstructable]
            point_plate_ids = point_plate_ids[points_are_reconstructable]
            point_ages = point_ages[points_are_reconstructable]
            points = [
                points[point_index]
                for point_index in range(num_points)
                if points_are_reconstructable[point_index]
            ]
            num_points = len(points)
            # All points are now reconstructable.
            points_are_reconstructable = np.full(num_points, True)

        # Create a feature for each point.
        #
        # Each feature has a point, a plate ID and a valid time range.
        #
        # Note: The valid time range always includes present day.
        point_features = []
        for point_index in range(num_points):
            point_feature = pygplates.Feature()
            # Set the geometry.
            point_feature.set_geometry(points[point_index])
            # Set the plate ID.
            point_feature.set_reconstruction_plate_id(point_plate_ids[point_index])  # type: ignore
            # Set the begin/end time.
            point_feature.set_valid_time(
                point_ages[point_index],  # begin (age)
                -np.inf,  # end (distant future; could also be zero for present day)
            )  # type: ignore
            point_features.append(point_feature)

        # If the points represent a snapshot at a *past* geological time then we need to reverse reconstruct them
        # such that their features contain present-day points.
        if time != 0:
            pygplates.reverse_reconstruct(  # type: ignore
                point_features,
                plate_reconstruction.rotation_model,
                time,
                anchor_plate_id=anchor_plate_id,
            )

        # Map each unique plate ID to indices of points assigned that plate ID.
        unique_plate_id_groups = {}
        unique_plate_ids = np.unique(point_plate_ids)
        for unique_plate_id in unique_plate_ids:
            # Determine which points have the current unique plate ID.
            unique_plate_id_point_indices = np.where(
                point_plate_ids == unique_plate_id
            )[
                0
            ]  # convert 1-tuple of 1D array to 1D array
            unique_plate_id_groups[unique_plate_id] = unique_plate_id_point_indices

        #
        # Assign data members.
        #

        # Note: These are documented attributes (in class docstring).
        #       And they cannot be changed later (they are properties with no setter).
        #       The other attributes probably should be readonly too (but at least they're not documented).
        self._plate_reconstruction = plate_reconstruction
        self._lons = lons
        self._lats = lats
        self._time = time
        self._plate_id = point_plate_ids
        self._age = point_ages
        self._anchor_plate_id = anchor_plate_id

        # get Cartesian coordinates
        self.x, self.y, self.z = _tools.lonlat2xyz(lons, lats, degrees=False)
        # scale by average radius of the Earth
        self.x *= _tools.EARTH_RADIUS
        self.y *= _tools.EARTH_RADIUS
        self.z *= _tools.EARTH_RADIUS
        # store concatenated arrays
        self.lonlat = np.c_[lons, lats]
        self.xyz = np.c_[self.x, self.y, self.z]

        self.points = points

        self.attributes = dict()

        self._reconstructable = points_are_reconstructable
        self._unique_plate_id_groups = unique_plate_id_groups

        self.features = point_features
        self.feature_collection = pygplates.FeatureCollection(point_features)

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
    def plate_reconstruction(self):
        """
        Object to provide the following essential components for reconstructing points.

        * :py:attr:`PlateReconstruction.rotation_model`
        * :py:attr:`PlateReconstruction.topology_featues`
        * :py:attr:`PlateReconstruction.static_polygons`

        :type: PlateReconstruction

        """
        return self._plate_reconstruction

    @property
    def lons(self):
        """
        Longitudes of the initial points at the initial ``time``.

        :type: float 1D array
        """
        return self._lons

    @property
    def lats(self):
        """
        Latitudes of the initial points at the initial ``time``.

        :type: float 1D array
        """
        return self._lats

    @property
    def plate_id(self):
        """
        A 1D array containing the plate IDs of the points.
        The length must match that of ``lons`` and ``lats``.

        :type: int 1D array
        """
        return self._plate_id

    @property
    def age(self):
        """
        A 1D array containing the ages (time of appearance) of the points.
        The length must match that of ``lons`` and ``lats``.
        For points on oceanic crust this is when they were created at a mid-ocean ridge.
        Any points existing for all time will have a value of ``numpy.inf`` (equivalent to ``float('inf')``).

        :type: float 1D array
        """
        return self._age

    @property
    def size(self):
        """
        Number of points.
        This is the size of ``lons``, ``lats``, ``plate_id`` and ``age``.

        :type: int
        """
        return len(self.points)

    @property
    def time(self):
        """
        The initial time (Ma) of the points.
        The initial ``lons`` and ``lats`` are the coordinates of the points at this time.

        :type: float
        """
        return self._time

    @property
    def anchor_plate_id(self):
        """
        Anchor plate that the initial ``lons`` and ``lats`` are relative to, at the initial ``time``.
        This is also used as the default anchor plate when reconstructing the points.
        It does not change, even if the anchor plate of ``plate_reconstruction`` subsequently changes.

        :type: int
        """
        return self._anchor_plate_id

    @staticmethod
    def _check_anchor_plate_id(id):
        id = int(id)
        if id < 0:
            raise ValueError("Invalid anchor plate ID: {}".format(id))
        return id

    def copy(self):
        """Returns a copy of the :py:class:`Points` object

        Returns
        -------
        Points
            A copy of the current :py:class:`Points` object
        """
        gpts = Points(
            self.plate_reconstruction,
            self.lons.copy(),
            self.lats.copy(),
            self.time,
            self.plate_id.copy(),
            self.age.copy(),
            anchor_plate_id=self.anchor_plate_id,
        )
        gpts.add_attributes(**self.attributes.copy())

    def add_attributes(self, **kwargs):
        """Adds the value of a feature attribute associated with a key.

        Parameters
        ----------
        **kwargs : sequence of key=item/s
            A single key=value pair, or a sequence of key=value pairs denoting the name and
            value of an attribute.

        Example
        -------
        .. code-block:: python
            :linenos:

            # Define latitudes and longitudes to set up a Points object
            pt_lons = np.array([140., 150., 160.])
            pt_lats = np.array([-30., -40., -50.])

            gpts = gplately.Points(model, pt_lons, pt_lats)

            # Add the attributes a, b and c to the points in the Points object
            gpts.add_attributes(
                a=[10,2,2],
                b=[2,3,3],
                c=[30,0,0],
            )

            print(gpts.attributes)

        The output would be:

        .. code:: console

            {'a': [10, 2, 2], 'b': [2, 3, 3], 'c': [30, 0, 0]}


        .. note::

            An **assertion** is raised if the number of points in the Points object is not equal
            to the number of values associated with an attribute key. For example, consider an instance
            of the Points object with 3 points. If the points are ascribed an attribute ``temperature``,
            there must be one ``temperature`` value per point, i.e. ``temperature = [20, 15, 17.5]``.

        """
        keys = kwargs.keys()

        for key in kwargs:
            attribute = kwargs[key]

            # make sure attribute is the same size as self.lons
            if type(attribute) is int or type(attribute) is float:
                array = np.full(self.lons.size, attribute)
                attribute = array
            elif isinstance(attribute, np.ndarray):
                if attribute.size == 1:
                    array = np.full(self.lons.size, attribute, dtype=attribute.dtype)
                    attribute = array

            assert (
                len(attribute) == self.lons.size
            ), "Size mismatch, ensure attributes have the same number of entries as Points"
            self.attributes[key] = attribute

        if any(kwargs):
            # add these to the FeatureCollection
            for f, feature in enumerate(self.feature_collection):
                for key in keys:
                    # extract value for each row in attribute
                    val = self.attributes[key][f]

                    # set this attribute on the feature
                    feature.set_shapefile_attribute(key, val)

    def get_geopandas_dataframe(self):
        """Return a ``geopandas.GeoDataFrame`` object for the points.

        Returns
        -------
        GeoDataFrame : geopandas.GeoDataFrame
            A ``geopandas.GeoDataFrame`` object with rows equal to the number of points in the :class:`Points` object,
            and an additional column containing a shapely ``geometry`` attribute.

        Example
        -------

        .. code-block:: python
            :linenos:

            pt_lons = np.array([140., 150., 160.])
            pt_lats = np.array([-30., -40., -50.])

            gpts = gplately.Points(model, pt_lons, pt_lats)

            # Add sample attributes a, b and c to the points in the Points object
            gpts.add_attributes(
                a=[10,2,2],
                b=[2,3,3],
                c=[30,0,0],
            )

            gpts.get_geopandas_dataframe()

        has the output:

        .. code:: console

                a  b   c                     geometry
            0  10  2  30  POINT (140.00000 -30.00000)
            1   2  3   0  POINT (150.00000 -40.00000)
            2   2  3   0  POINT (160.00000 -50.00000)


        """
        import geopandas as gpd
        from shapely import geometry

        # create shapely points
        points = []
        for lon, lat in zip(self.lons, self.lats):
            points.append(geometry.Point(lon, lat))

        attributes = self.attributes.copy()
        attributes["geometry"] = points

        return gpd.GeoDataFrame(attributes, geometry="geometry")

    def get_geodataframe(self):
        """The same as :meth:`get_geopandas_dataframe`."""
        return self.get_geopandas_dataframe()

    def reconstruct(
        self, time, anchor_plate_id=None, return_array=False, return_point_indices=False
    ):
        """Reconstruct points from the initial time (``self.time``) to the specified time (``time``).

        Only those points that are reconstructable (See :ref:`Notes <points-note>`) and that have ages greater than or equal
        to ``time`` (ie, at points that exist at ``time``) are reconstructed.

        Parameters
        ----------
        time : float
            The specific geological time (Ma) to reconstruct features to.

        anchor_plate_id : int, optional
            Reconstruct features with respect to a certain anchor plate.
            By default, reconstructions are made with respect to ``self.anchor_plate_id``
            (which is the anchor plate that the initial points at the initial time are relative to).

        return_array : bool, default=False
            Return a 2-tuple of ``numpy.ndarray``, rather than a :class:`Points` object.

        return_point_indices : bool, default=False
            Return the indices of the points that are reconstructed.
            Those points with an age less than ``time`` have not yet appeared at ``time``, and therefore are not reconstructed.
            These are indices into ``self.lons``, ``self.lats``, ``self.plate_id`` and ``self.age``.

        Returns
        -------
        reconstructed_points : :class:`Points`
            If the ``return_array`` is False, return the reconstructed points in a :class:`Points` object.
        rlons, rlats : ndarray
            If the ``return_array`` is True, return the longitude and latitude coordinate arrays of the reconstructed points.
        point_indices : ndarray
            If the ``return_point_indices`` is True, return the indices of the returned points (that are reconstructed).
            This array is the same size as ``rlons`` and ``rlats`` (or size of ``reconstructed_points``).
            These are indices into ``self.lons``, ``self.lats``, ``self.plate_id`` and ``self.age``.
        """
        if anchor_plate_id is None:
            anchor_plate_id = self.anchor_plate_id

        # Start with an empty array.
        lat_lon_points = np.empty((self.size, 2))

        # Determine which points are valid.
        #
        # These are those points that are reconstructable and have appeared before (or at) 'time'
        # (ie, have a time-of-appearance that's greater than or equal to 'time').
        valid_mask = self._reconstructable & (self.age >= time)

        # Iterate over groups of points with the same plate ID.
        for (
            plate_id,
            point_indices_with_plate_id,
        ) in self._unique_plate_id_groups.items():

            # Determine which points (indices) with the current unique plate ID are valid.
            point_indices_with_plate_id = point_indices_with_plate_id[
                valid_mask[point_indices_with_plate_id]
            ]
            # If none of the points (with the current unique plate ID) are valid then skip to next unique plate ID.
            if point_indices_with_plate_id.size == 0:
                continue

            # Get the reconstructed points with the current unique plate ID that have appeared before (or at) 'time'.
            reconstructed_points_with_plate_id = pygplates.MultiPointOnSphere(
                self.points[point_index] for point_index in point_indices_with_plate_id
            )

            # First reconstruct the internal points from the initial time ('self.time') to present day using
            # our internal anchor plate ID (the same anchor plate used in '__init__').
            # Then reconstruct from present day to 'time' using the *requested* anchor plate ID.
            #
            # Note 'self.points' (and hence 'reconstructed_points_with_plate_id') are the locations at 'self.time'
            #      (just like 'self.lons' and 'self.lats').
            reconstruct_rotation = (
                self.plate_reconstruction.rotation_model.get_rotation(
                    to_time=time,
                    moving_plate_id=plate_id,
                    from_time=0,
                    anchor_plate_id=anchor_plate_id,
                )
                * self.plate_reconstruction.rotation_model.get_rotation(
                    to_time=0,
                    moving_plate_id=plate_id,
                    from_time=self.time,
                    anchor_plate_id=self.anchor_plate_id,
                )
            )
            reconstructed_points_with_plate_id = (
                reconstruct_rotation * reconstructed_points_with_plate_id
            )

            # Write the reconstructed points.
            lat_lon_points[point_indices_with_plate_id] = [
                rpoint.to_lat_lon() for rpoint in reconstructed_points_with_plate_id
            ]

        rlonslats = lat_lon_points[valid_mask]  # remove invalid points
        rlons = rlonslats[:, 1]
        rlats = rlonslats[:, 0]

        return_tuple = ()

        if return_array:
            return_tuple += rlons, rlats
        else:
            reconstructed_points = Points(
                self.plate_reconstruction,
                rlons,
                rlats,
                time=time,
                plate_id=self.plate_id[valid_mask],  # remove invalid points
                age=self.age[valid_mask],  # remove invalid points
                anchor_plate_id=anchor_plate_id,
            )
            reconstructed_points.add_attributes(**self.attributes.copy())
            return_tuple += (reconstructed_points,)

        if return_point_indices:
            all_point_indices = np.arange(self.size, dtype=int)
            point_indices = all_point_indices[valid_mask]  # remove invalid points
            return_tuple += (point_indices,)

        # Return tuple of objects (unless only a single object, eg, just a 'Points' object).
        if len(return_tuple) == 1:
            return return_tuple[0]
        else:
            return return_tuple

    def reconstruct_to_birth_age(
        self, ages, anchor_plate_id=None, return_point_indices=False
    ):
        """Reconstruct points from the initial time (``self.time``) to a range of times.

        The number of supplied times must equal the number of points supplied to this :class:`Points` object (ie, ``self.size`` attribute).
        Only those points that are reconstructable (See :ref:`Notes <points-note>`) and that have ages greater than or equal
        to the respective supplied ages (ie, at points that exist at the supplied ages) are reconstructed.

        Parameters
        ----------
        ages : array
            Geological times to reconstruct points to. Must have the same length as the number of points (``self.size`` attribute).

        anchor_plate_id : int, optional
            Reconstruct points with respect to a certain anchor plate.
            By default, reconstructions are made with respect to ``self.anchor_plate_id``
            (which is the anchor plate that the initial points at the initial time are relative to).

        return_point_indices : bool, default=False
            Return the indices of the points that are reconstructed.
            Those points with an age less than their respective supplied age have not yet appeared, and therefore are not reconstructed.
            These are indices into ``self.lons``, ``self.lats``, ``self.plate_id`` and ``self.age``.

        Raises
        ------
        ValueError
            If the number of ages is not equal to the number of points supplied to this :class:`Points` object.

        Returns
        -------
        rlons, rlats : ndarray
            The longitude and latitude coordinate arrays of points reconstructed to the specified ages.
        point_indices : ndarray
            Only provided if ``return_point_indices`` is True.
            The indices of the returned points (that are reconstructed).
            This array is the same size as ``rlons`` and ``rlats``.
            These are indices into ``self.lons``, ``self.lats``, ``self.plate_id`` and ``self.age``.

        Examples
        --------
        To reconstruct n seed points to B Ma (for this example n=2, with (lon,lat) = (78,30) and (56,22) at time=0 Ma,
        and we reconstruct to B=10 Ma):

        .. code-block:: python
            :linenos:

            # Longitude and latitude of n=2 seed points
            pt_lon = np.array([78., 56])
            pt_lat = np.array([30., 22])

            # Call the Points object!
            gpts = gplately.Points(model, pt_lon, pt_lat)
            print(gpts.features[0].get_all_geometries())   # Confirms we have features represented as points on a sphere

            ages = numpy.linspace(10,10, len(pt_lon))
            rlons, rlats = gpts.reconstruct_to_birth_age(ages)

        """
        if anchor_plate_id is None:
            anchor_plate_id = self.anchor_plate_id

        # Call it 'reconstruct_ages' to avoid confusion with 'self.age' (which is time-of-appearance of points).
        reconstruct_ages = np.asarray(ages)

        if len(reconstruct_ages) != self.size:
            raise ValueError(
                "'ages' must be same length as number of points ({} != {})".format(
                    len(reconstruct_ages), self.size
                )
            )

        # Start with an empty array.
        lat_lon_points = np.empty((self.size, 2))

        # Determine which points are valid.
        #
        # These are those points that are reconstructable and have appeared before (or at) their respective reconstruct ages
        # (ie, have a time-of-appearance that's greater than or equal to the respective reconstruct age).
        valid_mask = self._reconstructable & (self.age >= reconstruct_ages)

        # Iterate over groups of points with the same plate ID.
        for (
            plate_id,
            point_indices_with_plate_id,
        ) in self._unique_plate_id_groups.items():

            # Determine which points (indices) with the current unique plate ID are valid.
            point_indices_with_plate_id = point_indices_with_plate_id[
                valid_mask[point_indices_with_plate_id]
            ]
            # If none of the points (with the current unique plate ID) are valid then skip to next unique plate ID.
            if point_indices_with_plate_id.size == 0:
                continue

            # Get all the unique reconstruct ages of all valid points with the current unique plate ID.
            point_reconstruct_ages_with_plate_id = reconstruct_ages[
                point_indices_with_plate_id
            ]
            unique_reconstruct_ages_with_plate_id = np.unique(
                point_reconstruct_ages_with_plate_id
            )
            for reconstruct_age in unique_reconstruct_ages_with_plate_id:
                # Indices of points with the current unique plate ID and the current unique reconstruct age.
                point_indices_with_plate_id_and_reconstruct_age = (
                    point_indices_with_plate_id[
                        point_reconstruct_ages_with_plate_id == reconstruct_age
                    ]
                )

                # Get the reconstructed points with the current unique plate ID and unique reconstruct age
                # (that exist at their respective reconstruct age).
                reconstructed_points_with_plate_id_and_reconstruct_age = pygplates.MultiPointOnSphere(
                    self.points[point_index]
                    for point_index in point_indices_with_plate_id_and_reconstruct_age
                )

                # First reconstruct the internal points from the initial time ('self.time') to present day using
                # our internal anchor plate ID (the same anchor plate used in '__init__').
                # Then reconstruct from present day to 'reconstruct_age' using the *requested* anchor plate ID.
                #
                # Note 'self.points' (and hence 'reconstructed_points_with_plate_id_and_reconstruct_age') are the locations at 'self.time'
                #      (just like 'self.lons' and 'self.lats').
                reconstruct_rotation = (
                    self.plate_reconstruction.rotation_model.get_rotation(
                        to_time=reconstruct_age,
                        moving_plate_id=plate_id,
                        from_time=0,
                        anchor_plate_id=anchor_plate_id,
                    )
                    * self.plate_reconstruction.rotation_model.get_rotation(
                        to_time=0,
                        moving_plate_id=plate_id,
                        from_time=self.time,
                        anchor_plate_id=self.anchor_plate_id,
                    )
                )
                reconstructed_points_with_plate_id_and_reconstruct_age = (
                    reconstruct_rotation
                    * reconstructed_points_with_plate_id_and_reconstruct_age
                )

                # Write the reconstructed points.
                lat_lon_points[point_indices_with_plate_id_and_reconstruct_age] = [
                    rpoint.to_lat_lon()
                    for rpoint in reconstructed_points_with_plate_id_and_reconstruct_age
                ]

        rlonslats = lat_lon_points[valid_mask]  # remove invalid points
        rlons = rlonslats[:, 1]
        rlats = rlonslats[:, 0]

        return_tuple = (rlons, rlats)

        if return_point_indices:
            all_point_indices = np.arange(self.size, dtype=int)
            point_indices = all_point_indices[valid_mask]  # remove invalid points
            return_tuple += (point_indices,)

        return return_tuple

    def plate_velocity(
        self,
        time,
        delta_time=1.0,
        *,
        velocity_delta_time_type=pygplates.VelocityDeltaTimeType.t_plus_delta_t_to_t,
        velocity_units=pygplates.VelocityUnits.cms_per_yr,
        earth_radius_in_kms=pygplates.Earth.mean_radius_in_kms,
        anchor_plate_id=None,
        return_reconstructed_points=False,
        return_point_indices=False,
    ):
        """Calculates the east and north components of the tectonic plate velocities of the internal points at a particular geological time.

        The point velocities are calculated using the plate IDs of the internal points and the rotation model of the internal :class:`PlateReconstruction` object.
        If the requested ``time`` differs from the initial time (``self.time``) then the internal points are first reconstructed to ``time`` before calculating velocities.
        Velocities are only calculated at points that are reconstructable (See :ref:`Notes <points-note>`) and that have ages greater than or equal to
        ``time`` (ie, at points that exist at ``time``).

        Parameters
        ----------
        time : float
            The specific geological time (Ma) at which to calculate plate velocities.

        delta_time : float, default=1.0
            The time interval used for velocity calculations. 1.0Ma by default.

        velocity_delta_time_type : pygplates.VelocityDeltaTimeType, default=pygplates.VelocityDeltaTimeType.t_plus_delta_t_to_t
            How the two velocity times are calculated relative to ``time`` (defaults to ``[time + velocity_delta_time, time]``).

        velocity_units : pygplates.VelocityUnits, default=pygplates.VelocityUnits.cms_per_yr
            Whether to return velocities in centimetres per year or kilometres per million years (defaults to centimetres per year).

        earth_radius_in_kms : float, default=pygplates.Earth.mean_radius_in_kms
            Radius of the Earth in kilometres.
            This is only used to calculate velocities (strain rates always use ``pygplates.Earth.equatorial_radius_in_kms``).

        anchor_plate_id : int, optional
            Anchor plate used to reconstruct the points and calculate velocities at their locations.
            By default, reconstructions are made with respect to ``self.anchor_plate_id``
            (which is the anchor plate that the initial points at the initial time are relative to).

        return_reconstructed_points : bool, default=False
            Return the reconstructed points (as longitude and latitude arrays) in addition to the velocities.

        return_point_indices : bool, default=False
            Return the indices of those internal points at which velocities are calculated.
            These are indices into ``self.lons``, ``self.lats``, ``self.plate_id`` and ``self.age``.
            Those points with an age less than ``time`` have not yet appeared at ``time``, and therefore will not have velocities returned.

        Returns
        -------
        velocity_lons, velocity_lats : ndarray
            The velocity arrays containing the **east** (longitude) and *north* (latitude) components of the velocity
            of each internal point that exists at ``time`` (ie, whose age greater than or equal to ``time``).
        rlons, rlats : ndarray
            Only provided if ``return_reconstructed_points`` is True.
            The longitude and latitude coordinate arrays of the reconstructed points (at which velocities are calculated).
            These arrays are the same size as ``velocity_lons`` and ``velocity_lats``.
        point_indices : ndarray
            Only provided if ``return_point_indices`` is True.
            The indices of the returned points (at which velocities are calculated).
            These are indices into ``self.lons``, ``self.lats``, ``self.plate_id`` and ``self.age``.
            This array is the same size as ``velocity_lons`` and ``velocity_lats``.


        .. note::

            The velocities are in **centimetres per year** by default (not **kilometres per million years**, the default
            in :meth:`PlateReconstruction.get_point_velocities`).
            This difference is maintained for backward compatibility.

            For each velocity, the *east** component is first followed by the *north* component.
            This is different to :meth:`PlateReconstruction.get_point_velocities` where the **north** component is first.
            This difference is maintained for backward compatibility.


        .. seealso::

            :meth:`PlateReconstruction.get_point_velocities` : Velocities of points calculated using topologies instead of plate IDs (assigned from static polygons).

        """
        if anchor_plate_id is None:
            anchor_plate_id = self.anchor_plate_id

        # Start with empty arrays.
        north_east_velocities = np.empty((self.size, 2))
        if return_reconstructed_points:
            lat_lon_points = np.empty((self.size, 2))
        else:
            lat_lon_points = np.array([])

        # Determine time interval for velocity calculation.
        if (
            velocity_delta_time_type
            == pygplates.VelocityDeltaTimeType.t_plus_delta_t_to_t
        ):
            from_time = time + delta_time
            to_time = time
        elif (
            velocity_delta_time_type
            == pygplates.VelocityDeltaTimeType.t_to_t_minus_delta_t
        ):
            from_time = time
            to_time = time - delta_time
        elif (
            velocity_delta_time_type
            == pygplates.VelocityDeltaTimeType.t_plus_minus_half_delta_t
        ):
            from_time = time + delta_time / 2
            to_time = time - delta_time / 2
        else:
            raise ValueError(
                "'velocity_delta_time_type' value not one of pygplates.VelocityDeltaTimeType enumerated values"
            )
        # Make sure time interval is non-negative.
        if to_time < 0:
            from_time -= to_time
            to_time = 0

        # Determine which points are valid.
        #
        # These are those points that are reconstructable and have appeared before (or at) 'time'
        # (ie, have a time-of-appearance that's greater than or equal to 'time').
        valid_mask = self._reconstructable & (self.age >= time)

        # Iterate over groups of points with the same plate ID.
        for (
            plate_id,
            point_indices_with_plate_id,
        ) in self._unique_plate_id_groups.items():

            # Determine which points (indices) with the current unique plate ID are valid.
            point_indices_with_plate_id = point_indices_with_plate_id[
                valid_mask[point_indices_with_plate_id]
            ]
            # If none of the points (with the current unique plate ID) are valid then skip to next unique plate ID.
            if point_indices_with_plate_id.size == 0:
                continue

            # Get the reconstructed points with the current unique plate ID that have appeared before (or at) 'time'.
            reconstructed_points_with_plate_id = pygplates.MultiPointOnSphere(
                self.points[point_index] for point_index in point_indices_with_plate_id
            )

            # Stage rotation for the current unique plate ID.
            velocity_equivalent_stage_rotation = (
                self.plate_reconstruction.rotation_model.get_rotation(
                    to_time, plate_id, from_time, anchor_plate_id=anchor_plate_id
                )
            )

            # First reconstruct the internal points from the initial time ('self.time') to present day using
            # our internal anchor plate ID (the same anchor plate used in '__init__').
            # Then reconstruct from present day to 'time' using the *requested* anchor plate ID.
            #
            # Note 'self.points' (and hence 'reconstructed_points_with_plate_id') are the locations at 'self.time'
            #      (just like 'self.lons' and 'self.lats').
            reconstruct_rotation = (
                self.plate_reconstruction.rotation_model.get_rotation(
                    to_time=time,
                    moving_plate_id=plate_id,
                    from_time=0,
                    anchor_plate_id=anchor_plate_id,
                )
                * self.plate_reconstruction.rotation_model.get_rotation(
                    to_time=0,
                    moving_plate_id=plate_id,
                    from_time=self.time,
                    anchor_plate_id=self.anchor_plate_id,
                )
            )
            reconstructed_points_with_plate_id = (
                reconstruct_rotation * reconstructed_points_with_plate_id
            )

            velocity_vectors_with_plate_id = pygplates.calculate_velocities(  # type: ignore
                reconstructed_points_with_plate_id,
                velocity_equivalent_stage_rotation,
                delta_time,
                velocity_units=velocity_units,
                earth_radius_in_kms=earth_radius_in_kms,
            )

            north_east_down_velocities_with_plate_id = (
                pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                    reconstructed_points_with_plate_id, velocity_vectors_with_plate_id
                )
            )

            # Write velocities of points with the current unique plate ID as (north, east) components.
            north_east_velocities[point_indices_with_plate_id] = [
                (ned.get_x(), ned.get_y())  # north, east
                for ned in north_east_down_velocities_with_plate_id
            ]

            # Also write the reconstructed points (if requested).
            if return_reconstructed_points:
                lat_lon_points[point_indices_with_plate_id] = [
                    rpoint.to_lat_lon() for rpoint in reconstructed_points_with_plate_id
                ]

        velocities = north_east_velocities[valid_mask]  # remove invalid points
        velocity_lons = velocities[:, 1]  # east
        velocity_lats = velocities[:, 0]  # north

        return_tuple = velocity_lons, velocity_lats

        if return_reconstructed_points:
            rlonslats = lat_lon_points[valid_mask]  # remove invalid points
            rlons = rlonslats[:, 1]
            rlats = rlonslats[:, 0]
            return_tuple += (rlons, rlats)

        if return_point_indices:
            all_point_indices = np.arange(self.size, dtype=int)
            point_indices = all_point_indices[valid_mask]  # remove invalid points
            return_tuple += (point_indices,)

        return return_tuple

    def motion_path(
        self, time_array, anchor_plate_id=None, return_rate_of_motion=False
    ):
        """Create a path of points to mark the trajectory of a plate's motion through geological time.

        Parameters
        ----------
        time_array : arr
            An array of reconstruction times at which to determine the trajectory of a point on a plate.

            For example,

            .. code-block:: python
                :linenos:

                import numpy as np
                min_time = 30
                max_time = 100
                time_step = 2.5
                time_array = np.arange(min_time, max_time + time_step, time_step)

        anchor_plate_id : int, optional
            Reconstruct features with respect to a certain anchor plate. By default, reconstructions are made
            with respect to the anchor plate ID specified in the :class:`PlateReconstruction` object.
        return_rate_of_motion : bool, default=False
            Choose whether to return the rate of plate motion through time for each

        Returns
        -------
        rlons : ndarray
            An n-dimensional array with columns containing the longitudes of
            the seed points at each timestep in ``time_array``. There are n columns for n seed points.
        rlats : ndarray
            An n-dimensional array with columns containing the latitudes of
            the seed points at each timestep in ``time_array``. There are n columns for n seed points.
        """
        time_array = np.atleast_1d(time_array)

        # ndarrays to fill with reconstructed points and
        # rates of motion (if requested)
        rlons = np.empty((len(time_array), len(self.lons)))
        rlats = np.empty((len(time_array), len(self.lons)))

        StepTimes = np.array([])
        StepRates = np.array([])

        for i, point_feature in enumerate(self.feature_collection):
            # Create the motion path feature
            motion_path_feature = pygplates.Feature.create_motion_path(
                point_feature.get_geometry(),
                time_array.tolist(),
                valid_time=(time_array.max(), time_array.min()),
                relative_plate=int(self.plate_id[i]),
                reconstruction_plate_id=(
                    anchor_plate_id  # if None then uses default anchor plate of 'self.plate_reconstruction'
                    if anchor_plate_id is not None
                    else self.plate_reconstruction.anchor_plate_id
                ),
            )

            reconstructed_motion_paths = self.plate_reconstruction.reconstruct(
                motion_path_feature,
                to_time=0,
                # from_time=0,
                reconstruct_type=pygplates.ReconstructType.motion_path,
                anchor_plate_id=anchor_plate_id,  # if None then uses default anchor plate of 'self.plate_reconstruction'
            )

            # Turn motion paths in to lat-lon coordinates
            trail = np.array([])
            for reconstructed_motion_path in reconstructed_motion_paths:
                # not sure about this. always set the "trail" to the last one in reconstructed_motion_paths?
                # or there is only one path in reconstructed_motion_paths? -- Michael Chin
                # again???
                trail = reconstructed_motion_path.get_motion_path().to_lat_lon_array()

            lon, lat = np.flipud(trail[:, 1]), np.flipud(trail[:, 0])

            rlons[:, i] = lon
            rlats[:, i] = lat

            # Obtain step-plot coordinates for rate of motion
            if return_rate_of_motion is True:
                StepTimes = np.empty(((len(time_array) - 1) * 2, len(self.lons)))
                StepRates = np.empty(((len(time_array) - 1) * 2, len(self.lons)))

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

        if return_rate_of_motion is True:
            return (
                np.squeeze(rlons),
                np.squeeze(rlats),
                np.squeeze(StepTimes),
                np.squeeze(StepRates),
            )
        else:
            return np.squeeze(rlons), np.squeeze(rlats)

    def flowline(
        self, time_array, left_plate_ID, right_plate_ID, return_rate_of_motion=False
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
            The plate ID of the polygon to the left of the spreading ridge.
        right_plate_ID : int
            The plate ID of the polygon to the right of the spreading ridge.
        return_rate_of_motion : bool, default False
            Choose whether to return a step time and step rate array for a step-plot of flowline motion.

        Returns
        -------
        left_lon : ndarray
            The longitudes of the **left** flowline for n seed points.
            There are n columns for n seed points, and m rows for m time steps in ``time_array``.
        left_lat : ndarray
            The latitudes of the **left** flowline of n seed points.
            There are n columns for n seed points, and m rows for m time steps in ``time_array``.
        right_lon : ndarray
            The longitudes of the **right** flowline of n seed points.
            There are n columns for n seed points, and m rows for m time steps in ``time_array``.
        right_lat : ndarray
            The latitudes of the **right** flowline of n seed points.
            There are n columns for n seed points, and m rows for m time steps in ``time_array``.

        Examples
        --------
        To access the i\\ :sup:`th`  seed point's left and right latitudes and longitudes:

        .. code-block:: python
            :linenos:

            for i in np.arange(0,len(seed_points)):
                left_flowline_longitudes = left_lon[:,i]
                left_flowline_latitudes = left_lat[:,i]
                right_flowline_longitudes = right_lon[:,i]
                right_flowline_latitudes = right_lat[:,i]
        """

        model = self.plate_reconstruction
        return model.create_flowline(
            self.lons,
            self.lats,
            time_array,
            left_plate_ID,
            right_plate_ID,
            return_rate_of_motion,
        )

    def _get_dataframe(self):
        import geopandas as gpd

        data = dict()
        data["Longitude"] = self.lons
        data["Latitude"] = self.lats
        data["Plate_ID"] = self.plate_id
        for key in self.attributes:
            data[key] = self.attributes[key]

        return gpd.GeoDataFrame(data)

    def save(self, filename):
        """Saves the feature collection used in the Points object under a given filename to the current directory.

        The file format is determined from the filename extension.

        Parameters
        ----------
        filename : str
            Can be provided as a string including the filename and the file format needed.

        Returns
        -------
        Feature collection saved under given filename to current directory.
        """
        filename = str(filename)

        if filename.endswith((".csv", ".txt", ".dat")):
            df = self._get_dataframe()
            df.to_csv(filename, index=False)

        elif filename.endswith((".xls", ".xlsx")):
            df = self._get_dataframe()
            df.to_excel(filename, index=False)

        elif filename.endswith("xml"):
            df = self._get_dataframe()
            df.to_xml(filename, index=False)

        elif (
            filename.endswith(".gpml")
            or filename.endswith(".gpmlz")
            or filename.endswith(".shp")
        ):
            self.feature_collection.write(filename)

        else:
            raise ValueError(
                "Cannot save to specified file type. Use csv, gpml, shp or xls file extension."
            )

    def rotate_reference_frames(
        self,
        reconstruction_time,
        from_rotation_features_or_model=None,  # filename(s), or pyGPlates feature(s)/collection(s) or a RotationModel
        to_rotation_features_or_model=None,  # filename(s), or pyGPlates feature(s)/collection(s) or a RotationModel
        from_rotation_reference_plate=0,
        to_rotation_reference_plate=0,
        non_reference_plate=701,
        output_name=None,
        return_array=False,
    ):
        """Rotate a grid defined in one plate model reference frame within a :class:`Raster` object to another plate
        reconstruction model reference frame.

        Parameters
        ----------
        reconstruction_time : float
            The time at which to rotate the reconstructed points.
        from_rotation_features_or_model : str/os.PathLike, list of str/os.PathLike, or instance of pygplates.RotationModel
            A filename, or a list of filenames, or a pyGPlates
            RotationModel object that defines the rotation model
            that the input grid is currently associated with.
            ``self.plate_reconstruction.rotation_model`` is default.
        to_rotation_features_or_model : str/os.PathLike, list of str/os.PathLike, or instance of pygplates.RotationModel
            A filename, or a list of filenames, or a pyGPlates
            RotationModel object that defines the rotation model
            that the input grid shall be rotated with.
            ``self.plate_reconstruction.rotation_model`` is default.
        from_rotation_reference_plate : int, default = 0
            The current reference plate for the plate model the points
            are defined in. Defaults to the anchor plate 0.
        to_rotation_reference_plate : int, default = 0
            The desired reference plate for the plate model the points
            to be rotated to. Defaults to the anchor plate 0.
        non_reference_plate : int, default = 701
            An arbitrary placeholder reference frame with which
            to define the "from" and "to" reference frames.
        output_name : str, default None
            If passed, the rotated points are saved as a gpml to this filename.

        Returns
        -------
        Points
            An instance of the :class:`Points` object containing the rotated points.
        """
        if output_name is not None:
            raise NotImplementedError("'output_name' parameter is not implemented")

        if from_rotation_features_or_model is None:
            from_rotation_features_or_model = self.plate_reconstruction.rotation_model
        if to_rotation_features_or_model is None:
            to_rotation_features_or_model = self.plate_reconstruction.rotation_model

        # Create the pygplates.FiniteRotation that rotates
        # between the two reference frames.
        from_rotation_model = pygplates.RotationModel(from_rotation_features_or_model)
        to_rotation_model = pygplates.RotationModel(to_rotation_features_or_model)
        from_rotation = from_rotation_model.get_rotation(
            reconstruction_time,
            non_reference_plate,
            anchor_plate_id=from_rotation_reference_plate,
        )
        to_rotation = to_rotation_model.get_rotation(
            reconstruction_time,
            non_reference_plate,
            anchor_plate_id=to_rotation_reference_plate,
        )
        reference_frame_conversion_rotation = to_rotation * from_rotation.get_inverse()

        # reconstruct points to reconstruction_time
        lons, lats = self.reconstruct(
            reconstruction_time,
            anchor_plate_id=from_rotation_reference_plate,
            return_array=True,
        )  # type: ignore

        # convert FeatureCollection to MultiPointOnSphere
        input_points = pygplates.MultiPointOnSphere(
            (lat, lon) for lon, lat in zip(lons, lats)  # type: ignore
        )

        # Rotate grid nodes to the other reference frame
        output_points = reference_frame_conversion_rotation * input_points

        # Assemble rotated points with grid values.
        out_lon = np.empty_like(self.lons)
        out_lat = np.empty_like(self.lats)
        for i, point in enumerate(output_points):
            out_lat[i], out_lon[i] = point.to_lat_lon()

        if return_array:
            return out_lon, out_lat
        else:
            return Points(
                self.plate_reconstruction,
                out_lon,
                out_lat,
                time=reconstruction_time,
                plate_id=self.plate_id.copy(),
                age=self.age.copy(),
                anchor_plate_id=to_rotation_reference_plate,
            )
