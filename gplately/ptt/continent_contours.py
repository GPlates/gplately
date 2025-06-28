#
#    Copyright (C) 2023-2025 The University of Sydney, Australia
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
Calculate continent contours, and fragmentation index (global perimeter-to-area ratio), at various times.
"""

import math
import os
import os.path
import sys
import time
from collections import deque
from inspect import signature

import numpy as np
import pygplates

from .utils import points_in_polygons, points_spatial_tree, proximity_query

# TODO: Replace the internal uniform lat/lon sampling with a sampling that's uniform across the sphere
#       (so that the contour resolution doesn't favour contours near the North/South poles).


class ContouredContinent(object):
    """
    Class to represent the contour around overlapping/abutting continental blocks.
    """

    def __init__(self):
        self._polygons_including_continent = []
        self._polygons_excluding_continent = []

    def _add_continent(self, continent_polygon):
        """
        Add a continent (landmass) polygon to this contoured continent.

        This is a polygon whose interior represents continental crust.
        And it can have interior rings (holes) which represent oceanic crust.

        Note that it's possible for a continent polygon to be inside the interior hole of another continent polygon.
        For example, a continental island inside an ocean basin that, in turn, is inside a larger continent.
        It's also possible for a continent polygon to be inside an ocean polygon.
        For example, a continental island inside an ocean basin (that itself has no continent containing it).
        """
        self._polygons_including_continent.append(continent_polygon)

    def _add_ocean(self, ocean_polygon):
        """
        Add an ocean polygon to this contoured continent.

        This is a polygon whose interior represents oceanic crust, but it has no continent that contains it.
        Unlike a continent polygon, an ocean polygon cannot have interior rings (holes). If there are any continent islands
        inside this ocean then they should be added as a separate (continent) polygons.

        Usually only a single continent polygon is needed to define a landmass.
        However, one or more ocean polygons can instead be needed if the landmass is actually a landmass that covers the entire globe
        except for a few oceanic holes. In which case the specified ocean polygon is actually one of those oceanic holes.
        This is a special case because you can't have a single global continent polygon with only interior holes (and no exterior ring).
        So intead we allow for multiple ocean polygons that represent these oceanic holes (and treat them specially).
        """
        if ocean_polygon.get_number_of_interior_rings() > 0:
            raise AssertionError("Ocean polygons cannot have interior rings")
        self._polygons_excluding_continent.append(ocean_polygon)

    def get_contours(self):
        """Return the *polyline* contours representing the boundaries of continental crust."""
        contours = []

        # Add each ring (of polygons *including* continent) as a polyline contour.
        for polygon in self._polygons_including_continent:
            exterior_polyline_points = list(polygon.get_exterior_ring_points())
            exterior_polyline_points.append(
                exterior_polyline_points[0]
            )  # polyline's last point should match first point
            contours.append(pygplates.PolylineOnSphere(exterior_polyline_points))
            for interior_ring_index in range(polygon.get_number_of_interior_rings()):
                interior_polyline_points = list(
                    polygon.get_interior_ring_points(interior_ring_index)
                )
                interior_polyline_points.append(
                    interior_polyline_points[0]
                )  # polyline's last point should match first point
                contours.append(pygplates.PolylineOnSphere(interior_polyline_points))

        # Add each *exterior* ring (of polygons *excluding* continent) as a polyline contour.
        # Note: Polygons that *exclude* continental crust have no interior rings.
        for polygon in self._polygons_excluding_continent:
            exterior_polyline_points = list(polygon.get_exterior_ring_points())
            exterior_polyline_points.append(
                exterior_polyline_points[0]
            )  # polyline's last point should match first point
            contours.append(pygplates.PolylineOnSphere(exterior_polyline_points))

        return contours

    def are_points_inside(self, points, points_spatial_tree=None):
        """Returns a numpy 1D boolean array with same length as 'points' (and in same order) containing True for each point inside this contoured continent."""

        # A special (unlikely) case is a single continent covering the entire globe (area 4*pi).
        # This happens when there are no inclusive polygons and no exclusive polygons.
        if (
            not self._polygons_including_continent
            and not self._polygons_excluding_continent
        ):
            # All points *include* continent.
            return np.full(len(points), True)

        # Improve efficiency by re-using spatial tree of points if caller provides it (otherwise create our own).
        if not points_spatial_tree:
            points_spatial_tree = points_spatial_tree.PointsSpatialTree(points)

        # If we have any polygons that *exclude* continent then it means the continent landmass covers the entire globe
        # except these excluding polygons.
        if self._polygons_excluding_continent:
            # By default all points are considered *inside* this contoured continent unless proven *outside*.
            points_inside = np.full(len(points), True)

            # See if the points are inside any of the exclusive polygons.
            exclusive_polygons_containing_points = (
                points_in_polygons.find_polygons_using_points_spatial_tree(
                    points, points_spatial_tree, self._polygons_excluding_continent
                )
            )

            # Any points *inside* an exclusive polygon are considered to be *outside* this contoured continent.
            #
            # Note: If there are any inclusive polygons (which must be inside these exclusive polygons) then later
            #       they will include some of these points that we just excluded.
            for point_index in range(len(points)):
                if exclusive_polygons_containing_points[point_index]:
                    points_inside[point_index] = False

        # else all polygons *include* continent (ie, none exclude continent)...
        else:
            # By default all points are considered *outside* this contoured continent unless proven *inside*.
            points_inside = np.full(len(points), False)

        if self._polygons_including_continent:
            # See if the points are inside any of the inclusive polygons.
            inclusive_polygons_containing_points = (
                points_in_polygons.find_polygons_using_points_spatial_tree(
                    points, points_spatial_tree, self._polygons_including_continent
                )
            )

            # Any points *inside* an inclusive polygon are considered to be *inside* this contoured continent.
            for point_index in range(len(points)):
                if inclusive_polygons_containing_points[point_index]:
                    points_inside[point_index] = True

        return points_inside

    def get_perimeter(self):
        """Sum of the length of the contour boundaries of this contoured continent (in radians)."""
        return math.fsum(
            polygon.get_arc_length()
            for polygon in (
                self._polygons_including_continent + self._polygons_excluding_continent
            )
        )

    def get_area(self):
        """The area of this contoured continent (in square radians, as known as steradians)."""

        # A special (unlikely) case is a single continent covering the entire globe (area 4*pi).
        # This happens when there are no inclusive polygons and no exclusive polygons.
        if (
            not self._polygons_including_continent
            and not self._polygons_excluding_continent
        ):
            return 4 * math.pi

        area = 0.0

        # Note that we can get one or more polygons that *exclude* continent.
        # This can happen when contouring a large landmass such that there is no contour that is an exterior ring.
        # In this case the landmass covers the entire globe except for a few oceanic holes.
        # However you can't have a single polygon with only interior holes (and no exterior ring).
        # So intead we allow for multiple exterior ring polygons to represent these holes (and treat this as a special case).
        if self._polygons_excluding_continent:
            # And since we can't have an exterior ring covering the entire globe we need to explicitly add the area of the entire globe.
            area += 4 * math.pi

            # Subtract the area of the exclusive holes.
            for polygon in self._polygons_excluding_continent:
                area -= polygon.get_area()

        # Add the areas of polygons that include continent.
        for polygon in self._polygons_including_continent:
            area += polygon.get_area()

        return area

    def get_perimeter_area_ratio(self):
        """The perimeter divided by the area (in units of 1/radians)."""
        return self.get_perimeter() / self.get_area()


# Default distance threshold (in radians) above which continents are separated.
#
# This is a small value to avoid numerical tolerance issues when two continent polygons abutt each other but actually
# have a very tiny non-zero distance between them (which would cause them to belong to separate continents, for a zero threshold,
# when they probably should belong to the same continent).
DEFAULT_CONTINENT_SEPARATION_DISTANCE_THRESHOLD_RADIANS = 1e-4  # ~1km


class ContinentContouring(object):
    """
    Class to calculate continent mask, contours, and fragmentation index (global perimeter-to-area ratio), at various times.
    """

    def __init__(
        self,
        rotaton_model_or_features,
        continent_features,  # regular features (not topologies)
        # The grid spacing (in degrees) between points in the grid used for contouring/aggregrating blocks of continental polygons.
        continent_contouring_point_spacing_degrees,
        # Optional parameter specifying a minimum area threshold (in square radians) for including contoured continents.
        #
        # Contoured continents with area smaller than this threshold will be excluded.
        # If this parameter is not specified then no area threshold is applied.
        #
        # Can also be a function (accepting time in Ma) and returning the area threshold.
        #
        # Note: Units here are for normalised sphere (ie, steradians or square radians) so full Earth area is 4*pi.
        #       So 0.1 covers an area of approximately 4,000,000 km^2 (ie, 0.1 * 6371^2, where Earth radius is 6371km).
        #       Conversely 4,000,000 km^2 is equivalent to (4,000,000 / 6371^2) steradians.
        continent_contouring_area_threshold_steradians=None,
        # Optional parameter specifying a distance (in radians) to expand contours ocean-ward - this also
        # ensures small gaps between continents are ignored during contouring.
        #
        # The continent(s) will be expanded by a buffer of this distance (in radians) when contouring/aggregrating blocks of continental polygons.
        # If this parameter is not specified then buffer expansion is not applied (to continent contours).
        #
        # This parameter can also be a function (that returns the distance).
        # The function can have a single function argument, accepting time (in Ma).
        # Or it can have two function arguments, with the second accepting the contoured continent (a 'ContouredContinent' object)
        # of the (unexpanded) contoured continent that the buffer/gap distance will apply to.
        # Or it can have three function arguments, with the third accepting a list of reconstructed polygons ('pygplates.ReconstructedFeatureGeometry' objects)
        # used to contour the (unexpanded) contoured continent that the buffer/gap distance will apply to.
        # Hence a function with *two* arguments means a different buffer/gap distance can be specified for each contoured continent (eg, based on its area).
        # And a function with *three* arguments can also use the feature properties (eg, plate ID) of the reconstructed polygons in the contoured continent.
        #
        # Note: Units here are for normalised sphere (ie, radians).
        #       So 1.0 radian is approximately 6371 km (where Earth radius is 6371 km).
        #       Also 1.0 degree is approximately 110 km.
        #
        # NOTE: This cannot be specified if 'continent_polygon_buffer_and_gap_distance_radians' is specified.
        #       You can only specify one or the other (or neither).
        continent_contouring_buffer_and_gap_distance_radians=None,
        # Optional parameter specifying a minimum area threshold (in square radians) for contours that exclude continental crust.
        #
        # Polygon contours that exclude continental crust and have an area smaller than this threshold will be excluded
        # (meaning they will now *include* continental crust, thus removing the contour).
        # This is useful for removing small holes inside continents.
        # If this parameter is not specified then no area threshold is applied.
        #
        # Can also be a function (accepting time in Ma) and returning the area threshold.
        #
        # Note: Units here are for normalised sphere (ie, steradians or square radians) so full Earth area is 4*pi.
        #       So 0.1 covers an area of approximately 4,000,000 km^2 (ie, 0.1 * 6371^2, where Earth radius is 6371km).
        #       Conversely 4,000,000 km^2 is equivalent to (4,000,000 / 6371^2) steradians.
        continent_exclusion_area_threshold_steradians=None,
        # Optional parameter specifying the distance threshold (in radians) above which continents are separated.
        #
        # Any continent polygons separated by a distance that is less than this threshold will become part of the same continent.
        #
        # If this parameter is not specified then it defaults to a small value to avoid numerical tolerance issues when two
        # continent polygons abutt each other but actually have a very tiny non-zero distance between them
        # (which would cause them to belong to separate continents, for a zero threshold, when they probably should belong to the same continent).
        #
        # Specifying None is the same as specifying a zero distance threshold.
        #
        # Can also be a function (accepting time in Ma) and returning the distance threshold.
        #
        # Note: Units here are for normalised sphere (ie, radians).
        #       So 1.0 radian is approximately 6371 km (where Earth radius is 6371 km).
        #       Also 1.0 degree is approximately 110 km.
        continent_separation_distance_threshold_radians=DEFAULT_CONTINENT_SEPARATION_DISTANCE_THRESHOLD_RADIANS,
        # Optional parameter specifying a distance (in radians) to expand each individual continental polygon ocean-ward - this also
        # ensures small gaps between continents are ignored during contouring.
        #
        # NOTE: This is similar to 'continent_contouring_buffer_and_gap_distance_radians' except it applies to each continental polygon
        #       (instead of applying to each aggregate block of continental polygons forming a continent contour).
        #
        # The continent polygons will be expanded by a buffer of this distance (in radians).
        # If this parameter is not specified then buffer expansion is not applied (to continental polygons).
        #
        # This parameter can also be a function (that returns the distance).
        # The function can have a single function argument, accepting time (in Ma).
        # Or it can have two function arguments, with the second accepting the reconstructed continental feature polygon
        # (a 'pygplates.ReconstructedFeatureGeometry' object) that the buffer/gap distance will apply to.
        # Hence a function with *two* arguments means a different buffer/gap distance can be specified for each continental polygon.
        # For example, you can use its feature properties (eg, plate ID), and/or its reconstructed polygon (eg, area).
        #
        # Note: Units here are for normalised sphere (ie, radians).
        #       So 1.0 radian is approximately 6371 km (where Earth radius is 6371 km).
        #       Also 1.0 degree is approximately 110 km.
        #
        # NOTE: This cannot be specified if 'continent_contouring_buffer_and_gap_distance_radians' is specified.
        #       You can only specify one or the other (or neither).
        continent_polygon_buffer_and_gap_distance_radians=None,
    ):

        self.rotation_model = pygplates.RotationModel(rotaton_model_or_features)
        self.continent_features = continent_features

        if continent_contouring_area_threshold_steradians:
            # Convert area threshold to a function of time, if not already a function.
            if callable(continent_contouring_area_threshold_steradians):
                # We can call the specified function directly.
                self.continent_contouring_area_threshold_steradians_function = (
                    continent_contouring_area_threshold_steradians
                )
            else:
                # Use a delegate function that returns the specified parameter.
                def continent_contouring_area_threshold_steradians_function(age):
                    return continent_contouring_area_threshold_steradians

                self.continent_contouring_area_threshold_steradians_function = (
                    continent_contouring_area_threshold_steradians_function
                )
        else:  # no area threshold (specified either None or zero)
            # Use a delegate function that returns zero.
            def continent_contouring_area_threshold_steradians_function(age):
                return 0

            self.continent_contouring_area_threshold_steradians_function = (
                continent_contouring_area_threshold_steradians_function
            )

        if (
            continent_contouring_buffer_and_gap_distance_radians is not None
            and continent_polygon_buffer_and_gap_distance_radians is not None
        ):
            raise RuntimeError(
                "You cannot specify both 'continent_contouring_buffer_and_gap_distance_radians' and "
                "'continent_polygon_buffer_and_gap_distance_radians'. You can only specify one or the other (or neither)."
            )

        if continent_contouring_buffer_and_gap_distance_radians is not None:
            # Convert continent contouring buffer/gap threshold to a function, if not already a function.
            if callable(continent_contouring_buffer_and_gap_distance_radians):
                callable_signature = signature(
                    continent_contouring_buffer_and_gap_distance_radians
                )
                callable_num_args = len(callable_signature.parameters)
                if not (
                    callable_num_args == 1
                    or callable_num_args == 2
                    or callable_num_args == 3
                ):
                    raise TypeError(
                        "Continent contouring buffer/gap distance is a callable but does not have 1 or 2 or 3 arguments"
                    )
                if callable_num_args == 3:
                    # We can call the specified function directly.
                    self.continent_contouring_buffer_and_gap_distance_radians_function = (
                        continent_contouring_buffer_and_gap_distance_radians
                    )
                elif callable_num_args == 2:
                    # The specified function only accepts age and contoured continent (not continent feature polygons).
                    # So use a delegate function that calls it and ignores continent feature polygons.
                    def continent_contouring_buffer_and_gap_distance_radians_function(
                        age, contoured_continent, continent_feature_polygons
                    ):
                        return continent_contouring_buffer_and_gap_distance_radians(
                            age, contoured_continent
                        )

                    self.continent_contouring_buffer_and_gap_distance_radians_function = (
                        continent_contouring_buffer_and_gap_distance_radians_function
                    )
                else:  # callable_num_args == 1
                    # The specified function only accepts age (not contoured continent or continent feature polygons).
                    # So use a delegate function that calls it and ignores contoured continent and continent feature polygons.
                    def continent_contouring_buffer_and_gap_distance_radians_function(
                        age, contoured_continent, continent_feature_polygons
                    ):
                        return continent_contouring_buffer_and_gap_distance_radians(age)

                    self.continent_contouring_buffer_and_gap_distance_radians_function = (
                        continent_contouring_buffer_and_gap_distance_radians_function
                    )
            else:
                # Use a delegate function that returns the specified parameter.
                def continent_contouring_buffer_and_gap_distance_radians_function(
                    age, contoured_continent, continent_feature_polygons
                ):
                    return continent_contouring_buffer_and_gap_distance_radians

                self.continent_contouring_buffer_and_gap_distance_radians_function = (
                    continent_contouring_buffer_and_gap_distance_radians_function
                )
        else:
            self.continent_contouring_buffer_and_gap_distance_radians_function = None

        if continent_polygon_buffer_and_gap_distance_radians is not None:
            # Convert continent polygon buffer/gap threshold to a function, if not already a function.
            if callable(continent_polygon_buffer_and_gap_distance_radians):
                callable_signature = signature(
                    continent_polygon_buffer_and_gap_distance_radians
                )
                callable_num_args = len(callable_signature.parameters)
                if not (callable_num_args == 1 or callable_num_args == 2):
                    raise TypeError(
                        "Continent polygon buffer/gap distance is a callable but does not have 1 or 2 arguments"
                    )
                if callable_num_args == 2:
                    # We can call the specified function directly.
                    self.continent_polygon_buffer_and_gap_distance_radians_function = (
                        continent_polygon_buffer_and_gap_distance_radians
                    )
                else:  # callable_num_args == 1
                    # The specified function only accepts age (not continent feature polygon).
                    # So use a delegate function that calls it and ignores the continent feature polygon.
                    def continent_polygon_buffer_and_gap_distance_radians_function(
                        age, continent_feature_polygon
                    ):
                        return continent_polygon_buffer_and_gap_distance_radians(age)

                    self.continent_polygon_buffer_and_gap_distance_radians_function = (
                        continent_polygon_buffer_and_gap_distance_radians_function
                    )
            else:
                # Use a delegate function that returns the specified parameter.
                def continent_polygon_buffer_and_gap_distance_radians_function(
                    age, continent_feature_polygon
                ):
                    return continent_polygon_buffer_and_gap_distance_radians

                self.continent_polygon_buffer_and_gap_distance_radians_function = (
                    continent_polygon_buffer_and_gap_distance_radians_function
                )
        else:
            self.continent_polygon_buffer_and_gap_distance_radians_function = None

        if continent_exclusion_area_threshold_steradians:
            # Convert area threshold to a function of time, if not already a function.
            if callable(continent_exclusion_area_threshold_steradians):
                # We can call the specified function directly.
                self.continent_exclusion_area_threshold_steradians_function = (
                    continent_exclusion_area_threshold_steradians
                )
            else:
                # Use a delegate function that returns the specified parameter.
                def continent_exclusion_area_threshold_steradians_function(age):
                    return continent_exclusion_area_threshold_steradians

                self.continent_exclusion_area_threshold_steradians_function = (
                    continent_exclusion_area_threshold_steradians_function
                )
        else:  # no area threshold (specified either None or zero)
            # Use a delegate function that returns zero.
            def continent_exclusion_area_threshold_steradians_function(age):
                return 0

            self.continent_exclusion_area_threshold_steradians_function = (
                continent_exclusion_area_threshold_steradians_function
            )

        if continent_separation_distance_threshold_radians:
            # Convert distance threshold to a function of time, if not already a function.
            if callable(continent_separation_distance_threshold_radians):
                # We can call the specified function directly.
                self.continent_separation_distance_threshold_radians_function = (
                    continent_separation_distance_threshold_radians
                )
            else:
                # Use a delegate function that returns the specified parameter.
                def continent_separation_distance_threshold_radians_function(age):
                    return continent_separation_distance_threshold_radians

                self.continent_separation_distance_threshold_radians_function = (
                    continent_separation_distance_threshold_radians_function
                )
        else:  # no distance threshold (specified either None or zero)
            # Use a delegate function that returns zero.
            def continent_separation_distance_threshold_radians_function(age):
                return 0

            self.continent_separation_distance_threshold_radians_function = (
                continent_separation_distance_threshold_radians_function
            )

        # The number of latitudes (including -90 and 90).
        self.contouring_grid_num_latitudes = (
            int(math.ceil(180.0 / continent_contouring_point_spacing_degrees)) + 1
        )
        # The number of longitudes (including -180 and 180).
        self.contouring_grid_num_longitudes = (
            2 * (self.contouring_grid_num_latitudes - 1) + 1
        )

        self.contouring_point_spacing_degrees = 180.0 / (
            self.contouring_grid_num_latitudes - 1
        )

        # A point grid to calculate contour polygons representing the boundary of reconstructed static polygons that overlap each other.
        #
        # NOTE: We must generate points on the dateline (ie, at both longitude -180 and 180) since the
        #       contouring alorithm depends on this. We also generate points at the North and South poles
        #       for the same reason.
        lats = np.linspace(-90.0, 90.0, self.contouring_grid_num_latitudes)
        lons = np.linspace(-180.0, 180.0, self.contouring_grid_num_longitudes)

        # Create a multipoint grid of points ordered by longitude first then latitude.
        contouring_longitude_array, contouring_latitude_array = np.meshgrid(lons, lats)
        self.contouring_points = pygplates.MultiPointOnSphere(
            zip(
                contouring_latitude_array.flatten(),
                contouring_longitude_array.flatten(),
            )
        )

        # Improve efficiency by re-using spatial tree of contouring points over time (when finding points in polygons and finding points near polygons).
        #
        # First calculate the subdivision depth to avoid doing too many point-in-polygon tests (for example) for each spatial tree node.
        # The lat/lon width of a root quad tree node in the spatial tree is 90 degrees which is 'n/2' points wide (where 'n' is 'self.contouring_grid_num_latitudes').
        # So a leaf node at subdivision depth 'D' is '(n/2) / 2^D' = 'n / 2^(D+1)'. The number of points is the square of that 'N = '(n / 2^(D+1)) ^ 2'.
        # Rearranging that gives the subdivision depth 'D' in terms of the number of points we'd like in a deepest (leaf) node N:
        #   D = log2(n / sqrt(N)) - 1
        max_num_points_per_spatial_tree_node = 64  # N
        points_spatial_tree_subdivision_depth = math.ceil(
            math.log(
                self.contouring_grid_num_latitudes
                / math.sqrt(max_num_points_per_spatial_tree_node),
                2,
            )
            - 1
        )  # D
        # We won't exceed 6 subdivision levels because it starts to use a lot of memory.
        # And 6 is still higher than the default of 4 in PointsSpatialTree.
        points_spatial_tree_subdivision_depth = min(
            6, points_spatial_tree_subdivision_depth
        )
        self.contouring_points_spatial_tree = points_spatial_tree.PointsSpatialTree(
            self.contouring_points, points_spatial_tree_subdivision_depth
        )

    def get_fragmentation(self, age):
        """
        Calculate the continental fragmentation index (global perimeter-to-area ratio) at the specified time.
        """

        # Get the contoured continents representing the boundary(s) of the reconstructed continent polygons that overlap each other.
        contoured_continents = self.get_contoured_continents(age)

        total_perimeter = 0.0
        total_area = 0.0

        # Update total perimeter and area.
        for contoured_continent in contoured_continents:
            total_perimeter += contoured_continent.get_perimeter()
            total_area += contoured_continent.get_area()

        # Avoid divide-by-zero.
        if total_area == 0.0:
            return 0.0

        # print(' global perimeter/area:', total_perimeter / total_area / pygplates.Earth.equatorial_radius_in_kms, 'km-1'); sys.stdout.flush()

        # print('age:', age, 'frag_index (1/km):', total_perimeter / total_area / 6371.0); sys.stdout.flush()
        return total_perimeter / total_area

    def get_continent_mask(self, age):
        """
        Reconstruct the continents (specified in constructor) and then find the latitude/longitude grid points that are inside them.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        Note that when writing to a NetCDF grid file you can convert to floating-point (with "continent_mask.astype('float')").
        """

        contoured_continents = self.get_contoured_continents(age)

        return self.calculate_continent_mask(contoured_continents)

    def get_contoured_continents(self, age):
        """
        Reconstruct the continents (specified in constructor) and then find their boundaries as contoured continents.

        Returns a list of 'ContouredContinent'.
        """

        reconstructed_continent_polygons = self.get_reconstructed_continent_polygons(
            age
        )

        return self.calculate_contoured_continents(
            reconstructed_continent_polygons, age
        )

    def get_continent_mask_and_contoured_continents(self, age):
        """
        Reconstruct the continents (specified in constructor) and then find both their boundaries as contoured continents and
        the latitude/longitude grid points that are inside them.

        Returns a 2-tuple of (a 2D boolean numpy array of shape (num_latitudes, num_longitudes), a list of 'ContouredContinent').
        """

        contoured_continents = self.get_contoured_continents(age)

        continent_mask = self.calculate_continent_mask(contoured_continents)

        return continent_mask, contoured_continents

    def get_reconstructed_continent_polygons(self, age):
        """
        Reconstruct the continents (specified in constructor).

        Note that these are just the original continent polygons (but reconstructed).
        They are NOT contoured, so they can still overlap/abutt each other.

        Returns a list of 2-tuple ('pygplates.PolygonOnSphere', 'pygplates.ReconstructedFeatureGeometry') where
        the polygon is obtained from the reconstructed feature geometry.
        The reconstructed feature geometry can be used to query the associated 'pygplates.Feature' and its properties.
        """

        # Reconstruct static continental polygons.
        reconstructed_feature_geometries = []
        pygplates.reconstruct(
            self.continent_features,
            self.rotation_model,
            reconstructed_feature_geometries,
            age,
        )

        # Return a list of 2-tuple ('pygplates.PolygonOnSphere', 'pygplates.ReconstructedFeatureGeometry').
        #
        # We should have polygons (not polylines) but turn into a polygon if happens to be a polyline
        # (but that actually only works if the polyline is a closed loop and not just part of a polygon's boundary).
        return [
            (
                pygplates.PolygonOnSphere(
                    reconstructed_feature_geometry.get_reconstructed_geometry()
                ),
                reconstructed_feature_geometry,
            )
            for reconstructed_feature_geometry in reconstructed_feature_geometries
        ]

    def calculate_continent_mask(self, contoured_continents):
        """
        Find the latitude/longitude grid points that are inside the specified contoured continents.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        Note that when writing to a NetCDF grid file you can convert to floating-point (with "continent_mask.astype('float')").
        """

        return self._find_grid_points_inside_contoured_continents(contoured_continents)

    def calculate_contoured_continents(self, continent_polygons, age=0):
        """
        Find the boundaries of the specified (potentially overlapping/abutting) continent polygons as contoured continents.

        Note that each continent polygon should be a 2-tuple ('pygplates.PolygonOnSphere', 'pygplates.ReconstructedFeatureGeometry').

        Note that small contoured continent islands with area less than the area threshold will NOT get returned.

        The 'age' is only used to look up the time-dependent thresholds (passed into constructor).
        If threshold does not vary with time then 'age' does not need to be specified (defaults to present day).

        Returns a list of 'ContouredContinent'.
        """

        # time1 = time.time()

        continent_separation_distance_threshold_radians = (
            self.continent_separation_distance_threshold_radians_function(age)
        )

        if self.continent_polygon_buffer_and_gap_distance_radians_function:

            # Convert 2-tuple of continent polygons to a 3-tuple where 2nd element is each continent polygon's buffer distance.
            continent_polygons = [
                (
                    polygon,
                    # Buffer distance for the current continent polygon...
                    self.continent_polygon_buffer_and_gap_distance_radians_function(
                        age, continent_feature_polygon
                    ),
                    continent_feature_polygon,
                )
                for polygon, continent_feature_polygon in continent_polygons
            ]

            # Find groups of continent polygons where each polygon in a group is within the specified distance of at least one other polygon in the group.
            continent_polygon_groups = self._find_continent_polygon_groups(
                continent_polygons, continent_separation_distance_threshold_radians
            )

            contoured_continents = []

            # Create the contoured continents, excluding those with area below the area threshold (if specified).
            for continent_polygons_in_group in continent_polygon_groups:

                # Find the grid points inside or near the current continent's polygons.
                #
                # Note: Each continental polygon may have a different buffer/gap distance (affecting which points are near each polygon).
                grid_points_inside_continent = (
                    self._find_grid_points_inside_or_near_continent_polygons(
                        continent_polygons_in_group
                    )
                )

                # Skip the current continent if its polygons (with buffer expansion) are too small such that they miss all the grid points.
                if not np.any(grid_points_inside_continent):
                    continue

                # Contour the grid points that are inside the current continent's polygons.
                contoured_continent = self._create_contoured_continent(
                    grid_points_inside_continent
                )

                # If the area threshold is non-zero then exclude the current contoured continents if its area is below the threshold.
                continent_contouring_area_threshold_steradians = (
                    self.continent_contouring_area_threshold_steradians_function(age)
                )
                if (
                    continent_contouring_area_threshold_steradians > 0
                    and contoured_continent.get_area()
                    < continent_contouring_area_threshold_steradians
                ):
                    continue

                contoured_continents.append(contoured_continent)

        else:  # not self.continent_polygon_buffer_and_gap_distance_radians_function ...

            # Convert 2-tuple of continent polygons to a 3-tuple where 2nd element is the polygon buffer distance of zero.
            continent_polygons = [
                (polygon, 0.0, continent_feature_polygon)
                for polygon, continent_feature_polygon in continent_polygons
            ]

            # Find groups of continent polygons where each polygon in a group is within the specified distance of at least one other polygon in the group.
            continent_polygon_groups = self._find_continent_polygon_groups(
                continent_polygons, continent_separation_distance_threshold_radians
            )

            continents = []

            # Create the initial contoured continents, only excluding those with area below the area threshold (if specified).
            for continent_polygons_in_group in continent_polygon_groups:
                # Find the grid points inside the current continent's polygons.
                #
                # Note: Each continental polygon has a zero buffer/gap distance
                #       (and so we don't need to consider points *near* each polygon).
                grid_points_inside_continent = (
                    self._find_grid_points_inside_continent_polygons(
                        continent_polygons_in_group
                    )
                )

                # Skip the current continent if its polygons are too small such that they miss all the grid points.
                if not np.any(grid_points_inside_continent):
                    continue

                # Contour the grid points that are inside the current continent's polygons.
                contoured_continent = self._create_contoured_continent(
                    grid_points_inside_continent
                )

                # If the area threshold is non-zero then exclude the current contoured continents if its area is below the threshold.
                continent_contouring_area_threshold_steradians = (
                    self.continent_contouring_area_threshold_steradians_function(age)
                )
                if (
                    continent_contouring_area_threshold_steradians > 0
                    and contoured_continent.get_area()
                    < continent_contouring_area_threshold_steradians
                ):
                    continue

                if self.continent_contouring_buffer_and_gap_distance_radians_function:
                    # The buffer distance for the current contoured continent.
                    #
                    # Note: Each continent polygon is actually an n-tuple with the third element being a pygplates.ReconstructedFeatureGeometry.
                    #       Passing 'pygplates.ReconstructedFeatureGeometry's to the buffer/gap distance function helps it decide the appropriate
                    #       buffer/gap for the contoured continent (that contours the associated polygons). For example, the function can look
                    #       at the plate IDs of the polygons (via their pygplates.Feature obtained from 'continent_feature_polygon.get_feature()').
                    continent_feature_polygons = [
                        continent_polygon[2]
                        for continent_polygon in continent_polygons_in_group
                    ]
                    contouring_buffer_and_gap_distance_radians = self.continent_contouring_buffer_and_gap_distance_radians_function(
                        age, contoured_continent, continent_feature_polygons
                    )

                else:
                    contouring_buffer_and_gap_distance_radians = 0

                # Add the current continent.
                continents.append(
                    self._Continent(
                        contoured_continent,
                        continent_polygons_in_group,
                        grid_points_inside_continent,
                        contouring_buffer_and_gap_distance_radians,
                    )
                )

            # time2 = time.time()
            # print(' contour continents({}): {:.2f}'.format(len(continents), time2 - time1))

            # If any continent has a non-zero buffer/gap distance expansion then this could cause it to join with nearby continents forming a single merged continent.
            merged_continents = self._find_merged_continents(
                continents, continent_separation_distance_threshold_radians
            )

            contoured_continents = []

            # Contour each merged continent.
            for merged_continent in merged_continents:
                # If any continents in the current merged continent have non-zero buffer/gap distances then we'll need to expand
                # those continents and re-contour the entire list of (merged) continents.
                if any(
                    continent.contouring_buffer_and_gap_distance_radians
                    for continent in merged_continent.continents
                ):

                    # The grids points inside the merged continent include the grid points inside all its (merged) continents.
                    grid_points_inside_merged_continent = np.full(
                        len(self.contouring_points), False
                    ).reshape(
                        (
                            self.contouring_grid_num_latitudes,
                            self.contouring_grid_num_longitudes,
                        )
                    )
                    for continent in merged_continent.continents:
                        grid_points_inside_merged_continent[
                            continent.grid_points_inside_continent
                        ] = True

                    # Find the grid points near the current merged continent's polygons.
                    #
                    # Note: Each continent (in the merged continent) may have a different buffer/gap distance.
                    grid_points_near_merged_continent = (
                        self._find_grid_points_near_merged_continent(merged_continent)
                    )
                    # Add these nearby grid points to those inside the merged continent.
                    grid_points_inside_merged_continent[
                        grid_points_near_merged_continent
                    ] = True

                    # Contour the grid points that are inside the merged continent's polygons.
                    contoured_continent = self._create_contoured_continent(
                        grid_points_inside_merged_continent
                    )

                elif len(merged_continent.continents) == 1:
                    # There's only one continent and it has no buffer/gap distance, so its contour will also be the merged continent's contour.
                    contoured_continent = merged_continent.continents[
                        0
                    ].contoured_continent
                else:
                    raise AssertionError(
                        "Shouldn't have multiple merged continents all with zero buffer/gap distances"
                    )

                contoured_continents.append(contoured_continent)

            # time3 = time.time()
            # print(' contour merged continents({}): {:.2f}'.format(len(merged_continents), time3 - time2))

        # Remove any ocean areas (regions which exclude continental crust) that are below the exclusion area threshold.
        continent_exclusion_area_threshold_steradians = (
            self.continent_exclusion_area_threshold_steradians_function(age)
        )
        if continent_exclusion_area_threshold_steradians > 0:
            self._remove_ocean_areas_below_exclusion_threshold(
                contoured_continents, continent_exclusion_area_threshold_steradians
            )

        # time_end = time.time()
        # print('calculate_contoured_continents({}): {:.2f}'.format(len(contoured_continents), time_end - time1))

        return contoured_continents

    class _Continent(object):
        """Private inner class containing information about a continent (before continents are merged due to non-zero buffer/gap distances)."""

        def __init__(
            self,
            contoured_continent,
            continent_polygons,
            grid_points_inside_continent,
            contouring_buffer_and_gap_distance_radians,
        ):
            self.contoured_continent = contoured_continent
            self.continent_polygons = continent_polygons
            self.grid_points_inside_continent = grid_points_inside_continent
            self.contouring_buffer_and_gap_distance_radians = (
                contouring_buffer_and_gap_distance_radians
            )

    class _MergedContinent(object):
        """Private inner class containing information about a merged continent (referencing several continents merged due to non-zero buffer/gap distances)."""

        def __init__(self, continents):
            self.continents = continents
            self.contoured_continent = None  # to be set later (after constructor)

    def _find_merged_continents(
        self, continents, continent_separation_distance_threshold_radians
    ):
        """
        Merge any continents that are within the continent separation distance of each other after expansion by their buffer/gap distances.
        """

        # time1 = time.time()

        def _are_continents_near_each_other(continent1_index, continent2_index):
            continent1 = continents[continent1_index]
            continent2 = continents[continent2_index]
            # If both continents have no buffer/gap expansion then we already know they are farther apart
            # than the continent separation distance (otherwise they wouldn't be separate continents).
            # Hence they are not near each other.
            if (
                continent1.contouring_buffer_and_gap_distance_radians == 0
                and continent2.contouring_buffer_and_gap_distance_radians == 0
            ):
                return False
            # The distance threshold is the continent separation distance plus the sum of the buffer/gap distance of both continents.
            # And we clamp to a maximum of PI.
            distance_threshold_radians = min(
                math.pi,
                continent_separation_distance_threshold_radians
                + continent1.contouring_buffer_and_gap_distance_radians
                + continent2.contouring_buffer_and_gap_distance_radians,
            )
            # Test all pairs of polygons between each continent.
            #
            # Note: Each continent polygon is actually a tuple where the first element is a pygplates.PolygonOnSphere.
            #       The second element of tuple is the buffer distance of each continent polygon.
            #       But since we're in this function then that should be zero
            #       (ie, we're only using buffer distances for *contoured* continents, not individual polygons).
            for polygon1, *_ in continent1.continent_polygons:
                for polygon2, *_ in continent2.continent_polygons:
                    # See if the two continent polygons are near each other (within the distance threshold).
                    if (
                        pygplates.GeometryOnSphere.distance(
                            polygon1,
                            polygon2,
                            distance_threshold_radians,
                            geometry1_is_solid=True,
                            geometry2_is_solid=True,
                        )
                        is not None
                    ):
                        return True
            return False

        # Each continent will have a list of other continents that are near it.
        nearby_continent_indices = [
            [] for _ in range(len(continents))
        ]  # initially a list of empty lists
        # Find the nearby continents.
        for continent1_index in range(len(continents) - 1):
            for continent2_index in range(continent1_index + 1, len(continents)):
                if _are_continents_near_each_other(continent1_index, continent2_index):
                    # Add links in both directions.
                    nearby_continent_indices[continent1_index].append(continent2_index)
                    nearby_continent_indices[continent2_index].append(continent1_index)

        # Whether a continent has been added to a group yet.
        have_added_continent_to_a_group = [False] * len(continents)

        def _add_nearby_continents_to_group(continent_index_group, continent_index):
            # Iterate over the continents near 'continent_index'.
            for nearby_continent_index in nearby_continent_indices[continent_index]:
                # Only add nearby continent to the group if it hasn't already been added to a group.
                if not have_added_continent_to_a_group[nearby_continent_index]:
                    # Add the nearby continent and mark it as having been added.
                    continent_index_group.append(nearby_continent_index)
                    have_added_continent_to_a_group[nearby_continent_index] = True
                    # Recursively add continents near the current nearby continent.
                    _add_nearby_continents_to_group(
                        continent_index_group, nearby_continent_index
                    )

        merged_continents = []

        # Create the merged continents.
        for continent_index in range(len(continents)):
            if not have_added_continent_to_a_group[continent_index]:
                # Add the current continent and mark it as having been added.
                continent_index_group = [continent_index]
                have_added_continent_to_a_group[continent_index] = True
                # Recursively add any nearby continents to the same group (if they haven't already been added to a group).
                _add_nearby_continents_to_group(continent_index_group, continent_index)
                # Create a merged continent from the group of nearby continents.
                merged_continents.append(
                    self._MergedContinent(
                        [continents[index] for index in continent_index_group]
                    )
                )

        if sum(
            len(merged_continent.continents) for merged_continent in merged_continents
        ) != len(continents):
            raise AssertionError("Sum of continents in merged continents is incorrect")

        # time2 = time.time()
        # print(' _find_merged_continents({}): {:.2f}'.format(len(merged_continents), time2 - time1))

        return merged_continents

    def _remove_ocean_areas_below_exclusion_threshold(
        self, contoured_continents, continent_exclusion_area_threshold_steradians
    ):
        """
        Remove any ocean polygon rings in contoured continents that are below the exclusion area threshold.

        And remove any continent landmasses inside those removed ocean polygon rings.
        """

        # Return early if there are no contoured continents.
        if not contoured_continents:
            return

        # time1 = time.time()

        #
        # Find those ocean polygon rings (ie, which exclude continental crust) that are below the exclusion area threshold.
        # These rings will now include (rather than exclude) continental crust.
        #
        # So we'll first remove the offending interior rings from their containing polygon (in the case of continent polygons) or
        # remove the polygon altogether (in the case of ocean polygons).
        #
        # Then we'll find all continent polygons that are inside these removed ocean polygon rings and remove them altogether
        # (since they are inside an ocean that no longer exists and hence they are no longer a separate continent).
        #
        # Doing all this will remove those ocean holes below the threshold, effectively turning them into continental crust.
        #

        # First find all ocean polygon rings below threshold and remove them.
        removed_ocean_ring_polygons = []
        for contoured_continent in contoured_continents:
            # For each polygon that *includes* continental crust we'll look at its *interior* rings (which exclude continental crust).
            for polygon_index, polygon in enumerate(
                contoured_continent._polygons_including_continent
            ):
                interior_rings_to_keep = []
                for interior_ring_index in range(
                    polygon.get_number_of_interior_rings()
                ):
                    interior_ring_points = polygon.get_interior_ring_points(
                        interior_ring_index
                    )
                    interior_ring_polygon = pygplates.PolygonOnSphere(
                        interior_ring_points
                    )
                    if (
                        interior_ring_polygon.get_area()
                        < continent_exclusion_area_threshold_steradians
                    ):
                        removed_ocean_ring_polygons.append(interior_ring_polygon)
                    else:
                        interior_rings_to_keep.append(interior_ring_points)

                # If some interior rings have area below the exclusion threshold then remove them by
                # creating a new polygon without them (need to do this because polygons are immutable).
                if len(interior_rings_to_keep) < polygon.get_number_of_interior_rings():
                    polygon = pygplates.PolygonOnSphere(
                        polygon.get_exterior_ring_points(), interior_rings_to_keep
                    )
                    contoured_continent._polygons_including_continent[polygon_index] = (
                        polygon
                    )

            # For each polygon that *excludes* continental crust we'll look at its *exterior* ring (which excludes continental crust).
            polygon_index = 0
            while polygon_index < len(
                contoured_continent._polygons_excluding_continent
            ):
                polygon = contoured_continent._polygons_excluding_continent[
                    polygon_index
                ]
                exterior_ring_polygon = pygplates.PolygonOnSphere(
                    polygon.get_exterior_ring_points()
                )
                if (
                    exterior_ring_polygon.get_area()
                    < continent_exclusion_area_threshold_steradians
                ):
                    removed_ocean_ring_polygons.append(exterior_ring_polygon)
                    # Remove the polygon.
                    # Note: Polygons that *exclude* continental crust have no interior rings.
                    del contoured_continent._polygons_excluding_continent[polygon_index]
                    continue
                polygon_index += 1

        def _is_polygon_inside_removed_ocean_ring_polygon(
            polygon, removed_ocean_ring_polygon
        ):
            # Find a vertex of the polygon that is not right *on* the outline of the removed ocean ring polygon
            # (it's possible the two polygons are abutting each other, ie, sharing part of their boundary outline).
            # If it's *on* the outline then we can't do a point-in-polygon test to determine whether polygon is inside or outside.
            for point_on_polygon in polygon.get_exterior_ring_points():
                if (
                    pygplates.GeometryOnSphere.distance(
                        point_on_polygon, removed_ocean_ring_polygon, 1e-5
                    )
                    is None
                ):
                    # Found a polygon vertex that is away from the outline of the removed ocean ring polygon.
                    # So, if the vertex is inside the removed ocean ring polygon then the entire polygon must also be inside it.
                    return removed_ocean_ring_polygon.is_point_in_polygon(
                        point_on_polygon
                    )

            # All points on the polygon's exterior ring are *on* the outline of the removed ocean ring polygon.
            # This is extremely unlikely, but if it happens then consider the polygon to be inside.
            return True

        # Next find all continent polygons that are inside the removed ocean polygon rings and remove them altogether.
        contoured_continent_index = 0
        while contoured_continent_index < len(contoured_continents):
            contoured_continent = contoured_continents[contoured_continent_index]

            polygon_index = 0
            while polygon_index < len(
                contoured_continent._polygons_including_continent
            ):
                polygon = contoured_continent._polygons_including_continent[
                    polygon_index
                ]
                for removed_ocean_ring_polygon in removed_ocean_ring_polygons:
                    if _is_polygon_inside_removed_ocean_ring_polygon(
                        polygon, removed_ocean_ring_polygon
                    ):
                        # Remove the current polygon contour from the contoured continent.
                        del contoured_continent._polygons_including_continent[
                            polygon_index
                        ]
                        polygon_index -= (
                            1  # undo the subsequent increment to next polygon
                        )
                        break  # skip to the next polygon
                polygon_index += 1

            # If the current contoured continent has no polygons then remove the contoured continent.
            if (
                not contoured_continent._polygons_including_continent
                and not contoured_continent._polygons_excluding_continent
            ):
                del contoured_continents[contoured_continent_index]
                continue

            contoured_continent_index += 1

        # If we are left with no contoured continents (after starting with at least one) then add a
        # single empty contoured continent (which implies a single landmass covering the entire globe).
        #
        # This can happen there was a very large landmass that only had polygons that *excluded* continental crust,
        # and all those *excluding* polygons were removed, leaving a landmass that covered the entire globe.
        if not contoured_continents:
            contoured_continents.append(ContouredContinent())

        # time2 = time.time()
        # print(' _remove_ocean_areas_below_exclusion_threshold({}): {:.2f}'.format(len(removed_ocean_ring_polygons), time2 - time1))

    def _find_continent_polygon_groups(
        self, continent_polygons, continent_separation_distance_threshold_radians
    ):
        """
        Find groups of polygons where each polygon in a group (when expanded by its individual polygon buffer distance) is within the
        specified separation distance of at least one other polygon in the group (also expanded by its individual polygon buffer distance).

        Note that each continent polygon should be an n-tuple where the first element is a pygplates.PolygonOnSphere and
        the second is its buffer/gap distance (in radians).
        Subsequent tuple elements (beyond the two) are optional, and the full tuple will be passed intact to the output groups.

        This is useful when creating an individual continent for each group.
        """

        # time1 = time.time()

        continent_polygon_groups = []
        for continent_polygon in continent_polygons:
            polygon, polygon_buffer_distance_radians, *_ = continent_polygon

            # See if the current continent polygon is near any polygon in any group.
            continent_polygon_group_index = None  # index of first group found (if any)

            # Iterate over the polygon groups.
            group_index = 0
            while group_index < len(continent_polygon_groups):
                # Iterate over polygons in the current group.
                for (
                    polygon_in_group,
                    polygon_in_group_buffer_distance_radians,
                    *_,
                ) in continent_polygon_groups[group_index]:

                    # The distance threshold is the continent separation distance plus the sum of each polygon's buffer distance, clamped to a maximum of PI.
                    distance_threshold_radians = min(
                        math.pi,
                        continent_separation_distance_threshold_radians
                        + polygon_buffer_distance_radians
                        + polygon_in_group_buffer_distance_radians,
                    )

                    # See if the current continent polygon is near the current polygon in the current group.
                    if (
                        pygplates.GeometryOnSphere.distance(
                            polygon,
                            polygon_in_group,
                            distance_threshold_radians,
                            geometry1_is_solid=True,
                            geometry2_is_solid=True,
                        )
                        is not None
                    ):

                        # If the current continent polygon hasn't been added to a group yet then add it now.
                        if continent_polygon_group_index is None:
                            continent_polygon_group_index = group_index
                            continent_polygon_groups[
                                continent_polygon_group_index
                            ].append(continent_polygon)
                        # Otherwise it is near another group, so merge that group into the current continent polygon's group.
                        else:
                            # Merge the current group into group that the current continent polygon belongs to.
                            continent_polygon_groups[
                                continent_polygon_group_index
                            ] += continent_polygon_groups[group_index]
                            # And then delete the current group.
                            del continent_polygon_groups[group_index]
                            group_index -= (
                                1  # undo the subsequent increment to next group
                            )

                        # Finished visiting polygons in the current group, so skip to the next group.
                        break

                # Next group.
                group_index += 1

            # If the current continent polygon is not near any group then add it to a new group.
            if continent_polygon_group_index is None:
                continent_polygon_groups.append([continent_polygon])

        # time2 = time.time()
        # print(' _find_continent_polygon_groups({}): {:.2f}'.format(len(continent_polygon_groups), time2 - time1))

        return continent_polygon_groups

    def _find_grid_points_inside_continent_polygons(self, continent_polygons):
        """
        Find the latitude/longitude grid points that are inside (one or more of) the specified polygons.

        Note that each continent polygon should be an n-tuple where the first element is a pygplates.PolygonOnSphere.
        Subsequent tuple elements (beyond the first) are optional and ignored.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        """

        # time1 = time.time()

        # Find the polygon (if any) containing each grid point.
        polygons = [continent_polygon[0] for continent_polygon in continent_polygons]
        polygons_containing_points = (
            points_in_polygons.find_polygons_using_points_spatial_tree(
                self.contouring_points, self.contouring_points_spatial_tree, polygons
            )
        )

        # time2 = time.time()

        # Determine which grid points are inside the polygons.
        points_inside_polygons = np.full(len(self.contouring_points), False)
        for contouring_point_index in range(len(self.contouring_points)):
            # If the current point is inside any polygon then mark it as such.
            if polygons_containing_points[contouring_point_index] is not None:
                points_inside_polygons[contouring_point_index] = True

        # time3 = time.time()
        # print('  _find_grid_points_inside_continent_polygons({}, {}): {:.2f} {:.2f}'.format(len(self.contouring_points), len(polygons), time2 - time1, time3 - time2))

        # Reshape 1D array as 2D array indexed by (latitude, longitude) - same order as the points.
        return points_inside_polygons.reshape(
            (self.contouring_grid_num_latitudes, self.contouring_grid_num_longitudes)
        )

    def _find_grid_points_inside_or_near_continent_polygons(self, continent_polygons):
        """
        Find the latitude/longitude grid points that are inside or near the specified continental polygons.

        Note that each continental polygon can have a different buffer/grap distance
        (affecting which points are near each polygon).

        Note that each continent polygon should be an n-tuple where the first element is a pygplates.PolygonOnSphere and
        the second is its buffer/gap distance (in radians).
        Subsequent tuple elements (beyond the two) are optional and ignored.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        """

        points_inside_or_near_polygons = np.full(len(self.contouring_points), False)

        # time1 = time.time()

        # Find the polygon (if any) containing each grid point.
        polygons = [continent_polygon[0] for continent_polygon in continent_polygons]
        polygons_containing_points = (
            points_in_polygons.find_polygons_using_points_spatial_tree(
                self.contouring_points, self.contouring_points_spatial_tree, polygons
            )
        )

        # Determine which grid points are inside the polygons.
        for contouring_point_index in range(len(self.contouring_points)):
            # If the current point is inside any polygon then mark it as such.
            if polygons_containing_points[contouring_point_index] is not None:
                points_inside_or_near_polygons[contouring_point_index] = True

        # time2 = time.time()

        # Group together polygons with the same polygon buffer distance.
        polygon_groups = {}
        for polygon, polygon_buffer_distance_radians, *_ in continent_polygons:
            # We can ignore polygons with a zero buffer distance (they don't get expanded).
            if polygon_buffer_distance_radians > 0:
                if polygon_buffer_distance_radians not in polygon_groups:
                    polygon_groups[polygon_buffer_distance_radians] = []
                polygon_groups[polygon_buffer_distance_radians].append(polygon)

        # Find nearest points to each group of polygons
        # (with all polygons in a group having the same polygon buffer distance).
        for (
            polygon_buffer_distance_radians,
            polygons_in_group,
        ) in polygon_groups.items():

            # The distance threshold is the polygon buffer distance clamped to a maximum of PI.
            distance_threshold_radians = min(math.pi, polygon_buffer_distance_radians)

            # Find the polygons in the current group (if any) near each point.
            points_near_polygons_in_group = proximity_query.find_closest_geometries_to_points_using_points_spatial_tree(
                self.contouring_points,
                self.contouring_points_spatial_tree,
                polygons_in_group,
                distance_threshold_radians=distance_threshold_radians,
            )

            for contouring_point_index in range(len(self.contouring_points)):
                if points_near_polygons_in_group[contouring_point_index] is not None:
                    points_inside_or_near_polygons[contouring_point_index] = True

        # time3 = time.time()
        # print('  _find_grid_points_inside_or_near_continent_polygons({}, {}): {:.2f} {:.2f}'.format(
        #    len(self.contouring_points), len(continent_polygons), time2 - time1, time3 - time2))

        # Reshape 1D array as 2D array indexed by (latitude, longitude) - same order as the points.
        return points_inside_or_near_polygons.reshape(
            (self.contouring_grid_num_latitudes, self.contouring_grid_num_longitudes)
        )

    def _find_grid_points_inside_contoured_continents(self, contoured_continents):
        """
        Find the latitude/longitude grid points that are inside (one or more of) the specified contoured continents.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        """

        # time1 = time.time()

        # Test all grid points against each contoured continent.
        #
        # Note that the original point-in-polygon boolean grid mask (calculated from the input continent polygons before contouring)
        # may not match the result of point-in-contoured-continent tests (done here) since any contoured continents with area below
        # the area threshold would have been removed. So we do our own point-in-contoured-continent tests here.
        points_inside_any_contoured_continent = np.full(
            len(self.contouring_points), False
        )
        for contoured_continent in contoured_continents:
            points_inside_contoured_continent = contoured_continent.are_points_inside(
                self.contouring_points, self.contouring_points_spatial_tree
            )

            # Combine the results of current contoured continent with previous contoured continents.
            #
            # Note that there is typically only a handful of contoured continents in general, so this should not be a bottleneck.
            points_inside_any_contoured_continent[points_inside_contoured_continent] = (
                True
            )

        # time2 = time.time()
        # print('  _find_grid_points_inside_contoured_continents({}): {:.2f}'.format(len(contoured_continents), time2 - time1))

        # Reshape 1D array as 2D array indexed by (latitude, longitude) - same order as the points.
        return points_inside_any_contoured_continent.reshape(
            (self.contouring_grid_num_latitudes, self.contouring_grid_num_longitudes)
        )

    def _find_grid_points_near_merged_continent(self, merged_continent):
        """
        Find the latitude/longitude grid points that are near the polygons of the continents in the specified merged continent.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        """

        # time1 = time.time()

        # Determine which grid points are near the polygons of the continents in the merged continent.
        points_near_merged_continent = np.full(len(self.contouring_points), False)
        for continent in merged_continent.continents:

            # The distance threshold is the continent buffer/gap distance clamped to a maximum of PI.
            distance_threshold_radians = min(
                math.pi, continent.contouring_buffer_and_gap_distance_radians
            )
            if distance_threshold_radians > 0:

                # Each continent polygon is actually a tuple where the first element is a pygplates.PolygonOnSphere.
                #
                # Note: The second element of tuple is the buffer distance of each continent polygon.
                #       But since we're in this function then that should be zero
                #       (ie, we're only using buffer distances for *contoured* continents, not individual polygons).
                polygons = [
                    continent_polygon[0]
                    for continent_polygon in continent.continent_polygons
                ]

                # Find the polygons (if any) near each point.
                points_near_continent = proximity_query.find_closest_geometries_to_points_using_points_spatial_tree(
                    self.contouring_points,
                    self.contouring_points_spatial_tree,
                    polygons,
                    distance_threshold_radians=distance_threshold_radians,
                )

                for contouring_point_index in range(len(self.contouring_points)):
                    if points_near_continent[contouring_point_index] is not None:
                        points_near_merged_continent[contouring_point_index] = True

        # time2 = time.time()
        # print('  _find_grid_points_near_merged_continent({}): {:.2f}'.format(len(merged_continent.continents), time2 - time1))

        # Reshape 1D array as 2D array indexed by (latitude, longitude) - same order as the points.
        return points_near_merged_continent.reshape(
            (self.contouring_grid_num_latitudes, self.contouring_grid_num_longitudes)
        )

    def _create_contoured_continent(self, points_inside_continent):
        """
        Create landmasses and their boundaries from the specified mask of grid points.

        Returns a single 'ContouredContinent' containing one or more landmasses and their contoured boundaries,
        or None if none of the grid points are inside continent.
        """

        # time1 = time.time()

        num_latitudes = self.contouring_grid_num_latitudes
        num_longitudes = self.contouring_grid_num_longitudes

        num_latitude_intervals = num_latitudes - 1
        num_longitude_intervals = num_longitudes - 1

        #
        # Use the Marching Squares algorithm.
        #
        # This is a 2D version (surface of the globe) of the 3D Marching Cubes algorithm.
        # However the main difference between this and using skimage.measure.find_contours(),
        # that we used previously and that also uses the Marching Squares algorithm, is we wrap across the
        # dateline and handle the poles. In this way we avoid contour polygons clamped to the dateline.
        #
        # The way we handle wrapping around the dateline is to have grid points on the dateline (ie, at both longitude -180 and 180).
        # This way lat/lon points on the left side of uniform lat/lon grid of points actually map to the same points on the globe
        # as the lat/lon points on the right side of the uniform lat/lon grid of points, and so they will generated the same
        # point-in-continent-polygon and point-near-continent-polygon results. This ensures the Marching Squares algorithm (below)
        # will produce continuous contour segments across the dateline (as we move from a square on one side of the dateline to the
        # adjacent square on the other side).
        #
        # We also handle the poles correctly by having the bottom row of lat/lon points map to the South pole and the top row
        # map to the North pole. Because all points in a (top or bottom) row map to the same point (pole) on the globe they will
        # generate the same point-in-continent-polygon and point-near-continent-polygon results. And because the entire row is either
        # inside or outside a contour the Marching Squares algorithm (below) cannot generate contour segments that penetrate the row.
        # This essentially avoids the problem at the poles.
        #

        #
        # First find those latitude/longitude squares (each square has 4 points from uniform lat/lon grid of points)
        # that have some of its 4 points inside continent and some outside. These are squares that will contain an edge (segment)
        # of a contour polygon. According to the Marching Squares algorithm there are 16 cases. Two of these have all 4 points
        # either inside or outside (and hence have no segments). Twelve cases have a single segment.
        # And two cases have two segments (because two diagonals points are inside and the other two diagonal points are outside).
        # Here we can choose to either join two separated contour islands or keep them separated.
        # We choose to join them because it makes the algorithm easier - if they were separated then a single square would contain
        # two contours belonging to two separate continents and we'd have to be careful that we visited only the contour belonging
        # to the continent we are currently filling. When we join them then the two contours belong to the same continent.
        #
        # Each segment starts at the middle of one side of the square and ends at the middle of another side.
        # Each side of the square is given an identifier...
        #
        #    ---2---
        #   |       |
        #   1       3
        #   |       |
        #    ---0---
        #
        # ...and each segment records a start and end identifier as a 2-tuple.
        #

        # Records the segments contained by all squares.
        marching_squares = [(None, None)] * (
            num_latitude_intervals * num_longitude_intervals
        )  # squares contain no segments by default
        # Records the lat/lon index of only those squares containing one (or two) segments.
        marching_squares_containing_segments = set()

        # Mark those points *inside* any landmass as requiring a visit.
        # We also need to remove them once they've been visited.
        #
        # Optimisation: The scattered code below (referencing this variable) is optimised for speed since it is a hotspot in the contouring algorithm.
        #               It's a replacement for the following simpler (but slower) code:
        #
        #               for latitude_index in range(num_latitudes):
        #                   for longitude_index in range(num_longitudes):
        #                       if points_inside_continent[latitude_index, longitude_index]:
        #                           points_inside_all_landmasses_to_visit.add((latitude_index, longitude_index))
        points_inside_all_landmasses_to_visit = set()

        bottom_squares_inside_continent = points_inside_continent[
            0, :
        ]  # whether first row of latitude points are inside continent
        for latitude_index in range(num_latitude_intervals):
            top_squares_inside_continent = points_inside_continent[
                latitude_index + 1, :
            ]  # whether next row of latitude points are inside continent

            # See if points on the left of the first square (at the current latitude) are inside continent.
            bottom_left_square_inside_continent = bottom_squares_inside_continent[0]
            top_left_square_inside_continent = top_squares_inside_continent[0]
            for longitude_index in range(num_longitude_intervals):

                # See if points on the right of the current square are inside continent.
                #
                # Note: These will become the left points in the next loop iteration.
                bottom_right_square_inside_continent = bottom_squares_inside_continent[
                    longitude_index + 1
                ]
                top_right_square_inside_continent = top_squares_inside_continent[
                    longitude_index + 1
                ]

                # Mark the point in the bottom-left of square (and that is *inside* any landmass) as requiring a visit.
                if bottom_left_square_inside_continent:
                    points_inside_all_landmasses_to_visit.add(
                        (latitude_index, longitude_index)
                    )

                # Handle the 16 cases of segments in a square.
                #
                # Store 2 segments (most of the time only 1 is needed).
                # Each segment stores a segment start and end identifier.
                if bottom_left_square_inside_continent:
                    if bottom_right_square_inside_continent:
                        if top_left_square_inside_continent:
                            if top_right_square_inside_continent:
                                # Current square doesn't contain a segment, so skip to the next square.
                                segments_in_square = None
                            else:
                                segments_in_square = (2, 3), None
                        else:
                            if top_right_square_inside_continent:
                                segments_in_square = (1, 2), None
                            else:
                                segments_in_square = (1, 3), None
                    else:
                        if top_left_square_inside_continent:
                            if top_right_square_inside_continent:
                                segments_in_square = (0, 3), None
                            else:
                                segments_in_square = (0, 2), None
                        else:
                            if top_right_square_inside_continent:
                                # Choose 2 segments that *do* join two islands.
                                segments_in_square = (0, 3), (1, 2)
                            else:
                                segments_in_square = (0, 1), None
                else:
                    if bottom_right_square_inside_continent:
                        if top_left_square_inside_continent:
                            if top_right_square_inside_continent:
                                segments_in_square = (0, 1), None
                            else:
                                # Choose 2 segments that *do* join two islands.
                                segments_in_square = (0, 1), (2, 3)
                        else:
                            if top_right_square_inside_continent:
                                segments_in_square = (0, 2), None
                            else:
                                segments_in_square = (0, 3), None
                    else:
                        if top_left_square_inside_continent:
                            if top_right_square_inside_continent:
                                segments_in_square = (1, 3), None
                            else:
                                segments_in_square = (1, 2), None
                        else:
                            if top_right_square_inside_continent:
                                segments_in_square = (2, 3), None
                            else:
                                # Current square doesn't contain a segment, so skip to the next square.
                                segments_in_square = None

                # If the current square contains a segment then record that.
                if segments_in_square is not None:
                    marching_squares[
                        latitude_index * num_longitude_intervals + longitude_index
                    ] = segments_in_square
                    marching_squares_containing_segments.add(
                        (latitude_index, longitude_index)
                    )

                # The next square is to the right, so the right side of current square becomes left side of next square.
                #
                # This is an optimisation since this loop is a hotspot in the contouring algorithm.
                # It halves the number of point-inside-contour lookups we need to do.
                bottom_left_square_inside_continent = (
                    bottom_right_square_inside_continent
                )
                top_left_square_inside_continent = top_right_square_inside_continent

            # Mark the point in the last longitude column (and that is *inside* any landmass) as requiring a visit.
            if bottom_right_square_inside_continent:
                points_inside_all_landmasses_to_visit.add(
                    (latitude_index, num_longitudes - 1)
                )

            # In the next loop iteration the bottom row of latitude points (inside continent) will be the current top row.
            bottom_squares_inside_continent = top_squares_inside_continent

        # Mark points in the last latitude row (and that are *inside* any landmass) as requiring a visit.
        for longitude_index in range(num_longitudes):
            if top_squares_inside_continent[longitude_index]:
                points_inside_all_landmasses_to_visit.add(
                    (num_latitudes - 1, longitude_index)
                )

        # time2 = time.time()

        # Return early if none of the grid points are inside continent.
        if not points_inside_all_landmasses_to_visit:
            return None

        #
        # Generate the sole contoured continent by adding one or more landmasses to it.
        #
        # Each landmass is found by picking an arbitrary point inside any contours and expanding around it until we've filled
        # the entire landmass. As we expand we detect when we reach a contour that has not yet been generated and generate it
        # for the current landmass. This expanding fill can detect more than one contour per landmass.
        #
        # This is repeated to find all landmasses (at which time we will have no more points left to visit inside contours).
        #
        contoured_continent = ContouredContinent()
        # Keep visting points *inside* any landmass until there are no points left to visit.
        while points_inside_all_landmasses_to_visit:
            # Contours of the current landmass.
            landmass_contours = []

            # Keep a queue of points inside the current landmass that we will search for contours.
            points_inside_landmass = deque()

            # Get any available point inside any landmass.
            lat_lon_indices_of_first_point_inside_landmass = (
                points_inside_all_landmasses_to_visit.pop()
            )
            point_index_of_first_point_inside_landmass = (
                lat_lon_indices_of_first_point_inside_landmass[0] * num_longitudes
                + lat_lon_indices_of_first_point_inside_landmass[1]
            )
            first_point_inside_landmass = self.contouring_points[
                point_index_of_first_point_inside_landmass
            ]
            # This will be the first point inside the current landmass.
            points_inside_landmass.append(
                lat_lon_indices_of_first_point_inside_landmass
            )

            # Find the remaining points inside the current landmass by recursively searching
            # nearbouring points until we reach a contour boundary of the current landmass.
            while points_inside_landmass:
                # Pop the current point to visit.
                latitude_index, longitude_index = points_inside_landmass.popleft()

                # Search the four squares, adjacent to the current point, for a contour.
                #
                # Note that, for an adjacent square containing a contour, we might already have generated
                # the contour in which case all segments of that contour will have been removed from
                # 'marching_squares_containing_segments' and hence we will be essentially searching for the
                # next contour (if any) of the current landmass (eg, an interior hole contour).
                # And note that all contours in the four adjacent squares belong to the current landmass because
                # we join continent islands (as opposed to separating them) as described above.
                #
                #  +--+--+
                #  |  |  |
                #  +--o--+
                #  |  |  |
                #  +--+--+
                #

                if latitude_index > 0:
                    if longitude_index > 0:
                        neighbour_square_location = (
                            latitude_index - 1,
                            longitude_index - 1,
                        )
                        if (
                            neighbour_square_location
                            in marching_squares_containing_segments
                        ):
                            landmass_contours.append(
                                self._extract_contour(
                                    neighbour_square_location,
                                    marching_squares,
                                    marching_squares_containing_segments,
                                )
                            )

                    if longitude_index < num_longitudes - 1:
                        neighbour_square_location = latitude_index - 1, longitude_index
                        if (
                            neighbour_square_location
                            in marching_squares_containing_segments
                        ):
                            landmass_contours.append(
                                self._extract_contour(
                                    neighbour_square_location,
                                    marching_squares,
                                    marching_squares_containing_segments,
                                )
                            )

                if latitude_index < num_latitude_intervals - 1:
                    if longitude_index > 0:
                        neighbour_square_location = latitude_index, longitude_index - 1
                        if (
                            neighbour_square_location
                            in marching_squares_containing_segments
                        ):
                            landmass_contours.append(
                                self._extract_contour(
                                    neighbour_square_location,
                                    marching_squares,
                                    marching_squares_containing_segments,
                                )
                            )

                    if longitude_index < num_longitudes - 1:
                        neighbour_square_location = latitude_index, longitude_index
                        if (
                            neighbour_square_location
                            in marching_squares_containing_segments
                        ):
                            landmass_contours.append(
                                self._extract_contour(
                                    neighbour_square_location,
                                    marching_squares,
                                    marching_squares_containing_segments,
                                )
                            )

                #
                # Propagate outwards from current point to progressively fill the inside of the current landmass.
                #
                # This requires visiting up to 8 neighbour points (the '+' in diagram below).
                # Only visit those points that are inside (the contour) and that have not yet been visited.
                #
                #  +--+--+
                #  |  |  |
                #  +--o--+
                #  |  |  |
                #  +--+--+
                #
                # Note that we need to wrap around the dateline (longitude) because we need to visit (and remove)
                # ALL points that are inside the *current* landmass (before we move onto the next landmass).
                # However we don't need to traverse beyond the poles (latitude) in the same way.
                #

                if latitude_index > 0:
                    neighbour_point_location = latitude_index - 1, longitude_index
                    if (
                        neighbour_point_location
                        in points_inside_all_landmasses_to_visit
                    ):
                        points_inside_landmass.append(neighbour_point_location)
                        points_inside_all_landmasses_to_visit.remove(
                            neighbour_point_location
                        )

                    if longitude_index > 0:
                        neighbour_point_location = (
                            latitude_index - 1,
                            longitude_index - 1,
                        )
                        if (
                            neighbour_point_location
                            in points_inside_all_landmasses_to_visit
                        ):
                            points_inside_landmass.append(neighbour_point_location)
                            points_inside_all_landmasses_to_visit.remove(
                                neighbour_point_location
                            )
                    else:
                        # Wrap around the dateline.
                        neighbour_point_location = (
                            latitude_index - 1,
                            num_longitudes - 1,
                        )
                        if (
                            neighbour_point_location
                            in points_inside_all_landmasses_to_visit
                        ):
                            points_inside_landmass.append(neighbour_point_location)
                            points_inside_all_landmasses_to_visit.remove(
                                neighbour_point_location
                            )

                    if longitude_index < num_longitudes - 1:
                        neighbour_point_location = (
                            latitude_index - 1,
                            longitude_index + 1,
                        )
                        if (
                            neighbour_point_location
                            in points_inside_all_landmasses_to_visit
                        ):
                            points_inside_landmass.append(neighbour_point_location)
                            points_inside_all_landmasses_to_visit.remove(
                                neighbour_point_location
                            )
                    else:
                        # Wrap around the dateline.
                        neighbour_point_location = latitude_index - 1, 0
                        if (
                            neighbour_point_location
                            in points_inside_all_landmasses_to_visit
                        ):
                            points_inside_landmass.append(neighbour_point_location)
                            points_inside_all_landmasses_to_visit.remove(
                                neighbour_point_location
                            )

                if latitude_index < num_latitudes - 1:
                    neighbour_point_location = latitude_index + 1, longitude_index
                    if (
                        neighbour_point_location
                        in points_inside_all_landmasses_to_visit
                    ):
                        points_inside_landmass.append(neighbour_point_location)
                        points_inside_all_landmasses_to_visit.remove(
                            neighbour_point_location
                        )

                    if longitude_index > 0:
                        neighbour_point_location = (
                            latitude_index + 1,
                            longitude_index - 1,
                        )
                        if (
                            neighbour_point_location
                            in points_inside_all_landmasses_to_visit
                        ):
                            points_inside_landmass.append(neighbour_point_location)
                            points_inside_all_landmasses_to_visit.remove(
                                neighbour_point_location
                            )
                    else:
                        # Wrap around the dateline.
                        neighbour_point_location = (
                            latitude_index + 1,
                            num_longitudes - 1,
                        )
                        if (
                            neighbour_point_location
                            in points_inside_all_landmasses_to_visit
                        ):
                            points_inside_landmass.append(neighbour_point_location)
                            points_inside_all_landmasses_to_visit.remove(
                                neighbour_point_location
                            )

                    if longitude_index < num_longitudes - 1:
                        neighbour_point_location = (
                            latitude_index + 1,
                            longitude_index + 1,
                        )
                        if (
                            neighbour_point_location
                            in points_inside_all_landmasses_to_visit
                        ):
                            points_inside_landmass.append(neighbour_point_location)
                            points_inside_all_landmasses_to_visit.remove(
                                neighbour_point_location
                            )
                    else:
                        # Wrap around the dateline.
                        neighbour_point_location = latitude_index + 1, 0
                        if (
                            neighbour_point_location
                            in points_inside_all_landmasses_to_visit
                        ):
                            points_inside_landmass.append(neighbour_point_location)
                            points_inside_all_landmasses_to_visit.remove(
                                neighbour_point_location
                            )

                if longitude_index > 0:
                    neighbour_point_location = latitude_index, longitude_index - 1
                    if (
                        neighbour_point_location
                        in points_inside_all_landmasses_to_visit
                    ):
                        points_inside_landmass.append(neighbour_point_location)
                        points_inside_all_landmasses_to_visit.remove(
                            neighbour_point_location
                        )
                else:
                    # Wrap around the dateline.
                    neighbour_point_location = latitude_index, num_longitudes - 1
                    if (
                        neighbour_point_location
                        in points_inside_all_landmasses_to_visit
                    ):
                        points_inside_landmass.append(neighbour_point_location)
                        points_inside_all_landmasses_to_visit.remove(
                            neighbour_point_location
                        )

                if longitude_index < num_longitudes - 1:
                    neighbour_point_location = latitude_index, longitude_index + 1
                    if (
                        neighbour_point_location
                        in points_inside_all_landmasses_to_visit
                    ):
                        points_inside_landmass.append(neighbour_point_location)
                        points_inside_all_landmasses_to_visit.remove(
                            neighbour_point_location
                        )
                else:
                    # Wrap around the dateline.
                    neighbour_point_location = latitude_index, 0
                    if (
                        neighbour_point_location
                        in points_inside_all_landmasses_to_visit
                    ):
                        points_inside_landmass.append(neighbour_point_location)
                        points_inside_all_landmasses_to_visit.remove(
                            neighbour_point_location
                        )

            # The current landmass should have encountered one or more contours since
            # it was filled until it reached a boundary (contour) between continent and ocean.
            if not landmass_contours:
                # However it is potentially possible for the landmass to cover the entire globe (ie, no contours).
                # In this case there must be only one landmass, which means all points must have been visited.
                if points_inside_all_landmasses_to_visit:
                    raise AssertionError(
                        "A single landmass covering entire globe must be the only landmass"
                    )

            # Add the current landmass (bounded by contours) to the ContouredContinent.
            # This uses an arbitrary point inside the landmass to determine which contours are its exterior rings and which are its interior rings.
            self._add_landmass_to_contoured_continent(
                contoured_continent, landmass_contours, first_point_inside_landmass
            )

        # time3 = time.time()
        # print('  _create_contoured_continent: {:.2f} {:.2f}'.format(time2 - time1, time3 - time2))

        return contoured_continent

    def _add_landmass_to_contoured_continent(
        self,
        contoured_continent,
        landmass_contours,
        any_point_inside_contoured_continent,
    ):
        """
        Add a landmass (bounded by the specified contours) to a ContouredContinent.
        """

        # For each contour create a ring (a polygon with only an exterior ring).
        contour_rings = [
            pygplates.PolygonOnSphere(contour) for contour in landmass_contours
        ]

        # Arrange the contour rings into those that *include* and those that *exclude* continent.
        #
        # Note that we can get multiple ocean polygons if there is no contour ring that includes continent.
        # In this case the continent covers the entire globe except for a few oceanic holes.
        # However you can't have a single polygon with only interior holes (and no exterior ring).
        # So intead we allow for multiple ocean polygons to represent these holes (and treat this as a special case).
        contour_rings_including_continent = []
        contour_rings_excluding_continent = []
        for contour_ring in contour_rings:
            # A point inside the continent might actually be outside the current contour ring.
            contour_ring_interior_contains_continent = contour_ring.is_point_in_polygon(
                any_point_inside_contoured_continent
            )
            if contour_ring_interior_contains_continent:
                contour_rings_including_continent.append(contour_ring)
            else:
                contour_rings_excluding_continent.append(contour_ring)

        # We should have either:
        #
        # 1) a single contour ring *including* continent and *zero or more* contour rings *excluding* continent, or
        # 2) no contour ring *including* continent and *one or more* contour rings *excluding* continent.
        #
        # For the case (1) we have a *single* polygon with zero or more interior rings. And it *includes* continent.
        # For the case (2) we have a *multiple* polygons, each with no interior rings. And they *exclude* continent with
        # the rest of the globe *including* continent.
        #
        # There's actually a third special case when there's no oceanic crust, just a single landmass covering the entire globe.
        # In that case there are no contours (and no continent or ocean polygons).
        if len(contour_rings_including_continent) == 1:
            # One continent polygon.
            continent_exterior_ring = contour_rings_including_continent[0]
            ocean_interior_rings = contour_rings_excluding_continent
            continent_polygon = pygplates.PolygonOnSphere(
                continent_exterior_ring, ocean_interior_rings
            )
            # Add the continent polygon.
            contoured_continent._add_continent(continent_polygon)
        elif len(contour_rings_including_continent) == 0:
            # Zero or more ocean polygons (each with no interior rings).
            # Note: We only have *zero* ocean polygons for the special case of a single landmass covering the entire globe.
            for ocean_exterior_ring in contour_rings_excluding_continent:
                ocean_polygon = pygplates.PolygonOnSphere(ocean_exterior_ring)
                # Add an ocean polygon.
                contoured_continent._add_ocean(ocean_polygon)
        else:  # len(contour_rings_including_continent) >= 2
            raise AssertionError(
                "A single landmass cannot have multiple polygons that include continent"
            )

    def _extract_contour(
        self,
        first_segment_lat_lon_indices,
        marching_squares,
        marching_squares_containing_segments,
    ):
        """
        Follow the segment in marching square at specified lat/lon index around contour back to that first segment.

        Note that a marching square can have two segments, in which case that square represents a thin connecting region between
        two larger islands of the contour (but still all just one contour). That square will get traversed twice (once in one direction
        through one segment and once in another the opposite direction through the second segment of that square).
        """

        contour_points = []

        interval_spacing = self.contouring_point_spacing_degrees
        num_latitude_intervals = self.contouring_grid_num_latitudes - 1
        num_longitude_intervals = self.contouring_grid_num_longitudes - 1

        #
        # When a contour is first encountered (during the caller's expanding fill) a contour ring is generated by starting at the first square
        # found that contains one (or two) segments, which represents the start of that contour. We then pick one of that square's segments
        # (in most cases it'll only have one segment) and generate the first contour point at that segment's start. Note that it doesn't matter
        # which segment we pick (if there's two segments) because the contour ring will traverse back to the second segment (since both segments
        # are part of the same contour because their containing square represents a thin connecting region between two larger islands of the ring).
        # We then find the adjacent square to the segment's end (since a segment ends in the middle of a side of the square we can find the adjacent square).
        # We then find the segment in the adjacent square that starts (or ends) at the that point (the previous segment end). The adjacent
        # square may contain two segments in which case we need to find the correct segment (that continues the previous segment).
        # We generate the next contour point at the segment start and continue this process following the contour through segments of squares
        # until we return to the first segment (thus closing the contour loop).
        #

        #
        # Starting at the first segment, follow the segments in a loop until they return to the first segment (thus forming a contour ring).
        #
        latitude_index, longitude_index = first_segment_lat_lon_indices
        prev_segment_end = None
        while True:

            # Get a segment from the current square.
            segment1, segment2 = marching_squares[
                latitude_index * num_longitude_intervals + longitude_index
            ]
            # If a square has only one segment then it will be in 'segment1' (not 'segment2').
            if segment1 is None:
                # Shouldn't be able to reach a square that doesn't have any segments (or previously had segments but now has none).
                raise AssertionError("Square has no segments")

            # Find a segment in the current square such that the segment start matches the
            # end of the previous segment (in the previous square).
            segment_start, segment_end = segment1
            if prev_segment_end is None:  # first segment of current contour...
                # Mark the start of the contour so that later we know when we've completed the contour.
                first_segment_start = segment_start
            else:
                # Continue from the side of the previous square that contains the end point of the previous segment to the side of
                # the current square that should contain the start point of the current segment (that continues previous segment).
                #
                # The right side of previous square continues to left side of current square.
                # The left side of previous square continues to right side of current square.
                # The top side of previous square continues to bottom side of current square.
                # The bottom side of previous square continues to top side of current square.
                #
                # The adjacency relation means 2->0, 0->2, 3->1 and 1->3...
                #
                #    ---2---
                #   |       |
                #   1       3
                #   |       |
                #    ---0---
                #
                # ...which is satisfied by 2^2->0, 0^2->2, 3^2->1 and 1^2->3 (where '^' is exclusive-or).
                #
                curr_segment_start = prev_segment_end ^ 0b10

                #
                # Find the right segment (if there's two segments) and reverse the segment if necessary
                # so that previous segment end matches current segment start.
                #
                if curr_segment_start == segment_start:
                    # We're traversing segment in the correct direction.
                    pass
                elif curr_segment_start == segment_end:
                    # Reverse segment direction (swap segment start and end points).
                    segment_start, segment_end = segment_end, segment_start
                else:
                    if segment2:
                        # Segment 1 didn't match so swap it with segment 2 (so it can be used later).
                        segment1, segment2 = segment2, segment1

                        segment_start, segment_end = segment1
                        if curr_segment_start == segment_start:
                            # We're traversing segment in the correct direction.
                            pass
                        elif curr_segment_start == segment_end:
                            # Reverse segment direction (swap segment start and end points).
                            segment_start, segment_end = segment_end, segment_start
                        else:
                            raise AssertionError("Unable to find connecting segment")
                    else:
                        raise AssertionError("Unable to find connecting segment")

            # The start position of 'segment'.
            # It will be at the midpoint of a side of the square.
            if segment_start == 0:
                segment_start_latitude = -90.0 + latitude_index * interval_spacing
                segment_start_longitude = (
                    -180.0 + (longitude_index + 0.5) * interval_spacing
                )
            elif segment_start == 1:
                segment_start_latitude = (
                    -90.0 + (latitude_index + 0.5) * interval_spacing
                )
                segment_start_longitude = -180.0 + longitude_index * interval_spacing
            elif segment_start == 2:
                segment_start_latitude = -90.0 + (latitude_index + 1) * interval_spacing
                segment_start_longitude = (
                    -180.0 + (longitude_index + 0.5) * interval_spacing
                )
            else:  # segment_start == 3
                segment_start_latitude = (
                    -90.0 + (latitude_index + 0.5) * interval_spacing
                )
                segment_start_longitude = (
                    -180.0 + (longitude_index + 1) * interval_spacing
                )

            # Generate a contour point at the start of the current segment.
            contour_point = pygplates.PointOnSphere(
                segment_start_latitude, segment_start_longitude
            )
            contour_points.append(contour_point)

            # We've just used 'segment1', so discard it by moving 'segment2' into its position to be used later.
            # And if 'segment2' is None then there are no more segments in current square so discard the entire square.
            marching_squares[
                latitude_index * num_longitude_intervals + longitude_index
            ] = (segment2, None)
            if segment2 is None:
                # There are no segments left in the current square, so we're finished with it.
                # Note: This will raise KeyError if not present in 'set'.
                marching_squares_containing_segments.remove(
                    (latitude_index, longitude_index)
                )

            # We're moving onto the next segment in the next square.
            prev_segment_end = segment_end

            # Move to the next square connected by the end of the current segment.
            #
            #    ---2---
            #   |       |
            #   1       3
            #   |       |
            #    ---0---
            #
            # As noted above, at each pole there is an entire row of lat/lon grid points that are all either inside or outside a contour.
            # This means the Marching Squares algorithm cannot generate contour segments that penetrate the row. So we should not be able
            # to move beyond the poles.
            #
            # Also as noted above, the both the leftmost and rightmost columns of the lat/lon grid of points will be on the dateline
            # (ie, at both longitude -180 and 180). This means the Marching Squares algorithm will produce continuous contour segments across
            # the dateline (as we move from a square on one side of the dateline to the adjacent square on the other side).
            if prev_segment_end == 0:
                latitude_index -= 1
                if latitude_index < 0:
                    raise AssertionError("Segment entered South Pole")
            elif prev_segment_end == 1:
                longitude_index -= 1
                if longitude_index < 0:
                    # Wrap around the dateline.
                    longitude_index += num_longitude_intervals
            elif prev_segment_end == 2:
                latitude_index += 1
                if latitude_index == num_latitude_intervals:
                    raise AssertionError("Segment entered North Pole")
            else:  # prev_segment_end == 3
                longitude_index += 1
                if longitude_index == num_longitude_intervals:
                    # Wrap around the dateline.
                    longitude_index -= num_longitude_intervals

            # See if we're returned to the first square (containing the first segment).
            if first_segment_lat_lon_indices == (latitude_index, longitude_index):
                # And make sure the end of the previous segment matches the start of the first segment.
                # See comment above about adjacency relation for explanatation of exclusive-or.
                if first_segment_start == (prev_segment_end ^ 0b10):
                    # Break out of current contour loop (we've completed the contour).
                    break

        # Return the ring of contour points.
        return contour_points
