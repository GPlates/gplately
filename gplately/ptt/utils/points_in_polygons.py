#
#    Copyright (C) 2017-2020 The University of Sydney, Australia
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
Efficient point-in-polygon testing when there are many relatively uniformly spaced points to be tested against polygons. 

For example, to find the plate ID of the polygon containing each point in a sequence of points:

```pyton
import points_in_polygons

# A list of 'pygplates.PointOnSphere' points.
points = [...]

# Some polygon features (eg, coastlines).
polygon_feature_collection = pygplates.FeatureCollection('polygons.gpml')

# Extract the polygons from the features.
polygons = []
polygon_features = []
for polygon_feature in polygon_feature_collection:
    polygons.append(polygon_feature.get_geometry())
    polygon_features.append(polygon_feature)

#
# Assuming the polygons are *non-overlapping*, find the single polygon (feature) containing each point.
#
polygon_features_containing_points = points_in_polygons.find_polygons(points, polygons, polygon_features)

# Assign a plate ID to each point (or 0 if point outside all polygons).
plate_ids = []
for polygon_feature in polygon_features_containing_points:
    if polygon_feature is not None:
        plate_id = polygon_feature.get_reconstruction_plate_id()
    else:
        plate_id = 0

    plate_ids.append(plate_id)

#
# Assuming the polygons are *overlapping*, find all polygons (features) containing each point.
#
polygon_features_containing_points = points_in_polygons.find_polygons(points, polygons, polygon_features, all_polygons=True)

# Assign multiple points to each plate ID.
plate_id_to_points_mapping = {} # Each plate ID has a list of points.
for point_index, polygon_feature_list in enumerate(polygon_features_containing_points):
    if polygon_feature_list:
        for polygon_feature in polygon_feature_list
            plate_id = polygon_feature.get_reconstruction_plate_id()
            points_with_plate_id = plate_id_to_points_mapping.setdefault(plate_id, [])
            points_with_plate_id.append(points[point_index])
```

"""


from __future__ import print_function

import math
import sys

import pygplates

from . import points_spatial_tree

# the following line would ensure more correct python2-3 compatibility,
# but requires a non-standard module
# from builtins import range


def find_polygons(
    points,
    polygons,
    polygon_proxies=None,
    all_polygons=False,
    subdivision_depth=points_spatial_tree.DEFAULT_SUBDIVISION_DEPTH,
):
    """
    Efficient point-in-polygon testing when there are many relatively uniformly spaced points to be tested against polygons.

    Parameters
    ----------
    points: a sequence of 'pygplates.PointOnSphere'
        a sequence of points

    polygons: a sequence of 'pygplates.PolygonOnSphere'
        a sequence of polygons

    polygon_proxies : Optional sequence of objects associated with 'polygons'
        If not specified then the proxies default to the polygons themselves.
        These can be any object (such as the 'pygplates.Feature' that the polygon came from).

    all_polygons: bool
        Whether to find all polygons containing each point or just the first one encountered.
        Set to True if polygons overlap each other, otherwise set to False (for non-overlapping polygons).
        Defaults to False (non-overlapping polygons).

    subdivision_depth: number
        The depth of the lat/lon quad tree used to speed up point-in-polygon queries.
        The lat/lon width of a leaf quad tree node is (90 / (2^subdivision_depth)) degrees.
        Generally the denser the 'points' the larger the depth should be.
        Setting this value too high causes unnecessary time to be spent generating a deep quad tree.
        Setting this value too low reduces the culling efficiency of the quad tree.
        However a value of 4 seems to work quite well for a uniform lat/lon spacing of 'points' of 1 degree and below
        without the cost of generating a deep quad tree.
        So most of the time the subdivision depth can be left at its default value.

    Returns
    -------
    A list of polygon proxies associated with 'points'
        The length of the returned list matches the length of 'points'.
        For each point in 'points', if the point is contained by a polygon then that polygon's proxy
        is stored (otherwise None is stored) at the same index (as the point) in the returned list.
        If 'all_polygons' is False then each item in returned list is a single polygon proxy (or a single None).
        If 'all_polygons' is True then each item in returned list is a *list* of polygon proxies (or a single None).

    Raises ValueError if the lengths of 'polygons' and 'polygon_proxies' (if specified) do not match.
    """

    spatial_tree_of_points = points_spatial_tree.PointsSpatialTree(
        points, subdivision_depth
    )
    return find_polygons_using_points_spatial_tree(
        points, spatial_tree_of_points, polygons, polygon_proxies, all_polygons
    )


def find_polygons_using_points_spatial_tree(
    points, spatial_tree_of_points, polygons, polygon_proxies=None, all_polygons=False
):
    """
    Same as 'find_polygons()' except 'spatial_tree_of_points' is a 'points_spatial_tree.PointsSpatialTree' of 'points'.

    This is useful when re-using a single 'points_spatial_tree.PointsSpatialTree'.
    For example, when using it both for point-in-polygon queries and minimum distance queries.

    Note that 'spatial_tree_of_points' should have been built from 'points' since it contains
    indices into the 'points' sequence.
    """

    # Use the polygons as proxies if no proxies have been specified.
    if polygon_proxies is None:
        polygon_proxies = polygons

    if len(polygons) != len(polygon_proxies):
        raise ValueError("Number of polygons must match number of proxies.")

    # Sort the polygons from largest to smallest area.
    # This makes searching for points/geometries more efficient.
    #
    # 'polygons_and_proxies' is a list of 2-tuples (polygon, polygon_proxy).
    polygons_and_proxies = sorted(
        ((polygons[index], polygon_proxies[index]) for index in range(len(polygons))),
        key=lambda polygon_and_proxy: polygon_and_proxy[0].get_area(),
        reverse=True,
    )

    # By default all points are outside all polygons.
    # If any are found to be inside then we'll set the relevant polygon proxy.
    polygon_proxies_containing_points = [None] * len(points)

    # Use a quad tree for efficiency - enables us to cull large groups of points that are either
    # outside all polygons or inside a polygon (avoids point-in-polygon tests for these points).
    for root_node in spatial_tree_of_points.get_root_nodes():
        _visit_spatial_tree_node(
            root_node,
            points,
            polygons_and_proxies,
            polygon_proxies_containing_points,
            all_polygons,
        )

    return polygon_proxies_containing_points


##################
# Implementation #
##################


def _visit_spatial_tree_node(
    node,
    points,
    parent_overlapping_polygons_and_proxies,
    polygon_proxies_containing_points,
    all_polygons,
):

    # See if the current quad tree node's bounding polygon overlaps any polygons.
    overlapping_polygons_and_proxies = []
    for polygon, polygon_proxy in parent_overlapping_polygons_and_proxies:

        # See if quad tree node and current polygon overlap.
        polygon_node_overlap_distance = pygplates.GeometryOnSphere.distance(
            node.get_bounding_polygon(),
            polygon,
            1e-4,  # Anything smaller than this is considered zero distance (intersection).
            geometry1_is_solid=True,
            geometry2_is_solid=True,
        )
        # Note: PyGPlates version 0.15 fixed an issue in GeometryOnSphere.distance() where very small distances (approx <= 1e-6)
        # should have been set to zero (to indicate an intersection) but were not. So we'll assume pyGPlates version 0.14 or less
        # test for small distances instead of zero (this also works for version 0.15 and above). Our small distance is the threshold
        # distance above. It's possible both polygons don't really overlap but are just very close, but that's OK since
        # we're being conservative in that we're looking for possibilities of points being in polygons.
        if polygon_node_overlap_distance is not None:

            # See if quad tree node is contained completely inside polygon.
            # We test this by ensuring the outlines of both polygons are not too close (ie, don't touch), and if they don't
            # touch then later testing if any point on the node bounding polygon outline is inside the current polygon
            # (and none of the current polygon's interior holes are inside the node bounding polygon).
            node_to_polygon_distance = pygplates.GeometryOnSphere.distance(
                node.get_bounding_polygon(), polygon, 1e-4
            )  # Anything smaller than this is considered zero distance (intersection).

            # See if the node bounding polygon is entirely inside the solid interior region of the current polygon.
            is_quad_tree_completely_inside_polygon = False
            if node_to_polygon_distance is None:
                # Test if any point on the node bounding polygon outline (we use first point) is inside the current polygon.
                if polygon.is_point_in_polygon(node.get_bounding_polygon()[0]):
                    is_quad_tree_completely_inside_polygon = True
                    # If the current polygon contains interior rings/holes then it's possible that the outline of the node bounding polygon
                    # lies entirely inside the current polygon *but* surrounds one of its interior rings/holes.
                    # And this means we cannot fill the node (like we could if the current polygon did not contain any holes).
                    # So we also detect if any interior holes of the current polygon are inside the node bounding polygon.
                    for interior_ring_index in range(
                        polygon.get_number_of_interior_rings()
                    ):
                        # Arbitrarily choose the first point on each interior ring.
                        if node.get_bounding_polygon().is_point_in_polygon(
                            polygon.get_interior_ring_points(interior_ring_index)[0]
                        ):
                            is_quad_tree_completely_inside_polygon = False
                            break

            if is_quad_tree_completely_inside_polygon:

                # Recursively fill the entire quad sub-tree as inside current polygon.
                _fill_spatial_tree_node_inside_polygon(
                    node, polygon_proxy, polygon_proxies_containing_points, all_polygons
                )

                if not all_polygons:
                    # Only storing first polygon proxy encountered, so skip remaining polygons.
                    return

                # Note: No need to add polygon to 'overlapping_polygons_and_proxies' since we've already taken care of it.

            else:
                overlapping_polygons_and_proxies.append((polygon, polygon_proxy))

    # If quad tree node is outside all polygons then nothing left to do since all points are marked as outside by default.
    if not overlapping_polygons_and_proxies:
        return

    # Visit child nodes (if internal node) or test each point (if leaf node).
    if node.is_internal_node():
        for child_node in node.get_child_nodes():
            _visit_spatial_tree_node(
                child_node,
                points,
                overlapping_polygons_and_proxies,
                polygon_proxies_containing_points,
                all_polygons,
            )
    else:
        for point_index in node.get_point_indices():
            point = points[point_index]
            for polygon, polygon_proxy in overlapping_polygons_and_proxies:
                if polygon.is_point_in_polygon(point):
                    # Point is inside a polygon.
                    if all_polygons:
                        # Each point has a *list* of polygon proxies (or None).
                        # Create list if first polygon proxy encountered for current point.
                        if polygon_proxies_containing_points[point_index] is None:
                            polygon_proxies_containing_points[point_index] = []
                        polygon_proxies_containing_points[point_index].append(
                            polygon_proxy
                        )
                    else:
                        # Each point has a *single* polygon proxy (or None).
                        polygon_proxies_containing_points[point_index] = polygon_proxy
                        # No need to visit remaining polygons for the current point.
                        break


def _fill_spatial_tree_node_inside_polygon(
    node, polygon_proxy, polygon_proxies_containing_points, all_polygons
):

    if node.is_internal_node():
        for child_node in node.get_child_nodes():
            _fill_spatial_tree_node_inside_polygon(
                child_node,
                polygon_proxy,
                polygon_proxies_containing_points,
                all_polygons,
            )
    else:
        for point_index in node.get_point_indices():
            # Point is inside a polygon.
            if all_polygons:
                # Each point has a *list* of polygon proxies (or None).
                # Create list if first polygon proxy encountered for current point.
                if polygon_proxies_containing_points[point_index] is None:
                    polygon_proxies_containing_points[point_index] = []
                polygon_proxies_containing_points[point_index].append(polygon_proxy)
            else:
                # Each point has a *single* polygon proxy (or None).
                polygon_proxies_containing_points[point_index] = polygon_proxy


# if __name__ == '__main__':
#
#    #
#    # Some testing/example code.
#    #
#
#    import time
#
#
#    print('Loading coastline polygons and rotation model...')
#    coastline_features = pygplates.FeatureCollection('../../../sample_data/2.0/SampleData/FeatureCollections/Coastlines/Matthews_etal_GPC_2016_Coastlines.gpmlz')
#    rotation_model = pygplates.RotationModel('../../../sample_data/2.0/SampleData/FeatureCollections/Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot')
#
#    print('Reconstructing coastline polygons...')
#    reconstruction_time = 200
#    coastline_reconstructed_feature_geometries = []
#    pygplates.reconstruct(coastline_features, rotation_model, coastline_reconstructed_feature_geometries, reconstruction_time)
#
#    polygons = []
#    polygon_features = []
#    for reconstructed_feature_geometry in coastline_reconstructed_feature_geometries:
#        polygons.append(reconstructed_feature_geometry.get_reconstructed_geometry())
#        polygon_features.append(reconstructed_feature_geometry.get_feature())
#
#    # Create uniform lat/lon distribution of points.
#    print('Creating lat/lon grid of points...')
#    num_latitudes = 180
#    num_longitudes = 360
#    lat_grid_spacing_degrees = 180.0 / num_latitudes
#    lon_grid_spacing_degrees = 360.0 / num_longitudes
#
#    points = []
#    for lat_index in xrange(num_latitudes):
#        # The 0.5 puts the point in the centre of the grid pixel.
#        # This also avoids sampling right on the poles.
#        lat = -90 + (lat_index + 0.5) * lat_grid_spacing_degrees
#
#        for lon_index in xrange(num_longitudes):
#            # The 0.5 puts the point in the centre of the grid pixel.
#            # This also avoids sampling right on the dateline where there might be
#            # age grid or static polygon artifacts.
#            lon = -180 + (lon_index + 0.5) * lon_grid_spacing_degrees
#
#            point = pygplates.PointOnSphere(lat, lon)
#            points.append(point)
#
#    print('Finding polygons containing points...')
#    time_begin = time.clock()
#
#    if True:
#        #
#        # The fast way (about 3 seconds).
#        #
#        polygon_features_containing_points = find_polygons(points, polygons, polygon_features, all_polygons=True)
#    else:
#        #
#        # The slow way (about 200 seconds).
#        #
#        # Similar to 'find_polygons()' except without using a quad tree.
#        #
#        polygons_and_features = sorted(
#                ((polygons[index], polygon_features[index]) for index in xrange(len(polygons))),
#                key=lambda polygon_and_feature: polygon_and_feature[0].get_area(),
#                reverse=True)
#        polygon_features_containing_points = [None] * len(points)
#        for point_index, point in enumerate(points):
#            for polygon, polygon_feature in polygons_and_features:
#                if polygon.is_point_in_polygon(point):
#                    if polygon_features_containing_points[point_index] is None:
#                        polygon_features_containing_points[point_index] = []
#                    polygon_features_containing_points[point_index].append(polygon_feature)
#
#    time_end = time.clock()
#    print('  {0} seconds'.format(time_end - time_begin))
#
#    print('Associate each point with zero or more polygon plate IDs...')
#
#    # Group points inside each polygon so can create one multi-point per polygon.
#    polygon_feature_to_points_mapping = {}
#    for point_index, polygon_feature_list in enumerate(polygon_features_containing_points):
#        if polygon_feature_list:
#            for polygon_feature in polygon_feature_list:
#                points_in_polygon = polygon_feature_to_points_mapping.setdefault(polygon_feature, [])
#                points_in_polygon.append(points[point_index])
#
#    # Create multi-point features.
#    multi_point_features = []
#    for polygon_feature, points_in_polygon in polygon_feature_to_points_mapping.iteritems():
#        multi_point_feature = pygplates.Feature()
#        multi_point_feature.set_geometry(
#                pygplates.MultiPointOnSphere(points_in_polygon))
#
#        begin_time, end_time = polygon_feature.get_valid_time()
#        multi_point_feature.set_valid_time(begin_time, end_time)
#
#        multi_point_feature.set_reconstruction_plate_id(
#                polygon_feature.get_reconstruction_plate_id())
#
#        multi_point_features.append(multi_point_feature)
#
#    print('Writing points feature collection...')
#    pygplates.FeatureCollection(multi_point_features).write('multi_point_features.gpml')
