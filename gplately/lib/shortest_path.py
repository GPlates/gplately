#
#    Copyright (C) 2017-2026 The University of Sydney, Australia
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

"""Shortest-distance path, around geometry obstacles, from source geometries to a target.

A port of EarthByte's `predicting-sediment-thickness
<https://github.com/EarthByte/predicting-sediment-thickness>`__ ``shortest_path.py``, used by
:mod:`gplately.grids.sediment_thickness` (`gplately#445
<https://github.com/GPlates/gplately/issues/445>`__) to route distances *around* continents
(rather than a great-circle straight line) when a straight line would otherwise cut through
land -- see that module's ``continent_obstacle_features`` argument, and its docstring for why
this matters.

Obstacles can be both polygons and polylines. Paths cannot go into the interior of polygons; if
a source (or target) geometry is entirely inside an obstacle polygon, the closest points on the
polygon's outline are used instead.

Example
-------
::

    grid = Grid(subdivision_depth=6)  # a spacing of (90 / 2^6) = 1.40625 degrees

    obstacle_grid = grid.create_obstacle_grid(obstacle_geometries)

    # Distances from every non-obstacle grid node to the nearest source geometry,
    # routed around the obstacles.
    distance_grid = obstacle_grid.create_distance_grid(source_geometries)

    distance_radians = distance_grid.shortest_distance(target_geometry)
    if distance_radians is not None:
        distance_km = distance_radians * pygplates.Earth.mean_radius_in_kms
    else:
        # Target is unreachable from every source (completely blocked by obstacles), or
        # exceeded an optional distance threshold.
        ...

Implementation
--------------
A :class:`Grid` is a fixed lon/lat mesh of points (with an associated quad tree, for fast
spatial queries) that doubles as a graph: each grid node is connected to its ~16 nearest
neighbours. An :class:`ObstacleGrid` removes nodes that fall inside an obstacle polygon, and
removes any node-to-neighbour edge that crosses an obstacle outline -- so paths through the
remaining graph can never cross an obstacle. A :class:`DistanceGrid` then runs Dijkstra's
algorithm from a set of source geometries across that graph, giving every remaining grid node
its shortest obstacle-respecting distance to the nearest source; :meth:`DistanceGrid.shortest_distance`
looks up (and locally interpolates) an arbitrary target point's distance from that precomputed
grid.
"""

import heapq
import math

import pygplates

__all__ = ["Grid", "ObstacleGrid", "DistanceGrid"]


class Grid(object):
    # Each quad tree node has:
    #   (1 << GRID_NODE_DEPTH_PER_QUAD_TREE_NODE) * (1 << GRID_NODE_DEPTH_PER_QUAD_TREE_NODE)
    # ...grid points.
    GRID_NODE_DEPTH_PER_QUAD_TREE_NODE = 2  # 16 grid points

    def __init__(self, subdivision_depth):
        if subdivision_depth < 0:
            raise ValueError("Subdivision depth must be a non-negative value.")

        self.num_latitudes = 2 * (1 << subdivision_depth)
        self.num_longitudes = 2 * self.num_latitudes
        self.grid_spacing_degrees = 180.0 / self.num_latitudes
        self.spacing_radians = math.radians(self.grid_spacing_degrees)
        self.subdivision_depth = subdivision_depth

        # In 'get_neighbour_grid_nodes()' we sample 16 nearest neighbours.
        # The 5x5 sampling pattern is such that the furthest neighbour is sqrt(2^2 + 1^2) pixels away.
        # We don't sample corners of 5x5 pattern (that would make it sqrt(2^2 + 2^2)).
        # This max spacing applies at the equator. Near the poles is gets smaller, but it only
        # needs to be conservative (too large is fine, too small is not).
        self.maximum_distance_radians_to_neighbour_grid_node = (
            math.sqrt(5.0) * self.spacing_radians
        )

        self._init_grid()
        self._init_quad_tree()

    def create_obstacle_grid(self, obstacle_geometries):
        return ObstacleGrid(self, obstacle_geometries)

    # Returns a list of 2-tuples (distance_radians, node_index) for up to 4 nearest grid nodes to 'point'.
    def get_nearest_grid_nodes(self, point):
        lat, lon = point.to_lat_lon()

        right_lat_index = int(0.5 + (lat - (-90)) / self.grid_spacing_degrees)
        right_lon_index = int(0.5 + (lon - (-180)) / self.grid_spacing_degrees)

        if right_lon_index == 0:
            left_lon_index = self.num_longitudes - 1
        elif right_lon_index == self.num_longitudes:
            left_lon_index = self.num_longitudes - 1
            right_lon_index = 0
        else:
            left_lon_index = right_lon_index - 1

        if right_lat_index == 0:
            lat_lon_indices = ((0, left_lon_index), (0, right_lon_index))
        elif right_lat_index == self.num_latitudes:
            lat_lon_indices = (
                (self.num_latitudes - 1, left_lon_index),
                (self.num_latitudes - 1, right_lon_index),
            )
        else:
            left_lat_index = right_lat_index - 1
            lat_lon_indices = (
                (left_lat_index, left_lon_index),
                (left_lat_index, right_lon_index),
                (right_lat_index, left_lon_index),
                (right_lat_index, right_lon_index),
            )

        nearest_grid_nodes = []
        for lat_index, lon_index in lat_lon_indices:
            node_index = lat_index * self.num_longitudes + lon_index
            distance_to_node = pygplates.GeometryOnSphere.distance(
                point, self.nodes[node_index].point
            )
            nearest_grid_nodes.append((distance_to_node, node_index))

        return nearest_grid_nodes

    # Returns a list of 2-tuples (distance_radians, node_index) for neighbour grid nodes of 'node_index'.
    def get_neighbour_grid_nodes(self, node_index):
        node = self.nodes[node_index]

        if node._neighbours is None:
            node._neighbours = []

            nodes = self.nodes
            num_latitudes = self.num_latitudes
            num_longitudes = self.num_longitudes

            lat_index, lon_index = divmod(node_index, num_longitudes)

            # Left neighbours.
            if lon_index > 0:
                node._add_neighbour(
                    lat_index * num_longitudes + (lon_index - 1), nodes
                )  # Nearest
                if lat_index > 0:
                    node._add_neighbour(
                        (lat_index - 1) * num_longitudes + (lon_index - 1), nodes
                    )  # Nearest
                    if lon_index > 1:
                        node._add_neighbour(
                            (lat_index - 1) * num_longitudes + (lon_index - 2), nodes
                        )  # Second nearest.
                    else:  # lon_index == 1
                        node._add_neighbour(
                            (lat_index - 1) * num_longitudes + (num_longitudes - 1),
                            nodes,
                        )  # Second nearest.
                    if lat_index > 1:
                        node._add_neighbour(
                            (lat_index - 2) * num_longitudes + (lon_index - 1), nodes
                        )  # Second nearest
                if lat_index < num_latitudes - 1:
                    node._add_neighbour(
                        (lat_index + 1) * num_longitudes + (lon_index - 1), nodes
                    )  # Nearest
                    if lon_index > 1:
                        node._add_neighbour(
                            (lat_index + 1) * num_longitudes + (lon_index - 2), nodes
                        )  # Second nearest.
                    else:  # lon_index == 1
                        node._add_neighbour(
                            (lat_index + 1) * num_longitudes + (num_longitudes - 1),
                            nodes,
                        )  # Second nearest.
                    if lat_index < num_latitudes - 2:
                        node._add_neighbour(
                            (lat_index + 2) * num_longitudes + (lon_index - 1), nodes
                        )  # Second nearest
            else:  # lon_index == 0
                node._add_neighbour(
                    lat_index * num_longitudes + (num_longitudes - 1), nodes
                )  # Nearest
                if lat_index > 0:
                    node._add_neighbour(
                        (lat_index - 1) * num_longitudes + (num_longitudes - 1), nodes
                    )  # Nearest
                    node._add_neighbour(
                        (lat_index - 1) * num_longitudes + (num_longitudes - 2), nodes
                    )  # Second nearest.
                    if lat_index > 1:
                        node._add_neighbour(
                            (lat_index - 2) * num_longitudes + (num_longitudes - 1),
                            nodes,
                        )  # Second nearest
                if lat_index < num_latitudes - 1:
                    node._add_neighbour(
                        (lat_index + 1) * num_longitudes + (num_longitudes - 1), nodes
                    )  # Nearest
                    node._add_neighbour(
                        (lat_index + 1) * num_longitudes + (num_longitudes - 2), nodes
                    )  # Second nearest.
                    if lat_index < num_latitudes - 2:
                        node._add_neighbour(
                            (lat_index + 2) * num_longitudes + (num_longitudes - 1),
                            nodes,
                        )  # Second nearest

            # Right neighbours.
            if lon_index < num_longitudes - 1:
                node._add_neighbour(
                    lat_index * num_longitudes + (lon_index + 1), nodes
                )  # Nearest
                if lat_index > 0:
                    node._add_neighbour(
                        (lat_index - 1) * num_longitudes + (lon_index + 1), nodes
                    )  # Nearest
                    if lon_index < num_longitudes - 2:
                        node._add_neighbour(
                            (lat_index - 1) * num_longitudes + (lon_index + 2), nodes
                        )  # Second nearest.
                    else:  # lon_index == num_longitudes - 2
                        node._add_neighbour(
                            (lat_index - 1) * num_longitudes + 0, nodes
                        )  # Second nearest.
                    if lat_index > 1:
                        node._add_neighbour(
                            (lat_index - 2) * num_longitudes + (lon_index + 1), nodes
                        )  # Second nearest
                if lat_index < num_latitudes - 1:
                    node._add_neighbour(
                        (lat_index + 1) * num_longitudes + (lon_index + 1), nodes
                    )  # Nearest
                    if lon_index < num_longitudes - 2:
                        node._add_neighbour(
                            (lat_index + 1) * num_longitudes + (lon_index + 2), nodes
                        )  # Second nearest.
                    else:  # lon_index == num_longitudes - 2
                        node._add_neighbour(
                            (lat_index + 1) * num_longitudes + 0, nodes
                        )  # Second nearest.
                    if lat_index < num_latitudes - 2:
                        node._add_neighbour(
                            (lat_index + 2) * num_longitudes + (lon_index + 1), nodes
                        )  # Second nearest
            else:  # lon_index == num_longitudes - 1
                node._add_neighbour(lat_index * num_longitudes, nodes)  # Nearest
                if lat_index > 0:
                    node._add_neighbour(
                        (lat_index - 1) * num_longitudes, nodes
                    )  # Nearest
                    node._add_neighbour(
                        (lat_index - 1) * num_longitudes + 1, nodes
                    )  # Second nearest.
                    if lat_index > 1:
                        node._add_neighbour(
                            (lat_index - 2) * num_longitudes, nodes
                        )  # Second nearest
                if lat_index < num_latitudes - 1:
                    node._add_neighbour(
                        (lat_index + 1) * num_longitudes, nodes
                    )  # Nearest
                    node._add_neighbour(
                        (lat_index + 1) * num_longitudes + 1, nodes
                    )  # Second nearest.
                    if lat_index < num_latitudes - 2:
                        node._add_neighbour(
                            (lat_index + 2) * num_longitudes, nodes
                        )  # Second nearest

            # Bottom neighbour.
            # Note: nodes in bottom row have no bottom neighbour.
            if lat_index > 0:
                node._add_neighbour(
                    (lat_index - 1) * num_longitudes + lon_index, nodes
                )  # Nearest

            # Top neighbour.
            # Note: nodes in top row have no top neighbour.
            if lat_index < num_latitudes - 1:
                node._add_neighbour(
                    (lat_index + 1) * num_longitudes + lon_index, nodes
                )  # Nearest

        return node._neighbours

    def _init_grid(self):

        num_latitudes = self.num_latitudes
        num_longitudes = self.num_longitudes
        grid_spacing_degrees = self.grid_spacing_degrees

        nodes = []

        # Create the nodes.
        for lat_index in range(num_latitudes):
            # The 0.5 puts the point in the centre of the grid pixel.
            # This also avoids sampling right on the poles.
            lat = -90 + (lat_index + 0.5) * grid_spacing_degrees

            for lon_index in range(num_longitudes):
                # The 0.5 puts the point in the centre of the grid pixel.
                # This also avoids sampling right on the dateline where there might be
                # age grid or static polygon artifacts.
                lon = -180 + (lon_index + 0.5) * grid_spacing_degrees

                node = GridNode(pygplates.PointOnSphere(lat, lon))
                nodes.append(node)

        self.nodes = nodes

    def _init_quad_tree(self):

        root_quad_tree_nodes = []

        # Each root quad tree node is quadrant of the globe (square in lat/lon space of size 90 x 90 degrees).
        for root_node_lat_index in range(2):
            for root_node_lon_index in range(4):
                root_node = self._create_quad_tree_node(
                    root_node_lon_index, root_node_lat_index, 0
                )
                root_quad_tree_nodes.append(root_node)

        self.root_quad_tree_nodes = root_quad_tree_nodes

    def _create_quad_tree_node(self, node_lon_index, node_lat_index, level):

        # Create the points of the polygon bounding the current quad tree node.
        bounding_polygon_points = []

        quad_tree_node_to_grid_node_factor = 1 << (self.subdivision_depth - level)

        start_lat_index = node_lat_index * quad_tree_node_to_grid_node_factor
        stop_lat_index = (node_lat_index + 1) * quad_tree_node_to_grid_node_factor
        start_lon_index = node_lon_index * quad_tree_node_to_grid_node_factor
        stop_lon_index = (node_lon_index + 1) * quad_tree_node_to_grid_node_factor

        left_lon = -180 + start_lon_index * self.grid_spacing_degrees
        right_lon = -180 + stop_lon_index * self.grid_spacing_degrees
        bottom_lat = -90 + start_lat_index * self.grid_spacing_degrees
        top_lat = -90 + stop_lat_index * self.grid_spacing_degrees

        # Northern and southern hemispheres handled separately.
        if start_lat_index >= self.num_latitudes / 2:
            # Northern hemisphere.
            left_boundary = pygplates.PolylineOnSphere([(0, left_lon), (90, left_lon)])
            right_boundary = pygplates.PolylineOnSphere(
                [(0, right_lon), (90, right_lon)]
            )

            # Midpoint of small circle arc bounding the bottom of quad tree node.
            bottom_mid_point = pygplates.PointOnSphere(
                bottom_lat, 0.5 * (left_lon + right_lon)
            )

            # Find the great circle (rotation) that passes through the bottom midpoint (and is oriented towards North pole).
            bottom_great_circle_rotation_axis = pygplates.Vector3D.cross(
                bottom_mid_point.to_xyz(),
                pygplates.Vector3D.cross(
                    pygplates.PointOnSphere.north_pole.to_xyz(),
                    bottom_mid_point.to_xyz(),
                ),
            ).to_normalised()
            bottom_great_circle_rotation = pygplates.FiniteRotation(
                bottom_great_circle_rotation_axis.to_xyz(), 0.5 * math.pi
            )

            # Intersect great circle bottom boundary with left and right boundaries to find bottom-left and bottom-right points.
            # The bottom boundary is actually a small circle (due to lat/lon grid), but since we need to use *great* circle arcs
            # in our geometries we need to be a bit loose with our bottom boundary otherwise it will go inside the quad tree node.
            bottom_boundary = pygplates.PolylineOnSphere(
                [
                    bottom_great_circle_rotation * bottom_mid_point,
                    bottom_mid_point,
                    bottom_great_circle_rotation.get_inverse() * bottom_mid_point,
                ]
            )
            _, _, bottom_left_point = pygplates.GeometryOnSphere.distance(
                bottom_boundary, left_boundary, return_closest_positions=True
            )
            _, _, bottom_right_point = pygplates.GeometryOnSphere.distance(
                bottom_boundary, right_boundary, return_closest_positions=True
            )

            bounding_polygon_points.append(bottom_left_point)
            bounding_polygon_points.append(bottom_right_point)

            bounding_polygon_points.append(pygplates.PointOnSphere(top_lat, right_lon))
            bounding_polygon_points.append(pygplates.PointOnSphere(top_lat, left_lon))
        else:
            # Southern hemisphere.
            left_boundary = pygplates.PolylineOnSphere([(0, left_lon), (-90, left_lon)])
            right_boundary = pygplates.PolylineOnSphere(
                [(0, right_lon), (-90, right_lon)]
            )

            # Midpoint of small circle arc bounding the top of quad tree node.
            top_mid_point = pygplates.PointOnSphere(
                top_lat, 0.5 * (left_lon + right_lon)
            )

            # Find the great circle (rotation) that passes through the top midpoint (and is oriented towards North pole).
            top_great_circle_rotation_axis = pygplates.Vector3D.cross(
                top_mid_point.to_xyz(),
                pygplates.Vector3D.cross(
                    pygplates.PointOnSphere.north_pole.to_xyz(), top_mid_point.to_xyz()
                ),
            ).to_normalised()
            top_great_circle_rotation = pygplates.FiniteRotation(
                top_great_circle_rotation_axis.to_xyz(), 0.5 * math.pi
            )

            # Intersect great circle top boundary with left and right boundaries to find top-left and top-right points.
            # The top boundary is actually a small circle (due to lat/lon grid), but since we need to use *great* circle arcs
            # in our geometries we need to be a bit loose with our top boundary otherwise it will go inside the quad tree node.
            top_boundary = pygplates.PolylineOnSphere(
                [
                    top_great_circle_rotation * top_mid_point,
                    top_mid_point,
                    top_great_circle_rotation.get_inverse() * top_mid_point,
                ]
            )
            _, _, top_left_point = pygplates.GeometryOnSphere.distance(
                top_boundary, left_boundary, return_closest_positions=True
            )
            _, _, top_right_point = pygplates.GeometryOnSphere.distance(
                top_boundary, right_boundary, return_closest_positions=True
            )

            bounding_polygon_points.append(top_left_point)
            bounding_polygon_points.append(top_right_point)

            bounding_polygon_points.append(
                pygplates.PointOnSphere(bottom_lat, right_lon)
            )
            bounding_polygon_points.append(
                pygplates.PointOnSphere(bottom_lat, left_lon)
            )

        bounding_polygon = pygplates.PolygonOnSphere(bounding_polygon_points)
        quad_tree_node = GridQuadTreeNode(bounding_polygon)

        if level + Grid.GRID_NODE_DEPTH_PER_QUAD_TREE_NODE >= self.subdivision_depth:
            # Reached leaf quad tree node, so add the grid point indices.
            quad_tree_node.grid_node_indices = []
            for lat_index in range(start_lat_index, stop_lat_index):
                for lon_index in range(start_lon_index, stop_lon_index):
                    node_index = lat_index * self.num_longitudes + lon_index
                    quad_tree_node.grid_node_indices.append(node_index)
        else:
            # Create four child quad tree nodes.
            quad_tree_node.child_quad_tree_nodes = []
            for child_node_lat_offset in range(2):
                for child_node_lon_offset in range(2):
                    quad_tree_node.child_quad_tree_nodes.append(
                        self._create_quad_tree_node(
                            2 * node_lon_index + child_node_lon_offset,
                            2 * node_lat_index + child_node_lat_offset,
                            level + 1,
                        )
                    )

        return quad_tree_node


class GridNode(object):
    def __init__(self, point):
        self.point = point
        self._neighbours = None

    def _add_neighbour(self, neighbour_node_index, nodes):
        distance_to_neighbour = pygplates.GeometryOnSphere.distance(
            self.point, nodes[neighbour_node_index].point
        )
        self._neighbours.append((distance_to_neighbour, neighbour_node_index))


class GridQuadTreeNode(object):
    def __init__(self, bounding_polygon):
        self.bounding_polygon = bounding_polygon
        # If an internal quad tree node then 'child_quad_tree_nodes' will be a list of
        # 4 child quad tree nodes and 'grid_node_indices' will be None.
        # Otherwise quad tree node is a leaf node where 'grid_node_indices' is a list of grid nodes and
        # 'child_quad_tree_nodes' will be None.
        self.child_quad_tree_nodes = None
        self.grid_node_indices = None


class ObstacleGrid(object):
    def __init__(self, grid, obstacle_geometries):
        self.grid = grid
        self.obstacle_geometries = obstacle_geometries

        # Separate obstacles into polygons and non-polygons.
        obstacle_polygons = []
        obstacle_non_polygons = []
        for obstacle_geometry in obstacle_geometries:
            if isinstance(obstacle_geometry, pygplates.PolygonOnSphere):
                obstacle_polygons.append(obstacle_geometry)
            else:
                obstacle_non_polygons.append(obstacle_geometry)

        # Sort the obstacle polygons from largest to smallest area.
        # This makes searching for points/geometries more efficient.
        obstacle_polygons = sorted(
            obstacle_polygons, key=lambda polygon: polygon.get_area(), reverse=True
        )
        self.obstacle_polygons = obstacle_polygons

        self._init_obstacle_grid()

    def create_distance_grid(self, source_geometries, distance_threshold_radians=None):
        return DistanceGrid(self, source_geometries, distance_threshold_radians)

    # Returns a list of 2-tuples (distance_radians, node_index) for neighbour grid nodes of 'node_index'.
    # This is similar to Grid.get_neighbour_grid_nodes() except neighbours that cross obstacle outlines are removed.
    def get_neighbour_grid_nodes(self, node_index):
        nodes = self.nodes
        node = nodes[node_index]

        if node._neighbours is None:
            node._neighbours = neighbours = []
            nearby_obstacle_geometries = node._nearby_obstacle_geometries

            # Get all neighbours from the grid.
            grid_node_neighbours = self.grid.get_neighbour_grid_nodes(node_index)
            for grid_node_neighbour in grid_node_neighbours:
                _, neighbour_node_index = grid_node_neighbour
                neighbour_node = nodes[neighbour_node_index]

                # Include neighbour node if it exists (ie, is outside all polygon obstacles) and
                # does not cross any nearby obstacle geometries.
                if neighbour_node is not None:
                    if nearby_obstacle_geometries is not None:
                        add_neighbour = True
                        grid_nodes = self.grid.nodes
                        # Single segment polyline from node point to neighbour node point.
                        node_to_neighbour_line = pygplates.PolylineOnSphere(
                            (
                                grid_nodes[node_index].point,
                                grid_nodes[neighbour_node_index].point,
                            )
                        )
                        # If intersects any obstacle outline then do not add as a neighbour.
                        for obstacle_geometry in nearby_obstacle_geometries:
                            if (
                                pygplates.GeometryOnSphere.distance(
                                    node_to_neighbour_line,
                                    obstacle_geometry,
                                    # Arbitrarily small threshold for efficiency since only interested in zero distance (intersection)...
                                    1e-4,
                                )
                                == 0
                            ):
                                add_neighbour = False
                                break

                        if add_neighbour:
                            neighbours.append(grid_node_neighbour)
                    else:
                        neighbours.append(grid_node_neighbour)

        return node._neighbours

    def _init_obstacle_grid(self):
        # Mark nodes that are inside obstacles.
        # By default all grid nodes are outside obstacles.
        # If any are found to be inside then we'll set the relevant grid nodes to False.
        self.node_is_outside_obstacle_polygons = [True] * len(self.grid.nodes)

        # Use a quad tree for efficiency - enables us to cull large groups of grid points that are either
        # outside all obstacles or inside an obstacle (avoids point-in-polygon tests for these points).
        for root_quad_tree_node in self.grid.root_quad_tree_nodes:
            self._init_nodes_outside_obstacle_polygons(
                root_quad_tree_node, self.obstacle_polygons
            )

        # Create a node for each grid node that is outside all obstacle polygons.
        self.nodes = [
            ObstacleGridNode() if outside_obstacle_polygons else None
            for outside_obstacle_polygons in self.node_is_outside_obstacle_polygons
        ]

        #
        # Find all obstacles (polygon and non-polygon) near each node.
        #
        # These will be used later when node neighbours are requested (only those neighbours that
        # don't cross an obstacle outline are considered neighbours).
        #

        # We want to find obstacles within the maximum distance between a node and one of its neighbours
        # used in the distance propagation because these obstacles can come between a node and its neighbour.
        nearby_distance_threshold = (
            self.grid.maximum_distance_radians_to_neighbour_grid_node
        )
        # Use a quad tree for efficiency - enables us to cull large groups of grid points that are not near
        # any obstacles (avoids neighbour intersection tests for these points).
        for root_quad_tree_node in self.grid.root_quad_tree_nodes:
            self._init_obstacle_geometries_near_nodes(
                root_quad_tree_node, self.obstacle_geometries, nearby_distance_threshold
            )

    def _init_nodes_outside_obstacle_polygons(
        self, quad_tree_node, parent_overlapping_obstacle_polygons
    ):
        # See if the current quad tree node's bounding polygon overlaps any obstacle polygons.
        overlapping_obstacle_polygons = []
        for obstacle_polygon in parent_overlapping_obstacle_polygons:

            # See if quad tree node and current obstacle polygon overlap.
            if (
                pygplates.GeometryOnSphere.distance(
                    quad_tree_node.bounding_polygon,
                    obstacle_polygon,
                    1e-4,  # Arbitrarily small threshold for efficiency since only interested in zero distance (intersection).
                    geometry1_is_solid=True,
                    geometry2_is_solid=True,
                )
                == 0
            ):

                overlapping_obstacle_polygons.append(obstacle_polygon)

                # See if quad tree node is contained completely inside obstacle polygon.
                # We test this by only considering the quad tree node polygon as solid (the obstacle polygon is an outline).
                if (
                    pygplates.GeometryOnSphere.distance(
                        quad_tree_node.bounding_polygon,
                        obstacle_polygon,
                        1e-4,  # Arbitrarily small threshold for efficiency since only interested in zero distance (intersection).
                        geometry1_is_solid=True,
                    )
                    != 0
                ):

                    # Recursively fill the entire quad sub-tree as inside obstacles.
                    self._fill_quad_tree_node_inside_obstacles(quad_tree_node)
                    return

        # If quad tree is outside all obstacles then nothing left to do since all grid nodes
        # are marked as outside by default.
        if not overlapping_obstacle_polygons:
            return

        # Visit child nodes (if internal node) or test each grid point (if leaf node).
        if quad_tree_node.child_quad_tree_nodes:
            for child_quad_tree_node in quad_tree_node.child_quad_tree_nodes:
                self._init_nodes_outside_obstacle_polygons(
                    child_quad_tree_node, overlapping_obstacle_polygons
                )
        else:
            grid_nodes = self.grid.nodes
            node_is_outside_obstacle_polygons = self.node_is_outside_obstacle_polygons
            for node_index in quad_tree_node.grid_node_indices:
                node_point = grid_nodes[node_index].point
                for polygon in overlapping_obstacle_polygons:
                    if polygon.is_point_in_polygon(node_point):
                        # Node is inside an obstacle.
                        node_is_outside_obstacle_polygons[node_index] = False
                        break

    def _fill_quad_tree_node_inside_obstacles(self, quad_tree_node):
        if quad_tree_node.child_quad_tree_nodes:
            for child_quad_tree_node in quad_tree_node.child_quad_tree_nodes:
                self._fill_quad_tree_node_inside_obstacles(child_quad_tree_node)
        else:
            node_is_outside_obstacle_polygons = self.node_is_outside_obstacle_polygons
            for node_index in quad_tree_node.grid_node_indices:
                # Node is inside an obstacle.
                node_is_outside_obstacle_polygons[node_index] = False

    def _init_obstacle_geometries_near_nodes(
        self,
        quad_tree_node,
        parent_nearby_obstacle_geometries,
        nearby_distance_threshold,
    ):
        # See if the current quad tree node's bounding polygon is near any obstacle geometries.
        nearby_obstacle_geometries = []
        for obstacle_geometry in parent_nearby_obstacle_geometries:

            # See if current obstacle geometry is near (or intersects or completely inside) quad tree node.
            if (
                pygplates.GeometryOnSphere.distance(
                    quad_tree_node.bounding_polygon,
                    obstacle_geometry,
                    nearby_distance_threshold,
                    geometry1_is_solid=True,
                )
                is not None
            ):

                nearby_obstacle_geometries.append(obstacle_geometry)

        # If quad tree is not near all obstacles then nothing left to do since all grid nodes
        # will have no nearby obstacle geometries by default.
        if not nearby_obstacle_geometries:
            return

        # Visit child nodes (if internal node) or assign each grid point the current list of nearby obstacles (if leaf node).
        if quad_tree_node.child_quad_tree_nodes:
            for child_quad_tree_node in quad_tree_node.child_quad_tree_nodes:
                self._init_obstacle_geometries_near_nodes(
                    child_quad_tree_node,
                    nearby_obstacle_geometries,
                    nearby_distance_threshold,
                )
        else:
            nodes = self.nodes
            for node_index in quad_tree_node.grid_node_indices:
                node = nodes[node_index]
                # We only have nodes at grid points outside all obstacle polygons.
                if node:
                    node._nearby_obstacle_geometries = nearby_obstacle_geometries


class ObstacleGridNode(object):
    def __init__(self):
        self._neighbours = None
        self._nearby_obstacle_geometries = None


class DistanceGrid(object):
    def __init__(
        self, obstacle_grid, source_geometries, distance_threshold_radians=None
    ):
        self.obstacle_grid = obstacle_grid
        self.grid = obstacle_grid.grid
        self.distance_threshold_radians = distance_threshold_radians

        # Distance around a source or target geometry that nodes must be within.
        # A grid spacing multiple of 0.5 is slightly too small and 1.0 includes a bit too many nodes.
        self._node_to_geometry_distance_threshold = 0.7 * self.grid.spacing_radians
        self._init_source_distances(source_geometries)

    def shortest_distance(self, target_geometry):
        # If target geometry is a point then use an optimised path.
        try:
            nearest_grid_nodes = self.grid.get_nearest_grid_nodes(target_geometry)

            # Find smoothed distance based with neighbour weights based on distance from target point.
            sum_weights = 0.0
            sum_weighted_distances = 0.0
            any_node_outside_obstacles = False
            for distance_node_to_target, node_index in nearest_grid_nodes:
                node = self.nodes[node_index]
                if node is not None:
                    # Node must be outside obstacles if we have a distance for it.
                    any_node_outside_obstacles = True
                    # If node right on target point then return it.
                    if distance_node_to_target == 0.0:
                        return node.distance_radians
                    # Weight node based on its distance from target geometry.
                    weight = 1.0 / distance_node_to_target
                    sum_weights += weight
                    sum_weighted_distances += weight * node.distance_radians
                elif self.obstacle_grid.node_is_outside_obstacle_polygons[node_index]:
                    # Node is unreachable (has no distance) but is still outside obstacles.
                    any_node_outside_obstacles = True

            if sum_weights != 0.0:
                return sum_weighted_distances / sum_weights
            elif any_node_outside_obstacles:
                # There is at least one nearest node outside obstacles but it is unreachable (has no distance)
                # so return None.
                return
            # else all nearest nodes are inside obstacles so we fall through and use code below
            # which searches wider to find the nearest nodes outside obstacles...

        except AttributeError:
            # Multipoints, polylines and polygons use the code below.
            pass

        # Use a quad tree for efficiency - enables to visit closer groups of nodes first and hence
        # reduce the distance threshold such that the further groups (visited afterwards) are culled.
        node_to_target_infos = []
        self._get_node_to_geometry_distances(
            self.grid.root_quad_tree_nodes,
            target_geometry,
            node_to_target_infos,
            # If threshold is None then set it to PI so we have something to compare against...
            (
                self.distance_threshold_radians
                if self.distance_threshold_radians is not None
                else math.pi
            ),
        )

        if not node_to_target_infos:
            return

        node_to_target_infos.sort()
        distance_closest_node_to_target, closest_node_index = node_to_target_infos[0]

        distance_threshold_node_to_target = (
            distance_closest_node_to_target + self._node_to_geometry_distance_threshold
        )
        if distance_threshold_node_to_target > math.pi:
            distance_threshold_node_to_target = math.pi

        min_distance_node_to_target = None
        for distance_node_to_target, node_index in node_to_target_infos:
            if distance_node_to_target > distance_threshold_node_to_target:
                # Skip all remaining nodes since their distance will be larger (since list sorted by distance).
                break

            node = self.nodes[node_index]
            if node is not None:
                # Find node with shortest distance from source geometry.
                if (
                    min_distance_node_to_target is None
                    or distance_node_to_target + node.distance_radians
                    < min_distance_node_to_target
                ):
                    min_distance_node_to_target = (
                        distance_node_to_target + node.distance_radians
                    )

        # Can be None if the target is unreachable from the source (ie, no path around obstacles) or
        # if the shortest path exceeds a user-specified threshold.
        return min_distance_node_to_target

    def _init_source_distances(self, source_geometries):
        distance_threshold_radians = self.distance_threshold_radians
        obstacle_grid = self.obstacle_grid
        grid = self.grid
        grid_nodes = grid.nodes
        num_nodes = len(grid_nodes)

        self.nodes = nodes = [None] * num_nodes

        # The minimum distance heap used by Dijkstra's algorithm below.
        distance_heap = []

        for source in source_geometries:

            # Find the grid nodes that we will initialise with our source geometries.
            #
            # Use a quad tree for efficiency - enables to visit closer groups of nodes first and hence
            # reduce the distance threshold such that the further groups (visited afterwards) are culled.
            node_to_source_infos = []
            self._get_node_to_geometry_distances(
                grid.root_quad_tree_nodes,
                source,
                node_to_source_infos,
                # If threshold is None then set it to PI so we have something to compare against...
                (
                    distance_threshold_radians
                    if distance_threshold_radians is not None
                    else math.pi
                ),
            )

            if node_to_source_infos:
                node_to_source_infos.sort()
                distance_closest_node_to_source, closest_node_index = (
                    node_to_source_infos[0]
                )

                # The source geometry might be deep inside the obstacle polygons so we need to expand
                # the minimum distance from source geometry to candidate grid nodes so that grid nodes
                # outside all obstacles can be found and initialised.
                # We choose the closest outside grid node and add the grid spacing to ensure it gets included.
                distance_threshold_node_to_source = (
                    distance_closest_node_to_source
                    + self._node_to_geometry_distance_threshold
                )
                if distance_threshold_node_to_source > math.pi:
                    distance_threshold_node_to_source = math.pi

                for distance_node_to_source, node_index in node_to_source_infos:
                    if distance_node_to_source > distance_threshold_node_to_source:
                        # Skip all remaining nodes since their distance will be larger (since list sorted by distance).
                        break

                    if nodes[node_index] is None:
                        # First time visiting this grid node.
                        nodes[node_index] = DistanceGridNode(distance_node_to_source)
                        heapq.heappush(
                            distance_heap, (distance_node_to_source, node_index)
                        )
                    elif distance_node_to_source < nodes[node_index].distance_radians:
                        # Already visited this grid node (from a previous source geometry).
                        nodes[node_index].distance_radians = distance_node_to_source
                        heapq.heappush(
                            distance_heap, (distance_node_to_source, node_index)
                        )

        # Keep track of which grid nodes have been processed by Dijkstra's algorithm below.
        processed_nodes = [False] * num_nodes

        # Propagate the source grid nodes initialised above to all outside grid nodes within
        # the threshold distance, or until can propagate no further (eg, if blocked by obstacles).
        while distance_heap:
            # Get the minimum distance.
            distance_node_to_source, node_index = heapq.heappop(distance_heap)

            # Skip grid node if it has already been processed.
            # We need this because we cannot update the priority of a grid node
            # (when visiting neighbours below) when it's in the middle of the heap.
            # Instead we add multiple entries into the heap for a single grid node and
            # only process the first entry encountered (it'll be the one with the minimum distance).
            if processed_nodes[node_index]:
                continue
            processed_nodes[node_index] = True

            # Visit the neighbours of the current grid node that are outside all obstacles
            # and have not yet been processed.
            grid_node_neighbours = obstacle_grid.get_neighbour_grid_nodes(node_index)
            node = nodes[node_index]
            for distance_to_neighbour, neighbour_node_index in grid_node_neighbours:
                if not processed_nodes[neighbour_node_index]:

                    neighbour_distance = node.distance_radians + distance_to_neighbour
                    # Skip update if distance exceeds threshold.
                    if (
                        distance_threshold_radians is not None
                        and neighbour_distance > distance_threshold_radians
                    ):
                        continue

                    neighbour_node = nodes[neighbour_node_index]
                    if neighbour_node is None:
                        # First time visiting this grid node.
                        nodes[neighbour_node_index] = neighbour_node = DistanceGridNode(
                            neighbour_distance
                        )
                        heapq.heappush(
                            distance_heap, (neighbour_distance, neighbour_node_index)
                        )
                    elif neighbour_distance < neighbour_node.distance_radians:
                        neighbour_node.distance_radians = neighbour_distance
                        # Already visited this grid node since it's already in the heap.
                        # Normally we'd adjust the priority of this heap entry but we cannot do that.
                        # Instead we add another entry into the heap for this grid node and only process
                        # the first entry encountered (it'll be the one with the minimum distance
                        # which will be this entry if no further updates to this node occur)
                        # and rejecting the second, third, etc, entries.
                        heapq.heappush(
                            distance_heap, (neighbour_distance, neighbour_node_index)
                        )

    def _get_node_to_geometry_distances(
        self,
        quad_tree_nodes,
        geometry,
        node_to_geometry_infos,
        current_distance_threshold_radians,
    ):
        # Sort quad tree nodes by distance.
        # This makes it very quick to converge on the minimum distance and means we can skip
        # processing a large number of subsequently visited nodes that are further away.
        distance_sorted_quad_tree_nodes = []
        for quad_tree_node in quad_tree_nodes:
            distance_geometry_to_quad_tree_node = pygplates.GeometryOnSphere.distance(
                quad_tree_node.bounding_polygon,
                geometry,
                current_distance_threshold_radians,
                geometry1_is_solid=True,
            )
            if distance_geometry_to_quad_tree_node is not None:
                distance_sorted_quad_tree_nodes.append(
                    (distance_geometry_to_quad_tree_node, quad_tree_node)
                )

        if not distance_sorted_quad_tree_nodes:
            return current_distance_threshold_radians

        # Visit the closest quad tree nodes first since they will reduce the distance threshold.
        distance_sorted_quad_tree_nodes.sort(key=lambda d: d[0])
        for (
            distance_geometry_to_quad_tree_node,
            quad_tree_node,
        ) in distance_sorted_quad_tree_nodes:
            if distance_geometry_to_quad_tree_node < current_distance_threshold_radians:
                # Visit child nodes (if internal node) or test each grid point (if leaf node).
                if quad_tree_node.child_quad_tree_nodes:
                    current_distance_threshold_radians = (
                        self._get_node_to_geometry_distances(
                            quad_tree_node.child_quad_tree_nodes,
                            geometry,
                            node_to_geometry_infos,
                            current_distance_threshold_radians,
                        )
                    )
                else:
                    grid_nodes = self.grid.nodes
                    node_to_geometry_distance_threshold = (
                        self._node_to_geometry_distance_threshold
                    )
                    node_is_outside_obstacle_polygons = (
                        self.obstacle_grid.node_is_outside_obstacle_polygons
                    )
                    for node_index in quad_tree_node.grid_node_indices:
                        if node_is_outside_obstacle_polygons[node_index]:
                            distance_node_to_geometry = (
                                pygplates.GeometryOnSphere.distance(
                                    grid_nodes[node_index].point,
                                    geometry,
                                    current_distance_threshold_radians,
                                )
                            )

                            if distance_node_to_geometry is not None:
                                # Reduce the distance threshold if possible, so we can converge quicker.
                                if (
                                    distance_node_to_geometry
                                    + node_to_geometry_distance_threshold
                                    < current_distance_threshold_radians
                                ):
                                    current_distance_threshold_radians = (
                                        distance_node_to_geometry
                                        + node_to_geometry_distance_threshold
                                    )

                                node_to_geometry_infos.append(
                                    (distance_node_to_geometry, node_index)
                                )

        return current_distance_threshold_radians


class DistanceGridNode(object):
    def __init__(self, distance_radians):
        self.distance_radians = distance_radians
