# -*- coding: utf-8 -*-

"""
    Copyright (C) 2017 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


#########################################################################
# Spatial tree for points (pygplates.PointOnSphere).                    #
#                                                                       #
# Can be used to speed up point-in-polygon and minimum distance queries #
# from points to arbitrary geometries.                                  #
#########################################################################


from __future__ import print_function
import math
import pygplates
import sys

# the following line would ensure more correct python2-3 compatibility,
# but requires a non-standard module
# from builtins import range

# The depth of the internal lat/lon quad tree.
# A value of 4 seems to work quite well for a uniform lat/lon spacing of 'points' of 1 degree
# and below without the cost of generating a deep quad tree.
DEFAULT_SUBDIVISION_DEPTH = 4


class PointsSpatialTree(object):
    
    def __init__(self, points, subdivision_depth = DEFAULT_SUBDIVISION_DEPTH):
        """
        Construct a spatial tree from a sequence of points up to a maximum tree depth.
        
        points: a sequence of 'pygplates.PointOnSphere'.
        
        subdivision_depth: The depth of the internal lat/lon quad tree.
                           The lat/lon width of a leaf quad tree node is (90 / (2^subdivision_depth)) degrees.
                           Generally the denser the 'points' the larger the depth should be.
                           Setting this value too high causes unnecessary time to be spent generating a deep quad tree.
                           Setting this value too low reduces the culling/querying efficiency of the quad tree.
                           However a value of 4 seems to work quite well for a uniform lat/lon spacing of 'points' of 1 degree
                           and below without the cost of generating a deep quad tree.
                           So most of the time the subdivision depth can be left at its default value.
        
        Raises ValueError if 'subdivision_depth' is not in the range [0, 100].
        """
        
        if subdivision_depth < 0:
            raise ValueError('Subdivision depth must be a non-negative value.')
        elif subdivision_depth > 100:
            raise ValueError('Subdivision depth is too large (should be 100 or less).')
        
        # Each root quad tree node is quadrant of the globe (square in lat/lon space of size 90 x 90 degrees).
        # So there are 8 of them.
        # We'll only create them as needed.
        self._root_nodes = [None] * 8
        
        # Place each point in a quad tree leaf node.
        for point_index, point in enumerate(points):
            point_lat, point_lon = point.to_lat_lon()
            
            # Get root node that current point is in.
            if point_lat < 0:
                root_node_lat_index = 0
            else:
                root_node_lat_index = 1
            
            if point_lon < 0:
                if point_lon < -90:
                    root_node_lon_index = 0
                else:
                    root_node_lon_index = 1
            else:
                if point_lon < 90:
                    root_node_lon_index = 2
                else:
                    root_node_lon_index = 3
            
            root_node_index = 4 * root_node_lat_index + root_node_lon_index
            root_node = self._root_nodes[root_node_index]
            
            is_north_hemisphere = (root_node_lat_index == 1)
            root_node_half_width_degrees = 45.0
            root_node_centre_lon = -180 + 90 * root_node_lon_index + root_node_half_width_degrees
            root_node_centre_lat = -90 + 90 * root_node_lat_index + root_node_half_width_degrees
            
            # Create root node if first time visited.
            if root_node is None:
                root_node = PointsSpatialTreeNode(
                        root_node_centre_lon,
                        root_node_centre_lat,
                        root_node_half_width_degrees,
                        is_north_hemisphere)
                self._root_nodes[root_node_index] = root_node
            
            # Iterate through the subdivision levels and place current point in the correct quad tree leaf node.
            node = root_node
            node_half_width_degrees = root_node_half_width_degrees
            node_centre_lon = root_node_centre_lon
            node_centre_lat = root_node_centre_lat
            for level in range(0, subdivision_depth):
                
                child_node_half_width_degrees = node_half_width_degrees / 2.0
                
                if point_lat < node_centre_lat:
                    child_node_lat_offset = 0
                    child_node_centre_lat = node_centre_lat - child_node_half_width_degrees
                else:
                    child_node_lat_offset = 1
                    child_node_centre_lat = node_centre_lat + child_node_half_width_degrees
                
                if point_lon < node_centre_lon:
                    child_node_lon_offset = 0
                    child_node_centre_lon = node_centre_lon - child_node_half_width_degrees
                else:
                    child_node_lon_offset = 1
                    child_node_centre_lon = node_centre_lon + child_node_half_width_degrees
                
                # The current node is an internal node (because it will have child nodes).
                # Create a list of child nodes if first time visiting node.
                if node._child_nodes is None:
                    # Only create each child node as needed.
                    node._child_nodes = [None] * 4
                
                child_node_index = 2 * child_node_lat_offset + child_node_lon_offset
                child_node = node._child_nodes[child_node_index]
                
                if child_node is None:
                    child_node = PointsSpatialTreeNode(
                            child_node_centre_lon,
                            child_node_centre_lat,
                            child_node_half_width_degrees,
                            is_north_hemisphere)
                    node._child_nodes[child_node_index] = child_node
                
                # Child node becomes parent node in next iteration.
                node = child_node
                node_half_width_degrees = child_node_half_width_degrees
                node_centre_lon = child_node_centre_lon
                node_centre_lat = child_node_centre_lat
            
            # Reached leaf node (end of subdivision).
            # Create a list of point indices if first time visiting node.
            if node._point_indices is None:
                node._point_indices = []
            
            # Add the current point (index) to the leaf node.
            node._point_indices.append(point_index)
    
    
    def get_root_nodes(self):
        """
        Return any root nodes that have points in their subtree.
        
        There are a maximum of 8 root nodes.
        
        Returns: A list of Node.
        """
        
        return [root_node for root_node in self._root_nodes if root_node is not None]


class PointsSpatialTreeNode(object):
    def __init__(self, centre_lon, centre_lat, half_width_degrees, is_north_hemisphere):
        # Parameters to describe location and extents of this node.
        self._centre_lon = centre_lon
        self._centre_lat = centre_lat
        self._half_width_degrees = half_width_degrees
        self._is_north_hemisphere = is_north_hemisphere
        
        # We'll create the bounding polygon/circle when/if they are requested.
        self._bounding_polygon = None
        self._bounding_circle = None
        
        # If an internal quad tree node then '_child_nodes' will be a list of
        # 4 child quad tree nodes and '_point_indices' will be None.
        # Otherwise quad tree node is a leaf node where 'point_indices' is a list of points and
        # 'child_nodes' will be None.
        self._child_nodes = None
        self._point_indices = None
    
    
    def get_bounding_polygon(self):
        """
        Returns a polygon that bounds the current node.
        
        The returned polygon is guaranteed to contain all points in this node (and any child nodes, etc).
        However it is not guaranteed to contain the bounding circle.
        
        Returns: pygplates.PolygonOnSphere
        """
        
        if self._bounding_polygon is None:
            self._create_bounding_polygon()
        
        return self._bounding_polygon
    
    
    def get_bounding_circle(self):
        """
        Returns a small circle that bounds the current node.
        
        The returned small circle is guaranteed to contain all points in this node (and any child nodes, etc).
        However it is not guaranteed to contain the bounding polygon.
        
        Returns: The centre of small circle and its radius (in radians) as the 2-tuple of type (pygplates.PointOnSphere, float).
        """
        
        if self._bounding_circle is None:
            self._create_bounding_circle()
        
        return self._bounding_circle
    
    
    def is_leaf_node(self):
        """
        Returns True if this node is a leaf node.
        
        A leaf node has indices into the sequence of points passed into
        spatial tree constructor - see 'get_point_indices()'.
        But it does not have any child nodes.
        """
        
        return self._point_indices is not None
    
    
    def is_internal_node(self):
        """
        Returns True if this node is an internal node.
        
        An internal node has child nodes - see 'get_child_nodes()'.
        But it does not have any point indices.
        """
        
        return self._child_nodes is not None
    
    
    def get_child_nodes(self):
        """
        Return any child nodes that have points in their subtree.
        
        There are a maximum of 4 child nodes.
        
        Should only be called if 'is_internal_node()' returns True.
        
        Returns: A list of Node.
        """
        
        return [child_node for child_node in self._child_nodes if child_node is not None]
    
    
    def get_point_indices(self):
        """
        Return indices of points that exist in this leaf node.
        
        The indices refer to the sequence of points passed into spatial tree constructor.
        
        Should only be called if 'is_leaf_node()' returns True.
        
        Returns: A list of int.
        """
        
        return self._point_indices
    
    
    def _create_bounding_polygon(self):
        
        # Create the points of the polygon bounding the current quad tree node.
        bounding_polygon_points = []
        
        left_lon = self._centre_lon - self._half_width_degrees
        right_lon = self._centre_lon + self._half_width_degrees
        bottom_lat = self._centre_lat - self._half_width_degrees
        top_lat = self._centre_lat + self._half_width_degrees
        
        # Northern and southern hemispheres handled separately.
        if self._is_north_hemisphere:
            # Northern hemisphere.
            left_boundary = pygplates.PolylineOnSphere([(0, left_lon), (90, left_lon)])
            right_boundary = pygplates.PolylineOnSphere([(0, right_lon), (90, right_lon)])
            
            # Midpoint of small circle arc bounding the bottom of quad tree node.
            bottom_mid_point = pygplates.PointOnSphere(bottom_lat, 0.5 * (left_lon + right_lon))
            
            # Find the great circle (rotation) that passes through the bottom midpoint (and is oriented towards North pole).
            bottom_great_circle_rotation_axis = pygplates.Vector3D.cross(
                    bottom_mid_point.to_xyz(),
                    pygplates.Vector3D.cross(pygplates.PointOnSphere.north_pole.to_xyz(), bottom_mid_point.to_xyz())
                            ).to_normalised()
            bottom_great_circle_rotation = pygplates.FiniteRotation(bottom_great_circle_rotation_axis.to_xyz(), 0.5 * math.pi)
            
            # Intersect great circle bottom boundary with left and right boundaries to find bottom-left and bottom-right points.
            # The bottom boundary is actually a small circle (due to lat/lon grid), but since we need to use *great* circle arcs
            # in our geometries we need to be a bit loose with our bottom boundary otherwise it will go inside the quad tree node.
            bottom_boundary = pygplates.PolylineOnSphere(
                    [bottom_great_circle_rotation * bottom_mid_point, bottom_mid_point, bottom_great_circle_rotation.get_inverse() * bottom_mid_point])
            _, _, bottom_left_point = pygplates.GeometryOnSphere.distance(bottom_boundary, left_boundary, return_closest_positions = True)
            _, _, bottom_right_point = pygplates.GeometryOnSphere.distance(bottom_boundary, right_boundary, return_closest_positions = True)
            
            bounding_polygon_points.append(bottom_left_point)
            bounding_polygon_points.append(bottom_right_point)
            
            bounding_polygon_points.append(pygplates.PointOnSphere(top_lat, right_lon))
            bounding_polygon_points.append(pygplates.PointOnSphere(top_lat, left_lon))
        else:
            # Southern hemisphere.
            left_boundary = pygplates.PolylineOnSphere([(0, left_lon), (-90, left_lon)])
            right_boundary = pygplates.PolylineOnSphere([(0, right_lon), (-90, right_lon)])
            
            # Midpoint of small circle arc bounding the top of quad tree node.
            top_mid_point = pygplates.PointOnSphere(top_lat, 0.5 * (left_lon + right_lon))
            
            # Find the great circle (rotation) that passes through the top midpoint (and is oriented towards North pole).
            top_great_circle_rotation_axis = pygplates.Vector3D.cross(
                    top_mid_point.to_xyz(),
                    pygplates.Vector3D.cross(pygplates.PointOnSphere.north_pole.to_xyz(), top_mid_point.to_xyz())
                            ).to_normalised()
            top_great_circle_rotation = pygplates.FiniteRotation(top_great_circle_rotation_axis.to_xyz(), 0.5 * math.pi)
            
            # Intersect great circle top boundary with left and right boundaries to find top-left and top-right points.
            # The top boundary is actually a small circle (due to lat/lon grid), but since we need to use *great* circle arcs
            # in our geometries we need to be a bit loose with our top boundary otherwise it will go inside the quad tree node.
            top_boundary = pygplates.PolylineOnSphere(
                    [top_great_circle_rotation * top_mid_point, top_mid_point, top_great_circle_rotation.get_inverse() * top_mid_point])
            _, _, top_left_point = pygplates.GeometryOnSphere.distance(top_boundary, left_boundary, return_closest_positions = True)
            _, _, top_right_point = pygplates.GeometryOnSphere.distance(top_boundary, right_boundary, return_closest_positions = True)
            
            bounding_polygon_points.append(top_left_point)
            bounding_polygon_points.append(top_right_point)
            
            bounding_polygon_points.append(pygplates.PointOnSphere(bottom_lat, right_lon))
            bounding_polygon_points.append(pygplates.PointOnSphere(bottom_lat, left_lon))
        
        self._bounding_polygon = pygplates.PolygonOnSphere(bounding_polygon_points)
    
    
    def _create_bounding_circle(self):
        
        bound_circle_centre = pygplates.PointOnSphere(self._centre_lat, self._centre_lon)
        
        left_lon = self._centre_lon - self._half_width_degrees
        right_lon = self._centre_lon + self._half_width_degrees
        bottom_lat = self._centre_lat - self._half_width_degrees
        top_lat = self._centre_lat + self._half_width_degrees
        
        # The small circle that bound the four corner points will also bound the associated lat/lon box.
        # This is because the small circle centre will be inside the box and the small circle curvature
        # is greater than all four sides of the box (two great circle arcs and two small circle arcs).
        dist_to_bottom_left = pygplates.GeometryOnSphere.distance(
                bound_circle_centre,
                pygplates.PointOnSphere(bottom_lat, left_lon))
        dist_to_bottom_right = pygplates.GeometryOnSphere.distance(
                bound_circle_centre,
                pygplates.PointOnSphere(bottom_lat, right_lon))
        dist_to_top_left = pygplates.GeometryOnSphere.distance(
                bound_circle_centre,
                pygplates.PointOnSphere(top_lat, left_lon))
        dist_to_top_right = pygplates.GeometryOnSphere.distance(
                bound_circle_centre,
                pygplates.PointOnSphere(top_lat, right_lon))
        
        bounding_circle_radius_radians = max(
                dist_to_bottom_left,
                dist_to_bottom_right,
                dist_to_top_left,
                dist_to_top_right)
        
        self._bounding_circle = (bound_circle_centre, bounding_circle_radius_radians)
