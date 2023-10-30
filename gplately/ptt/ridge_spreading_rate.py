#!/usr/bin/env python
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


#############################################################################################################################
# Calculate mid-ocean ridge spreading rates. #
#############################################################################################################################


from __future__ import print_function
import math
import pygplates
from . import separate_ridge_transform_segments
import sys



def spreading_rates(
        rotation_features_or_model,
        topology_features,
        time,
        threshold_sampling_distance_radians,
        spreading_feature_types = None,
        transform_segment_deviation_in_radians = separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS,
        velocity_delta_time = 1.0,
        anchor_plate_id = 0):
    """
    Calculates spreading rate and length of ridge segments of spreading features (mid-ocean ridges) of resolved topologies at specified time.
    
    The transform segments of spreading features are ignored.
    
    Resolves topologies at 'time', tessellates all resolved spreading features to within 'threshold_sampling_distance_radians' radians and
    returns a list of tuples where each tuple represents a tessellated point and contains the following parameters:
    
    - point longitude
    - point latitude
    - spreading velocity magnitude (in cm/yr)
    - length of arc segment (in degrees) that current point is on
    
    
    rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).
    
    topology_features: Topology feature collection(s), or list of features, or filename(s) or any combination of those.
    
    time: Reconstruction time to resolved topologies.
    
    threshold_sampling_distance_radians: Threshold sampling distance along spreading features (in radians).
    
    spreading_feature_types: Only spreading features with a feature type contained in this list are considered.
                             If None then all spreading features are considered.
    
    transform_segment_deviation_in_radians: How much a segment can deviate from the stage pole before
                                            it's considered a transform segment (in radians).
    
    velocity_delta_time: Delta time interval used to calculate spreading velocity.
    
    Returns: List of the tuples described above.
    """
    time = float(time)
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn topology data into a list of features (if not already).
    topology_features = pygplates.FeaturesFunctionArgument(topology_features)
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(topology_features.get_features(), rotation_model, resolved_topologies, time, shared_boundary_sections, anchor_plate_id)
    
    # List of tesselated spreading points and associated spreading parameters for the current 'time'.
    output_data = []
    
    # Iterate over the shared boundary sections of all resolved topologies.
    for shared_boundary_section in shared_boundary_sections:
        
        spreading_feature = shared_boundary_section.get_feature()
        
        # Skip sections that are not spreading features (if requested).
        if (spreading_feature_types and
            spreading_feature.get_feature_type() not in spreading_feature_types):
            continue
        
        # Find the stage rotation of the spreading feature in the frame of reference of its reconstructed geometry at the current 'time'.
        # The stage pole can then be directly geometrically compared to the reconstructed spreading geometry.
        spreading_stage_rotation = separate_ridge_transform_segments.get_stage_rotation_for_reconstructed_geometry(
                spreading_feature, rotation_model, time)
        if not spreading_stage_rotation:
            # Skip current feature - it's not a spreading feature.
            continue
        
        # Iterate over the shared sub-segments of the current line.
        # These are the parts of the line that actually contribute to topological boundaries.
        for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
            
            # Split into ridge and transform segments.
            ridge_and_transform_segment_geometries = separate_ridge_transform_segments.separate_geometry_into_ridges_and_transforms(
                    spreading_stage_rotation,
                    shared_sub_segment.get_resolved_geometry(),
                    transform_segment_deviation_in_radians)
            if not ridge_and_transform_segment_geometries:
                # Skip shared sub segment - it's not a polyline (or polygon).
                continue
            
            # Only interested in ridge segments.
            ridge_sub_segment_geometries, _ = ridge_and_transform_segment_geometries
            
            # Ensure the ridge sub-segments are tessellated to within the threshold sampling distance.
            tessellated_shared_sub_segment_polylines = [
                ridge_sub_segment_geometry.to_tessellated(threshold_sampling_distance_radians)
                    for ridge_sub_segment_geometry in ridge_sub_segment_geometries]
            
            # Iterate over the great circle arcs of the tessellated polylines to get the arc midpoints and lengths.
            # There is an arc between each adjacent pair of points in the polyline.
            arc_midpoints = []
            arc_lengths = []
            for tessellated_shared_sub_segment_polyline in tessellated_shared_sub_segment_polylines:
                for arc in tessellated_shared_sub_segment_polyline.get_segments():
                    if not arc.is_zero_length():
                        arc_midpoints.append(arc.get_arc_point(0.5))
                        arc_lengths.append(arc.get_arc_length())
            
            # Shouldn't happen, but just in case ridge sub-segment polylines coincide with points.
            if not arc_midpoints:
                continue
            
            # Calculate the spreading velocities at the arc midpoints.
            #
            # Note that the stage rotation can be used directly on the reconstructed geometries because
            # it is already in the frame of reference of the reconstructed geometries.
            spreading_velocity_vectors = pygplates.calculate_velocities(
                    arc_midpoints,
                    spreading_stage_rotation,
                    velocity_delta_time,
                    pygplates.VelocityUnits.cms_per_yr)
            
            for arc_index in range(len(arc_midpoints)):
                
                arc_midpoint = arc_midpoints[arc_index]
                arc_length = arc_lengths[arc_index]
                lat, lon = arc_midpoint.to_lat_lon()
                
                spreading_velocity_magnitude = spreading_velocity_vectors[arc_index].get_magnitude()
                
                # The data will be output in GMT format (ie, lon first, then lat, etc).
                output_data.append((
                        lon,
                        lat,
                        spreading_velocity_magnitude,
                        math.degrees(arc_length)))
    
    return output_data

def spreading_rates_dense(
        rotation_features_or_model,
        topology_features,
        time,
        threshold_sampling_distance_radians,
        spreading_feature_types = None,
        transform_segment_deviation_in_radians = separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS,
        velocity_delta_time = 1.0,
        anchor_plate_id = 0):
    """
    Calculates spreading rate and length of ridge segments of spreading features (mid-ocean ridges) of resolved topologies at specified time.
    
    The transform segments of spreading features are ignored.
    
    Resolves topologies at 'time', tessellates all resolved spreading features to within 'threshold_sampling_distance_radians' radians and
    returns a list of tuples where each tuple represents a tessellated point and contains the following parameters:
    
    - point longitude
    - point latitude
    - spreading velocity magnitude (in cm/yr)
    - spreading obliquity in degrees
    - length of arc segment (in degrees) that current point is on
    - azimuth of vector normal to the arc segment in degrees
    - left plate ID
    - right plate ID
    
    rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).
    
    topology_features: Topology feature collection(s), or list of features, or filename(s) or any combination of those.
    
    time: Reconstruction time to resolved topologies.
    
    threshold_sampling_distance_radians: Threshold sampling distance along spreading features (in radians).
    
    spreading_feature_types: Only spreading features with a feature type contained in this list are considered.
                             If None then all spreading features are considered.
    
    transform_segment_deviation_in_radians: How much a segment can deviate from the stage pole before
                                            it's considered a transform segment (in radians).
    
    velocity_delta_time: Delta time interval used to calculate spreading velocity.
    
    Returns: List of the tuples described above.
    """
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn topology data into a list of features (if not already).
    topology_features = pygplates.FeaturesFunctionArgument(topology_features)
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(topology_features.get_features(), rotation_model, resolved_topologies, time, shared_boundary_sections, anchor_plate_id)
    
    # List of tesselated spreading points and associated spreading parameters for the current 'time'.
    output_data = []
    
    # Iterate over the shared boundary sections of all resolved topologies.
    for shared_boundary_section in shared_boundary_sections:
        
        spreading_feature = shared_boundary_section.get_feature()
        
        # Skip sections that are not spreading features (if requested).
        if (spreading_feature_types and
            spreading_feature.get_feature_type() not in spreading_feature_types):
            continue
        
        # Find the stage rotation of the spreading feature in the frame of reference of its reconstructed geometry at the current 'time'.
        # The stage pole can then be directly geometrically compared to the reconstructed spreading geometry.
        spreading_stage_rotation = separate_ridge_transform_segments.get_stage_rotation_for_reconstructed_geometry(
                spreading_feature, rotation_model, time)
        if not spreading_stage_rotation:
            # Skip current feature - it's not a spreading feature.
            continue

        pID_left = spreading_feature.get_left_plate()
        pID_right = spreading_feature.get_right_plate()
        
        # Iterate over the shared sub-segments of the current line.
        # These are the parts of the line that actually contribute to topological boundaries.
        for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
            
            # Split into ridge and transform segments.
            ridge_and_transform_segment_geometries = separate_ridge_transform_segments.separate_geometry_into_ridges_and_transforms(
                    spreading_stage_rotation,
                    shared_sub_segment.get_resolved_geometry(),
                    transform_segment_deviation_in_radians)
            if not ridge_and_transform_segment_geometries:
                # Skip shared sub segment - it's not a polyline (or polygon).
                continue
            
            # Only interested in ridge segments.
            ridge_sub_segment_geometries, _ = ridge_and_transform_segment_geometries
            
            # Ensure the ridge sub-segments are tessellated to within the threshold sampling distance.
            tessellated_shared_sub_segment_polylines = [
                ridge_sub_segment_geometry.to_tessellated(threshold_sampling_distance_radians)
                    for ridge_sub_segment_geometry in ridge_sub_segment_geometries]
            
            # Iterate over the great circle arcs of the tessellated polylines to get the arc midpoints and lengths.
            # There is an arc between each adjacent pair of points in the polyline.
            arc_midpoints = []
            arc_lengths = []
            spreading_arc_normals = []
            for tessellated_shared_sub_segment_polyline in tessellated_shared_sub_segment_polylines:
                for arc in tessellated_shared_sub_segment_polyline.get_segments():
                    if not arc.is_zero_length():
                        arc_midpoints.append(arc.get_arc_point(0.5))
                        arc_lengths.append(arc.get_arc_length())
                        # The normal to the spreading ridge.
                        spreading_arc_normals.append(arc.get_great_circle_normal())
            
            # Shouldn't happen, but just in case ridge sub-segment polylines coincide with points.
            if not arc_midpoints:
                continue

            # The spreading arc normals relative to North (azimuth).
            # Convert global 3D normal vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
            spreading_arc_local_normals = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                    arc_midpoints, spreading_arc_normals)
            
            # Calculate the spreading velocities at the arc midpoints.
            #
            # Note that the stage rotation can be used directly on the reconstructed geometries because
            # it is already in the frame of reference of the reconstructed geometries.
            spreading_velocity_vectors = pygplates.calculate_velocities(
                    arc_midpoints,
                    spreading_stage_rotation,
                    velocity_delta_time,
                    pygplates.VelocityUnits.cms_per_yr)
            
            for arc_index in range(len(arc_midpoints)):
                arc_midpoint = arc_midpoints[arc_index]
                arc_length = arc_lengths[arc_index]
                
                arc_midpoint = arc_midpoints[arc_index]
                arc_length = arc_lengths[arc_index]
                spreading_arc_normal = spreading_arc_normals[arc_index]
                spreading_arc_normal_azimuth = spreading_arc_local_normals[arc_index][1]
                lat, lon = arc_midpoint.to_lat_lon()

                # The direction towards which we rotate from the subducting normal in a clockwise fashion.
                clockwise_direction = pygplates.Vector3D.cross(spreading_arc_normal, arc_midpoint.to_xyz())

                # calculate the spreading rate parameters
                spreading_velocity_vector = spreading_velocity_vectors[arc_index]
                spreading_velocity_magnitude = spreading_velocity_vector.get_magnitude()


                # Angle range [0, 180]
                spreading_obliquity_degrees = math.degrees(pygplates.Vector3D.angle_between(
                            spreading_velocity_vector, spreading_arc_normal))
                # Minimum deviation from 'spreading_arc_normal' and '-spreading_arc_normal'.
                # Angle range [0, 90]
                if spreading_obliquity_degrees > 90:
                    spreading_obliquity_degrees = 180 - spreading_obliquity_degrees

                # The data will be output in GMT format (ie, lon first, then lat, etc).
                output_data.append((
                        lon,
                        lat,
                        spreading_velocity_magnitude,
                        spreading_obliquity_degrees,
                        math.degrees(arc_length),
                        math.degrees(spreading_arc_normal_azimuth),
                        pID_left,
                        pID_right))
    
    return output_data
