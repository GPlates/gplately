#
#    Copyright (C) 2017-2025 The University of Sydney, Australia
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
Calculate mid-ocean ridge spreading rates.
"""
from __future__ import print_function

import math
import sys

import pygplates

from . import separate_ridge_transform_segments


def spreading_rates(
    rotation_features_or_model,
    topology_features,
    time,
    threshold_sampling_distance_radians,
    spreading_feature_types=[pygplates.FeatureType.gpml_mid_ocean_ridge],
    transform_segment_deviation_in_radians=separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS,
    velocity_delta_time=1.0,
    anchor_plate_id=None,
    include_network_boundaries=False,
    output_obliquity_and_normal_and_left_right_plates=False,
):
    """Calculate spreading rate of ridge segments (and other statistics like spreading obqliquity, and ridge length/normals and left/right plates).

    The transform segments of spreading features are ignored (unless ``transform_segment_deviation_in_radians`` is ``None``).

    Resolves topologies at ``time``, tessellates all resolved spreading features to within ``threshold_sampling_distance_radians`` radians and
    returns a list of tuples where each tuple represents a tessellated point.

    Parameters
    ----------
    rotation_features_or_model : str/os.PathLike, or `pygplates.FeatureCollection`_, or `pygplates.Feature`_, or sequence of `pygplates.Feature`_, or a sequence (eg, a list or tuple) of any combination of those four types, or `pygplates.RotationModel`_
        A rotation model to calculate spreading rate.
        Can be provided as a rotation filename, or rotation feature collection, or rotation feature, or
        sequence of rotation features, or a sequence (eg, a list or tuple) of any combination of those four types.
    topology_features : str/os.PathLike, or `pygplates.FeatureCollection`_, or `pygplates.Feature`_, or a sequence `pygplates.Feature`_, or a sequence (eg, a list or tuple) of any combination of those four types
        Reconstructable topological features representing the spreading features.
    time : float
        The reconstruction time (Ma) at which to calculate spreading rate.
    threshold_sampling_distance_radians : float
        Threshold sampling distance along spreading features (in radians).
    spreading_feature_types : `pygplates.FeatureType`_ or sequence of `pygplates.FeatureType`_, default=pygplates.FeatureType.gpml_mid_ocean_ridge
        Only spreading features with a feature type contained in this list are considered.
        If ``None`` then all spreading features are considered.
    transform_segment_deviation_in_radians : float, default=<implementation-defined>
        How much a segment can deviate from the stage pole before it's considered a transform segment (in radians).
        The default value has been empirically determined to give the best results for typical models.
        If ``None`` then the full feature geometry is used (ie, it is not split into ridge and transform segments, with the transform segments getting ignored).
    velocity_delta_time : float, default=1.0
        Delta time interval used to calculate spreading velocity. Defaults to 1 Myr.
    anchor_plate_id : int, optional
        Anchor plate ID for reconstruction.
        If not specified then uses the default anchor plate of ``rotation_model`` if it's a `pygplates.RotationModel`_ (otherwise uses zero).
    include_network_boundaries : bool, default=False
        Whether to calculate spreading rate along network boundaries that are not also plate boundaries (defaults to False).
        If a deforming network shares a boundary with a plate then it'll get included regardless of this option.
    output_obliquity_and_normal_and_left_right_plates : bool, default=False
        Whether to also return spreading obliquity, normal azimuth and left/right plates.

    Returns
    -------
    output_data : a list of tuples
        The results for all tessellated points sampled along spreading features.
        The size of the returned list is equal to the number of tessellated points.
        Each tuple in the list corresponds to a tessellated point and has the following tuple items
        (depending on ``output_obliquity_and_normal_and_left_right_plates``):

        If ``output_obliquity_and_normal_and_left_right_plates`` is ``False`` (the default):

        * longitude of tessellated point
        * latitude of tessellated point
        * spreading velocity magnitude (in cm/yr)
        * length of arc segment (in degrees) that tessellated point is on

        If ``output_obliquity_and_normal_and_left_right_plates`` is ``True``:

        * longitude of tessellated point
        * latitude of tessellated point
        * spreading velocity magnitude (in cm/yr)
        * spreading obliquity in degrees (deviation from normal line in range 0 to 90 degrees)
        * length of arc segment (in degrees) that tessellated point is on
        * azimuth of vector normal to the arc segment in degrees (clockwise starting at North, ie, 0 to 360 degrees)
        * left plate ID
        * right plate ID


    .. _pygplates.RotationModel: https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel
    .. _pygplates.Feature: https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature
    .. _pygplates.FeatureCollection: https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection
    .. _pygplates.FeatureType: https://www.gplates.org/docs/pygplates/generated/pygplates.featuretype
    """
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)

    # Turn topology data into a list of features (if not already).
    topology_features = pygplates.FeaturesFunctionArgument(topology_features)

    # If requested, exclude resolved topological *networks*.
    # Caller may only want mid-ocean ridges (and transforms) along *plate* boundaries.
    resolve_topology_types = pygplates.ResolveTopologyType.boundary
    if include_network_boundaries:
        resolve_topology_types = (
            resolve_topology_types | pygplates.ResolveTopologyType.network
        )

    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(
        topology_features.get_features(),
        rotation_model,
        resolved_topologies,
        time,
        shared_boundary_sections,
        anchor_plate_id,
        resolve_topology_types=resolve_topology_types,
    )

    # List of tesselated spreading points and associated spreading parameters for the current 'time'.
    output_data = []

    # If spreading feature types is a single feature type then convert to a list.
    # Othewise it should already be a list (or None).
    if isinstance(spreading_feature_types, pygplates.FeatureType):
        spreading_feature_types = [spreading_feature_types]

    # Iterate over the shared boundary sections of all resolved topologies.
    for shared_boundary_section in shared_boundary_sections:

        spreading_feature = shared_boundary_section.get_feature()

        # Skip sections that are not spreading features (if requested).
        if (
            spreading_feature_types is not None
            and spreading_feature.get_feature_type() not in spreading_feature_types
        ):
            continue

        # Find the stage rotation of the spreading feature in the frame of reference of its reconstructed geometry at the current 'time'.
        # The stage pole can then be directly geometrically compared to the reconstructed spreading geometry.
        spreading_stage_rotation = separate_ridge_transform_segments.get_stage_rotation_for_reconstructed_geometry(
            spreading_feature,
            rotation_model,
            time,
            return_left_right_plates=output_obliquity_and_normal_and_left_right_plates,
        )
        if not spreading_stage_rotation:
            # Skip current feature - it's not a spreading feature.
            continue

        if output_obliquity_and_normal_and_left_right_plates:
            spreading_stage_rotation, left_plate_id, right_plate_id = (
                spreading_stage_rotation
            )

        # Iterate over the shared sub-segments of the current line.
        # These are the parts of the line that actually contribute to topological boundaries.
        for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
            spreading_feature_resolved_geometry = (
                shared_sub_segment.get_resolved_geometry()
            )

            if transform_segment_deviation_in_radians is not None:
                # Split into ridge and transform segments.
                ridge_and_transform_segment_geometries = separate_ridge_transform_segments.separate_geometry_into_ridges_and_transforms(
                    spreading_stage_rotation,
                    spreading_feature_resolved_geometry,
                    transform_segment_deviation_in_radians,
                )
                if not ridge_and_transform_segment_geometries:
                    # Skip shared sub segment - it's not a polyline (or polygon).
                    continue

                # Only interested in ridge segments.
                ridge_sub_segment_geometries, _ = ridge_and_transform_segment_geometries

            else:
                # Not splitting into ridge and transform segments, so just the full unsplit geometry.
                ridge_sub_segment_geometries = [spreading_feature_resolved_geometry]

            # Ensure the ridge sub-segments are tessellated to within the threshold sampling distance.
            tessellated_shared_sub_segment_polylines = [
                ridge_sub_segment_geometry.to_tessellated(
                    threshold_sampling_distance_radians
                )
                for ridge_sub_segment_geometry in ridge_sub_segment_geometries
            ]

            # Iterate over the great circle arcs of the tessellated polylines to get the arc midpoints and lengths.
            # There is an arc between each adjacent pair of points in the polyline.
            arc_midpoints = []
            arc_lengths = []
            if output_obliquity_and_normal_and_left_right_plates:
                arc_normals = []
            for (
                tessellated_shared_sub_segment_polyline
            ) in tessellated_shared_sub_segment_polylines:
                for arc in tessellated_shared_sub_segment_polyline.get_segments():
                    if not arc.is_zero_length():
                        arc_midpoints.append(arc.get_arc_point(0.5))
                        arc_lengths.append(arc.get_arc_length())
                        if output_obliquity_and_normal_and_left_right_plates:
                            arc_normals.append(arc.get_great_circle_normal())

            # Shouldn't happen, but just in case ridge sub-segment polylines coincide with points.
            if not arc_midpoints:
                continue

            if output_obliquity_and_normal_and_left_right_plates:
                # The arc normals relative to North (azimuth).
                # Convert global 3D normal vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
                arc_local_normals = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                    arc_midpoints, arc_normals
                )

            # Calculate the spreading velocities at the arc midpoints.
            #
            # Note that the stage rotation can be used directly on the reconstructed geometries because
            # it is already in the frame of reference of the reconstructed geometries.
            spreading_velocity_vectors = pygplates.calculate_velocities(
                arc_midpoints,
                spreading_stage_rotation,
                velocity_delta_time,
                pygplates.VelocityUnits.cms_per_yr,
            )

            for arc_index in range(len(arc_midpoints)):

                arc_midpoint = arc_midpoints[arc_index]
                arc_length = arc_lengths[arc_index]
                lat, lon = arc_midpoint.to_lat_lon()

                spreading_velocity_vector = spreading_velocity_vectors[arc_index]
                spreading_velocity_magnitude = spreading_velocity_vector.get_magnitude()

                if output_obliquity_and_normal_and_left_right_plates:
                    arc_normal = arc_normals[arc_index]
                    arc_normal_azimuth = arc_local_normals[arc_index][1]

                    if spreading_velocity_vector.is_zero_magnitude():
                        spreading_obliquity_degrees = 0
                    else:
                        # Angle range [0, 180].
                        spreading_obliquity_degrees = math.degrees(
                            pygplates.Vector3D.angle_between(
                                spreading_velocity_vector, arc_normal
                            )
                        )
                        # Minimum deviation from 'arc_normal' and '-arc_normal'.
                        # Angle range [0, 90].
                        if spreading_obliquity_degrees > 90:
                            spreading_obliquity_degrees = (
                                180 - spreading_obliquity_degrees
                            )

                    # The data will be output in GMT format (ie, lon first, then lat, etc).
                    output_data.append(
                        (
                            lon,
                            lat,
                            spreading_velocity_magnitude,
                            spreading_obliquity_degrees,
                            math.degrees(arc_length),
                            math.degrees(arc_normal_azimuth),
                            left_plate_id,
                            right_plate_id,
                        )
                    )
                else:
                    # The data will be output in GMT format (ie, lon first, then lat, etc).
                    output_data.append(
                        (
                            lon,
                            lat,
                            spreading_velocity_magnitude,
                            math.degrees(arc_length),
                        )
                    )

    return output_data


def spreading_rates_dense(
    rotation_features_or_model,
    topology_features,
    time,
    threshold_sampling_distance_radians,
    spreading_feature_types=[pygplates.FeatureType.gpml_mid_ocean_ridge],
    transform_segment_deviation_in_radians=separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS,
    velocity_delta_time=1.0,
    anchor_plate_id=None,
    include_network_boundaries=False,
):
    """
    Equivalent to `spreading_rates(..., output_obliquity_and_normal_and_left_right_plates=True)`.

    See Also
    --------
    spreading_rates
    """

    return spreading_rates(
        rotation_features_or_model,
        topology_features,
        time,
        threshold_sampling_distance_radians,
        spreading_feature_types=spreading_feature_types,
        transform_segment_deviation_in_radians=transform_segment_deviation_in_radians,
        velocity_delta_time=velocity_delta_time,
        anchor_plate_id=anchor_plate_id,
        include_network_boundaries=include_network_boundaries,
        output_obliquity_and_normal_and_left_right_plates=True,
    )
