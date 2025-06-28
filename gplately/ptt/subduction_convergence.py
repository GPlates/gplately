#
#    Copyright (C) 2016-2025 The University of Sydney, Australia
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

"""Find the convergence rate of trenches (subduction zones)."""

from __future__ import print_function

import argparse
import math
import os
import sys
import warnings

import pygplates

# The default threshold sampling distance along trenches (subduction zones).
DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES = 0.5
DEFAULT_THRESHOLD_SAMPLING_DISTANCE_KMS = (
    math.radians(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES)
    * pygplates.Earth.equatorial_radius_in_kms
)

DEFAULT_TIME_RANGE_YOUNG_TIME = 0
DEFAULT_TIME_RANGE_OLD_TIME = 200
DEFAULT_TIME_INCREMENT = 1

DEFAULT_VELOCITY_DELTA_TIME = 1


# Determine the subducting plate of the subduction shared sub-segment.
#
# Note: There is now a similar method in PyGPlates version 30 called pygplates.ResolvedTopologicalSharedSubSegment.get_subducting_plate().
def find_subducting_plate(
    subduction_shared_sub_segment,
    include_slab_topologies=False,
):
    """Determine the subducting plate of the subduction shared sub-segment."""
    # Get the subduction polarity of the subducting line.
    subduction_polarity = subduction_shared_sub_segment.get_feature().get_enumeration(
        pygplates.PropertyName.gpml_subduction_polarity
    )
    if (not subduction_polarity) or (subduction_polarity == "Unknown"):
        return None

    subducting_plate = None

    # Iterate over the resolved topologies sharing the subduction sub-segment.
    # We are looking for exactly one subducting plate, or one boundary and
    # one network.
    #
    # There can be zero, one or more overriding plates but that does not affect us (since only looking for subducting plate).
    # This actually makes things more robust because it's possible the topologies were built in such a way that a subduction line
    # is inadvertently duplicated such that the subducting plate uses one of the subduction lines as its boundary and the overriding plate
    # uses the other. In this case the subducting line attached to the subducting plate will not also be attached to the overriding plate
    # and hence there will be zero overriding plates here.
    # Another example is having two overriding plates (or at least there will be two topologies on the overriding side of the subduction line).
    #
    # So all these overriding cases do not affect us, which means we will actually get a more accurate total subduction zone length in this script
    # because we are not forced to ignore these overriding cases (normally we would be forced to find *one* overriding plate if we were looking for
    # both the subducting and overriding plates). And also we're not counting duplicate subduction lines because we only count one of the duplicate
    # subduction lines (the one attached to the subducting plate). However we will still have a problem if too many subducting plates are found.
    sharing_resolved_topologies = (
        subduction_shared_sub_segment.get_sharing_resolved_topologies()
    )
    geometry_reversal_flags = (
        subduction_shared_sub_segment.get_sharing_resolved_topology_geometry_reversal_flags()
    )

    def key(t):
        """Boundaries first, then networks - boundaries will be replaced by networks."""
        if isinstance(t[0], pygplates.ResolvedTopologicalBoundary):
            return 0
        if isinstance(t[0], pygplates.ResolvedTopologicalNetwork):
            return 1
        return 2

    zipped = list(zip(sharing_resolved_topologies, geometry_reversal_flags))
    zipped.sort(key=key)
    sharing_resolved_topologies, geometry_reversal_flags = zip(*zipped)

    n_subducting_plates = 0
    subducting_topology_types = []
    for (
        sharing_resolved_topology,
        geometry_reversal_flag,
    ) in zip(sharing_resolved_topologies, geometry_reversal_flags):
        if include_slab_topologies and (
            sharing_resolved_topology.get_feature()
            .get_feature_type()
            .to_qualified_string()
            == "gpml:TopologicalSlabBoundary"
        ):
            # Ignore slab topologies (e.g. flat slabs)
            continue
        if (
            sharing_resolved_topology.get_resolved_boundary().get_orientation()
            == pygplates.PolygonOnSphere.Orientation.clockwise
        ):
            # The current topology sharing the subducting line has clockwise orientation (when viewed from above the Earth).
            # If the overriding plate (subduction polarity) is to the 'left' of the subducting line (when following its vertices in order)
            # and the subducting line is not reversed when contributing to the topology then that topology is the subducting plate.
            # A similar test applies to the 'right' but with the subducting line reversed in the topology.
            if (subduction_polarity == "Left" and not geometry_reversal_flag) or (
                subduction_polarity == "Right" and geometry_reversal_flag
            ):
                n_subducting_plates += 1
                subducting_topology_types.append(type(sharing_resolved_topology))
                # Make sure this topology actually has a plate ID
                if (
                    sharing_resolved_topology.get_feature().get_reconstruction_plate_id()
                    is None
                ):
                    continue
                subducting_plate = sharing_resolved_topology
        else:
            # The current topology sharing the subducting line has counter-clockwise orientation (when viewed from above the Earth).
            # If the overriding plate (subduction polarity) is to the 'left' of the subducting line (when following its vertices in order)
            # and the subducting line is reversed when contributing to the topology then that topology is the subducting plate.
            # A similar test applies to the 'right' but with the subducting line not reversed in the topology.
            if (subduction_polarity == "Left" and geometry_reversal_flag) or (
                subduction_polarity == "Right" and not geometry_reversal_flag
            ):
                n_subducting_plates += 1
                subducting_topology_types.append(type(sharing_resolved_topology))
                # Make sure this topology actually has a plate ID
                if (
                    sharing_resolved_topology.get_feature().get_reconstruction_plate_id()
                    is None
                ):
                    continue
                subducting_plate = sharing_resolved_topology

    if subducting_plate is None or n_subducting_plates > 2:
        # Unable to find subducting plate, so return None.
        return None
    if len(subducting_topology_types) != len(set(subducting_topology_types)):
        # More than one rigid plate or more than one deforming network
        return None

    return (subducting_plate, subduction_polarity)


def subduction_convergence(
    rotation_features_or_model,
    topology_features,
    threshold_sampling_distance_radians,
    time,
    velocity_delta_time=1.0,
    anchor_plate_id=None,
    include_slab_topologies=False,
    include_network_boundaries=False,
    **kwargs,
):
    """Find the convergence and absolute velocities sampled along trenches (subduction zones) at a particular geological time.

    Each sampled point along trench returns the following information:

    - 0 longitude of sample point
    - 1 latitude of sample point
    - 2 subducting convergence (relative to trench) velocity magnitude (in cm/yr)
    - 3 subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
    - 4 trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
    - 5 trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
    - 6 length of arc segment (in degrees) that current point is on
    - 7 trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
    - 8 subducting plate ID
    - 9 trench plate ID

    The optional ``output_*`` parameters can be used to append extra data to the above output for each sampled trench point.
    The order of any extra data is the same order in which the parameters are listed below.

    The obliquity angles are in the range (-180, 180). The range (0, 180) goes clockwise (when viewed from above the Earth) from the
    trench normal direction to the velocity vector. The range (0, -180) goes counter-clockwise.
    You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.
    The trench normal is perpendicular to the trench and pointing toward the overriding plate.

    Note that the convergence velocity magnitude is negative if the plates are diverging (if convergence obliquity angle
    is greater than 90 or less than -90). And note that the absolute velocity magnitude is negative if the trench (subduction zone)
    is moving towards the overriding plate (if absolute obliquity angle is less than 90 and greater than -90) - note that this
    ignores the kinematics of the subducting plate. Similiarly for the subducting plate absolute velocity magnitude
    (if keyword argument ``output_subducting_absolute_velocity`` is ``True``).

    Parameters
    ----------
    rotation_features_or_model : `pygplates.RotationModel`_, or any combination of str, `pygplates.FeatureCollection`_, `pygplates.Feature`_
        The rotation model can be specified as a `pygplates.RotationModel`_. Or it can be specified as a rotation feature collection,
        or rotation filename, or rotation feature, or sequence of rotation features, or a sequence (eg, list or tuple) of any combination
        of those four types.
    topology_features: any combination of str, `pygplates.FeatureCollection`_, `pygplates.Feature`_
        The topological boundary and network features and the topological section features they reference (regular and topological lines).
        Can be specified as a feature collection, or filename, or feature, or sequence of features, or a sequence (eg, list or tuple)
        of any combination of those four types.
    threshold_sampling_distance_radians: float
        Threshold sampling distance along trench (in radians).
    time: float
        The reconstruction time at which to query subduction convergence.
    velocity_delta_time: float, optional
        The delta time interval used for velocity calculations. Defaults to 1My.
    anchor_plate_id: int, optional
        Anchor plate ID for reconstruction.
        If not specified then uses the default anchor plate of ``rotation_model`` if it's a `pygplates.RotationModel`_ (otherwise uses zero).
    include_slab_topologies : bool, default False
        Include slab topologies (``gpml:TopologicalSlabBoundary``) in analysis.
    include_network_boundaries : bool, default False
        Whether to calculate subduction convergence along network boundaries that are not also plate boundaries (defaults to False).
        If a deforming network shares a boundary with a plate then it'll get included regardless of this option.
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
    list of tuples
        The results for all points sampled along trench.
        The size of the returned list is equal to the number of sampled points.
        Each tuple in the list corresponds to a point and has the following tuple items:

        * longitude of sample point
        * latitude of sample point
        * subducting convergence (relative to trench) velocity magnitude (in cm/yr)
        * subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
        * trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
        * trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
        * length of arc segment (in degrees) that current point is on
        * trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
        * subducting plate ID
        * trench plate ID

        The optional ``output_*`` parameters can be used to append extra data to the tuple of each sampled trench point.
        The order of any extra data is the same order in which the parameters are listed in this function.


    .. note::

        Each point in the output is the midpoint of a great circle arc between two adjacent points in the trench polyline.
        The trench normal vector used in the obliquity calculations is perpendicular to the great circle arc of each point (arc midpoint)
        and pointing towards the overriding plate (rather than away from it).

        Each trench is sampled at approximately uniform intervals along its length (specified via a threshold sampling distance).
        The sampling along the entire length of a trench is not exactly uniform. Each segment along a trench is sampled
        such that the samples have a uniform spacing that is less than or equal to the threshold sampling distance. However each segment
        in a trench might have a slightly different spacing distance (since segment lengths are not integer multiples of the threshold sampling distance).

        The trench normal (at each arc segment mid-point) always points **towards** the overriding plate.


    .. _pygplates.RotationModel: https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel
    .. _pygplates.Feature: https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature
    .. _pygplates.FeatureCollection: https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection
    .. _pygplates.FeatureType: https://www.gplates.org/docs/pygplates/generated/pygplates.featuretype
    """
    time = float(time)

    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)

    # Turn topology data into a list of features (if not already).
    topology_features = pygplates.FeaturesFunctionArgument(
        topology_features
    ).get_features()
    if not include_slab_topologies:
        # Ignore slab topologies (usually flat slabs)
        topology_features = [
            i
            for i in topology_features
            if i.get_feature_type().to_qualified_string()
            != "gpml:TopologicalSlabBoundary"
        ]

    # If requested, exclude resolved topological *networks*.
    # Caller may only want subduction zones along *plate* boundaries.
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
        topology_features,
        rotation_model,
        resolved_topologies,
        time,
        shared_boundary_sections,
        anchor_plate_id,
        resolve_topology_types=resolve_topology_types,
    )

    # List of tesselated subduction zone (trench) shared subsegment points and associated convergence parameters
    # for the current 'time'.
    output_data = []

    # Iterate over the shared boundary sections of all resolved topologies.
    for shared_boundary_section in shared_boundary_sections:
        # Skip sections that are not subduction zones (trenches).
        if (
            shared_boundary_section.get_feature().get_feature_type()
            != pygplates.FeatureType.gpml_subduction_zone
        ):
            continue

        # The trench length is the sum of the lengths of the shared sub-segments.
        # The shared sub-segments represent portions of the trench with different subducting plates.
        # We need the trench length to determine whether a point along the trench is closer to the start or end of trench.
        trench_length_radians = math.fsum(
            shared_sub_segment.get_resolved_geometry().get_arc_length()
            for shared_sub_segment in shared_boundary_section.get_shared_sub_segments()
        )

        # The distance-along-trench will accumulate as we traverse the shared sub-segments.
        distance_along_trench_radians = 0.0

        # Iterate over the shared sub-segments of the current subducting line.
        # These are the parts of the subducting line that actually contribute to topological boundaries.
        for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
            # Find the subducting plate of the shared sub-segment.
            #
            # Note that prior to pyGPlates 0.22 we also looked for the overriding plate since it couldn't extract the individual trench plate IDs
            # from a "resolved topological line" trench (now we use the trench plate IDs since our minimum requirement is above version 0.28).
            # Not having to find the overriding plate means we actually get a more accurate total subduction zone length in this script.
            # This is because we are not forced to ignore trench sections where there's not exactly one overriding plate
            # (and optionally a deforming network overlapping it). And also we're not counting duplicate subduction lines
            # (where one duplicate is attached only to the overriding plate and the other attached only to the subducting plate)
            # because we only count the subduction line attached to the subducting plate.
            subducting_plate_and_polarity = find_subducting_plate(
                shared_sub_segment,
                include_slab_topologies=include_slab_topologies,
            )
            if not subducting_plate_and_polarity:
                warnings.warn(
                    'Unable to find the subducting plate of the subducting sub-segment "{0}" at {1}Ma.\n'
                    "    Either the subduction polarity is not properly set or there"
                    " are too many subducting plates sharing the sub-segment.\n"
                    "    Ignoring current sub-segment.".format(
                        shared_sub_segment.get_feature().get_name(), time
                    ),
                    RuntimeWarning,
                )
                continue
            subducting_plate, subduction_polarity = subducting_plate_and_polarity
            subducting_plate_id = (
                subducting_plate.get_feature().get_reconstruction_plate_id()
            )

            # We need to reverse the trench normal direction if overriding plate is to
            # the right of the subducting line since great circle arc normal is always to the left.
            if subduction_polarity == "Left":
                trench_normal_reversal = 1
            else:
                trench_normal_reversal = -1

            # The plate ID of the trench line (as opposed to the subducting plate).
            #
            # Update: The plate IDs of the trench line and overriding plate can differ
            # even in a non-deforming model due to smaller plates, not modelled by topologies, moving
            # differently than the larger topological plate being modelled - and the trench line
            # having plate IDs of the smaller plates near them. For that reason we use the plate ID
            # of the trench line whenever we can.
            #
            # If the current shared sub-segment is part of a topological line then we obtain
            # its sub-sub-segments. This is because trench lines that are topological lines might
            # actually be deforming (or intended to be deforming) and hence their
            # plate ID is not meaningful or at least we can't be sure whether it will be zero or the
            # overriding plate (or something else). In this case we look at the plate IDs of the
            # sub-sub-segments.
            sub_segments_of_topological_line_sub_segment = (
                shared_sub_segment.get_sub_segments()
            )
            if sub_segments_of_topological_line_sub_segment:
                # Iterate over the sub-sub-segments associated with the topological line shared sub-segment.
                for sub_sub_segment in sub_segments_of_topological_line_sub_segment:
                    trench_plate_id = (
                        sub_sub_segment.get_feature().get_reconstruction_plate_id()
                    )

                    sub_sub_segment_geometry = sub_sub_segment.get_resolved_geometry()
                    # If sub-sub-segment was reversed when it contributed to the topological line shared sub-segment then
                    # we need to use that reversed geometry so that it has the same order of points as the topological line.
                    if sub_sub_segment.was_geometry_reversed_in_topology():
                        # Create a new sub-sub-segment polyline with points in reverse order.
                        sub_sub_segment_geometry = pygplates.PolylineOnSphere(
                            sub_sub_segment_geometry[::-1]
                        )

                    _sub_segment_subduction_convergence(
                        output_data,
                        time,
                        sub_sub_segment_geometry,
                        trench_plate_id,
                        subducting_plate_id,
                        trench_normal_reversal,
                        trench_length_radians,
                        distance_along_trench_radians,
                        threshold_sampling_distance_radians,
                        velocity_delta_time,
                        rotation_model,
                        anchor_plate_id,
                        **kwargs,
                    )

                    # Accumulate distance-along-trench.
                    distance_along_trench_radians += (
                        sub_sub_segment_geometry.get_arc_length()
                    )

            else:  # It's not a topological line...
                trench_plate_id = (
                    shared_sub_segment.get_feature().get_reconstruction_plate_id()
                )
                sub_segment_geometry = shared_sub_segment.get_resolved_geometry()
                _sub_segment_subduction_convergence(
                    output_data,
                    time,
                    sub_segment_geometry,
                    trench_plate_id,
                    subducting_plate_id,
                    trench_normal_reversal,
                    trench_length_radians,
                    distance_along_trench_radians,
                    threshold_sampling_distance_radians,
                    velocity_delta_time,
                    rotation_model,
                    anchor_plate_id,
                    **kwargs,
                )

                # Accumulate distance-along-trench by length of sub-segment geometry.
                distance_along_trench_radians += sub_segment_geometry.get_arc_length()

    return output_data


def _sub_segment_subduction_convergence(
    output_data,
    time,
    sub_segment_geometry,
    trench_plate_id,
    subducting_plate_id,
    trench_normal_reversal,
    trench_length_radians,
    distance_along_trench_radians,
    threshold_sampling_distance_radians,
    velocity_delta_time,
    rotation_model,
    anchor_plate_id,
    **kwargs,
):
    #
    # Process keyword arguments.
    #
    output_distance_to_nearest_edge_of_trench = kwargs.get(
        "output_distance_to_nearest_edge_of_trench", False
    )
    output_distance_to_start_edge_of_trench = kwargs.get(
        "output_distance_to_start_edge_of_trench", False
    )
    output_convergence_velocity_components = kwargs.get(
        "output_convergence_velocity_components", False
    )
    output_trench_absolute_velocity_components = kwargs.get(
        "output_trench_absolute_velocity_components", False
    )
    output_subducting_absolute_velocity = kwargs.get(
        "output_subducting_absolute_velocity", False
    )
    output_subducting_absolute_velocity_components = kwargs.get(
        "output_subducting_absolute_velocity_components", False
    )
    output_trench_normal = kwargs.get("output_trench_normal", False)

    # Get the rotation of the subducting plate relative to the trench line
    # from 'time + velocity_delta_time' to 'time'.
    convergence_relative_stage_rotation = rotation_model.get_rotation(
        time,
        subducting_plate_id,
        time + velocity_delta_time,
        trench_plate_id,
        anchor_plate_id=anchor_plate_id,
    )
    #
    # In the following:
    #   * T is for Trench (subduction zone line)
    #   * S is subducting plate
    #   * A is anchor plate
    #
    # The trenches have been reconstructed using the rotation "R(0->t,A->T)":
    #
    #   reconstructed_geometry = R(0->t,A->T) * present_day_geometry
    #
    # We can write "R(0->t,A->T)" in terms of the convergence stage rotation "R(t+dt->t,T->S)" as:
    #
    #   R(0->t,A->T)  = R(0->t,A->S) * R(0->t,S->T)
    #                 = R(0->t,A->S) * inverse[R(0->t,T->S)]
    #                 = R(0->t,A->S) * inverse[R(t+dt->t,T->S) * R(0->t+dt,T->S)]
    #                 = R(0->t,A->S) * inverse[stage_rotation * R(0->t+dt,T->S)]
    #                 = R(0->t,A->S) * inverse[R(0->t+dt,T->S)] * inverse[stage_rotation]
    #                 = R(0->t,A->S) * R(0->t+dt,S->T) * inverse[stage_rotation]
    #
    # So to get the *reconstructed* subduction line geometry into the stage rotation reference frame
    # we need to rotate it by "inverse[R(0->t,A->S) * R(0->t+dt,S->T)]":
    #
    #   reconstructed_geometry = R(0->t,A->T) * present_day_geometry
    #                          = R(0->t,A->S) * R(0->t+dt,S->T) * inverse[stage_rotation] * present_day_geometry
    #   inverse[R(0->t,A->S) * R(0->t+dt,S->T)] * reconstructed_geometry = inverse[stage_rotation] * present_day_geometry
    #
    # Once we've done that we can calculate the velocities of those geometry points
    # using the stage rotation. Then the velocities need to be rotated back from the
    # stage rotation reference frame using the rotation "R(0->t,A->S) * R(0->t+dt,S->T)".
    #
    from_convergence_stage_frame = rotation_model.get_rotation(
        time, subducting_plate_id, anchor_plate_id=anchor_plate_id
    ) * rotation_model.get_rotation(
        time + velocity_delta_time,
        trench_plate_id,
        fixed_plate_id=subducting_plate_id,
        anchor_plate_id=anchor_plate_id,
    )
    to_convergence_stage_frame = from_convergence_stage_frame.get_inverse()

    # Get the rotation of the trench relative to the anchor plate
    # from 'time + velocity_delta_time' to 'time'.
    #
    # Note: We don't need to convert to and from the stage rotation reference frame
    # like the above convergence because...
    #
    #   R(0->t,A->T)  = R(t+dt->t,A->T) * R(0->t+dt,A->T)
    #
    #   reconstructed_geometry = R(0->t,A->T) * present_day_geometry
    #                          = R(t+dt->t,A->T) * R(0->t+dt,A->T) * present_day_geometry
    #
    # ...where *reconstructed* subduction line geometry is already in the frame of the stage rotation "R(t+dt->t,A->T)".
    trench_equivalent_stage_rotation = rotation_model.get_rotation(
        time,
        trench_plate_id,
        time + velocity_delta_time,
        anchor_plate_id=anchor_plate_id,
    )

    if (
        output_subducting_absolute_velocity
        or output_subducting_absolute_velocity_components
    ):
        # Get the rotation of the subducting plate relative to the anchor plate
        # from 'time + velocity_delta_time' to 'time'.
        #
        # Note: We don't need to convert to and from the stage rotation reference frame
        # like the above convergence because...
        #
        #   R(0->t,A->T)  = R(0->t,A->S) * R(0->t,S->T)
        #                 = R(t+dt->t,A->S) * R(0->t+dt,A->S) * R(0->t,S->T)
        #
        #   reconstructed_geometry = R(0->t,A->T) * present_day_geometry
        #                          = R(t+dt->t,A->S) * R(0->t+dt,A->S) * R(0->t,S->T) * present_day_geometry
        #
        # ...where *reconstructed* subduction line geometry is already in the frame of the stage rotation "R(t+dt->t,A->S)".
        subducting_equivalent_stage_rotation = rotation_model.get_rotation(
            time,
            subducting_plate_id,
            time + velocity_delta_time,
            anchor_plate_id=anchor_plate_id,
        )

    # Ensure the shared sub-segment is tessellated to within the threshold sampling distance.
    tessellated_shared_sub_segment_polyline = sub_segment_geometry.to_tessellated(
        threshold_sampling_distance_radians
    )

    # Iterate over the great circle arcs of the tessellated polyline to get the
    # arc midpoints, lengths and trench normals.
    # There is an arc between each adjacent pair of points in the polyline.
    arc_midpoints = []
    arc_lengths = []
    trench_normals = []
    for arc in tessellated_shared_sub_segment_polyline.get_segments():
        if not arc.is_zero_length():
            arc_midpoints.append(arc.get_arc_point(0.5))
            arc_lengths.append(arc.get_arc_length())
            # The normal to the trench in the direction of subduction (towards overriding plate).
            trench_normals.append(
                trench_normal_reversal * arc.get_great_circle_normal()
            )

    # Shouldn't happen, but just in case the shared sub-segment polyline coincides with a point.
    if not arc_midpoints:
        return

    # The trench normals relative to North (azimuth).
    # Convert global 3D normal vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
    trench_local_normals = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
        arc_midpoints, trench_normals
    )

    # Calculate the convergence velocities at the arc midpoints.
    #
    # Note; We need to convert the reconstructed geometry points into the convergence stage rotation
    # reference frame to calculate velocities and then convert the velocities using the
    # reverse transform as mentioned above.
    arc_midpoints_in_convergence_stage_frame = [
        to_convergence_stage_frame * arc_midpoint for arc_midpoint in arc_midpoints
    ]
    convergence_velocity_vectors_in_convergence_stage_frame = (
        pygplates.calculate_velocities(
            arc_midpoints_in_convergence_stage_frame,
            convergence_relative_stage_rotation,
            velocity_delta_time,
            pygplates.VelocityUnits.cms_per_yr,
        )
    )
    convergence_velocity_vectors = [
        from_convergence_stage_frame * velocity
        for velocity in convergence_velocity_vectors_in_convergence_stage_frame
    ]

    # Calculate the trench absolute velocities at the arc midpoints.
    trench_absolute_velocity_vectors = pygplates.calculate_velocities(
        arc_midpoints,
        trench_equivalent_stage_rotation,
        velocity_delta_time,
        pygplates.VelocityUnits.cms_per_yr,
    )

    if (
        output_subducting_absolute_velocity
        or output_subducting_absolute_velocity_components
    ):
        # Calculate the subducting absolute velocities at the arc midpoints.
        subducting_absolute_velocity_vectors = pygplates.calculate_velocities(
            arc_midpoints,
            subducting_equivalent_stage_rotation,
            velocity_delta_time,
            pygplates.VelocityUnits.cms_per_yr,
        )

    for arc_index in range(len(arc_midpoints)):
        arc_midpoint = arc_midpoints[arc_index]
        arc_length = arc_lengths[arc_index]
        trench_normal = trench_normals[arc_index]
        trench_normal_azimuth = trench_local_normals[arc_index][1]
        lat, lon = arc_midpoint.to_lat_lon()

        # The direction towards which we rotate from the trench normal in a clockwise fashion.
        clockwise_direction = pygplates.Vector3D.cross(
            trench_normal, arc_midpoint.to_xyz()
        )

        # Calculate the convergence rate parameters.
        convergence_velocity_vector = convergence_velocity_vectors[arc_index]
        if convergence_velocity_vector.is_zero_magnitude():
            convergence_velocity_magnitude = 0
            convergence_obliquity_degrees = 0
        else:
            convergence_velocity_magnitude = convergence_velocity_vector.get_magnitude()
            convergence_obliquity_degrees = math.degrees(
                pygplates.Vector3D.angle_between(
                    convergence_velocity_vector, trench_normal
                )
            )
            # Anti-clockwise direction has range (0, -180) instead of (0, 180).
            if (
                pygplates.Vector3D.dot(convergence_velocity_vector, clockwise_direction)
                < 0
            ):
                convergence_obliquity_degrees = -convergence_obliquity_degrees

            # See if plates are diverging (moving away from each other).
            # If plates are diverging (moving away from each other) then make the
            # velocity magnitude negative to indicate this. This could be inferred from
            # the obliquity but it seems this is the standard way to output convergence rate.
            if math.fabs(convergence_obliquity_degrees) > 90:
                convergence_velocity_magnitude = -convergence_velocity_magnitude

        # Calculate the trench absolute velocity magnitude and obliquity.
        trench_absolute_velocity_vector = trench_absolute_velocity_vectors[arc_index]
        if trench_absolute_velocity_vector.is_zero_magnitude():
            trench_absolute_velocity_magnitude = 0
            trench_absolute_obliquity_degrees = 0
        else:
            trench_absolute_velocity_magnitude = (
                trench_absolute_velocity_vector.get_magnitude()
            )
            trench_absolute_obliquity_degrees = math.degrees(
                pygplates.Vector3D.angle_between(
                    trench_absolute_velocity_vector, trench_normal
                )
            )
            # Anti-clockwise direction has range (0, -180) instead of (0, 180).
            if (
                pygplates.Vector3D.dot(
                    trench_absolute_velocity_vector, clockwise_direction
                )
                < 0
            ):
                trench_absolute_obliquity_degrees = -trench_absolute_obliquity_degrees

            # See if the trench absolute motion is heading in the direction of the
            # overriding plate. If it is then make the velocity magnitude negative to
            # indicate this. This could be inferred from the obliquity but it seems this
            # is the standard way to output trench velocity magnitude.
            #
            # Note that we are not calculating the motion of the trench
            # relative to the overriding plate - they are usually attached to each other
            # and hence wouldn't move relative to each other.
            if math.fabs(trench_absolute_obliquity_degrees) < 90:
                trench_absolute_velocity_magnitude = -trench_absolute_velocity_magnitude

        # Start with the standard tuple, and add extra data later (if requested).
        #
        # The data will be output in GMT format (ie, lon first, then lat, etc).
        output_tuple = (
            lon,
            lat,
            convergence_velocity_magnitude,
            convergence_obliquity_degrees,
            trench_absolute_velocity_magnitude,
            trench_absolute_obliquity_degrees,
            math.degrees(arc_length),
            math.degrees(trench_normal_azimuth),
            subducting_plate_id,
            trench_plate_id,
        )

        if (
            output_distance_to_nearest_edge_of_trench
            or output_distance_to_start_edge_of_trench
        ):
            # Increase by distance from previous segment mid-point to current segment mid-point.
            # Which is half previous segment length and half current segment length.
            if arc_index > 0:
                prev_arc_length = arc_lengths[arc_index - 1]
                distance_along_trench_radians += 0.5 * prev_arc_length
            distance_along_trench_radians += 0.5 * arc_length

            # Distance to nearest edge of the trench.
            if output_distance_to_nearest_edge_of_trench:
                if distance_along_trench_radians < 0.5 * trench_length_radians:
                    distance_to_nearest_edge_of_trench_radians = (
                        distance_along_trench_radians
                    )
                else:
                    distance_to_nearest_edge_of_trench_radians = (
                        trench_length_radians - distance_along_trench_radians
                    )

                output_tuple += (
                    math.degrees(distance_to_nearest_edge_of_trench_radians),
                )

            # Distance to start edge of the trench.
            if output_distance_to_start_edge_of_trench:
                # We want the distance to be along the clockwise direction around the overriding plate.
                if trench_normal_reversal < 0:
                    # The overriding plate is on the right of the trench.
                    # So the clockwise direction starts at the beginning of the trench.
                    distance_to_start_edge_of_trench_radians = (
                        distance_along_trench_radians
                    )
                else:
                    # The overriding plate is on the left of the trench.
                    # So the clockwise direction starts at the end of the trench.
                    distance_to_start_edge_of_trench_radians = (
                        trench_length_radians - distance_along_trench_radians
                    )

                output_tuple += (
                    math.degrees(distance_to_start_edge_of_trench_radians),
                )

        if output_convergence_velocity_components:
            # The orthogonal and parallel components are just magnitude multiplied by cosine and sine.
            convergence_velocity_orthogonal = math.cos(
                math.radians(convergence_obliquity_degrees)
            ) * math.fabs(convergence_velocity_magnitude)
            convergence_velocity_parallel = math.sin(
                math.radians(convergence_obliquity_degrees)
            ) * math.fabs(convergence_velocity_magnitude)
            output_tuple += (
                convergence_velocity_orthogonal,
                convergence_velocity_parallel,
            )

        if output_trench_absolute_velocity_components:
            # The orthogonal and parallel components are just magnitude multiplied by cosine and sine.
            trench_absolute_velocity_orthogonal = math.cos(
                math.radians(trench_absolute_obliquity_degrees)
            ) * math.fabs(trench_absolute_velocity_magnitude)
            trench_absolute_velocity_parallel = math.sin(
                math.radians(trench_absolute_obliquity_degrees)
            ) * math.fabs(trench_absolute_velocity_magnitude)
            output_tuple += (
                trench_absolute_velocity_orthogonal,
                trench_absolute_velocity_parallel,
            )

        if (
            output_subducting_absolute_velocity
            or output_subducting_absolute_velocity_components
        ):
            # Calculate the subducting absolute velocity magnitude and obliquity.
            subducting_absolute_velocity_vector = subducting_absolute_velocity_vectors[
                arc_index
            ]
            if subducting_absolute_velocity_vector.is_zero_magnitude():
                if output_subducting_absolute_velocity:
                    output_tuple += (0.0, 0.0)
                if output_subducting_absolute_velocity_components:
                    output_tuple += (0.0, 0.0)
            else:
                subducting_absolute_velocity_magnitude = (
                    subducting_absolute_velocity_vector.get_magnitude()
                )
                subducting_absolute_obliquity_degrees = math.degrees(
                    pygplates.Vector3D.angle_between(
                        subducting_absolute_velocity_vector, trench_normal
                    )
                )
                # Anti-clockwise direction has range (0, -180) instead of (0, 180).
                if (
                    pygplates.Vector3D.dot(
                        subducting_absolute_velocity_vector, clockwise_direction
                    )
                    < 0
                ):
                    subducting_absolute_obliquity_degrees = (
                        -subducting_absolute_obliquity_degrees
                    )

                # See if the subducting absolute motion is heading in the direction of the overriding plate.
                # If it is then make the velocity magnitude negative to indicate this.
                if math.fabs(subducting_absolute_obliquity_degrees) < 90:
                    subducting_absolute_velocity_magnitude = (
                        -subducting_absolute_velocity_magnitude
                    )

                if output_subducting_absolute_velocity:
                    output_tuple += (
                        subducting_absolute_velocity_magnitude,
                        subducting_absolute_obliquity_degrees,
                    )
                if output_subducting_absolute_velocity_components:
                    # The orthogonal and parallel components are just magnitude multiplied by cosine and sine.
                    subducting_absolute_velocity_orthogonal = math.cos(
                        math.radians(subducting_absolute_obliquity_degrees)
                    ) * math.fabs(subducting_absolute_velocity_magnitude)
                    subducting_absolute_velocity_parallel = math.sin(
                        math.radians(subducting_absolute_obliquity_degrees)
                    ) * math.fabs(subducting_absolute_velocity_magnitude)
                    output_tuple += (
                        subducting_absolute_velocity_orthogonal,
                        subducting_absolute_velocity_parallel,
                    )

        if output_trench_normal:
            output_tuple += trench_normal.to_xyz()

        output_data.append(output_tuple)


def write_output_file(output_filename, output_data):
    with open(output_filename, "w") as output_file:
        for output_line in output_data:
            output_file.write(" ".join(str(item) for item in output_line) + "\n")


def create_coverage_feature_from_convergence_data(
    subduction_convergence_data, time, **kwargs
):
    """Create a feature with a coverage geometry containing the calculated convergence and absolute velocity data.

    Parameters
    ----------
    subduction_convergence_data : list of tuples
        The subduction convergence data calculated by :func:`subduction_convergence`.
        Each tuple in the list contains the calculated data for a single sample point on a trench line.
    time: float
        The reconstruction time associated with the subduction convergence data.
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
    pygplates.Feature
        The feature with a coverage geometry containing the calculated convergence and absolute velocity data.

    Notes
    -----
    The optional 'output_*' parameters are the same as that of the `subduction_convergence` function , which is used to append extra data
    to the output of each sample point. Here it is used to interpret the same output data and create scalar coverages.
    """

    #
    # Process keyword arguments.
    #
    output_distance_to_nearest_edge_of_trench = kwargs.get(
        "output_distance_to_nearest_edge_of_trench", False
    )
    output_distance_to_start_edge_of_trench = kwargs.get(
        "output_distance_to_start_edge_of_trench", False
    )
    output_convergence_velocity_components = kwargs.get(
        "output_convergence_velocity_components", False
    )
    output_trench_absolute_velocity_components = kwargs.get(
        "output_trench_absolute_velocity_components", False
    )
    output_subducting_absolute_velocity = kwargs.get(
        "output_subducting_absolute_velocity", False
    )
    output_subducting_absolute_velocity_components = kwargs.get(
        "output_subducting_absolute_velocity_components", False
    )
    output_trench_normal = kwargs.get("output_trench_normal", False)

    # Convert the list of tuples (one tuple per sample point) into a tuple of lists (one list per data parameter).
    parameter_lists = list(zip(*subduction_convergence_data))

    # Put all convergence data for the current reconstruction time into a single feature.
    coverage_feature = pygplates.Feature()

    # Make it only appear at 'time'.
    coverage_feature.set_valid_time(time + 0.5, time - 0.5)

    # Extract the non-optional parameters.
    all_lon = parameter_lists[0]
    all_lat = parameter_lists[1]
    all_convergence_velocity_magnitude = parameter_lists[2]
    all_convergence_obliquity_degrees = parameter_lists[3]
    all_trench_absolute_velocity_magnitude = parameter_lists[4]
    all_trench_absolute_obliquity_degrees = parameter_lists[5]
    all_subducting_length_degrees = parameter_lists[6]
    all_trench_normal_azimuth_degrees = parameter_lists[7]
    all_subducting_plate_id = parameter_lists[8]
    all_trench_plate_id = parameter_lists[9]

    # Add each data parameter as a separate scalar coverage.
    coverage_geometry = pygplates.MultiPointOnSphere(zip(all_lat, all_lon))
    coverage_scalars = {
        pygplates.ScalarType.create_gpml(
            "ConvergenceVelocityMagnitude"
        ): all_convergence_velocity_magnitude,
        pygplates.ScalarType.create_gpml(
            "ConvergenceObliquityDegrees"
        ): all_convergence_obliquity_degrees,
        pygplates.ScalarType.create_gpml(
            "TrenchAbsoluteVelocityMagnitude"
        ): all_trench_absolute_velocity_magnitude,
        pygplates.ScalarType.create_gpml(
            "TrenchAbsoluteObliquityDegrees"
        ): all_trench_absolute_obliquity_degrees,
        pygplates.ScalarType.create_gpml(
            "SubductingLengthDegrees"
        ): all_subducting_length_degrees,
        pygplates.ScalarType.create_gpml(
            "TrenchNormalAzimuthDegrees"
        ): all_trench_normal_azimuth_degrees,
        pygplates.ScalarType.create_gpml("SubductingPlateId"): all_subducting_plate_id,
        pygplates.ScalarType.create_gpml("TrenchPlateId"): all_trench_plate_id,
    }

    # Extract the optional parameters.
    optional_parameter_index = 10

    if output_distance_to_nearest_edge_of_trench:
        all_distance_to_nearest_edge_of_trench_degrees = parameter_lists[
            optional_parameter_index
        ]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("DistanceToNearestEdgeOfTrenchDegrees")
        ] = all_distance_to_nearest_edge_of_trench_degrees
        optional_parameter_index += 1
    if output_distance_to_start_edge_of_trench:
        all_distance_to_start_edge_of_trench_degrees = parameter_lists[
            optional_parameter_index
        ]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("DistanceToStartEdgeOfTrenchDegrees")
        ] = all_distance_to_start_edge_of_trench_degrees
        optional_parameter_index += 1
    if output_convergence_velocity_components:
        all_convergence_velocity_orthogonal = parameter_lists[optional_parameter_index]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("ConvergenceVelocityOrthogonal")
        ] = all_convergence_velocity_orthogonal
        optional_parameter_index += 1
        all_convergence_velocity_parallel = parameter_lists[optional_parameter_index]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("ConvergenceVelocityParallel")
        ] = all_convergence_velocity_parallel
        optional_parameter_index += 1
    if output_trench_absolute_velocity_components:
        all_trench_absolute_velocity_orthogonal = parameter_lists[
            optional_parameter_index
        ]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("TrenchAbsoluteVelocityOrthogonal")
        ] = all_trench_absolute_velocity_orthogonal
        optional_parameter_index += 1
        all_trench_absolute_velocity_parallel = parameter_lists[
            optional_parameter_index
        ]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("TrenchAbsoluteVelocityParallel")
        ] = all_trench_absolute_velocity_parallel
        optional_parameter_index += 1
    if output_subducting_absolute_velocity:
        all_subducting_absolute_velocity_magnitude = parameter_lists[
            optional_parameter_index
        ]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("SubductingAbsoluteVelocityMagnitude")
        ] = all_subducting_absolute_velocity_magnitude
        optional_parameter_index += 1
        all_subducting_absolute_obliquity_degrees = parameter_lists[
            optional_parameter_index
        ]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("SubductingAbsoluteObliquityDegrees")
        ] = all_subducting_absolute_obliquity_degrees
        optional_parameter_index += 1
    if output_subducting_absolute_velocity_components:
        all_subducting_absolute_velocity_orthogonal = parameter_lists[
            optional_parameter_index
        ]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("SubductingAbsoluteVelocityOrthogonal")
        ] = all_subducting_absolute_velocity_orthogonal
        optional_parameter_index += 1
        all_subducting_absolute_velocity_parallel = parameter_lists[
            optional_parameter_index
        ]
        coverage_scalars[
            pygplates.ScalarType.create_gpml("SubductingAbsoluteVelocityParallel")
        ] = all_subducting_absolute_velocity_parallel
        optional_parameter_index += 1
    if output_trench_normal:
        all_trench_normal_x = parameter_lists[optional_parameter_index]
        coverage_scalars[pygplates.ScalarType.create_gpml("TrenchNormalX")] = (
            all_trench_normal_x
        )
        optional_parameter_index += 1
        all_trench_normal_y = parameter_lists[optional_parameter_index]
        coverage_scalars[pygplates.ScalarType.create_gpml("TrenchNormalY")] = (
            all_trench_normal_y
        )
        optional_parameter_index += 1
        all_trench_normal_z = parameter_lists[optional_parameter_index]
        coverage_scalars[pygplates.ScalarType.create_gpml("TrenchNormalZ")] = (
            all_trench_normal_z
        )
        optional_parameter_index += 1

    coverage_feature.set_geometry((coverage_geometry, coverage_scalars))

    return coverage_feature


def subduction_convergence_over_time(
    output_filename_prefix,
    output_filename_extension,
    rotation_filenames,
    topology_filenames,
    threshold_sampling_distance_radians,
    time_young,
    time_old,
    time_increment,
    velocity_delta_time=1.0,
    anchor_plate_id=None,
    output_gpml_filename=None,
    include_slab_topologies=False,
    **kwargs,
):
    """calculate subduction convergence over time"""
    if time_increment <= 0:
        raise ValueError(
            'The time increment "{0}" is not positive and non-zero.'.format(
                time_increment
            )
        )

    if time_young > time_old:
        raise ValueError(
            "The young time {0} is older (larger) than the old time {1}.".format(
                time_young, time_old
            )
        )

    rotation_model = pygplates.RotationModel(rotation_filenames)

    # Read/parse the topological features once so we're not doing at each time iteration.
    topology_features = [
        pygplates.FeatureCollection(topology_filename)
        for topology_filename in topology_filenames
    ]

    if output_gpml_filename:
        coverage_features = []

    # Iterate over the time rage.
    time = time_young
    while time <= pygplates.GeoTimeInstant(time_old):
        # print('Time {0}'.format(time))

        # Returns a list of tesselated trench points and associated convergence parameters
        # to write to the output file for the current 'time'.
        output_data = subduction_convergence(
            rotation_model,
            topology_features,
            threshold_sampling_distance_radians,
            time,
            velocity_delta_time,
            anchor_plate_id,
            include_slab_topologies=include_slab_topologies,
            **kwargs,
        )

        if output_data:
            output_filename = "{0}_{1:0.2f}.{2}".format(
                output_filename_prefix, time, output_filename_extension
            )
            write_output_file(output_filename, output_data)

            # Also keep track of convergence data if we need to write out a GPML file.
            if output_gpml_filename:
                coverage_feature = create_coverage_feature_from_convergence_data(
                    output_data, time, **kwargs
                )
                coverage_features.append(coverage_feature)

        # Increment the time further into the past.
        time += time_increment

    if output_gpml_filename:
        # Write out all coverage features to a single GPML file.
        pygplates.FeatureCollection(coverage_features).write(output_gpml_filename)

    return 0  # Success


# Convert mm/yr to cm/yr.
_MM_YR_TO_CM_YR = 0.1

# Earth mean radius in kms.
_EARTH_RADIUS_KMS = 6371.009


def convert_old_convergence_output(old_convergence_data, **kwargs):
    """Converts the output of the old "convergence.py" script to the current output format of the `subduction_convergence` function.

    This function's purpose is to help transition workflows that use "convergence.py" to use this script instead.
    At the time this function was implemented, both scripts have essentially the same output with just minor format differences (as seen below).
    However in the future this script will generate additional outputs (such as deformation strain rates) that the old script will not.

    Parameters
    ----------
    old_convergence_data : list of tuples
        The subduction convergence data output from the old "convergence.py" script.
        Each tuple in the list contains the data for a single sample point on a trench line.
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
    list of tuples
        The converted results for all input points from the old "convergence.py" script.
        The size of the returned list is equal to the number of input points.
        Each tuple in the list corresponds to a point and has the following tuple items:

        * longitude of sample point
        * latitude of sample point
        * subducting convergence (relative to trench) velocity magnitude (in cm/yr)
        * subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
        * trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
        * trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
        * length of arc segment (in degrees) that current point is on
        * trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
        * subducting plate ID
        * trench plate ID

        The optional 'output_*' parameters can be used to append extra data to the tuple of each sampled trench point.
        The order of any extra data is the same order in which the parameters are listed in this function.

    Notes
    -----
    The optional 'output_*' parameters are the same as that of the `subduction_convergence` function , which is used to append extra data
    to the output of each sample point.
    """

    #
    # Process keyword arguments.
    #
    output_distance_to_nearest_edge_of_trench = kwargs.get(
        "output_distance_to_nearest_edge_of_trench", False
    )
    output_distance_to_start_edge_of_trench = kwargs.get(
        "output_distance_to_start_edge_of_trench", False
    )
    output_convergence_velocity_components = kwargs.get(
        "output_convergence_velocity_components", False
    )
    output_trench_absolute_velocity_components = kwargs.get(
        "output_trench_absolute_velocity_components", False
    )
    output_subducting_absolute_velocity = kwargs.get(
        "output_subducting_absolute_velocity", False
    )
    output_subducting_absolute_velocity_components = kwargs.get(
        "output_subducting_absolute_velocity_components", False
    )
    output_trench_normal = kwargs.get("output_trench_normal", False)

    new_convergence_data = []
    for old_convergence_data_sample in old_convergence_data:
        (
            # Longitude of mid-point of segment (degrees).
            old_longitude,
            # Latitude of mid-point of segment (degrees).
            old_latitude,
            # Convergence velocity magnitude (mm/yr).
            old_convRate,
            # Length of segment containing the sample mid-point (km).
            old_segmentLength,
            # Subducting plate velocity orthogonal to the trench in direction of overriding plate (mm/yr).
            old_orthAbs,
            # Overriding plate velocity orthogonal to the trench in direction of overriding plate (mm/yr).
            old_orthOP,
            # Trench velocity orthogonal to the trench in direction of overriding plate (mm/yr).
            old_orthTrench,
            # Convergence obliquity (degrees) varies from 0 to 90 degrees (regardless of diverging/converging).
            old_subObliquity,
            # Subduction line polarity (degrees).
            # Angle from North (-180 to 180 degrees) along subduction line
            # in order of its points (ie, first to last).
            # Note that the subduction line's points are reversed (if necessary)
            # so that the overriding plate is always on its right.
            old_subPolarity,
            # Distance (km) to start of subduction line, where start is the first point along
            # subduction line where subducting/overriding plate IDs are same as current location.
            # Note that the subduction line's points are reversed (if necessary)
            # so that the overriding plate is always on its right.
            old_distEdge,
            # Convergence velocity orthogonal to the trench in direction of overriding plate (mm/yr).
            old_convPerp,
            # Convergence velocity parallel to the trench (mm/yr) along subduction line
            # in order of its points (ie, first to last).
            # Note that the subduction line's points are reversed (if necessary)
            # so that the overriding plate is always on its right.
            old_convPar,
            # Subducting plate velocity parallel to the trench (mm/yr) along subduction line
            # in order of its points (ie, first to last).
            # Note that the subduction line's points are reversed (if necessary)
            # so that the overriding plate is always on its right.
            old_parAbs,
            # Overriding plate velocity parallel to the trench (mm/yr) along subduction line
            # in order of its points (ie, first to last).
            # Note that the subduction line's points are reversed (if necessary)
            # so that the overriding plate is always on its right.
            old_parOP,
            # Trench velocity parallel to the trench (mm/yr) along subduction line
            # in order of its points (ie, first to last).
            # Note that the subduction line's points are reversed (if necessary)
            # so that the overriding plate is always on its right.
            old_parTrench,
            # Distance (km) to edge of subduction line, where edge is nearest location along
            # subduction line where subducting/overriding plate IDs differ from current location.
            old_distEdgeTotal,
            # Subducting plate ID.
            old_SPid,
            # Trench plate ID.
            old_TRENCHid,
            # Overriding plate ID.
            old_OPid,
        ) = old_convergence_data_sample

        new_longitude = old_longitude
        new_latitude = old_latitude

        # Convergence velocity magnitude.
        new_convergence_velocity_magnitude = _MM_YR_TO_CM_YR * old_convRate
        # The convergence velocity magnitude is negative if the plates are diverging.
        if old_convPerp < 0:
            new_convergence_velocity_magnitude = -new_convergence_velocity_magnitude

        #
        # The vertical line is subduction line in old convergence (normal is always to the right).
        # The new obliquity goes clockwise (0 to 180) and counter-clockwise (0 to -180) from the normal
        # (which points in the direction of the overriding plate).
        # This is atan2(-y, x) where x is component of convergence velocity in the normal direction, and
        # y is component of convergence along old subduction line.
        #
        #                 --> normal
        #                ^\    positive
        #                | \  obliquity
        #                |  \
        #
        #             <-- normal
        #    negative   /^
        #   obliquity  / |
        #             /  |
        #
        new_convergence_obliquity_degrees = math.degrees(
            math.atan2(-old_convPar, old_convPerp)
        )

        # Trench absolute velocity magnitude obtained parallel and orthogonal components.
        new_trench_absolute_velocity_magnitude = _MM_YR_TO_CM_YR * math.sqrt(
            old_parTrench * old_parTrench + old_orthTrench * old_orthTrench
        )
        # The trench absolute velocity magnitude is negative if trench is moving toward the overriding plate.
        if old_orthTrench > 0:
            new_trench_absolute_velocity_magnitude = (
                -new_trench_absolute_velocity_magnitude
            )

        # Same reasoning as for convergence obliquity above.
        new_trench_absolute_obliquity_degrees = math.degrees(
            math.atan2(-old_parTrench, old_orthTrench)
        )

        # Convert distance in kms to degrees.
        new_arc_length_degrees = math.degrees(old_segmentLength / _EARTH_RADIUS_KMS)

        # The old subduction polarity points along the subduction line, and since the overriding plate
        # (and hence normal) is always to the right, we simply add 90 degrees to point along normal.
        new_subducting_arc_normal_azimuth_degrees = old_subPolarity + 90.0
        # Ensure range is [0,360].
        new_subducting_arc_normal_azimuth_degrees = (
            new_subducting_arc_normal_azimuth_degrees + 360.0
        ) % 360.0

        new_subducting_plate_id = old_SPid
        new_subduction_zone_plate_id = old_TRENCHid

        new_convergence_data_sample = (
            new_longitude,
            new_latitude,
            new_convergence_velocity_magnitude,
            new_convergence_obliquity_degrees,
            new_trench_absolute_velocity_magnitude,
            new_trench_absolute_obliquity_degrees,
            new_arc_length_degrees,
            new_subducting_arc_normal_azimuth_degrees,
            new_subducting_plate_id,
            new_subduction_zone_plate_id,
        )

        if output_distance_to_nearest_edge_of_trench:
            # Convert distance in kms to degrees.
            new_distance_to_nearest_edge_of_trench_degrees = math.degrees(
                old_distEdgeTotal / _EARTH_RADIUS_KMS
            )

            new_convergence_data_sample += (
                new_distance_to_nearest_edge_of_trench_degrees,
            )

        if output_distance_to_start_edge_of_trench:
            # Convert distance in kms to degrees.
            new_distance_to_start_edge_of_trench_degrees = math.degrees(
                old_distEdge / _EARTH_RADIUS_KMS
            )

            new_convergence_data_sample += (
                new_distance_to_start_edge_of_trench_degrees,
            )

        if output_convergence_velocity_components:
            new_convergence_velocity_orthogonal = _MM_YR_TO_CM_YR * old_convPerp
            new_convergence_velocity_parallel = _MM_YR_TO_CM_YR * -old_convPar

            # Should give same values as the following...
            #
            # # The orthogonal and parallel components are just magnitude multiplied by cosine and sine.
            # new_convergence_velocity_orthogonal = (
            #     math.cos(math.radians(new_convergence_obliquity_degrees)) * math.fabs(new_convergence_velocity_magnitude))
            # new_convergence_velocity_parallel = (
            #     math.sin(math.radians(new_convergence_obliquity_degrees)) * math.fabs(new_convergence_velocity_magnitude))

            new_convergence_data_sample += (
                new_convergence_velocity_orthogonal,
                new_convergence_velocity_parallel,
            )

        if output_trench_absolute_velocity_components:
            new_trench_absolute_velocity_orthogonal = _MM_YR_TO_CM_YR * old_orthTrench
            new_trench_absolute_velocity_parallel = _MM_YR_TO_CM_YR * -old_parTrench

            # Should give same values as the following...
            #
            # # The orthogonal and parallel components are just magnitude multiplied by cosine and sine.
            # new_trench_absolute_velocity_orthogonal = (
            #     math.cos(math.radians(new_trench_absolute_obliquity_degrees)) * math.fabs(new_trench_absolute_velocity_magnitude))
            # new_trench_absolute_velocity_parallel = (
            #     math.sin(math.radians(new_trench_absolute_obliquity_degrees)) * math.fabs(new_trench_absolute_velocity_magnitude))

            new_convergence_data_sample += (
                new_trench_absolute_velocity_orthogonal,
                new_trench_absolute_velocity_parallel,
            )

        if (
            output_subducting_absolute_velocity
            or output_subducting_absolute_velocity_components
        ):
            # Subducting absolute velocity magnitude obtained parallel and orthogonal components.
            new_subducting_absolute_velocity_magnitude = _MM_YR_TO_CM_YR * math.sqrt(
                old_parAbs * old_parAbs + old_orthAbs * old_orthAbs
            )
            # The subducting absolute velocity magnitude is negative if subducting plate is moving toward the overriding plate.
            if old_orthAbs > 0:
                new_subducting_absolute_velocity_magnitude = (
                    -new_subducting_absolute_velocity_magnitude
                )

            # Same reasoning as for convergence obliquity above.
            new_subducting_absolute_obliquity_degrees = math.degrees(
                math.atan2(-old_parAbs, old_orthAbs)
            )

            if output_subducting_absolute_velocity:
                new_convergence_data_sample += (
                    new_subducting_absolute_velocity_magnitude,
                    new_subducting_absolute_obliquity_degrees,
                )
            if output_subducting_absolute_velocity_components:
                new_subducting_absolute_velocity_orthogonal = (
                    _MM_YR_TO_CM_YR * old_orthAbs
                )
                new_subducting_absolute_velocity_parallel = (
                    _MM_YR_TO_CM_YR * -old_parAbs
                )

                # Should give same values as the following...
                #
                # # The orthogonal and parallel components are just magnitude multiplied by cosine and sine.
                # new_subducting_absolute_velocity_orthogonal = (
                #     math.cos(math.radians(new_subducting_absolute_obliquity_degrees)) * math.fabs(new_subducting_absolute_velocity_magnitude))
                # new_subducting_absolute_velocity_parallel = (
                #     math.sin(math.radians(new_subducting_absolute_obliquity_degrees)) * math.fabs(new_subducting_absolute_velocity_magnitude))

                new_convergence_data_sample += (
                    new_subducting_absolute_velocity_orthogonal,
                    new_subducting_absolute_velocity_parallel,
                )

        if output_trench_normal:
            # Convert trench normal from local (magnitude, azimuth, inclination) to a global 3D vector.
            local_cartesian = pygplates.LocalCartesian((new_latitude, new_longitude))
            trench_normal = (
                local_cartesian.from_magnitude_azimuth_inclination_to_geocentric(
                    (1.0, math.radians(new_subducting_arc_normal_azimuth_degrees), 0.0)
                )
            )
            new_convergence_data_sample += trench_normal.to_xyz()

        new_convergence_data.append(new_convergence_data_sample)

    return new_convergence_data


# List of optional output parameters available to append to the output of each sample point.
#
# NOTE: This should be in the same order as output parameters are actually added.
_OUTPUT_PARAMETER_NAME_LIST = [
    "output_distance_to_nearest_edge_of_trench",
    "output_distance_to_start_edge_of_trench",
    "output_convergence_velocity_components",
    "output_trench_absolute_velocity_components",
    "output_subducting_absolute_velocity",
    "output_subducting_absolute_velocity_components",
    "output_trench_normal",
]


########################
# Command-line parsing #
########################


def warning_format(message, category, filename, lineno, file=None, line=None):
    # return '{0}:{1}: {1}:{1}\n'.format(filename, lineno, category.__name__, message)
    return "{0}: {1}\n".format(category.__name__, message)


# Print the warnings without the filename and line number.
# Users are not going to want to see that.
warnings.formatwarning = warning_format

__description__ = """Find the convergence rates along trenches (subduction zones) over time.

For each time (over a range of times) an output xy file is generated containing the resolved trenches
(with point locations as the first two columns x and y) and the following convergence rate parameters in subsequent columns:

    - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
    - subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
    - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
    - trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
    - length of arc segment (in degrees) that current point is on
    - trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
    - subducting plate ID
    - trench plate ID
    
    - extra data can be appended by specifying optional extra output parameters using the '-x' command-line option.
    
    The order of any extra data is the same order in which the parameters are listed below.
    But note that all parameters are optional, you only need to specify the ones you want.
    
    +------------------------------------------------+-----------------------------------------------------------------------------------+
    | Name                                           | Description                                                                       |
    +------------------------------------------------+-----------------------------------------------------------------------------------+
    | output_distance_to_nearest_edge_of_trench      | Append the distance (in degrees) along the trench line to the nearest             |
    |                                                | trench edge to each returned sample point. A trench edge is the farthermost       |
    |                                                | location on the current trench feature that contributes to a plate boundary.      |
    +------------------------------------------------+-----------------------------------------------------------------------------------+
    | output_distance_to_start_edge_of_trench        | Append the distance (in degrees) along the trench line from the start edge of     |
    |                                                | the trench to each returned sample point. The start of the trench is along the    |
    |                                                | clockwise direction around the overriding plate.                                  |
    +------------------------------------------------+-----------------------------------------------------------------------------------+
    | output_convergence_velocity_components         | Append the convergence velocity orthogonal and parallel                           |
    |                                                | components (in cm/yr) to each returned sample point.                              |
    |                                                | Orthogonal is normal to trench (in direction of overriding plate when positive).  |
    |                                                | Parallel is along trench (90 degrees clockwise from trench normal when positive). |
    +------------------------------------------------+-----------------------------------------------------------------------------------+
    | output_trench_absolute_velocity_components     | Append the trench absolute velocity orthogonal and parallel                       |
    |                                                | components (in cm/yr) to each returned sample point.                              |
    |                                                | Orthogonal is normal to trench (in direction of overriding plate when positive).  |
    |                                                | Parallel is along trench (90 degrees clockwise from trench normal when positive). |
    +------------------------------------------------+-----------------------------------------------------------------------------------+
    | output_subducting_absolute_velocity            | Append the subducting plate absolute velocity magnitude (in cm/yr) and            |
    |                                                | obliquity angle (in degrees) to each returned sample point.                       |
    +------------------------------------------------+-----------------------------------------------------------------------------------+
    | output_subducting_absolute_velocity_components | Append the subducting plate absolute velocity orthogonal and parallel             |
    |                                                | components (in cm/yr) to each returned sample point.                              |
    |                                                | Orthogonal is normal to trench (in direction of overriding plate when positive).  |
    |                                                | Parallel is along trench (90 degrees clockwise from trench normal when positive). |
    +------------------------------------------------+-----------------------------------------------------------------------------------+
    | output_trench_normal                           | Append the x, y and z components of the trench normal unit-length 3D vectors.     |
    |                                                | These vectors are normal to the trench in the direction of subduction             |
    |                                                | (towards overriding plate). These are global 3D vectors which differ from         |
    |                                                | trench normal azimuth angles (ie, angles relative to North).                      |
    +------------------------------------------------+-----------------------------------------------------------------------------------+

The obliquity angles are in the range (-180 180). The range (0, 180) goes clockwise (when viewed from above the Earth) from the
trench normal direction to the velocity vector. The range (0, -180) goes counter-clockwise.
You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.

Note that the convergence velocity magnitude is negative if the plates are diverging (if convergence obliquity angle
is greater than 90 or less than -90). And note that the absolute velocity magnitude is negative if the trench
is moving towards the overriding plate (if absolute obliquity angle is less than 90 or greater than -90) - note that this
ignores the kinematics of the subducting plate. Similiarly for the subducting plate absolute velocity magnitude
(if extra output parameter 'output_subducting_absolute_velocity' is specified).

Each point in the output is the midpoint of a great circle arc between two adjacent points in the trench polyline.
The trench normal vector used in the obliquity calculations is perpendicular to the great circle arc of each point (arc midpoint)
and pointing towards the overriding plate (rather than away from it).

Each trench is sampled at approximately uniform intervals along its length (specified via a threshold sampling distance).
The sampling along the entire length of a trench is not exactly uniform. Each segment along a trench is sampled
such that the samples have a uniform spacing that is less than or equal to the threshold sampling distance. However each segment
in a trench might have a slightly different spacing distance (since segment lengths are not integer multiples of
the threshold sampling distance).

The trench normal (at each arc segment mid-point) always points *towards* the overriding plate.

NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
For example...

%(prog)s -r rotations.rot -m topologies.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- convergence
    """


def add_arguments(parser: argparse.ArgumentParser):
    """add command line argument parser"""

    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.description = __description__

    parser.set_defaults(func=main)

    parser.add_argument(
        "-r",
        "--rotation_filenames",
        type=str,
        nargs="+",
        required=True,
        metavar="rotation_filename",
        help="One or more rotation files.",
    )
    parser.add_argument(
        "-m",
        "--topology_filenames",
        type=str,
        nargs="+",
        required=True,
        metavar="topology_filename",
        help="One or more topology files to generate resolved subducting lines.",
    )
    parser.add_argument(
        "-a",
        "--anchor",
        type=int,
        default=0,
        dest="anchor_plate_id",
        help="Anchor plate id used for reconstructing. Defaults to zero.",
    )

    # Can specify only one of '-i', '-l' or '-t'.
    threshold_sampling_distance_group = parser.add_mutually_exclusive_group()
    threshold_sampling_distance_group.add_argument(
        "-d",
        "--threshold_sampling_distance_degrees",
        type=float,
        help="Threshold sampling distance along trenches (in degrees). "
        "Defaults to {0} degrees.".format(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES),
    )
    threshold_sampling_distance_group.add_argument(
        "-k",
        "--threshold_sampling_distance_kms",
        type=float,
        help="Threshold sampling distance along trenches (in Kms). "
        "Defaults to {0:.2f} Kms (which is equivalent to {1} degrees).".format(
            DEFAULT_THRESHOLD_SAMPLING_DISTANCE_KMS,
            DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES,
        ),
    )

    parser.add_argument(
        "-t",
        "--time_range",
        type=float,
        nargs=2,
        metavar=("young_time", "old_time"),
        default=[DEFAULT_TIME_RANGE_YOUNG_TIME, DEFAULT_TIME_RANGE_OLD_TIME],
        help="The time range (in Ma) from young time to old time. "
        "Defaults to {0} -> {1} Ma.".format(
            DEFAULT_TIME_RANGE_YOUNG_TIME, DEFAULT_TIME_RANGE_OLD_TIME
        ),
    )

    def parse_positive_number(value_string):
        try:
            value = float(value_string)
        except ValueError:
            raise argparse.ArgumentTypeError("%s is not a number" % value_string)

        if value <= 0:
            raise argparse.ArgumentTypeError("%g is not a positive number" % value)

        return value

    parser.add_argument(
        "-i",
        "--time_increment",
        type=parse_positive_number,
        default=DEFAULT_TIME_INCREMENT,
        help="The time increment in My. Defaults to {0} My.".format(
            DEFAULT_TIME_INCREMENT
        ),
    )

    parser.add_argument(
        "-v",
        "--velocity_delta_time",
        type=parse_positive_number,
        default=DEFAULT_VELOCITY_DELTA_TIME,
        help="The delta time interval used to calculate velocities in My. "
        "Defaults to {0} My.".format(DEFAULT_VELOCITY_DELTA_TIME),
    )

    parser.add_argument(
        "-x",
        "--extra_output_parameters",
        type=str,
        nargs="+",
        metavar="output_parameter",
        help="Optional extra output parameters to be appended to the usual (non-optional) "
        "output parameters for each sample point. "
        "Choices include {0}. "
        "The order in which they are output is in the same order as those choices.".format(
            ", ".join(_OUTPUT_PARAMETER_NAME_LIST)
        ),
    )

    parser.add_argument(
        "-g",
        "--output_gpml_filename",
        type=str,
        help="Optional GPML output filename to contain the subduction convergence data for all specified times. "
        "This can then be loaded into GPlates to display the data as scalar coverages.",
    )

    parser.add_argument(
        "output_filename_prefix",
        type=str,
        help="The output filename prefix. An output file is created for each geological time in the sequence where "
        "the filename suffix contains the time and the filename extension.",
    )
    parser.add_argument(
        "-e",
        "--output_filename_extension",
        type=str,
        default="xy",
        help='The output xy filename extension. Defaults to "xy".',
    )
    parser.add_argument(
        "-w",
        "--ignore_topology_warnings",
        action="store_true",
        help="If specified then topology warnings are ignored (not output). "
        "These are the warnings about not finding the overriding and subducting plates.",
    )


def main(args):
    # Topology warnings correspond to calls to "warnings.warn(... , RuntimeWarning)".
    # Ignore them if requested.
    if args.ignore_topology_warnings:
        warnings.simplefilter("ignore", RuntimeWarning)

    if args.time_range[0] > args.time_range[1]:
        raise argparse.ArgumentTypeError(
            "First (young) value in time range is greater than second (old) value"
        )

    # Determine threshold sampling distance.
    if args.threshold_sampling_distance_degrees:
        threshold_sampling_distance_radians = math.radians(
            args.threshold_sampling_distance_degrees
        )
    elif args.threshold_sampling_distance_kms:
        threshold_sampling_distance_radians = (
            args.threshold_sampling_distance_kms
            / pygplates.Earth.equatorial_radius_in_kms
        )
    else:  # default...
        threshold_sampling_distance_radians = math.radians(
            DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES
        )

    # Specify any optional output parameters to append for each sample point.
    kwargs = {}
    if (
        args.extra_output_parameters
    ):  # Default argument is None (not an empty list), so check first before iterating.
        for extra_output_parameter in args.extra_output_parameters:
            if extra_output_parameter in _OUTPUT_PARAMETER_NAME_LIST:
                kwargs[extra_output_parameter] = True
            else:
                raise ValueError(
                    "{0} is not a valid output parameter".format(extra_output_parameter)
                )

    return_code = subduction_convergence_over_time(
        args.output_filename_prefix,
        args.output_filename_extension,
        args.rotation_filenames,
        args.topology_filenames,
        threshold_sampling_distance_radians,
        args.time_range[0],
        args.time_range[1],
        args.time_increment,
        args.velocity_delta_time,
        args.anchor_plate_id,
        args.output_gpml_filename,
        **kwargs,
    )
    if return_code is None:
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    # The command-line parser.
    parser = argparse.ArgumentParser(
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # add arguments
    add_arguments(parser)

    # Parse command-line options.
    args = parser.parse_args()

    # call main function
    try:
        main(args)
    except Exception as exc:
        print(
            "ERROR: {0}: {1}".format(os.path.basename(__file__), exc), file=sys.stderr
        )
        # Uncomment this to print traceback to location of raised exception.
        # import traceback
        # traceback.print_exc()
        sys.exit(1)
