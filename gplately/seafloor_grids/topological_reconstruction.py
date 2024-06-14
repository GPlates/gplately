import math

import pygplates

from ..ptt.utils import points_in_polygons


def reconstruct(points, rotation_files, topology_files, from_time, to_time):
    """reconstruct points by using topological polygons"""
    rotation_model = pygplates.RotationModel(rotation_files)
    resolved_topological_boundaries = _get_resolved_topological_boundary_for_points(
        points, rotation_model, topology_files, from_time
    )
    pids = [
        b.get_feature().get_reconstruction_plate_id() if b else None
        for b in resolved_topological_boundaries
    ]
    # Cache stage rotations by plate ID.
    reconstruct_stage_rotation_dict = {}

    points_r = []
    for point, pid in zip(points, pids):
        if not pid:
            points_r.append(None)
            continue

        # Speed up by caching stage rotations in a dict.
        stage_rotation = reconstruct_stage_rotation_dict.get(pid)
        if not stage_rotation:
            stage_rotation = rotation_model.get_rotation(
                to_time,
                pid,
                from_time,
            )
            reconstruct_stage_rotation_dict[pid] = stage_rotation

        # Use the stage rotation to reconstruct the point
        points_r.append(stage_rotation * point)

    new_pids = [
        b.get_feature().get_reconstruction_plate_id() if b else None
        for b in _get_resolved_topological_boundary_for_points(
            points, rotation_model, topology_files, to_time
        )
    ]
    for idx in range(len(points)):
        if _detect_collision(
            rotation_model,
            to_time,
            from_time - to_time,
            points[idx],
            points_r[idx],
            pids[idx],
            resolved_topological_boundaries[idx],
            new_pids[idx],
        ):
            points_r[idx] = None

    return points_r


def _get_resolved_topological_boundary_for_points(
    points, rotation_model, topology_files, time
):
    """return a list of pygplates.ResolvedTopologicalBoundary.
    one pygplates.ResolvedTopologicalBoundary object or None for each point
    """
    # Resolve the plate polygons for the current time.
    resolved_topologies = []
    pygplates.resolve_topologies(
        topology_files,
        rotation_model,
        resolved_topologies,
        time,
    )

    # For each valid current point find the resolved topology containing it.
    return points_in_polygons.find_polygons(
        points,
        [
            resolved_topology.get_resolved_boundary()
            for resolved_topology in resolved_topologies
        ],
        resolved_topologies,
    )


def _detect_collision(
    rotation_model,
    time,
    reconstruction_time_interval,
    prev_point,
    curr_point,
    prev_topology_plate_id,
    prev_resolved_plate_boundary,
    curr_topology_plate_id,
):
    """
    Returns True if a collision was detected.

    If transitioning from a rigid plate to another rigid plate with a different plate ID then
    calculate the difference in velocities and continue testing as follows
    (otherwise, if there's no transition, then the point is still active)...

    If the velocity difference is below a threshold then we assume the previous plate was split,
    or two plates joined. In this case the point has not subducted (forward in time) or
    been consumed by a mid-ocean (backward in time) and hence is still active.

    If the velocity difference is large enough then we see if the distance of the *previous* position
    to the polygon boundary (of rigid plate containing it) exceeds a threshold.
    If the distance exceeds the threshold then the point is far enough away from the boundary that it
    cannot be subducted or consumed by it and hence the point is still active.
    However if the point is close enough then we assume the point was subducted/consumed
    (remember that the point switched plate IDs).
    Also note that the threshold distance increases according to the velocity difference to account for fast
    moving points (that would otherwise tunnel through the boundary and accrete onto the other plate).
    The reason for testing the distance from the *previous* point, and not from the *current* point, is:

      (i)  A topological boundary may *appear* near the current point (such as a plate split at the current time)
           and we don't want that split to consume the current point regardless of the velocity difference.
           It won't get consumed because the *previous* point was not near a boundary (because before split happened).
           If the velocity difference is large enough then it might cause the current point to transition to the
           adjacent split plate in the *next* time step (and that's when it should get consumed, not in the current time step).
           An example of this is a mid-ocean ridge suddenly appearing (going forward in time).

      (ii) A topological boundary may *disappear* near the current point (such as a plate merge at the current time)
           and we want that merge to consume the current point if the velocity difference is large enough.
           In this case the *previous* point is near a boundary (because before plate merged) and hence can be
           consumed (provided velocity difference is large enough). And since the boundary existed in the previous
           time step, it will affect position of the current point (and whether it gets consumed or not).
           An example of this is a mid-ocean ridge suddenly disappearing (going backward in time).

    ...note that items (i) and (ii) above apply both going forward and backward in time.
    """
    if not (
        curr_topology_plate_id != prev_topology_plate_id
        and prev_topology_plate_id is not None
        and curr_topology_plate_id is not None
    ):
        return False

    prev_location_velocity_stage_rotation = rotation_model.get_rotation(
        time + 1, prev_topology_plate_id, time
    )

    curr_location_velocity_stage_rotation = rotation_model.get_rotation(
        time + 1, curr_topology_plate_id, time
    )

    # Note that even though the current point is not inside the previous boundary (because different plate ID), we can still
    # calculate a velocity using its plate ID (because we really should use the same point in our velocity comparison).
    prev_location_velocity = pygplates.calculate_velocities(
        (curr_point,),
        prev_location_velocity_stage_rotation,
        1,
        pygplates.VelocityUnits.kms_per_my,
    )[0]
    curr_location_velocity = pygplates.calculate_velocities(
        (curr_point,),
        curr_location_velocity_stage_rotation,
        1,
        pygplates.VelocityUnits.kms_per_my,
    )[0]

    delta_velocity = curr_location_velocity - prev_location_velocity
    delta_velocity_magnitude = delta_velocity.get_magnitude()

    threshold_velocity_delta = 5
    threshold_distance_to_boundary_per_my = 10

    for (
        prev_boundary_sub_segment
    ) in prev_resolved_plate_boundary.get_boundary_sub_segments():
        if _detect_collision_using_collision_parameters(
            reconstruction_time_interval,
            delta_velocity_magnitude,
            prev_point,
            prev_boundary_sub_segment.get_resolved_geometry(),
            threshold_velocity_delta,
            threshold_distance_to_boundary_per_my,
        ):
            return True
    return False


def _detect_collision_using_collision_parameters(
    reconstruction_time_interval,
    delta_velocity_magnitude,
    prev_point,
    prev_boundary_geometry,
    threshold_velocity_delta,
    threshold_distance_to_boundary_per_my,
):
    if delta_velocity_magnitude > threshold_velocity_delta:
        # Add the minimum distance threshold to the delta velocity threshold.
        # The delta velocity threshold only allows those points that are close enough to the boundary to reach
        # it given their current relative velocity.
        # The minimum distance threshold accounts for sudden changes in the shape of a plate boundary
        # which are no supposed to represent a new or shifted boundary but are just a result of the topology
        # builder/user digitising a new boundary line that differs noticeably from that of the previous time period.
        distance_threshold_radians = (
            (threshold_distance_to_boundary_per_my + delta_velocity_magnitude)
            * math.fabs(reconstruction_time_interval)
            / pygplates.Earth.equatorial_radius_in_kms
        )
        distance_threshold_radians = min(distance_threshold_radians, math.pi)
        distance_threshold_radians = max(distance_threshold_radians, 0.0)

        distance = pygplates.GeometryOnSphere.distance(
            prev_point,
            prev_boundary_geometry,
            distance_threshold_radians=float(distance_threshold_radians),
        )

        if distance is not None:
            # Detected a collision.
            return True

    return False
