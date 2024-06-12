import pygplates

from ..ptt.utils import points_in_polygons


def reconstruct(points, rotation_files, topology_files, from_time, to_time):
    """reconstruct points by using topological polygons"""
    rotation_model = pygplates.RotationModel(rotation_files)
    pids = get_topological_plate_ids(points, rotation_model, topology_files, from_time)
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
        point_r = stage_rotation * point
        points_r.append(point_r)

    return points_r


def get_topological_plate_ids(points, rotation_model, topology_files, time):
    """get plate ids from topological polygons"""
    # Resolve the plate polygons for the current time.
    resolved_topologies = []
    pygplates.resolve_topologies(
        topology_files,
        rotation_model,
        resolved_topologies,
        time,
    )

    # For each valid current point find the resolved topology containing it.
    pip_result = points_in_polygons.find_polygons(
        points,
        [
            resolved_topology.get_resolved_boundary()
            for resolved_topology in resolved_topologies
        ],
        resolved_topologies,
    )

    pids = []
    for resolved_topology in pip_result:
        if resolved_topology:
            pids.append(resolved_topology.get_feature().get_reconstruction_plate_id())
        else:
            pids.append(None)
    return pids
