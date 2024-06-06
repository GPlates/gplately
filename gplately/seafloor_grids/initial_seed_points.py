import logging
import os
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
import pygplates
import stripy

from ..ptt.utils import points_in_polygons
from ..utils.log_utils import get_debug_level

logger = logging.getLogger("gplately")


def _generate_initial_ocean_seed_points(
    time: float,
    rotation_files,
    continent_files,
    level=5,
):
    """generate initial ocean seed points at `time`
    return the seed points as a list of pygplates.PointOnSphere
    """
    # reconstruct the continental polygons to `time`
    reconstructed_feature_geometries = []
    pygplates.reconstruct(
        continent_files,
        rotation_files,
        reconstructed_feature_geometries,
        reconstruction_time=time,
        anchor_plate_id=0,
    )

    polygons = []
    for rfg in reconstructed_feature_geometries:
        geom = rfg.get_reconstructed_geometry()
        begin_time, end_time = rfg.get_feature().get_valid_time()
        if isinstance(geom, pygplates.PolygonOnSphere) and begin_time >= end_time:
            polygons.append(geom)

    # generate the ocean basin points using Stripy's icosahedral spherical mesh
    icosahedral_global_mesh = stripy.spherical_meshes.icosahedral_mesh(
        level, include_face_points=False, trisection=False, tree=False
    )
    lats = np.rad2deg(icosahedral_global_mesh.lats)
    lons = np.rad2deg(icosahedral_global_mesh.lons)
    points = [pygplates.PointOnSphere(lat, lon) for lat, lon in zip(lats, lons)]

    pip_result = points_in_polygons.find_polygons(points, polygons, all_polygons=False)

    sea_points = []
    for idx, polygon in enumerate(pip_result):
        if not polygon:
            sea_points.append(points[idx])

    return sea_points


def _find_distance_to_nearest_ridge(
    resolved_topologies, shared_boundary_sections, points, default_value=5000.0
):
    """for each seed point, calculate the nearest distance to middle ocean ridges
    return a list of distances(floating point numbers)
    """
    # first, let's get all "MidOceanRidge" features and create a dict with plate ids as keys
    # later we will need these "MidOceanRidge" features to calculate the distance
    mid_ocean_ridges_on_plate = {}
    for shared_boundary_section in shared_boundary_sections:
        if (
            shared_boundary_section.get_feature().get_feature_type()
            == pygplates.FeatureType.create_gpml("MidOceanRidge")
        ):
            for shared_subsegment in shared_boundary_section.get_shared_sub_segments():
                sharing_resolved_topologies = (
                    shared_subsegment.get_sharing_resolved_topologies()
                )
                for resolved_polygon in sharing_resolved_topologies:
                    plate_id = (
                        resolved_polygon.get_feature().get_reconstruction_plate_id()
                    )
                    if plate_id in mid_ocean_ridges_on_plate:
                        mid_ocean_ridges_on_plate[plate_id].append(
                            shared_subsegment.get_resolved_geometry()
                        )
                    else:
                        mid_ocean_ridges_on_plate[plate_id] = [
                            shared_subsegment.get_resolved_geometry()
                        ]

    # iterate through all the points to calculate the distance
    all_distance_to_ridge = []
    for point in points:
        min_distance_to_ridge = None
        for topology in resolved_topologies:
            if topology.get_resolved_geometry().is_point_in_polygon(point):
                pid = topology.get_resolved_feature().get_reconstruction_plate_id()
                if pid in mid_ocean_ridges_on_plate:
                    for ridge in mid_ocean_ridges_on_plate[pid]:
                        distance_to_ridge = pygplates.GeometryOnSphere.distance(
                            point, ridge, min_distance_to_ridge
                        )
                        if distance_to_ridge is not None:
                            min_distance_to_ridge = distance_to_ridge
                break

        all_distance_to_ridge.append(min_distance_to_ridge)

    return [
        d * pygplates.Earth.mean_radius_in_kms if d is not None else default_value
        for d in all_distance_to_ridge
    ]


def _get_initial_active_points_df(
    time: float,
    rotation_files,
    topology_files,
    continent_files,
    level: int = 5,
    initial_ocean_mean_spreading_rate=75.0,
):
    """return the initial seed points as Pandas dataframe"""
    ocean_seed_points = _generate_initial_ocean_seed_points(
        time=time,
        rotation_files=rotation_files,
        continent_files=continent_files,
        level=level,
    )

    # calculate the seafloor age from "spreading rate" and "distance to ridge"
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(
        topology_files,
        rotation_files,
        resolved_topologies,
        time,
        shared_boundary_sections,
    )

    distances = _find_distance_to_nearest_ridge(
        resolved_topologies, shared_boundary_sections, ocean_seed_points
    )
    # Divide spreading rate by 2 to use half the mean spreading rate
    pAge = np.array(distances) / (initial_ocean_mean_spreading_rate / 2.0)

    # prepare DataFrame to return
    lats_lons = [p.to_lat_lon() for p in ocean_seed_points]
    data = {
        "lon": [i[1] for i in lats_lons],
        "lat": [i[0] for i in lats_lons],
        "begin_time": pAge + time,
        "end_time": 0,
        "SPREADING_RATE": initial_ocean_mean_spreading_rate,
    }

    # the code below is for debug purpose only
    if get_debug_level() > 100:
        initial_ocean_point_features = []
        for point, age in zip(ocean_seed_points, pAge):
            point_feature = pygplates.Feature()
            point_feature.set_geometry(point)
            point_feature.set_valid_time(age + time, 0)
            initial_ocean_point_features.append(point_feature)

        pygplates.FeatureCollection(initial_ocean_point_features).write(
            os.path.join(
                ".",
                "ocean_basin_seed_points_{}_RLs_{}Ma.gpmlz".format(
                    level,
                    time,
                ),
            )
        )
    logger.info("Finished building initial_ocean_seed_points!")
    return pd.DataFrame(data=data)
