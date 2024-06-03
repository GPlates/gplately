import logging
import os
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
import pygplates
import stripy
from plate_model_manager import PlateModel, PlateModelManager

from ..ptt.utils import points_in_polygons
from ..utils.log_utils import get_debug_level

logger = logging.getLogger("gplately")

initial_ocean_mean_spreading_rate = 75.0


def make_seafloor_grids(
    model_name: str,
    times: List[float],
    rotation_files: List[str] = [],
    topology_files: List[str] = [],
    continent_files: List[str] = [],
):
    """make seafloor grids, such as SEAFLOOR_AGE, SPREADING_RATE, etc"""
    if model_name:
        try:
            plate_model = PlateModelManager().get_model(
                model_name, data_dir="plate-model-repo"
            )
        except:
            plate_model = PlateModel(
                model_name, data_dir="plate-model-repo", readonly=True
            )

    if not rotation_files:
        rotation_files = plate_model.get_rotation_model()
    if not topology_files:
        topology_files = plate_model.get_topologies()
    if not continent_files:
        continent_files = plate_model.get_layer("ContinentalPolygons")
        if "Cratons" in plate_model.get_avail_layers():
            continent_files += plate_model.get_layer("Cratons")

    level = 5
    active_points_df = prepare_initial_active_points(
        times[0],
        rotation_files=rotation_files,
        topology_files=topology_files,
        continent_files=continent_files,
        level=level,
    )
    print(active_points_df)

    """
    self.build_all_continental_masks()
    self.build_all_MOR_seedpoints()

    # not necessary, but put here for readability purpose only
    self.current_active_points_df = self.initial_ocean_point_df

    time = int(self._max_time)
    while True:
        self.current_active_points_df.to_pickle(
            self.sample_points_file_path.format(time)
        )
        self._save_gridding_input_data(time)
        # save debug file
        if get_debug_level() > 100:
            _save_age_grid_sample_points_to_gpml(
                self.current_active_points_df["lon"],
                self.current_active_points_df["lat"],
                self.current_active_points_df["begin_time"] - time,
                time,
                self.sample_points_dir,
            )
        next_time = time - int(self._ridge_time_step)
        if next_time >= int(self._min_time):
            points = [
                pygplates.PointOnSphere(row.lat, row.lon)
                for index, row in self.current_active_points_df.iterrows()
            ]
            # reconstruct_geometry() needs time to be integral value
            # https://www.gplates.org/docs/pygplates/generated/pygplates.topologicalmodel#pygplates.TopologicalModel.reconstruct_geometry
            reconstructed_time_span = self.topological_model.reconstruct_geometry(
                points,
                initial_time=time,
                youngest_time=next_time,
                time_increment=int(self._ridge_time_step),
                deactivate_points=pygplates.ReconstructedGeometryTimeSpan.DefaultDeactivatePoints(
                    threshold_velocity_delta=self.subduction_collision_parameters[0]
                    / 10,  # cms/yr
                    threshold_distance_to_boundary=self.subduction_collision_parameters[
                        1
                    ],  # kms/myr
                    deactivate_points_that_fall_outside_a_network=True,
                ),
            )

            reconstructed_points = reconstructed_time_span.get_geometry_points(
                next_time, return_inactive_points=True
            )
            logger.info(
                f"Finished topological reconstruction of {len(self.current_active_points_df)} points from {time} to {next_time} Ma."
            )
            # update the current activate points to prepare for the reconstruction to "next time"
            self._update_current_active_points_coordinates(reconstructed_points)
            self._remove_continental_points(next_time)
            self._load_middle_ocean_ridge_points(next_time)
            time = next_time
        else:
            break
    """


def generate_initial_ocean_seed_points(
    time: float,
    rotation_files,
    continent_files,
    level=5,
):
    """generate initial ocean seed points at `time`"""
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


def find_distance_to_nearest_ridge(
    resolved_topologies, shared_boundary_sections, points, default_value=5000.0
):
    """calculate the nearest distance to middle ocean ridges"""
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


def prepare_initial_active_points(
    time: float, rotation_files, topology_files, continent_files, level: int = 5
):
    ocean_seed_points = generate_initial_ocean_seed_points(
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

    distances = find_distance_to_nearest_ridge(
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
