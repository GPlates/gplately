import logging
import os
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
import pygplates
import stripy
from plate_model_manager import PlateModel, PlateModelManager

from ..ptt.utils import points_in_polygons
from ..utils.log_utils import get_debug_level
from .continental_mask import _build_continental_masks
from .initial_seed_points import _get_initial_active_points_df
from .mid_ocean_ridges import _generate_mid_ocean_ridge_points

logger = logging.getLogger("gplately")


def make_seafloor_grids(
    model_name: str,
    times: List[float],
    rotation_files: List[str] = [],
    topology_files: List[str] = [],
    continent_files: List[str] = [],
    initial_ocean_mean_spreading_rate=75.0,
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
    active_points_df = _get_initial_active_points_df(
        times[0],
        rotation_files=rotation_files,
        topology_files=topology_files,
        continent_files=continent_files,
        level=level,
        initial_ocean_mean_spreading_rate=initial_ocean_mean_spreading_rate,
    )
    print(active_points_df)

    continental_mask_dir = "./continent_masks"
    Path(continental_mask_dir).mkdir(parents=True, exist_ok=True)
    continental_mask_file_path = os.path.join(
        continental_mask_dir, "continent_mask_{:0.2f}Ma.nc"
    )
    _build_continental_masks(
        output_file_path=continental_mask_file_path,
        times=times,
        rotation_files=rotation_files,
        continent_files=continent_files,
    )

    mid_ocean_ridge_dir = "./mid_ocean_ridges"
    Path(mid_ocean_ridge_dir).mkdir(parents=True, exist_ok=True)
    mid_ocean_ridge_file_path = os.path.join(
        mid_ocean_ridge_dir, "MOR_df_{:0.2f}_Ma.pkl"
    )
    _generate_mid_ocean_ridge_points(
        times,
        delta_time=1,
        output_file_path=mid_ocean_ridge_file_path,
        rotation_files=rotation_files,
        topology_files=topology_files,
        ridge_sampling=0.5,
    )

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
