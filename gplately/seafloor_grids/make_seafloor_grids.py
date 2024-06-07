import logging
import os
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
import pygplates
import stripy
from plate_model_manager import PlateModel, PlateModelManager

from ..grids import read_netcdf_grid
from ..ptt.utils import points_in_polygons
from ..utils.log_utils import get_debug_level
from .continental_mask import _build_all_continental_masks
from .create_netcdf import _create_netcdf_files
from .initial_seed_points import _get_initial_active_points_df
from .mid_ocean_ridges import _generate_all_mid_ocean_ridge_points

logger = logging.getLogger("gplately")


def make_seafloor_grids(
    model_name: str,
    initial_time: int,
    youngest_time: int,
    time_increment: int = 1,
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

    times = list(range(initial_time, youngest_time - 1, -time_increment))
    assert len(times) > 1

    level = 5
    current_active_points_df = _get_initial_active_points_df(
        initial_time,
        rotation_files=rotation_files,
        topology_files=topology_files,
        continent_files=continent_files,
        level=level,
        initial_ocean_mean_spreading_rate=initial_ocean_mean_spreading_rate,
    )

    continental_mask_dir = "./continent_masks"
    Path(continental_mask_dir).mkdir(parents=True, exist_ok=True)
    continental_mask_file_path = os.path.join(
        continental_mask_dir, "continent_mask_{:0.2f}Ma.nc"
    )
    _build_all_continental_masks(
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
    _generate_all_mid_ocean_ridge_points(
        times,
        delta_time=1,
        output_file_path=mid_ocean_ridge_file_path,
        rotation_files=rotation_files,
        topology_files=topology_files,
        ridge_sampling=0.5,
    )

    # seed points file path
    seed_points_dir = "./seed_points"
    Path(seed_points_dir).mkdir(parents=True, exist_ok=True)
    seed_points_file_path = os.path.join(seed_points_dir, "seed_points_{:0.2f}_Ma.pkl")

    current_time = initial_time
    next_time = initial_time - time_increment
    topological_model = pygplates.TopologicalModel(topology_files, rotation_files)
    threshold_velocity_delta = 0.5
    threshold_distance_to_boundary = 10

    # save the initial active seed points
    current_active_points_df.to_pickle(seed_points_file_path.format(initial_time))
    # save debug file
    if get_debug_level() > 100:
        _save_age_grid_seed_points_to_gpml(
            current_active_points_df["lon"],
            current_active_points_df["lat"],
            current_active_points_df["begin_time"] - initial_time,
            initial_time,
            seed_points_dir,
        )

    while next_time >= youngest_time:
        points = [
            pygplates.PointOnSphere(row.lat, row.lon)
            for index, row in current_active_points_df.iterrows()
        ]
        # reconstruct_geometry() needs time to be integral value
        # https://www.gplates.org/docs/pygplates/generated/pygplates.topologicalmodel#pygplates.TopologicalModel.reconstruct_geometry
        reconstructed_time_span = topological_model.reconstruct_geometry(
            points,
            initial_time=current_time,
            youngest_time=next_time,
            time_increment=time_increment,
            deactivate_points=pygplates.ReconstructedGeometryTimeSpan.DefaultDeactivatePoints(
                threshold_velocity_delta=threshold_velocity_delta,  # cms/yr
                threshold_distance_to_boundary=threshold_distance_to_boundary,  # kms/myr
                deactivate_points_that_fall_outside_a_network=True,
            ),
        )

        reconstructed_points = reconstructed_time_span.get_geometry_points(
            next_time, return_inactive_points=True
        )
        logger.info(
            f"Finished topological reconstruction of {len(current_active_points_df)} points from {current_time} to {next_time} Ma."
        )

        # update the current activate points to prepare for the reconstruction to "next time"
        current_active_points_df = _update_current_active_points_coordinates(
            current_active_points_df, reconstructed_points
        )
        current_active_points_df = _remove_continental_points(
            current_active_points_df, next_time, continental_mask_file_path
        )
        current_active_points_df = _load_middle_ocean_ridge_points(
            current_active_points_df, next_time, mid_ocean_ridge_file_path
        )

        current_time = next_time
        next_time = current_time - time_increment

        # save the reconstructed active seed points
        current_active_points_df.to_pickle(seed_points_file_path.format(current_time))
        # save debug file
        if get_debug_level() > 100:
            _save_age_grid_seed_points_to_gpml(
                current_active_points_df["lon"],
                current_active_points_df["lat"],
                current_active_points_df["begin_time"] - current_time,
                current_time,
                seed_points_dir,
            )

    seafloor_age_dir = "./seafloor_age"
    Path(seafloor_age_dir).mkdir(parents=True, exist_ok=True)
    seafloor_age_file_path = os.path.join(
        seafloor_age_dir, "SEAFLOOR_AGE_grid_{:0.2f}_Ma.nc"
    )
    _create_netcdf_files(
        variable_name="SEAFLOOR_AGE",
        times=times,
        seed_points_file_path=seed_points_file_path,
        continent_mask_file_path=continental_mask_file_path,
        output_file_path=seafloor_age_file_path,
    )


def _update_current_active_points_coordinates(
    current_active_points_df, reconstructed_points: List[pygplates.PointOnSphere]
):
    """Update the current active points with the reconstructed coordinates.
    The length of `reconstructed_points` must be the same with the length of self.current_active_points_df
    """
    # TODO: use a class to manage the seed points data
    assert len(reconstructed_points) == len(current_active_points_df)
    lons = []
    lats = []
    begin_times = []
    end_times = []
    spread_rates = []
    for i in range(len(reconstructed_points)):
        if reconstructed_points[i]:
            lat_lon = reconstructed_points[i].to_lat_lon()
            lons.append(lat_lon[1])
            lats.append(lat_lon[0])
            begin_times.append(current_active_points_df.loc[i, "begin_time"])
            end_times.append(current_active_points_df.loc[i, "end_time"])
            spread_rates.append(current_active_points_df.loc[i, "SPREADING_RATE"])
    data = {
        "lon": lons,
        "lat": lats,
        "begin_time": begin_times,
        "end_time": end_times,
        "SPREADING_RATE": spread_rates,
    }
    return pd.DataFrame(data=data)


def _remove_continental_points(current_active_points_df, time, continent_mask_filepath):
    """remove all the points which are inside continents at `time` from current_active_points_df"""
    gridZ, gridX, gridY = read_netcdf_grid(
        continent_mask_filepath.format(time), return_grids=True
    )
    ni, nj = gridZ.shape
    xmin = np.nanmin(gridX)
    xmax = np.nanmax(gridX)
    ymin = np.nanmin(gridY)
    ymax = np.nanmax(gridY)

    # TODO: take this function out somewhere so that other functions can use it
    def remove_points_on_continents(row):
        i = int(round((ni - 1) * ((row.lat - ymin) / (ymax - ymin))))
        j = int(round((nj - 1) * ((row.lon - xmin) / (xmax - xmin))))
        i = 0 if i < 0 else i
        j = 0 if j < 0 else j
        i = ni - 1 if i > ni - 1 else i
        j = nj - 1 if j > nj - 1 else j

        if gridZ[i, j] > 0:
            return False
        else:
            return True

    m = current_active_points_df.apply(remove_points_on_continents, axis=1)
    return current_active_points_df[m]


def _load_middle_ocean_ridge_points(
    current_active_points_df, time, mid_ocean_ridges_file_path
):
    """add middle ocean ridge points at `time` to current_active_points_df"""
    df = pd.read_pickle(mid_ocean_ridges_file_path.format(time))
    df["begin_time"] = time
    df["end_time"] = 0
    return pd.concat(
        [
            current_active_points_df,
            df,
        ],
        ignore_index=True,
    )


def _save_age_grid_seed_points_to_gpml(
    lons, lats, seafloor_ages, paleo_time, output_file_dir
):
    """save sample points to .gpmlz for debug purpose"""
    logger.debug(f"saving age grid sample points to gpml file -- {paleo_time} Ma")
    features = []
    for lon, lat, age in zip(lons, lats, seafloor_ages):
        f = pygplates.Feature()
        p = pygplates.PointOnSphere(lat, lon)
        f.set_geometry(p)
        f.set_valid_time(age + paleo_time, paleo_time)
        features.append(f)
    pygplates.FeatureCollection(features).write(
        os.path.join(
            output_file_dir, "MOR_plus_one_points_{:0.2f}.gpmlz".format(paleo_time)
        )
    )
