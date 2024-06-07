import logging
import math
import multiprocessing
import os
from functools import partial
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
import pygplates

from ..grids import Raster, write_netcdf_grid
from ..tools import griddata_sphere
from ..utils.log_utils import get_debug_level

logger = logging.getLogger("gplately")


def _create_netcdf_files(
    variable_name,
    times,
    seed_points_file_path,
    continent_mask_file_path,
    output_file_path,
    num_cpus=None,
):
    if num_cpus is None:
        try:
            num_cpus = multiprocessing.cpu_count() - 1
        except NotImplementedError:
            num_cpus = 1

    if num_cpus > 1:
        with multiprocessing.Pool(num_cpus) as pool:
            pool.map(
                partial(
                    _create_netcdf_file,
                    variable_name=variable_name,
                    output_file_path=output_file_path,
                    continent_mask_file_path=continent_mask_file_path,
                    seed_points_file_path=seed_points_file_path,
                ),
                times,
            )
    else:
        for time in times:
            _create_netcdf_file(
                time=time,
                variable_name=variable_name,
                output_file_path=output_file_path,
                continent_mask_file_path=continent_mask_file_path,
                seed_points_file_path=seed_points_file_path,
            )


def _create_netcdf_file(
    time,
    variable_name,
    output_file_path,
    continent_mask_file_path,
    seed_points_file_path,
):
    df = pd.read_pickle(seed_points_file_path.format(time))
    # drop invalid data
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    # Drop duplicate latitudes and longitudes
    unique_data = df.drop_duplicates(subset=["lon", "lat"])

    # Acquire lons, lats and zvalues for each time
    lons = unique_data["lon"].to_list()
    lats = unique_data["lat"].to_list()
    if variable_name == "SEAFLOOR_AGE":
        zdata = (unique_data["begin_time"] - time).to_numpy()
    elif variable_name == "SPREADING_RATE":
        zdata = unique_data["SPREADING_RATE"].to_numpy()
    else:
        raise Exception(f"Unknown variable name: {variable_name}")

    # Create a regular grid on which to interpolate lats, lons and zdata
    grid_lon = np.linspace(-180, 180, 3601)
    grid_lat = np.linspace(-90, 90, 1801)
    X, Y = np.meshgrid(grid_lon, grid_lat)

    # Interpolate lons, lats and zvals over a regular grid using nearest neighbour interpolation
    Z = griddata_sphere((lons, lats), zdata, (X, Y), method="nearest")

    # Identify regions in the grid in the continental mask
    cont_mask = Raster(data=continent_mask_file_path.format(time))

    # We need the continental mask to match the number of nodes
    # in the uniform grid defined above. This is important if we
    # pass our own continental mask to SeafloorGrid
    if cont_mask.shape[1] != 3601:
        cont_mask.resize(3601, 1801, inplace=True)

    # Use the continental mask
    Z = np.ma.array(
        Raster(data=Z.astype("float")).data.data,
        mask=cont_mask.data.data,
        fill_value=np.nan,
    )

    write_netcdf_grid(
        output_file_path.format(time),
        Z,
    )
    logger.info(f"Save {variable_name} netCDF grid at {time:0.2f} Ma completed!")
