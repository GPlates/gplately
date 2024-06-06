import logging
import math
import multiprocessing
import os
from functools import partial

import pygplates

from ..grids import write_netcdf_grid
from ..ptt import continent_contours

logger = logging.getLogger("gplately")


def _continent_contouring_buffer_and_gap_distance_radians(time, contoured_continent):
    # One distance for time interval [1000, 300] and another for time interval [250, 0].
    # And linearly interpolate between them over the time interval [300, 250].
    pre_pangea_distance_radians = math.radians(2.5)  # convert degrees to radians
    post_pangea_distance_radians = math.radians(0.0)  # convert degrees to radians
    if time > 300:
        buffer_and_gap_distance_radians = pre_pangea_distance_radians
    elif time < 250:
        buffer_and_gap_distance_radians = post_pangea_distance_radians
    else:
        # Linearly interpolate between 250 and 300 Ma.
        interp = float(time - 250) / (300 - 250)
        buffer_and_gap_distance_radians = (
            interp * pre_pangea_distance_radians
            + (1 - interp) * post_pangea_distance_radians
        )

    # Area of the contoured continent.
    area_steradians = contoured_continent.get_area()

    # Linearly reduce the buffer/gap distance for contoured continents with area smaller than 1 million km^2.
    area_threshold_square_kms = 500000
    area_threshold_steradians = area_threshold_square_kms / (
        pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms
    )
    if area_steradians < area_threshold_steradians:
        buffer_and_gap_distance_radians *= area_steradians / area_threshold_steradians

    return buffer_and_gap_distance_radians


def _build_continental_mask_with_contouring(
    time: float,
    continent_mask_filepath,
    rotation_files,
    continent_files,
    overwrite=False,
):
    """build the continent mask for a given time with continent contouring method"""
    mask_fn = continent_mask_filepath.format(time)
    if os.path.isfile(mask_fn) and not overwrite:
        logger.info(f"Continent mask file exists and will not create again.\n{mask_fn}")
        return

    continent_contouring_point_spacing_degrees = 0.25
    continent_contouring_area_threshold_square_kms = 0
    continent_contouring_area_threshold_steradians = (
        continent_contouring_area_threshold_square_kms
        / (pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms)
    )
    continent_exclusion_area_threshold_square_kms = 800000
    continent_exclusion_area_threshold_steradians = (
        continent_exclusion_area_threshold_square_kms
        / (pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms)
    )
    continent_separation_distance_threshold_radians = (
        continent_contours.DEFAULT_CONTINENT_SEPARATION_DISTANCE_THRESHOLD_RADIANS
    )

    continent_features = [pygplates.FeatureCollection(f) for f in continent_files]
    continent_contouring = continent_contours.ContinentContouring(
        pygplates.RotationModel(rotation_files),
        continent_features,
        continent_contouring_point_spacing_degrees,
        continent_contouring_area_threshold_steradians,
        _continent_contouring_buffer_and_gap_distance_radians,
        continent_exclusion_area_threshold_steradians,
        continent_separation_distance_threshold_radians,
    )
    continent_mask, _ = (
        continent_contouring.get_continent_mask_and_contoured_continents(time)
    )
    write_netcdf_grid(
        continent_mask_filepath.format(time),
        continent_mask.astype("float"),
        extent=[-180, 180, -90, 90],
    )
    logger.info(
        f"Finished building a continental mask at {time} Ma! (continent_contouring)"
    )


def _build_continental_masks(
    output_file_path,
    times,
    rotation_files,
    continent_files,
    num_cpus=None,
    overwrite=False,
):
    """Create a continental mask to define the ocean basin for all times between `min_time` and `max_time`.
    The continental masks will be saved to 'output_file_path' as compressed netCDF4 files.
    The 'output_file_path' should be something like "./continent_masks/continent_mask_{:0.2f}Ma.nc"
    """
    if not num_cpus:
        try:
            num_cpus = multiprocessing.cpu_count() - 1
        except NotImplementedError:
            num_cpus = 1

        if num_cpus > 1:
            with multiprocessing.Pool(num_cpus) as pool:
                pool.map(
                    partial(
                        _build_continental_mask_with_contouring,
                        continent_mask_filepath=output_file_path,
                        rotation_files=rotation_files,
                        continent_files=continent_files,
                        overwrite=overwrite,
                    ),
                    times,
                )
        else:
            for time in times:
                _build_continental_mask_with_contouring(
                    time,
                    continent_mask_filepath=output_file_path,
                    rotation_files=rotation_files,
                    continent_files=continent_files,
                    overwrite=overwrite,
                )
