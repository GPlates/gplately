import logging
import math
import multiprocessing
import os
from functools import partial
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
import pygplates

from ..ptt import separate_ridge_transform_segments
from ..tools import calculate_spreading_rates
from ..utils.log_utils import get_debug_level

logger = logging.getLogger("gplately")
MOR_GPMLZ_FILE_NAME = "MOR_plus_one_points_{:0.2f}.gpmlz"


def _generate_all_mid_ocean_ridge_points(
    times,
    delta_time,
    output_file_path,
    rotation_files,
    topology_files,
    ridge_sampling,
    num_cpus=None,
    overwrite=False,
):
    if not num_cpus:
        try:
            num_cpus = multiprocessing.cpu_count() - 1
        except NotImplementedError:
            num_cpus = 1

        if num_cpus > 1:
            with multiprocessing.Pool(num_cpus) as pool:
                pool.map(
                    partial(
                        _generate_mid_ocean_ridge_points,
                        delta_time=delta_time,
                        output_file_path=output_file_path,
                        rotation_files=rotation_files,
                        topology_files=topology_files,
                        ridge_sampling=ridge_sampling,
                        overwrite=overwrite,
                    ),
                    times,
                )
        else:
            for time in times:
                _generate_mid_ocean_ridge_points(
                    time,
                    delta_time=delta_time,
                    output_file_path=output_file_path,
                    rotation_files=rotation_files,
                    topology_files=topology_files,
                    ridge_sampling=ridge_sampling,
                    overwrite=overwrite,
                )


def _generate_mid_ocean_ridge_points(
    time: float,
    delta_time: float,
    output_file_path: str,
    rotation_files: List[str],
    topology_files: List[str],
    ridge_sampling,
    overwrite=True,
):
    """generate middle ocean ridge seed points at a given time"""
    mor_fn = output_file_path.format(time)
    if os.path.isfile(mor_fn) and not overwrite:
        logger.info(
            f"Middle ocean ridge file exists and will not create again.\n{mor_fn}"
        )
        return

    topology_features_extracted = pygplates.FeaturesFunctionArgument(topology_files)
    rotation_model = pygplates.RotationModel(rotation_files)

    # Points and their z values that emerge from MORs at this time.
    shifted_mor_points = []
    point_spreading_rates = []

    # Resolve topologies to the current time.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(
        topology_features_extracted.get_features(),
        rotation_model,
        resolved_topologies,
        time,
        shared_boundary_sections,
    )

    # pygplates.ResolvedTopologicalSection objects.
    for shared_boundary_section in shared_boundary_sections:
        if (
            shared_boundary_section.get_feature().get_feature_type()
            == pygplates.FeatureType.create_gpml("MidOceanRidge")
        ):
            spreading_feature = shared_boundary_section.get_feature()

            # Find the stage rotation of the spreading feature in the
            # frame of reference of its geometry at the current
            # reconstruction time (the MOR is currently actively spreading).
            # The stage pole can then be directly geometrically compared
            # to the *reconstructed* spreading geometry.
            stage_rotation = separate_ridge_transform_segments.get_stage_rotation_for_reconstructed_geometry(
                spreading_feature, rotation_model, time
            )
            if not stage_rotation:
                # Skip current feature - it's not a spreading feature.
                continue

            # Get the stage pole of the stage rotation.
            # Note that the stage rotation is already in frame of
            # reference of the *reconstructed* geometry at the spreading time.
            stage_pole, _ = stage_rotation.get_euler_pole_and_angle()

            # One way rotates left and the other right, but don't know
            # which - doesn't matter in our example though.
            rotate_slightly_off_mor_one_way = pygplates.FiniteRotation(
                stage_pole, np.radians(0.01)
            )
            rotate_slightly_off_mor_opposite_way = (
                rotate_slightly_off_mor_one_way.get_inverse()
            )

            # Iterate over the shared sub-segments.
            for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                # Tessellate MOR section.
                mor_points = pygplates.MultiPointOnSphere(
                    shared_sub_segment.get_resolved_geometry().to_tessellated(
                        np.radians(ridge_sampling)
                    )
                )

                coords = mor_points.to_lat_lon_list()
                lats = [i[0] for i in coords]
                lons = [i[1] for i in coords]
                boundary_feature = shared_boundary_section.get_feature()
                left_plate = boundary_feature.get_left_plate(None)
                right_plate = boundary_feature.get_right_plate(None)
                if left_plate is not None and right_plate is not None:
                    # Get the spreading rates for all points in this sub segment
                    (
                        spreading_rates,
                        _,
                    ) = calculate_spreading_rates(
                        time=time,
                        lons=lons,
                        lats=lats,
                        left_plates=[left_plate] * len(lons),
                        right_plates=[right_plate] * len(lons),
                        rotation_model=rotation_model,
                        delta_time=delta_time,
                    )

                else:
                    spreading_rates = [np.nan] * len(lons)

                # Loop through all but the 1st and last points in the current sub segment
                for point, rate in zip(
                    mor_points.get_points()[1:-1],
                    spreading_rates[1:-1],
                ):
                    # Add the point "twice" to the main shifted_mor_points list; once for a L-side
                    # spread, another for a R-side spread. Then add the same spreading rate twice
                    # to the list - this therefore assumes spreading rate is symmetric.
                    shifted_mor_points.append(rotate_slightly_off_mor_one_way * point)
                    shifted_mor_points.append(
                        rotate_slightly_off_mor_opposite_way * point
                    )
                    point_spreading_rates.extend([rate] * 2)

    # save the middle ocean ridges points to .pkl file
    lats_lons = [list(point.to_lat_lon()) for point in shifted_mor_points]
    df = pd.DataFrame(lats_lons, columns=["lat", "lon"])
    df["SPREADING_RATE"] = point_spreading_rates
    df.to_pickle(mor_fn)

    if get_debug_level() > 100:
        # Summarising get_isochrons_for_ridge_snapshot;
        # Write out the ridge point born at 'ridge_time' but their position at 'ridge_time - time_step'.
        mor_point_features = []
        for curr_point in shifted_mor_points:
            feature = pygplates.Feature()
            feature.set_geometry(curr_point)
            feature.set_valid_time(time, -999)  # delete - time_step
            mor_point_features.append(feature)

        mor_points = pygplates.FeatureCollection(mor_point_features)

        # Write MOR points at `time` to gpmlz
        mor_points.write(
            os.path.join(os.path.dirname(mor_fn), MOR_GPMLZ_FILE_NAME.format(time))
        )

    logger.info(f"Finished building MOR seedpoints at {time} Ma!")
