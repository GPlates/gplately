# SPDX-License-Identifier: GPL-2.0-only
"""CLI method selection without changing the public PTT or reconstruction APIs."""

import argparse
import math
import warnings

import pygplates

from ..ptt import subduction_convergence as ptt
from ..reconstruction import PlateReconstruction


def add_arguments(parser):
    """Reuse the existing CLI contract and add explicit method selection."""
    ptt.add_arguments(parser)
    parser.add_argument(
        "--method",
        choices=("ptt", "pygplates"),
        default="ptt",
        help="Calculation method. Defaults to ptt for compatibility; pygplates "
        "uses PlateReconstruction.tessellate_subduction_zones(). "
        "The methods can differ in sampling and calculated values.",
    )
    parser.set_defaults(func=main)


def main(args):
    if args.method == "ptt":
        return ptt.main(args)

    young, old = args.time_range
    if not all(math.isfinite(value) for value in (young, old)):
        raise ValueError("Time range must be finite.")
    if young > old:
        raise argparse.ArgumentTypeError(
            "First (young) value in time range is greater than second (old) value"
        )
    if not math.isfinite(args.time_increment) or args.time_increment <= 0:
        raise ValueError("Time increment must be finite and positive.")
    if not math.isfinite(args.velocity_delta_time) or args.velocity_delta_time <= 0:
        raise ValueError("Velocity delta time must be finite and positive.")

    if args.threshold_sampling_distance_degrees is not None:
        threshold = math.radians(args.threshold_sampling_distance_degrees)
    elif args.threshold_sampling_distance_kms is not None:
        threshold = (
            args.threshold_sampling_distance_kms
            / pygplates.Earth.equatorial_radius_in_kms
        )
    else:
        threshold = math.radians(ptt.DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES)
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError("Sampling distance must be finite and positive.")

    output_flags = {}
    for name in args.extra_output_parameters or ():
        if name not in ptt._OUTPUT_PARAMETER_NAME_LIST:
            raise ValueError(f"{name} is not a valid output parameter")
        output_flags[name] = True

    reconstruction = PlateReconstruction(
        args.rotation_filenames, args.topology_filenames
    )
    coverage_features = []
    time = young
    with warnings.catch_warnings():
        if args.ignore_topology_warnings:
            warnings.simplefilter("ignore", RuntimeWarning)
        # Match the inclusive legacy comparison, including its time tolerance.
        while time <= pygplates.GeoTimeInstant(old):
            next_time = time + args.time_increment
            if next_time == time:
                raise ValueError("Time increment cannot advance reconstruction time.")
            data = reconstruction.tessellate_subduction_zones(
                time,
                tessellation_threshold_radians=threshold,
                use_ptt=False,
                return_geodataframe=False,
                anchor_plate_id=args.anchor_plate_id,
                velocity_delta_time=args.velocity_delta_time,
                **output_flags,
            )
            if len(data):
                # The API ndarray stores IDs as floats. Preserve the legacy text
                # field types without changing any calculated numeric value.
                rows = []
                for sample in data:
                    row = [float(value) for value in sample]
                    for index in (8, 9):
                        if not row[index].is_integer():
                            raise ValueError("Plate IDs must be integers.")
                        row[index] = int(row[index])
                    rows.append(row)
                filename = (
                    f"{args.output_filename_prefix}_{time:.2f}."
                    f"{args.output_filename_extension}"
                )
                ptt.write_output_file(filename, rows)
                if args.output_gpml_filename:
                    coverage_features.append(
                        ptt.create_coverage_feature_from_convergence_data(
                            rows, time, **output_flags
                        )
                    )
            time = next_time

    if args.output_gpml_filename:
        pygplates.FeatureCollection(coverage_features).write(args.output_gpml_filename)
    return 0
