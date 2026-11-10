#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

import argparse
import math

import pygplates

from ..ptt import subduction_convergence as ptt_subduction_convergence
from ..reconstruction import PlateReconstruction

help_str = "Find the convergence rates along trenches (subduction zones) over time."

__description__ = f"""{help_str}

For each time in the time range, an output file is written containing points sampled along the
resolved trenches, with one point per line and the following columns:

    - longitude and latitude of the sample point
    - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
    - subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
    - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
    - trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
    - length of trench (in degrees) represented by the sample point
    - trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at the sample point
    - subducting plate ID
    - trench plate ID

Extra columns can be appended using '-x' (see its help below). They're output in the order they're listed
there, regardless of the order specified. A time is skipped (no output file) if no trench points are found.

The trench normal always points towards the overriding plate. The obliquity angles are in the range (-180, 180),
where (0, 180) is clockwise (when viewed from above the Earth) from the trench normal to the velocity vector.
The convergence velocity magnitude is negative if the plates are diverging, and the trench absolute velocity
magnitude is negative if the trench is moving towards the overriding plate.

By default the convergence is calculated from the velocities of the plates on either side of each sample point
(the same as gplately.PlateReconstruction.tessellate_subduction_zones), and points are sampled at exactly the
threshold sampling distance along each trench. Specify '--use-ptt' to instead use the older PlateTectonicTools
calculation (which uses the stage rotation between the subducting and trench plate IDs), where each sample is
the midpoint of a trench segment, and segments are no longer than the threshold sampling distance.

NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
For example...

    gplately subduction-convergence -r rotations.rot -m topologies.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- convergence
"""


def add_parser(parser):
    """Add command line argument parser."""

    cmd = parser.add_parser(
        "subduction-convergence",
        aliases=("sc",),
        help=help_str,
        add_help=True,
    )
    # keep the old underscore name working (issue #450) without showing it in --help
    parser._name_parser_map["subduction_convergence"] = cmd

    # The input, time range, sampling and output arguments are the same as PlateTectonicTools.
    ptt_subduction_convergence.add_arguments(cmd)
    # Replace the PTT description (that also sets the raw-description formatter).
    cmd.description = __description__
    cmd.set_defaults(func=_run_subduction_convergence)

    # Only the PTT calculation produces topology warnings.
    # Look up the action by its option string so that a rename in PTT fails loudly.
    cmd._option_string_actions["--ignore-topology-warnings"].help = (
        "If specified then topology warnings are ignored (not output). "
        "These are the warnings about not finding the overriding and subducting plates, "
        "and are only produced with '--use-ptt' (the default calculation uses plate boundary "
        "statistics, which is better at finding the overriding and subducting plates)."
    )

    cmd.add_argument(
        "--use-ptt",
        action="store_true",
        help="Use the older PlateTectonicTools calculation of subduction convergence.",
    )
    cmd.add_argument(
        "--convergence-threshold",
        metavar="CM_PER_YR",
        type=float,
        default=None,
        dest="convergence_threshold_in_cm_per_yr",
        help="Only output sample points whose orthogonal convergence velocity (in cm/yr) is at least this value. "
        "For example, 0 removes diverging points. Can be negative. Cannot be used with '--use-ptt'. "
        "By default all (converging and diverging) points are output.",
    )
    cmd.add_argument(
        "--include-all-subducting-boundary-types",
        action="store_true",
        help="Also sample features with a subduction polarity that aren't subduction zones "
        "(such as orogenic belts). By default only subduction zone features are sampled.",
    )
    cmd.add_argument(
        "--include-network-boundaries",
        action="store_true",
        help="Also sample deforming network boundaries that aren't plate boundaries. "
        "By default only plate boundaries are sampled.",
    )


def _run_subduction_convergence(args):
    young_time, old_time = args.time_range
    if young_time > old_time:
        raise ValueError(
            f"The young time {young_time} is older (larger) than the old time {old_time}."
        )

    if args.use_ptt and args.convergence_threshold_in_cm_per_yr is not None:
        raise ValueError("'--convergence-threshold' cannot be used with '--use-ptt'.")

    if args.threshold_sampling_distance_degrees is not None:
        threshold_sampling_distance_radians = math.radians(
            args.threshold_sampling_distance_degrees
        )
    elif args.threshold_sampling_distance_kms is not None:
        threshold_sampling_distance_radians = (
            args.threshold_sampling_distance_kms
            / pygplates.Earth.equatorial_radius_in_kms
        )
    else:
        threshold_sampling_distance_radians = math.radians(
            ptt_subduction_convergence.DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES
        )
    if threshold_sampling_distance_radians <= 0:
        raise ValueError("The threshold sampling distance must be positive.")

    # Optional output parameters to append to each sample point.
    output_parameters = {}
    for name in args.extra_output_parameters or []:
        if name not in ptt_subduction_convergence._OUTPUT_PARAMETER_NAME_LIST:
            raise ValueError(f"{name} is not a valid output parameter")
        output_parameters[name] = True

    plate_reconstruction = PlateReconstruction(
        args.rotation_filenames,
        args.topology_filenames,
        anchor_plate_id=args.anchor_plate_id,
    )

    coverage_features = []

    time = young_time
    while time <= pygplates.GeoTimeInstant(old_time):
        subduction_data = plate_reconstruction.tessellate_subduction_zones(
            time,
            threshold_sampling_distance_radians,
            ignore_warnings=args.ignore_topology_warnings,
            use_ptt=args.use_ptt,
            include_network_boundaries=args.include_network_boundaries,
            include_all_subducting_boundary_types=args.include_all_subducting_boundary_types,
            convergence_threshold_in_cm_per_yr=args.convergence_threshold_in_cm_per_yr,
            velocity_delta_time=args.velocity_delta_time,
            **output_parameters,
        )

        if len(subduction_data):
            # The returned array is all floats, so write the plate IDs (columns 8 and 9) as integers.
            output_data = [
                sample[:8].tolist()
                + [int(sample[8]), int(sample[9])]
                + sample[10:].tolist()
                for sample in subduction_data
            ]

            ptt_subduction_convergence.write_output_file(
                f"{args.output_filename_prefix}_{time:0.2f}.{args.output_filename_extension}",
                output_data,
            )

            if args.output_gpml_filename:
                coverage_features.append(
                    ptt_subduction_convergence.create_coverage_feature_from_convergence_data(
                        output_data, time, **output_parameters
                    )
                )

        time += args.time_increment

    if args.output_gpml_filename:
        pygplates.FeatureCollection(coverage_features).write(args.output_gpml_filename)
