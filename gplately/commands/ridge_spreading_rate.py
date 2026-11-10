#
#    Copyright (C) 2026 The University of Sydney, Australia
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
import functools
import math

import pygplates

from ..ptt import separate_ridge_transform_segments
from ..reconstruction import PlateReconstruction

DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES = 0.5
DEFAULT_TIME_RANGE_YOUNG_TIME = 0
DEFAULT_TIME_RANGE_OLD_TIME = 200
DEFAULT_TIME_INCREMENT = 1
DEFAULT_VELOCITY_DELTA_TIME = 1
DEFAULT_SPREADING_FEATURE_TYPES = [pygplates.FeatureType.gpml_mid_ocean_ridge]

# Names of the scalar coverages written to the GPML file (one per output column after lon/lat).
_COVERAGE_SCALAR_NAMES = [
    "SpreadingVelocityMagnitude",
    "SpreadingObliquityDegrees",
    "SpreadingLengthDegrees",
    "SpreadingNormalAzimuthDegrees",
    "LeftPlateId",
    "RightPlateId",
]

help_str = "Find the spreading rates along mid-ocean ridges over time."

__description__ = f"""{help_str}

For each time in the time range, an output file is written containing points sampled along the
resolved mid-ocean ridges, with one point per line and the following columns:

    - longitude and latitude of the sample point
    - spreading velocity magnitude (in cm/yr)
    - spreading obliquity (in degrees), the deviation of the spreading direction from the ridge normal (0 to 90 degrees)
    - length of ridge (in degrees) represented by the sample point
    - ridge normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at the sample point
    - left plate ID
    - right plate ID

A time is skipped (no output file) if no ridge points are found.

Ridge segments are separated from transform segments using the spreading direction, and transform
segments are skipped (see '--transform-segment-deviation' and '--include-transform-segments').

By default the spreading rate is calculated from the velocities of the plates on either side of each sample point
(the same as gplately.PlateReconstruction.tessellate_mid_ocean_ridges), and points are sampled at exactly the
threshold sampling distance along each ridge segment. Specify '--use-ptt' to instead use the older PlateTectonicTools
calculation (which uses the stage rotation between the left and right plate IDs), where each sample is the midpoint
of a ridge segment, and segments are no longer than the threshold sampling distance.

NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
For example...

    gplately ridge-spreading-rate -r rotations.rot -m topologies.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- spreading_rate
"""


def _positive_float(value_string):
    try:
        value = float(value_string)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value_string} is not a number")
    if not (value > 0 and math.isfinite(value)):
        raise argparse.ArgumentTypeError(
            f"{value_string} is not a positive (finite) number"
        )
    return value


@functools.cache
def _known_feature_types():
    """Qualified names of all feature types known to pyGPlates (eg, "gpml:MidOceanRidge")."""
    return {
        getattr(pygplates.FeatureType, name).to_qualified_string()
        for name in dir(pygplates.FeatureType)
        if name.startswith("gpml_")
    }


def _feature_type(value_string):
    # Unqualified names (eg, "MidOceanRidge") default to the "gpml" namespace.
    # Reject other prefixes explicitly, since pyGPlates maps unknown prefixes to "gpml".
    namespace, _, name = value_string.rpartition(":")
    feature_type = (
        pygplates.FeatureType.create_gpml(name) if namespace in ("", "gpml") else None
    )
    if (
        feature_type is None
        or feature_type.to_qualified_string() not in _known_feature_types()
    ):
        raise argparse.ArgumentTypeError(
            f"{value_string} is not a GPlates feature type (eg, gpml:MidOceanRidge)"
        )
    return feature_type


def add_parser(parser):
    """Add command line argument parser."""

    cmd = parser.add_parser(
        "ridge-spreading-rate",
        aliases=("rsr",),
        help=help_str,
        add_help=True,
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    cmd.set_defaults(func=_run_ridge_spreading_rate)

    cmd.add_argument(
        "-r",
        "--rotation-filenames",
        dest="rotation_filenames",
        nargs="+",
        required=True,
        metavar="rotation_filename",
        help="One or more rotation files.",
    )
    cmd.add_argument(
        "-m",
        "--topology-filenames",
        dest="topology_filenames",
        nargs="+",
        required=True,
        metavar="topology_filename",
        help="One or more topology files.",
    )
    cmd.add_argument(
        "-a",
        "--anchor",
        type=int,
        default=0,
        dest="anchor_plate_id",
        help="Anchor plate id used for reconstructing. Defaults to zero.",
    )

    threshold_sampling_distance_group = cmd.add_mutually_exclusive_group()
    threshold_sampling_distance_group.add_argument(
        "-d",
        "--threshold-sampling-distance-degrees",
        dest="threshold_sampling_distance_degrees",
        type=_positive_float,
        help="Threshold sampling distance along ridges (in degrees). "
        f"Defaults to {DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES} degrees.",
    )
    threshold_sampling_distance_group.add_argument(
        "-k",
        "--threshold-sampling-distance-kms",
        dest="threshold_sampling_distance_kms",
        type=_positive_float,
        help="Threshold sampling distance along ridges (in Kms). "
        f"Defaults to {math.radians(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES) * pygplates.Earth.equatorial_radius_in_kms:.2f} Kms "
        f"(which is equivalent to {DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES} degrees).",
    )

    cmd.add_argument(
        "-t",
        "--time-range",
        dest="time_range",
        type=float,
        nargs=2,
        metavar=("young_time", "old_time"),
        default=[DEFAULT_TIME_RANGE_YOUNG_TIME, DEFAULT_TIME_RANGE_OLD_TIME],
        help="The time range (in Ma) from young time to old time. "
        f"Defaults to {DEFAULT_TIME_RANGE_YOUNG_TIME} -> {DEFAULT_TIME_RANGE_OLD_TIME} Ma.",
    )
    cmd.add_argument(
        "-i",
        "--time-increment",
        dest="time_increment",
        type=_positive_float,
        default=DEFAULT_TIME_INCREMENT,
        help=f"The time increment in My. Defaults to {DEFAULT_TIME_INCREMENT} My.",
    )
    cmd.add_argument(
        "-v",
        "--velocity-delta-time",
        dest="velocity_delta_time",
        type=_positive_float,
        default=DEFAULT_VELOCITY_DELTA_TIME,
        help="The delta time interval used to calculate velocities in My. "
        f"Defaults to {DEFAULT_VELOCITY_DELTA_TIME} My.",
    )

    cmd.add_argument(
        "-g",
        "--output-gpml-filename",
        dest="output_gpml_filename",
        help="Optional GPML output filename to contain the spreading rate data for all specified times. "
        "This can then be loaded into GPlates to display the data as scalar coverages.",
    )
    cmd.add_argument(
        "output_filename_prefix",
        metavar="output-filename-prefix",
        help="The output filename prefix. An output file is created for each geological time in the sequence where "
        "the filename suffix contains the time and the filename extension.",
    )
    cmd.add_argument(
        "-e",
        "--output-filename-extension",
        dest="output_filename_extension",
        default="xy",
        help='The output xy filename extension. Defaults to "xy".',
    )

    cmd.add_argument(
        "--use-ptt",
        action="store_true",
        help="Use the older PlateTectonicTools calculation of spreading rates.",
    )
    cmd.add_argument(
        "--divergence-threshold",
        metavar="CM_PER_YR",
        type=float,
        default=None,
        dest="divergence_threshold_in_cm_per_yr",
        help="Only output sample points whose orthogonal diverging velocity (in cm/yr) is at least this value. "
        "For example, 0 removes converging points. Can be negative. Cannot be used with '--use-ptt'. "
        "By default all (diverging and converging) points are output.",
    )

    transform_group = cmd.add_mutually_exclusive_group()
    transform_group.add_argument(
        "--transform-segment-deviation",
        metavar="DEGREES",
        type=_positive_float,
        default=separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_DEGREES,
        dest="transform_segment_deviation_degrees",
        help="How much the spreading direction can deviate from the segment normal (in degrees) before the segment "
        "is considered a transform segment (and skipped). Defaults to "
        f"{separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_DEGREES} degrees.",
    )
    transform_group.add_argument(
        "--include-transform-segments",
        action="store_true",
        help="Don't separate ridges into ridge and transform segments (sample the entire ridge).",
    )

    feature_types_group = cmd.add_mutually_exclusive_group()
    feature_types_group.add_argument(
        "--spreading-feature-types",
        metavar="FEATURE_TYPE",
        nargs="+",
        type=_feature_type,
        default=DEFAULT_SPREADING_FEATURE_TYPES,
        help="Only sample plate boundaries of these feature types (eg, gpml:MidOceanRidge). "
        "Defaults to "
        + " ".join(t.to_qualified_string() for t in DEFAULT_SPREADING_FEATURE_TYPES)
        + ".",
    )
    feature_types_group.add_argument(
        "--all-plate-boundaries",
        action="store_true",
        help="Sample all plate boundaries, regardless of feature type. "
        "For example, combine with '--include-transform-segments' and '--divergence-threshold 0' "
        "to sample all diverging plate boundaries. "
        "Note that '--use-ptt' still only samples spreading feature types (such as mid-ocean ridges).",
    )

    cmd.add_argument(
        "--include-network-boundaries",
        action="store_true",
        help="Also sample deforming network boundaries that aren't plate boundaries. "
        "By default only plate boundaries are sampled.",
    )


def _create_coverage_feature(ridge_data, time):
    """Create a feature with a coverage geometry containing the ridge spreading data at 'time'."""
    coverage_feature = pygplates.Feature()

    # Make it only appear at 'time'.
    coverage_feature.set_valid_time(time + 0.5, time - 0.5)

    coverage_geometry = pygplates.MultiPointOnSphere(
        (lat, lon) for lon, lat, *_ in ridge_data
    )
    coverage_scalars = {
        pygplates.ScalarType.create_gpml(name): [row[column] for row in ridge_data]
        for column, name in enumerate(_COVERAGE_SCALAR_NAMES, start=2)
    }
    coverage_feature.set_geometry((coverage_geometry, coverage_scalars))

    return coverage_feature


def _run_ridge_spreading_rate(args):
    young_time, old_time = args.time_range
    if young_time > old_time:
        raise ValueError(
            f"The young time {young_time} is older (larger) than the old time {old_time}."
        )

    if args.use_ptt and args.divergence_threshold_in_cm_per_yr is not None:
        raise ValueError("'--divergence-threshold' cannot be used with '--use-ptt'.")

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
            DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES
        )

    if args.include_transform_segments:
        transform_segment_deviation_radians = None
    else:
        transform_segment_deviation_radians = math.radians(
            args.transform_segment_deviation_degrees
        )

    if args.all_plate_boundaries:
        spreading_feature_types = None
    else:
        spreading_feature_types = args.spreading_feature_types

    plate_reconstruction = PlateReconstruction(
        args.rotation_filenames,
        args.topology_filenames,
        anchor_plate_id=args.anchor_plate_id,
    )

    coverage_features = []

    time = young_time
    while time <= pygplates.GeoTimeInstant(old_time):
        ridge_data = plate_reconstruction.tessellate_mid_ocean_ridges(
            time,
            threshold_sampling_distance_radians,
            use_ptt=args.use_ptt,
            spreading_feature_types=spreading_feature_types,
            transform_segment_deviation_in_radians=transform_segment_deviation_radians,
            include_network_boundaries=args.include_network_boundaries,
            divergence_threshold_in_cm_per_yr=args.divergence_threshold_in_cm_per_yr,
            output_obliquity_and_normal_and_left_right_plates=True,
            velocity_delta_time=args.velocity_delta_time,
        )

        if len(ridge_data):
            # The returned array is all floats, so write the plate IDs (columns 6 and 7) as integers.
            output_data = [
                sample[:6].tolist() + [int(sample[6]), int(sample[7])]
                for sample in ridge_data
            ]

            output_filename = f"{args.output_filename_prefix}_{time:0.2f}.{args.output_filename_extension}"
            with open(output_filename, "w") as output_file:
                for row in output_data:
                    output_file.write(" ".join(str(item) for item in row) + "\n")

            if args.output_gpml_filename:
                coverage_features.append(_create_coverage_feature(output_data, time))

        time += args.time_increment

    if args.output_gpml_filename:
        pygplates.FeatureCollection(coverage_features).write(args.output_gpml_filename)
