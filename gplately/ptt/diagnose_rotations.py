#
#    Copyright (C) 2020-2025 The University of Sydney, Australia
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

from __future__ import print_function

import argparse
import os.path
import sys

import pygplates

# The default is useful when the rotations were loaded from a PLATES rotation file
# that stored rotation lat/lon/angle to 2 decimal places of accuracy.
DEFAULT_ROTATION_THRESHOLD_DEGREES = 0.01


def diagnose_rotations(
    rotation_features, rotation_threshold_degrees=DEFAULT_ROTATION_THRESHOLD_DEGREES
):
    """
    Diagnose one or more rotation files to check for inconsistencies.

    'rotation_features' can be a rotation feature collection, or rotation filename, or rotation feature, or
    sequence of rotation features, or a sequence (eg, list or tuple) of any combination of those four types.

    'rotation_threshold_degrees' allows two rotations to compare equal if their rotation latitude, longitude or angle
    differ by less than the specified amount (in degrees). The default value is useful for some PLATES rotation files
    that are typically accurate to 2 decimal places (or threshold of 0.01).
    Specifing 'None' is equivalent to having no threshold (see pygplates.FiniteRotation.are_equal).
    """

    # Use helper class to convert 'rotation_features' argument to a list of features.
    rotation_features = pygplates.FeaturesFunctionArgument(rotation_features)
    rotation_feature_sequence = rotation_features.get_features()

    # Make sure threshold is a number.
    if rotation_threshold_degrees is not None:
        rotation_threshold_degrees = float(rotation_threshold_degrees)

    # A 'dict' to map each moving plate to a list of total reconstruction poles
    # (one per moving/fixed plate pair)
    total_reconstruction_poles_by_moving_plate = {}

    # Get the moving/fixed total reconstruction poles.
    for rotation_feature in rotation_feature_sequence:
        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        # If the current feature is a valid rotation feature...
        if total_reconstruction_pole:
            (
                fixed_plate_id,
                moving_plate_id,
                rotation_sequence,
            ) = total_reconstruction_pole
            # Each moving plate has a list of total reconstruction poles.
            total_reconstruction_poles = (
                total_reconstruction_poles_by_moving_plate.setdefault(
                    moving_plate_id, []
                )
            )
            total_reconstruction_poles.append(total_reconstruction_pole)

    # Iterate over the moving plates.
    for (
        moving_plate_id,
        total_reconstruction_poles,
    ) in total_reconstruction_poles_by_moving_plate.items():
        moving_plate_time_samples = []

        # Iterate over the total reconstruction poles associated with the current moving plate.
        for total_reconstruction_pole in total_reconstruction_poles:
            fixed_plate_id = total_reconstruction_pole[0]
            rotation_sequence = total_reconstruction_pole[2]

            # Skip moving/fixed plates 999 since they are used as comments in PLATES rotation file.
            if moving_plate_id == 999 and fixed_plate_id == 999:
                continue

            time_samples = rotation_sequence.get_enabled_time_samples()
            if not time_samples:
                print(
                    "Moving_pid({0}), fixed_pid({1}): no enabled poles in sequence.".format(
                        moving_plate_id, fixed_plate_id
                    )
                )
                continue
            if len(time_samples) == 1:
                print(
                    "Moving_pid({0}), fixed_pid({1}): only one enabled pole in sequence.".format(
                        moving_plate_id, fixed_plate_id
                    )
                )
                continue
            times = [time_sample.get_time() for time_sample in time_samples]
            if sorted(times) != times:
                print(
                    "Moving_pid({0}), fixed_pid({1}): times of enabled samples not monotonically increasing between {2} and {3}.".format(
                        moving_plate_id, fixed_plate_id, times[0], times[-1]
                    )
                )
                continue

            moving_plate_time_samples.append((fixed_plate_id, time_samples))

        # Go through all rotation sequences with the same moving plate (and potentially different
        # fixed plates) to see if there's any overlap in their time ranges and to see if two
        # adjacent rotation sequences have matching rotations (or crossovers if fixed plates differ)
        # at the joint time.
        for index1 in range(len(moving_plate_time_samples) - 1):
            for index2 in range(index1 + 1, len(moving_plate_time_samples)):
                fixed_plate_id1 = moving_plate_time_samples[index1][0]
                fixed_plate_id2 = moving_plate_time_samples[index2][0]
                fixed_plate_time_samples1 = moving_plate_time_samples[index1][1]
                fixed_plate_time_samples2 = moving_plate_time_samples[index2][1]
                time_range1 = (
                    pygplates.GeoTimeInstant(fixed_plate_time_samples1[0].get_time()),
                    pygplates.GeoTimeInstant(fixed_plate_time_samples1[-1].get_time()),
                )
                time_range2 = (
                    pygplates.GeoTimeInstant(fixed_plate_time_samples2[0].get_time()),
                    pygplates.GeoTimeInstant(fixed_plate_time_samples2[-1].get_time()),
                )

                if time_range1[0] > time_range2[0]:
                    if time_range1[0] < time_range2[1]:
                        print(
                            "Moving_pid({0}), fixed_pid1({1}), fixed_pid2({2}): sequences overlap between {3} and {4}.".format(
                                moving_plate_id,
                                fixed_plate_id1,
                                fixed_plate_id2,
                                time_range1[0],
                                time_range2[1],
                            )
                        )
                    elif time_range1[0] == time_range2[1]:
                        # Sequence at 'index1' is *older* than sequence at 'index2' and abuts it.
                        if fixed_plate_id1 == fixed_plate_id2:
                            # Check that both rotations are equal.
                            if not pygplates.FiniteRotation.are_equal(
                                fixed_plate_time_samples1[0]
                                .get_value()
                                .get_finite_rotation(),
                                fixed_plate_time_samples2[-1]
                                .get_value()
                                .get_finite_rotation(),
                                rotation_threshold_degrees,
                            ):
                                print(
                                    "Moving_pid({0}), fixed_pid({1}): two different rotations at time {2}.".format(
                                        moving_plate_id, fixed_plate_id1, time_range1[0]
                                    )
                                )
                        else:
                            # We have a crossover (different fixed plates).
                            # Check that the crossover is synchronised.
                            # Implementation follows that in pygplates.synchronise_crossovers():
                            #  1. Disable crossover in older rotation sequence,
                            #  2. Calculate rotation relative to older crossover's fixed plate,
                            #  3. Compare calculated rotation with previous rotation,
                            #  4. Re-enable crossover in older rotation sequence.
                            fixed_plate_time_samples1[0].set_disabled()
                            rotation_model = pygplates.RotationModel(
                                rotation_feature_sequence
                            )
                            if not pygplates.FiniteRotation.are_equal(
                                rotation_model.get_rotation(
                                    time_range1[0],
                                    moving_plate_id,
                                    anchor_plate_id=fixed_plate_id1,
                                ),
                                fixed_plate_time_samples1[0]
                                .get_value()
                                .get_finite_rotation(),
                                rotation_threshold_degrees,
                            ):
                                print(
                                    "Moving_pid({0}), young_fixed_pid({1}), old_fixed_pid({2}): "
                                    "crossover at time {3} needs fixing.".format(
                                        moving_plate_id,
                                        fixed_plate_id2,
                                        fixed_plate_id1,
                                        time_range1[0],
                                    )
                                )
                            fixed_plate_time_samples1[0].set_enabled()
                else:  # time_range1[0] <= time_range2[0] ...
                    if time_range2[0] < time_range1[1]:
                        print(
                            "Moving_pid({0}), fixed_pid1({1}), fixed_pid2({2}): sequences overlap between {3} and {4}.".format(
                                moving_plate_id,
                                fixed_plate_id1,
                                fixed_plate_id2,
                                time_range2[0],
                                time_range1[1],
                            )
                        )
                    elif time_range2[0] == time_range1[1]:
                        # Sequence at 'index1' is *younger* than sequence at 'index2' and abuts it.
                        if fixed_plate_id1 == fixed_plate_id2:
                            # Check that both rotations are equal.
                            if not pygplates.FiniteRotation.are_equal(
                                fixed_plate_time_samples1[-1]
                                .get_value()
                                .get_finite_rotation(),
                                fixed_plate_time_samples2[0]
                                .get_value()
                                .get_finite_rotation(),
                                rotation_threshold_degrees,
                            ):
                                print(
                                    "Moving_pid({0}), fixed_pid({1}): two different rotations at time {2}.".format(
                                        moving_plate_id, fixed_plate_id1, time_range2[0]
                                    )
                                )
                        else:
                            # We have a crossover (different fixed plates).
                            # Check that the crossover is synchronised.
                            # Implementation follows that in pygplates.synchronise_crossovers():
                            #  1. Disable crossover in older rotation sequence,
                            #  2. Calculate rotation relative to older crossover's fixed plate,
                            #  3. Compare calculated rotation with previous rotation,
                            #  4. Re-enable crossover in older rotation sequence.
                            fixed_plate_time_samples2[0].set_disabled()
                            rotation_model = pygplates.RotationModel(
                                rotation_feature_sequence
                            )
                            if not pygplates.FiniteRotation.are_equal(
                                rotation_model.get_rotation(
                                    time_range2[0],
                                    moving_plate_id,
                                    anchor_plate_id=fixed_plate_id2,
                                ),
                                fixed_plate_time_samples2[0]
                                .get_value()
                                .get_finite_rotation(),
                                rotation_threshold_degrees,
                            ):
                                print(
                                    "Moving_pid({0}), young_fixed_pid({1}), old_fixed_pid({2}): "
                                    "crossover at time {3} needs fixing.".format(
                                        moving_plate_id,
                                        fixed_plate_id1,
                                        fixed_plate_id2,
                                        time_range2[0],
                                    )
                                )
                            fixed_plate_time_samples2[0].set_enabled()


__description__ = """Diagnose one or more rotation files to check for inconsistencies.

NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
For example...

%(prog)s input_rotations1.rot input_rotations2.rot
"""


def add_arguments(parser: argparse.ArgumentParser):
    """add command line argument parser"""

    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.description = __description__

    parser.set_defaults(func=main)

    def parse_positive_number(value_string):
        try:
            value = float(value_string)
        except ValueError:
            raise argparse.ArgumentTypeError("%s is not a number" % value_string)

        if value < 0:
            raise argparse.ArgumentTypeError("%g is not a positive number" % value)

        return value

    parser.add_argument(
        "-t",
        "--rotation_threshold_degrees",
        type=parse_positive_number,
        default=DEFAULT_ROTATION_THRESHOLD_DEGREES,
        help="Two rotations differ if either the rotation latitude, longitude or angle differ by "
        "the specified amount (in degrees). The default (0.01 degrees) is useful for some "
        "PLATES rotation files that are typically accurate to 2 decimal places.",
    )

    parser.add_argument(
        "rotation_filenames",
        type=str,
        nargs="+",
        metavar="rotation_filename",
        help="One or more rotation filenames.",
    )


def main(args):
    diagnose_rotations(args.rotation_filenames, args.rotation_threshold_degrees)


if __name__ == "__main__":

    # The command-line parser.
    parser = argparse.ArgumentParser(
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # add arguments
    add_arguments(parser)

    # Parse command-line options.
    args = parser.parse_args()

    # call main function
    main(args)
