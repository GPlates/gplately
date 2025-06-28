#
#    Copyright (C) 2014-2025 The University of Sydney, Australia
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

DEFAULT_OUTPUT_FILENAME_SUFFIX = "_fixed_crossovers"


def fix_crossovers(
    rotation_feature_collections,
    crossover_threshold_degrees,
    crossover_type_function,
    ignore_moving_plates,
    debug,
):
    """
    Fix crossovers.
    If any errors occurred we will still write to the output files.
    """
    crossover_filter = None
    # Ignore a subset of moving plates if requested.
    if ignore_moving_plates:
        crossover_filter = (
            lambda crossover: crossover.moving_plate_id not in ignore_moving_plates
        )

    # Synchronise crossovers.
    crossover_results = []
    synchronise_crossovers_success = pygplates.synchronise_crossovers(
        rotation_feature_collections,
        crossover_filter,
        crossover_threshold_degrees,
        crossover_type_function,
        crossover_results,
    )

    # Print debug output if requested.
    if debug:
        num_crossovers_error = 0
        num_crossovers_ignored = 0
        num_crossovers_passed = 0
        num_crossovers_corrected = 0

        for crossover_result in crossover_results:
            crossover = crossover_result[0]
            result = crossover_result[1]

            if crossover.type == pygplates.CrossoverType.synch_old_crossover_and_stages:
                type_str = "synch old crossover and stages"
            elif crossover.type == pygplates.CrossoverType.synch_old_crossover_only:
                type_str = "synch old crossover only"
            elif (
                crossover.type
                == pygplates.CrossoverType.synch_young_crossover_and_stages
            ):
                type_str = "synch young crossover and stages"
            elif crossover.type == pygplates.CrossoverType.synch_young_crossover_only:
                type_str = "synch young crossover only"
            elif crossover.type == pygplates.CrossoverType.ignore:
                type_str = "ignore"
            else:
                type_str = "unknown"

            if result == pygplates.CrossoverResult.not_synchronised:
                num_crossovers_passed += 1
                result_str = "passed"
            elif result == pygplates.CrossoverResult.synchronised:
                num_crossovers_corrected += 1
                result_str = "corrected"
            elif result == pygplates.CrossoverResult.ignored:
                num_crossovers_ignored += 1
                result_str = "ignored"
            else:
                num_crossovers_error += 1
                result_str = "error"

            print(
                "Time({0}), moving_pid({1}), young_fixed_pid({2}), old_fixed_pid({3}), type({4}): {5}".format(
                    crossover.time,
                    crossover.moving_plate_id,
                    crossover.young_crossover_fixed_plate_id,
                    crossover.old_crossover_fixed_plate_id,
                    type_str,
                    result_str,
                )
            )

        print("Results:")
        print("  Total number of crossovers = {0}".format(len(crossover_results)))
        print("  Total errors = {0}".format(num_crossovers_error))
        print("  Total ignored = {0}".format(num_crossovers_ignored))
        print("  Total passed = {0}".format(num_crossovers_passed))
        print("  Total corrected = {0}".format(num_crossovers_corrected))

    return synchronise_crossovers_success


def parse_positive_number(value_string):
    """parse and return a positive number"""
    try:
        value = float(value_string)
    except ValueError:
        raise argparse.ArgumentTypeError("%s is not a number" % value_string)

    if value < 0:
        raise argparse.ArgumentTypeError("%g is not a positive number" % value)

    return value


def add_arguments(parser: argparse.ArgumentParser):
    """add command line argument parser"""

    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.description = __description__

    parser.set_defaults(func=main)

    parser.add_argument(
        "-d", "--debug", action="store_true", help="Print debug output."
    )
    parser.add_argument(
        "-c",
        "--crossover_threshold_degrees",
        type=parse_positive_number,
        help="If specified then crossovers are fixed only if post-crossover rotation latitude, "
        "longitude or angle differ from those in pre-crossover rotation by the specified amount "
        "(in degrees). This is useful for some PLATES rotation files that are typically accurate "
        "to 2 decimal places (or threshold of 0.01).",
    )

    # Can specify only one of '-x' or '-g'.
    crossover_type_group = parser.add_mutually_exclusive_group()
    crossover_type_group.add_argument(
        "-x",
        "--default_xo_ys",
        action="store_true",
        dest="crossover_type_default_xo_ys",
        help="If specified, then if a crossover's type is unknown it will default to "
        '"synch old crossover and stages", which is equivalent to the "@xo_ys" comment tag.',
    )
    crossover_type_group.add_argument(
        "-g",
        "--default_xo_ig",
        action="store_true",
        dest="crossover_type_default_xo_ig",
        help="If specified, then if a crossover's type is unknown it will default to "
        'ignoring the crossover, which is equivalent to the "@xo_ig" comment tag.',
    )

    parser.add_argument(
        "-i",
        "--ignore_moving_plates",
        type=parse_positive_number,
        nargs="+",
        metavar="MOVING_PLATE_ID",
        help="If specified then is a list of moving plate ids to ignore when fixing crossovers.",
    )
    parser.add_argument(
        "-s",
        "--output_filename_suffix",
        type=str,
        default="{0}".format(DEFAULT_OUTPUT_FILENAME_SUFFIX),
        help="The suffix to append to each input rotation filename to get each output rotation "
        "filename - the default suffix is '{0}'".format(DEFAULT_OUTPUT_FILENAME_SUFFIX),
    )

    parser.add_argument(
        "input_rotation_filenames",
        type=str,
        nargs="+",
        metavar="input_rotation_filename",
        help="One or more input rotation filenames (original files).",
    )


__description__ = """Loads one or more input rotation files, fixes any crossovers and saves the rotations to \
output rotation files.

    The name of each output file is the associated input filename with a suffix appended. For example:
       'rotations/input_rotations.rot' -> 'rotations/input_rotations{0}.rot'

    The method used to synchronise a crossover depends on following strings found in each 'young'
    crossover pole:
      * @xo_ig : Ignore the crossover. Has the same effect as not specifying any tag, except
                 it avoids a warning/error message. All finite rotations in the young and old
                 crossover sequences are preserved.
      * @xo_ys : All finite rotations in the old crossover sequence will be synchronised
                 (such that old stage rotations are preserved). All finite rotations in the
                 young crossover sequence are preserved.
      * @xo_yf : Only the crossover finite rotation in the old crossover sequence will be
                 synchronised (such that the older finite rotations are preserved).
                 All finite rotations in the young crossover sequence are preserved.
      * @xo_os : All finite rotations in the young crossover sequence will be synchronised
                 (such that young stage rotations are preserved). All finite rotations in
                 old crossover sequence are preserved. Note: This can result in non-zero
                 finite rotations at present day if the younger sequence includes present day.
      * @xo_of : Only the crossover finite rotation in the young crossover sequence will be
                 synchronised (such that the younger finite rotations are preserved).
                 All finite rotations in the old crossover sequence are preserved.
    
    ...if any of the above tags are missing in a crossover then it will not be processed and
    a warning/error message will be printed.
    However the '-x' or '-g' option can optionally be used to default to @xo_ys or @xo_ig behaviour
    (respectively) for each crossover that does not have any of the above text strings (see below).

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    %(prog)s -d -c 0.01 -i 201 701 -- input_rotations1.rot input_rotations2.rot
    """.format(
    DEFAULT_OUTPUT_FILENAME_SUFFIX
)


def main(args):

    file_registry = pygplates.FeatureCollectionFileFormatRegistry()

    # Read/parse the input rotation feature collections.
    rotation_feature_collections = [
        file_registry.read(input_rotation_filename)
        for input_rotation_filename in args.input_rotation_filenames
    ]

    # Whether to crossover types should default to '@xo_ys' or '@xo_ig' if type not found...
    if args.crossover_type_default_xo_ys:
        crossover_type_function = (
            pygplates.CrossoverTypeFunction.type_from_xo_tags_in_comment_default_xo_ys
        )
    elif args.crossover_type_default_xo_ig:
        crossover_type_function = (
            pygplates.CrossoverTypeFunction.type_from_xo_tags_in_comment_default_xo_ig
        )
    else:
        crossover_type_function = (
            pygplates.CrossoverTypeFunction.type_from_xo_tags_in_comment
        )

    # Fix crossovers.
    # If anfix_crossoversed we will still write to the output files.
    if not fix_crossovers(
        rotation_feature_collections,
        args.crossover_threshold_degrees,
        crossover_type_function,
        args.ignore_moving_plates,
        args.debug,
    ):
        print(
            "Warning: One or more crossovers were not processed since unable to determine crossover "
            "correction method or infinite cycle detected.",
            file=sys.stderr,
        )

    # Write the modified rotation feature collections to disk.
    for rotation_feature_collection_index in range(len(rotation_feature_collections)):
        rotation_feature_collection = rotation_feature_collections[
            rotation_feature_collection_index
        ]

        # Each output filename is the input filename with a suffix appended.
        input_rotation_filename = args.input_rotation_filenames[
            rotation_feature_collection_index
        ]
        filename_root, filename_ext = os.path.splitext(input_rotation_filename)
        output_rotation_filename = "".join(
            (filename_root, args.output_filename_suffix, filename_ext)
        )

        file_registry.write(rotation_feature_collection, output_rotation_filename)


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
