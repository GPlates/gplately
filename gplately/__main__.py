#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
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
import os
import sys
from typing import List

import pygplates

from gplately import __version__

from .commands import (
    create_age_grids,
    feature_filter,
    list_models,
    regrid,
    reset_feature_type,
)
from .ptt import (
    cleanup_topologies,
    convert_xy_to_gplates,
    diagnose_rotations,
    fix_crossovers,
    gpmdb,
    remove_plate_rotations,
    resolve_topologies,
    rotation_tools,
    separate_ridge_transform_segments,
    subduction_convergence,
)


def combine_feature_collections(input_files: List[str], output_file: str):
    """Combine multiple feature collections into one.

    Usage example: gplately combine input_file_1 input_file_2 input_file_3 output_file
    """
    feature_collection = pygplates.FeatureCollection()
    for file in input_files:
        if not os.path.isfile(file):
            raise Exception(f"{file} is not a file.")
        feature_collection.add(pygplates.FeatureCollection(file))

    feature_collection.write(output_file)

    print(f"Done! The combined feature collection has been saved to {output_file}.")


def _run_combine_feature_collections(args):
    combine_feature_collections(
        [args.combine_first_input_file] + args.combine_other_input_files,
        args.combine_output_file,
    )


class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"error: {message}\n")
        self.print_help()
        sys.exit(1)


def main():
    parser = ArgParser()

    parser.add_argument("-v", "--version", action="store_true")

    # sub-commands
    subparser = parser.add_subparsers(
        dest="command",
        title="subcommands",
        description="valid subcommands",
    )
    # add "list models" sub-command
    list_models.add_parser(subparser)

    # add "combine feature" sub-command
    combine_cmd = subparser.add_parser(
        "combine",
        help="Combine multiple feature collections into one.",
        description=combine_feature_collections.__doc__,
    )
    combine_cmd.formatter_class = argparse.RawDescriptionHelpFormatter

    # add "feature filter" sub-command
    feature_filter.add_parser(subparser)

    # add "reset feature type" sub-command
    reset_feature_type.add_parser(subparser)

    # add "create age grids" sub-command
    create_age_grids.add_parser(subparser)

    # add "regrid" sub-command
    regrid.add_parser(subparser)

    # add "fix crossovers" sub-command
    fix_crossovers_cmd = subparser.add_parser(
        "fix_crossovers",
        help="Loads one or more input rotation files, fixes any crossovers and saves the rotations to output rotation files.",
        add_help=True,
    )
    fix_crossovers.add_arguments(fix_crossovers_cmd)

    # add "remove plate rotations" sub-command
    remove_plate_rotations_cmd = subparser.add_parser(
        "remove_rotations",
        help="Remove one or more plate IDs from a rotation model (consisting of one or more rotation files).",
        add_help=True,
    )
    remove_plate_rotations.add_arguments(remove_plate_rotations_cmd)

    # add "cleanup topologies" sub-command
    cleanup_topologies_cmd = subparser.add_parser(
        "cleanup_topologies",
        help="Remove any regular features not referenced by topological features.",
        add_help=True,
    )
    cleanup_topologies.add_arguments(cleanup_topologies_cmd)

    # add convert_xy_to_gplates sub-command
    convert_xy_to_gplates_cmd = subparser.add_parser(
        "convert_xy_to_gplates",
        help="Converts geometry in one or more input ascii files (such as '.xy' files) to output files suitable for loading into GPlates.",
        add_help=True,
    )
    convert_xy_to_gplates.add_arguments(convert_xy_to_gplates_cmd)

    # add diagnose_rotations sub-command
    diagnose_rotations_cmd = subparser.add_parser(
        "diagnose_rotations",
        help="Diagnose one or more rotation files to check for inconsistencies.",
        add_help=True,
    )
    diagnose_rotations.add_arguments(diagnose_rotations_cmd)

    # add resolve_topologies sub-command
    resolve_topologies_cmd = subparser.add_parser(
        "resolve_topologies",
        help="Resolve topological plate polygons (and deforming networks) and saves (to separate files) the resolved topologies, and their boundary sections as subduction zones, mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges).",
        add_help=True,
    )
    resolve_topologies.add_arguments(resolve_topologies_cmd)

    # add rotation_tools sub-command
    rotation_tools_cmd = subparser.add_parser(
        "rotation_tools",
        help="Calculate stage rotations between consecutive finite rotations in plate pairs.",
        add_help=True,
    )
    rotation_tools.add_arguments(rotation_tools_cmd)

    # add separate_ridge_transform_segments sub-command
    separate_ridge_transform_segments_cmd = subparser.add_parser(
        "separate_ridge_transform_segments",
        help="Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments.",
        add_help=True,
    )
    separate_ridge_transform_segments.add_arguments(
        separate_ridge_transform_segments_cmd
    )

    # add subduction_convergence sub-command
    subduction_convergence_cmd = subparser.add_parser(
        "subduction_convergence",
        help="Find the convergence rates along trenches (subduction zones) over time.",
        add_help=True,
    )
    subduction_convergence.add_arguments(subduction_convergence_cmd)

    # add gpmdb sub-command
    gpmdb_cmd = subparser.add_parser(
        "gpmdb",
        help="Retrieve paleomagnetic data from https://www.gpmdb.net, create GPlates-compatible VGP features and save the VGP features in a .gpmlz file.",
        add_help=True,
    )
    gpmdb.add_arguments(gpmdb_cmd)

    # combine command arguments
    combine_cmd.set_defaults(func=_run_combine_feature_collections)
    combine_cmd.add_argument("combine_first_input_file", type=str)
    combine_cmd.add_argument("combine_other_input_files", nargs="+", type=str)
    combine_cmd.add_argument("combine_output_file", type=str)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.version:
        print(__version__)
        sys.exit(0)

    args.func(args)


if __name__ == "__main__":
    main()
