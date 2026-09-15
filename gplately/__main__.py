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
from importlib.resources import files
import os
import sys
from typing import List

import pygplates

from gplately import __version__

from .commands import (
    seafloor_grids,
    feature_filter_cmd,
    list_models,
    regrid,
    reset_feature_type,
    rotate_grid,
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


def _print_cli_config_example():
    print(files("gplately").joinpath("data", "gplately-cli-config.toml").read_text())


class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"error: {message}\n")
        self.print_help()
        sys.exit(1)


def _hide_alias(subparser, alias: str, cmd: argparse.ArgumentParser):
    """Register *alias* as a working subcommand name without showing it in --help.

    argparse's `aliases=` shows every alias in the subcommand list and usage
    line; the old underscore-separated names (issue #450) should keep working
    but stay out of --help, so they're added directly to the subparsers
    action's name->parser map instead.
    """
    subparser._name_parser_map[alias] = cmd


def main():
    parser = ArgParser()

    parser.add_argument("-v", "--version", action="store_true")

    # sub-commands
    subparser = parser.add_subparsers(
        dest="command",
        title="subcommands",
        description="valid subcommands",
        metavar="<subcommand>",
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
    combine_cmd.set_defaults(func=_run_combine_feature_collections)
    combine_cmd.add_argument("combine_first_input_file", type=str)
    combine_cmd.add_argument("combine_other_input_files", nargs="+", type=str)
    combine_cmd.add_argument("combine_output_file", type=str)

    # add "feature filter" sub-command
    feature_filter_cmd.add_parser(subparser)

    # add "reset feature type" sub-command
    reset_feature_type.add_parser(subparser)

    # add "create age grids" sub-command
    seafloor_grids.add_parser(subparser)

    # add "regrid" sub-command
    regrid.add_parser(subparser)

    # add "rotate_grid" sub-command
    rotate_grid.add_parser(subparser)

    # add "fix crossovers" sub-command
    fix_crossovers_cmd = subparser.add_parser(
        "fix-crossovers",
        aliases=("fc",),
        help="Loads one or more input rotation files, fixes any crossovers and saves the rotations to output rotation files.",
        add_help=True,
    )
    _hide_alias(subparser, "fix_crossovers", fix_crossovers_cmd)
    fix_crossovers.add_arguments(fix_crossovers_cmd)

    # add "remove plate rotations" sub-command
    remove_plate_rotations_cmd = subparser.add_parser(
        "remove-rotations",
        aliases=("rr",),
        help="Remove one or more plate IDs from a rotation model (consisting of one or more rotation files).",
        add_help=True,
    )
    _hide_alias(subparser, "remove_rotations", remove_plate_rotations_cmd)
    remove_plate_rotations.add_arguments(remove_plate_rotations_cmd)

    # add "cleanup topologies" sub-command
    cleanup_topologies_cmd = subparser.add_parser(
        "cleanup-topologies",
        aliases=("ct",),
        help="Remove any regular features not referenced by topological features.",
        add_help=True,
    )
    _hide_alias(subparser, "cleanup_topologies", cleanup_topologies_cmd)
    cleanup_topologies.add_arguments(cleanup_topologies_cmd)

    # add "convert_xy_to_gplates" sub-command
    convert_xy_to_gplates_cmd = subparser.add_parser(
        "convert-xy-to-gplates",
        aliases=("cxg",),
        help="Converts geometry in one or more input ascii files (such as '.xy' files) to output files suitable for loading into GPlates.",
        add_help=True,
    )
    _hide_alias(subparser, "convert_xy_to_gplates", convert_xy_to_gplates_cmd)
    convert_xy_to_gplates.add_arguments(convert_xy_to_gplates_cmd)

    # add "diagnose_rotations" sub-command
    diagnose_rotations_cmd = subparser.add_parser(
        "diagnose-rotations",
        aliases=("dr",),
        help="Diagnose one or more rotation files to check for inconsistencies.",
        add_help=True,
    )
    _hide_alias(subparser, "diagnose_rotations", diagnose_rotations_cmd)
    diagnose_rotations.add_arguments(diagnose_rotations_cmd)

    # add "resolve_topologies" sub-command
    resolve_topologies_cmd = subparser.add_parser(
        "resolve-topologies",
        aliases=("rt",),
        help="Resolve topological plate polygons (and deforming networks) and saves (to separate files) the resolved topologies, and their boundary sections as subduction zones, mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges).",
        add_help=True,
    )
    _hide_alias(subparser, "resolve_topologies", resolve_topologies_cmd)
    resolve_topologies.add_arguments(resolve_topologies_cmd)

    # add "rotation_tools" sub-command
    rotation_tools_cmd = subparser.add_parser(
        "rotation-tools",
        aliases=("rots",),
        help="Calculate stage rotations between consecutive finite rotations in plate pairs.",
        add_help=True,
    )
    _hide_alias(subparser, "rotation_tools", rotation_tools_cmd)
    rotation_tools.add_arguments(rotation_tools_cmd)

    # add "separate_ridge_transform_segments" sub-command
    separate_ridge_transform_segments_cmd = subparser.add_parser(
        "separate-ridge-transform-segments",
        aliases=("srts",),
        help="Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments.",
        add_help=True,
    )
    _hide_alias(
        subparser,
        "separate_ridge_transform_segments",
        separate_ridge_transform_segments_cmd,
    )
    separate_ridge_transform_segments.add_arguments(
        separate_ridge_transform_segments_cmd
    )

    # add "subduction_convergence" sub-command
    subduction_convergence_cmd = subparser.add_parser(
        "subduction-convergence",
        aliases=("sc",),
        help="Find the convergence rates along trenches (subduction zones) over time.",
        add_help=True,
    )
    _hide_alias(subparser, "subduction_convergence", subduction_convergence_cmd)
    subduction_convergence.add_arguments(subduction_convergence_cmd)

    # add "gpmdb" sub-command
    gpmdb_cmd = subparser.add_parser(
        "gpmdb",
        help="Retrieve paleomagnetic data from https://www.gpmdb.net, create GPlates-compatible VGP features and save the VGP features in a .gpmlz file.",
        add_help=True,
    )
    gpmdb.add_arguments(gpmdb_cmd)

    # add "get_cli_config_example" sub-command
    get_cli_config_example_cmd = subparser.add_parser(
        "get-cli-config-example",
        aliases=("gcce",),
        help="Print an example CLI configuration in TOML format to stdout (for use with --config); "
        "redirect it to save, e.g. 'gplately get-cli-config-example > my-gplately-cli-config.toml'",
        add_help=True,
    )
    _hide_alias(subparser, "get_cli_config_example", get_cli_config_example_cmd)
    get_cli_config_example_cmd.set_defaults(
        func=lambda args: _print_cli_config_example()
    )

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
