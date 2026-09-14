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
import logging

from ..grids.paleobathymetry import AGE_DEPTH_MODELS, simple_paleobathymetry
from .sediment_thickness import (
    _add_common_arguments,
    _add_distance_arguments,
    _resolve_age_grid_filenames_and_times,
    _resolve_rotation_topology_proximity_files,
)

_logger = logging.getLogger("gplately")


def _run_paleobathymetry(args):
    age_grid_filenames_and_times, plate_model = _resolve_age_grid_filenames_and_times(
        args
    )
    rotation_files, topology_files, proximity_files = (
        _resolve_rotation_topology_proximity_files(args, plate_model)
    )

    simple_paleobathymetry(
        rotation_model=rotation_files,
        proximity_features=proximity_files,
        topological_features=topology_files,
        age_grid_filenames_and_times=age_grid_filenames_and_times,
        age_depth_model=args.age_depth_model,
        grid_spacing=args.grid_spacing,
        time_increment=args.time_step,
        max_reconstruction_time=args.max_reconstruction_time,
        anchor_plate_id=args.anchor_plate_id or 0,
        clamp_distance_km=args.clamp_distance_km,
        richards_table_filename=args.richards_table,
        output_directory=args.output_dir,
    )
    _logger.info(f"Paleobathymetry grids written to {args.output_dir}")


def add_parser(parser):
    """add command line argument parser for 'paleobathymetry'"""

    cmd = parser.add_parser(
        "paleobathymetry",
        aliases=("pb",),
        help="Run the simple_paleobathymetry workflow (Steps 1-4) end to end.",
        add_help=True,
        description=(
            "Reconstruct paleobathymetry of ocean crust from a seafloor-age grid: age -> "
            "basement depth (Step 1), distance to the nearest passive continental margin "
            "(Step 2), predicted sediment thickness (Step 3), and isostatically-compensated "
            "paleobathymetry (Step 4). A port of EarthByte's simple_paleobathymetry workflow; "
            "see gplately.grids.paleobathymetry for the Python API and its docstring for what "
            "is not yet included (continent-obstacle routing, pyBacktrack merge).\n\n"
            "Example usage:\n"
            "    gplately pb output_dir -m muller2025 --proximity-features cobs.gpml -e 0 -s 10\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_arguments(cmd)
    _add_distance_arguments(cmd)
    cmd.add_argument(
        "--age-depth-model",
        metavar="age_depth_model",
        choices=AGE_DEPTH_MODELS,
        default="gdh1",
        dest="age_depth_model",
        help="thermal-subsidence (age -> basement depth) model; default: gdh1",
    )
    cmd.add_argument(
        "--richards-table",
        metavar="richards_table",
        default=None,
        dest="richards_table",
        help="age-depth lookup table for --age-depth-model rhcw18; "
        "default: the table shipped with gplately",
    )
    cmd.set_defaults(func=_run_paleobathymetry)
