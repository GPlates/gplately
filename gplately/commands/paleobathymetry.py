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
import tempfile

import pygplates

from ..grids.paleobathymetry import AGE_DEPTH_MODELS, simple_paleobathymetry
from .sediment_thickness import (
    _add_common_arguments,
    _add_distance_arguments,
    _resolve_age_grid_filenames_and_times,
    _resolve_distance_grid_kwargs,
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

    kwargs = _resolve_distance_grid_kwargs(args, plate_model)
    if args.pybacktrack:
        static_polygon_filename = args.static_polygons or (
            plate_model.get_layer("StaticPolygons") if plate_model else None
        )
        present_day_age_grid_filename = args.present_day_age_grid or (
            plate_model.get_age_grid(0) if plate_model else None
        )
        if not static_polygon_filename:
            raise Exception(
                "--pybacktrack requires --static-polygons (unless -m/--model's plate "
                "model provides a StaticPolygons layer)."
            )
        if not present_day_age_grid_filename:
            raise Exception(
                "--pybacktrack requires --present-day-age-grid (unless -m/--model is used, "
                "which fetches the age grid at 0 Ma automatically)."
            )
        if isinstance(static_polygon_filename, (list, tuple)):
            if len(static_polygon_filename) > 1:
                # pyBacktrack's static_polygon_filename takes exactly one file; merge
                # rather than silently dropping every file but the first.
                merged = pygplates.FeatureCollection()
                for filename in static_polygon_filename:
                    merged.add(pygplates.FeatureCollection(filename))
                merged_file = tempfile.NamedTemporaryFile(suffix=".gpmlz", delete=False)
                merged_file.close()
                merged.write(merged_file.name)
                _logger.info(
                    f"Merged {len(static_polygon_filename)} StaticPolygons files into "
                    f"{merged_file.name} for --pybacktrack"
                )
                static_polygon_filename = merged_file.name
            else:
                static_polygon_filename = static_polygon_filename[0]
        kwargs.update(
            pybacktrack=True,
            static_polygon_filename=static_polygon_filename,
            present_day_age_grid_filename=present_day_age_grid_filename,
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
        decimal_places_in_time=args.decimal_places_in_time,
        **kwargs,
    )
    _logger.info(f"Paleobathymetry grids written to {args.output_dir}")


def add_parser(parser):
    """add command line argument parser for 'paleobathymetry'"""

    cmd = parser.add_parser(
        "paleobathymetry",
        aliases=("pb",),
        help="Run the simple_paleobathymetry workflow (Steps 1-4, optionally 5) end to end.",
        add_help=True,
        description=(
            "Reconstruct paleobathymetry of ocean crust from a seafloor-age grid: age -> "
            "basement depth (Step 1), distance to the nearest passive continental margin "
            "(Step 2), predicted sediment thickness (Step 3), and isostatically-compensated "
            "paleobathymetry (Step 4). Optionally (--pybacktrack) also merge in pyBacktrack's "
            "present-day paleobathymetry (Step 5) to also cover submerged continental crust "
            "and crust that has since subducted. A port of EarthByte's simple_paleobathymetry "
            "workflow; see gplately.grids.paleobathymetry for the Python API and its docstring "
            "for what is not yet included (continent-obstacle routing in Step 2).\n\n"
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
    cmd.add_argument(
        "--pybacktrack",
        action="store_true",
        dest="pybacktrack",
        help="also run Step 5: merge in pyBacktrack's present-day paleobathymetry, to also "
        "cover submerged continental crust and crust that has since subducted. Requires the "
        "optional 'pybacktrack' package (pip install pybacktrack, or gplately[paleobathymetry]).",
    )
    cmd.add_argument(
        "--static-polygons",
        metavar="static_polygon_filename",
        default=None,
        dest="static_polygons",
        help="static polygons, for --pybacktrack (pyBacktrack uses these to assign plate IDs); "
        "required unless -m/--model's plate model provides a StaticPolygons layer",
    )
    cmd.add_argument(
        "--present-day-age-grid",
        metavar="present_day_age_grid_filename",
        default=None,
        dest="present_day_age_grid",
        help="the seafloor-age grid at 0 Ma, for --pybacktrack (regardless of -e/-s, since "
        "pyBacktrack backtracks from the present day); default: fetched automatically "
        "when -m/--model is used",
    )
    cmd.set_defaults(func=_run_paleobathymetry)
