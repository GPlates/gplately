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

from plate_model_manager import PlateModelManager

from ..grids.continent_contouring import generate_passive_margins

_logger = logging.getLogger("gplately")


def _run_generate_passive_margins(args):
    times = list(range(int(args.min_time), int(args.max_time) + 1, int(args.time_step)))

    plate_model = None
    if args.model_name:
        plate_model = PlateModelManager().get_model(
            args.model_name, data_dir=args.plate_model_repo
        )
        if not plate_model:
            raise Exception(
                f"Unable to create PlateModel object for model {args.model_name}."
            )

    rotation_files = args.rotation_filenames or (
        plate_model.get_rotation_model() if plate_model else None
    )
    topology_files = args.topology_filenames or (
        plate_model.get_layer("Topologies") if plate_model else None
    )
    continent_files = args.continent_filenames or (
        plate_model.get_layer("ContinentalPolygons") if plate_model else None
    )
    if not rotation_files or not topology_files:
        raise Exception(
            "No rotation/topology files found: use -m/--model, or --rotations/--topologies."
        )
    if not continent_files:
        raise Exception(
            "No continental polygon files found: use -m/--model, or --continents."
        )

    _logger.info(f"Using rotation files: {rotation_files}")
    _logger.info(f"Using topology files: {topology_files}")
    _logger.info(f"Using continent files: {continent_files}")

    generate_passive_margins(
        rotation_model=rotation_files,
        continent_features=continent_files,
        topological_features=topology_files,
        times=times,
        point_spacing_degrees=args.point_spacing,
        area_threshold_square_kms=args.area_threshold_km2,
        buffer_and_gap_distance_kms=args.buffer_and_gap_km,
        exclusion_area_threshold_square_kms=args.exclusion_area_threshold_km2,
        max_distance_of_subduction_from_active_margin_kms=args.max_distance_active_margin_km,
        anchor_plate_id=args.anchor_plate_id or 0,
        time_step=args.time_step,
        output_directory=args.output_dir,
    )
    _logger.info(f"Passive margins written to {args.output_dir}")


def add_parser(parser):
    """add command line argument parser for 'generate-passive-margins'"""

    cmd = parser.add_parser(
        "generate-passive-margins",
        aliases=("gpm",),
        help="Dynamically contour continents through time and split each contour into "
        "passive margins.",
        add_help=True,
        description=(
            "For each time, contour continental polygons into continents (gplately's "
            "ContinentContouring engine), then split each contour into passive-margin "
            "segments by removing the parts close to a subduction zone. A port of "
            "EarthByte's continent-contouring create_passive_margins.py.\n\n"
            "Example usage:\n"
            "    gplately gpm output_dir -m muller2025 -e 0 -s 100\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    cmd.set_defaults(func=_run_generate_passive_margins)
    cmd.add_argument(
        metavar="output_dir",
        help="(required) output directory",
        dest="output_dir",
    )
    cmd.add_argument(
        "-m",
        "--model",
        metavar="model_name",
        dest="model_name",
        default=None,
        help="reconstruction model name (fetched via the Plate Model Manager); "
        "supplies rotations/topologies/continental polygons unless overridden below",
    )
    cmd.add_argument(
        "-f",
        "--plate-model-repo",
        metavar="plate_model_repo",
        dest="plate_model_repo",
        default="plate-model-repo",
        help="local cache directory for -m/--model",
    )
    cmd.add_argument(
        "--rotations",
        metavar="rotation_filenames",
        nargs="+",
        dest="rotation_filenames",
        default=[],
        help="alternative to -m/--model",
    )
    cmd.add_argument(
        "--topologies",
        metavar="topology_filenames",
        nargs="+",
        dest="topology_filenames",
        default=[],
        help="alternative to -m/--model",
    )
    cmd.add_argument(
        "--continents",
        metavar="continent_filenames",
        nargs="+",
        dest="continent_filenames",
        default=[],
        help="continental polygon (or craton) file(s) to contour; alternative to -m/--model",
    )
    cmd.add_argument(
        "-e",
        "--min-time",
        metavar="min_time",
        type=float,
        default=0,
        dest="min_time",
        help="minimum time (Ma); default: 0",
    )
    cmd.add_argument(
        "-s",
        "--max-time",
        metavar="max_time",
        type=float,
        default=100,
        dest="max_time",
        help="maximum time (Ma); default: 100",
    )
    cmd.add_argument(
        "--time-step",
        metavar="time_step",
        type=float,
        default=1,
        dest="time_step",
        help="time increment (Myr); default: 1",
    )
    cmd.add_argument(
        "-r",
        "--point-spacing",
        metavar="point_spacing_degrees",
        type=float,
        default=0.25,
        dest="point_spacing",
        help="grid spacing (degrees) used to contour/aggregate continental polygons; "
        "default: 0.25",
    )
    cmd.add_argument(
        "--area-threshold-km2",
        metavar="area_threshold_km2",
        type=float,
        default=0.0,
        dest="area_threshold_km2",
        help="exclude contoured continents smaller than this (km^2); default: 0",
    )
    cmd.add_argument(
        "--buffer-and-gap-km",
        metavar="buffer_and_gap_km",
        type=float,
        default=0.0,
        dest="buffer_and_gap_km",
        help="expand continents ocean-ward by this distance (km) before contouring; "
        "default: 0",
    )
    cmd.add_argument(
        "--exclusion-area-threshold-km2",
        metavar="exclusion_area_threshold_km2",
        type=float,
        default=800000.0,
        dest="exclusion_area_threshold_km2",
        help="drop enclosed interior gaps (e.g. lakes) smaller than this (km^2); "
        "default: 800000",
    )
    cmd.add_argument(
        "--max-distance-active-margin-km",
        metavar="max_distance_active_margin_km",
        type=float,
        default=500.0,
        dest="max_distance_active_margin_km",
        help="a contour segment within this distance (km) of a subduction zone is an "
        "active margin; default: 500",
    )
    cmd.add_argument(
        "-a",
        "--anchor-plate-id",
        metavar="anchor_plate_id",
        type=int,
        default=None,
        dest="anchor_plate_id",
        help="anchor plate ID; default: 0",
    )
