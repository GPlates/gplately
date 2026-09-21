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
import os

from plate_model_manager import PlateModelManager

from ..grids import read_netcdf_grid
from ..grids._utils import (
    DEFAULT_DECIMAL_PLACES_IN_TIME,
    DEFAULT_DISTANCE_GRID_DECIMAL_PLACES_IN_TIME,
    check_times_are_distinct_in_filenames,
    distance_grid_filename,
    resolve_decimal_places_in_time,
)
from ..grids.sediment_thickness import (
    generate_distance_grids,
    generate_sediment_thickness_grids,
)

_logger = logging.getLogger("gplately")


def _time_range(min_time, max_time, time_step):
    """min_time, max_time, time_step -> a list of times, min_time to max_time inclusive.

    Unlike ``range(int(min_time), int(max_time) + 1, int(time_step))``, this works for
    fractional (sub-Myr) values -- --min-time/--max-time/--time-step are all declared as
    floats on the CLI, so truncating them to int here would silently corrupt (or, for a
    time_step below 1, crash outright: ``range()`` rejects a step of 0) a perfectly valid
    request like ``--time-step 0.5``.
    """
    if time_step <= 0:
        raise ValueError("time_step must be positive.")
    num_steps = int(round((max_time - min_time) / time_step)) + 1
    return [min_time + i * time_step for i in range(num_steps)]


def _resolve_age_grid_filenames_and_times(args):
    times = _time_range(args.min_time, args.max_time, args.time_step)
    if args.model_name:
        plate_model = PlateModelManager().get_model(
            args.model_name, data_dir=args.plate_model_repo
        )
        if not plate_model:
            raise Exception(
                f"Unable to create PlateModel object for model {args.model_name}."
            )
        age_grid_filenames_and_times = [
            (plate_model.get_age_grid(t), float(t)) for t in times
        ]
    elif args.age_grid_template:
        age_grid_filenames_and_times = [
            (args.age_grid_template.format(time=t), float(t)) for t in times
        ]
    else:
        raise Exception(
            "No age grid source given: use -m/--model or --age-grid-template."
        )
    return age_grid_filenames_and_times, plate_model if args.model_name else None


def _resolve_rotation_topology_proximity_files(args, plate_model):
    rotation_files = args.rotation_filenames or (
        plate_model.get_rotation_model() if plate_model else None
    )
    topology_files = args.topology_filenames or (
        plate_model.get_layer("Topologies", return_none_if_not_exist=True)
        if plate_model
        else None
    )
    proximity_files = args.proximity_filenames or (
        plate_model.get_layer("COBs", return_none_if_not_exist=True)
        if plate_model
        else None
    )
    if not rotation_files or not topology_files:
        raise Exception(
            "No rotation/topology files found: use -m/--model, or --rotations/--topologies."
        )
    if not proximity_files:
        raise Exception(
            "No proximity feature files found: use --proximity-features to supply "
            "passive-margin continent-ocean-boundary line segments (the Plate Model "
            "Manager does not deliver these for most models)."
        )

    _logger.info(f"Using rotation files: {rotation_files}")
    _logger.info(f"Using topology files: {topology_files}")
    _logger.info(f"Using proximity feature files: {proximity_files}")
    return rotation_files, topology_files, proximity_files


def _resolve_distance_grid_kwargs(args, plate_model):
    """Resolve the continent-obstacle-routing kwargs shared by 'generate-distance-grids' and
    'paleobathymetry'. Returns a dict suitable for **-splatting into generate_distance_grids()
    (or gplately.grids.paleobathymetry.simple_paleobathymetry()) -- empty if obstacle routing
    was not requested.
    """
    if not (args.route_around_continents or args.continent_obstacle_filenames):
        return {}

    continent_obstacle_files = args.continent_obstacle_filenames or (
        plate_model.get_layer("Coastlines", return_none_if_not_exist=True)
        if plate_model
        else None
    )
    if not continent_obstacle_files:
        raise Exception(
            "--route-around-continents (or --continent-obstacles) requires continent/coastline "
            "files: use --continent-obstacles, or -m/--model's plate model must provide a "
            "Coastlines layer."
        )
    _logger.info(f"Using continent obstacle files: {continent_obstacle_files}")

    return dict(
        continent_obstacle_features=continent_obstacle_files,
        shortest_path_grid_subdivision_depth=args.shortest_path_grid_depth,
    )


def _run_generate_distance_grids(args):
    age_grid_filenames_and_times, plate_model = _resolve_age_grid_filenames_and_times(
        args
    )
    rotation_files, topology_files, proximity_files = (
        _resolve_rotation_topology_proximity_files(args, plate_model)
    )
    distance_grid_kwargs = _resolve_distance_grid_kwargs(args, plate_model)

    generate_distance_grids(
        rotation_model=rotation_files,
        proximity_features=proximity_files,
        topological_features=topology_files,
        age_grid_filenames_and_times=age_grid_filenames_and_times,
        grid_spacing=args.grid_spacing,
        time_increment=args.time_increment,
        max_reconstruction_time=args.max_reconstruction_time,
        anchor_plate_id=args.anchor_plate_id or 0,
        clamp_distance_km=args.clamp_distance_km,
        output_directory=args.output_dir,
        decimal_places_in_time=args.decimal_places_in_time,
        **distance_grid_kwargs,
    )
    _logger.info(f"Distance grids written to {args.output_dir}")


def _run_generate_sediment_grids(args):
    age_grid_filenames_and_times, _ = _resolve_age_grid_filenames_and_times(args)

    # Check the output names before reading anything: the library raises on a collision,
    # but only after this function has already read every distance grid off disk.
    check_times_are_distinct_in_filenames(
        [time for _, time in age_grid_filenames_and_times],
        resolve_decimal_places_in_time(
            args.decimal_places_in_time, DEFAULT_DECIMAL_PLACES_IN_TIME
        ),
        "sediment_thickness_{}Ma.nc",
    )

    # Must match the names generate_distance_grids() wrote, which default to one decimal
    # place of time rather than the zero used by the paleobathymetry grids.
    distance_grid_decimal_places = resolve_decimal_places_in_time(
        args.decimal_places_in_time, DEFAULT_DISTANCE_GRID_DECIMAL_PLACES_IN_TIME
    )
    distance_grids = {}
    for _, t in age_grid_filenames_and_times:
        distance_path = os.path.join(
            args.distance_grids_dir,
            distance_grid_filename(args.grid_spacing, t, distance_grid_decimal_places),
        )
        grid, lon, lat = read_netcdf_grid(distance_path, return_grids=True)
        distance_grids[t] = (lon, lat, grid)

    generate_sediment_thickness_grids(
        age_grid_filenames_and_times,
        distance_grids,
        output_directory=args.output_dir,
        decimal_places_in_time=args.decimal_places_in_time,
    )
    _logger.info(f"Sediment thickness grids written to {args.output_dir}")


def _add_common_arguments(cmd):
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
        "supplies rotations/topologies/age-grids unless overridden below",
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
        "--age-grid-template",
        metavar="age_grid_template",
        dest="age_grid_template",
        default=None,
        help="alternative to -m/--model: a filename template for local age grids, "
        "using '{time}' for the reconstruction age in Ma, e.g. 'agegrids/age_{time:.0f}Ma.nc'",
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
        help="spacing (Myr) of the times to generate output for, between --min-time and "
        "--max-time; default: 1. This selects which times get output; how finely each "
        "point's lifetime is sampled is --time-increment, which every output time must be "
        "a multiple of (so a fractional step below 1 Myr needs --time-increment set to "
        "match)",
    )
    cmd.add_argument(
        "--decimal-places-in-time",
        metavar="decimal_places_in_time",
        type=int,
        default=None,
        dest="decimal_places_in_time",
        help="decimal places of the reconstruction time in output filenames; by default 1 "
        "for the mean-distance grids and 0 for all the others, which reproduces the "
        "original workflows' names. Raise it when using a fractional time step, otherwise "
        "consecutive times share a filename and only the last is kept",
    )
    cmd.add_argument(
        "-r",
        "--grid-spacing",
        metavar="grid_spacing",
        type=float,
        default=0.5,
        dest="grid_spacing",
        help="grid spacing (degrees); default: 0.5",
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


def _add_distance_arguments(cmd):
    cmd.add_argument(
        "--time-increment",
        metavar="time_increment",
        type=float,
        default=1,
        dest="time_increment",
        help="increment (Myr) used to step the backward reconstruction and sample distance "
        "along each ocean point's lifetime; default: 1, as in the original workflow. This is "
        "independent of --time-step: coarser output does not mean coarser sampling",
    )
    cmd.add_argument(
        "--proximity-features",
        metavar="proximity_filenames",
        nargs="+",
        dest="proximity_filenames",
        default=[],
        help="passive-margin continent-ocean-boundary line-segment file(s); "
        "required unless -m/--model's plate model provides a COBs layer",
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
        "--max-reconstruction-time",
        metavar="max_reconstruction_time",
        type=float,
        default=None,
        dest="max_reconstruction_time",
        help="do not reconstruct ocean points older than this (Ma); default: unlimited",
    )
    cmd.add_argument(
        "--clamp-distance-km",
        metavar="clamp_distance_km",
        type=float,
        default=3000.0,
        dest="clamp_distance_km",
        help="clamp mean distances above this (km); default: 3000",
    )
    cmd.add_argument(
        "--route-around-continents",
        action="store_true",
        dest="route_around_continents",
        help="route distances around continents instead of a great-circle straight line "
        "(auto-resolves --continent-obstacles from -m/--model's Coastlines layer if not "
        "given explicitly)",
    )
    cmd.add_argument(
        "--continent-obstacles",
        metavar="continent_obstacle_filenames",
        nargs="+",
        dest="continent_obstacle_filenames",
        default=[],
        help="continent/coastline file(s) to route around; implies --route-around-continents",
    )
    cmd.add_argument(
        "--shortest-path-grid-depth",
        metavar="shortest_path_grid_subdivision_depth",
        type=int,
        default=6,
        dest="shortest_path_grid_depth",
        help="subdivision depth of the grid used for continent-obstacle routing "
        "(spacing = 90/2^depth degrees); default: 6",
    )


def add_parser(parser):
    """add command line argument parsers for 'generate-distance-grids' and 'generate-sediment-grids'"""

    distance_cmd = parser.add_parser(
        "generate-distance-grids",
        aliases=("gdg",),
        help="Generate distance-to-passive-margin grids (for predicting sediment thickness).",
        add_help=True,
        description=(
            "For each ocean point in a seafloor-age grid, reconstruct it backward through "
            "time and compute its lifetime-mean distance to the nearest passive-margin "
            "continent-ocean-boundary line segment.\n\n"
            "Example usage:\n"
            "    gplately gdg output_dir -m muller2025 --proximity-features cobs.gpml -e 0 -s 10\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_arguments(distance_cmd)
    _add_distance_arguments(distance_cmd)
    distance_cmd.set_defaults(func=_run_generate_distance_grids)

    sediment_cmd = parser.add_parser(
        "generate-sediment-grids",
        aliases=("gsg",),
        help="Predict sediment-thickness grids from seafloor age and distance-to-margin grids.",
        add_help=True,
        description=(
            "Combine a seafloor-age grid with the distance-to-passive-margin grids from "
            "'generate-distance-grids' into predicted compacted sediment-thickness grids "
            "(Dutkiewicz et al., 2017).\n\n"
            "Example usage:\n"
            "    gplately gsg output_dir -m muller2025 --distance-grids-dir distances/ -e 0 -s 10\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_arguments(sediment_cmd)
    sediment_cmd.add_argument(
        "--distance-grids-dir",
        metavar="distance_grids_dir",
        required=True,
        dest="distance_grids_dir",
        help="directory of mean_distance_<spacing>d_<time>.nc grids from 'generate-distance-grids'",
    )
    sediment_cmd.set_defaults(func=_run_generate_sediment_grids)
