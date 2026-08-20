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
import datetime
import logging
import multiprocessing
import time
import warnings
from pathlib import Path
from typing import List, Union

from plate_model_manager import PlateModel, PlateModelManager, plate_model
from pygplates import FeaturesFunctionArgument

from ..reconstruction import PlateReconstruction
from ..grids import (
    SeafloorGrid,
    IsochronSeafloorGrid,
    TopologySeafloorGrid,
    OutputScalarType,
)
from ..auxiliary import get_plate_reconstruction

_logger = logging.getLogger("gplately")

# Single source of truth for default values of every option that can be set
# via the command line and/or a --config file. output_dir is deliberately
# excluded: it's a required positional argument (see add_parser()) because
# argparse cannot unambiguously split two consecutive variable-length
# positional arguments (INPUT_FILES and OUTPUT_DIR) if OUTPUT_DIR were made
# optional too - so OUTPUT_DIR must always be given on the command line, even
# when a --config file is used for everything else.
_DEFAULTS = {
    "model_name": None,
    "input_filenames": [],
    "rotation_filenames": [],
    "topology_filenames": [],
    "continents_filenames": [],
    "ridges_filenames": [],
    "isochron_filenames": [],
    "isochron_cob_filenames": [],
    "grid_spacing": 0.1,
    "interval_spacing": 0.1,
    "refinement_levels": 6,
    "ridge_sampling": 0.5,
    "initial_spreadrate": 75,
    "times": [],
    "min_time": 0,
    "max_time": 140,
    "time_interval": 1,
    "n_jobs": None,
    "file_collection": None,
    "unmasked": False,
    "resume_from_checkpoints": False,
    "use_continent_contouring": False,
    "use_isochron_interp": False,
    "anchor_plate_id": None,
}


def _default_values() -> dict:
    """Return a fresh copy of _DEFAULTS (so callers never share mutable list defaults)."""
    return {
        key: (list(value) if isinstance(value, list) else value)
        for key, value in _DEFAULTS.items()
    }


def _load_config_file(path: str) -> dict:
    """Load argument values from a TOML config file.

    The file may either be flat (keys matching the CLI's dest names directly
    at the top level) or namespace the agegrid settings under an [agegrid]
    table, so one config file can be shared across multiple gplately
    subcommands, e.g.:

        [seafloorgrid]
        model_name = "merdith2021"
        min_time = 0
        max_time = 10
        grid_spacing = 0.5
    """
    try:
        import tomllib  # Python 3.11+
    except ImportError:
        try:
            import tomli as tomllib  # backport for Python < 3.11
        except ImportError as exc:
            raise Exception(
                "Reading a --config file requires the 'tomli' package on "
                "Python < 3.11. Install it with 'pip install tomli'."
            ) from exc

    config_path = Path(path)
    if not config_path.is_file():
        raise Exception(f"Config file not found: {path}")

    with open(config_path, "rb") as f:
        try:
            data = tomllib.load(f)
        except Exception as exc:
            raise Exception(f"Failed to parse config file {path}: {exc}") from exc

    if "seafloorgrid" in data and isinstance(data["seafloorgrid"], dict):
        return data["seafloorgrid"]
    return data


def _merge_config_and_args(args: argparse.Namespace) -> dict:
    """Merge built-in defaults, an optional --config file, and CLI arguments.

    Precedence (highest wins): command line flags > --config file > built-in
    defaults. Arguments not explicitly given on the command line are absent
    from `vars(args)`, because every optional argparse argument below is
    defined with default=argparse.SUPPRESS - that's what lets us tell "user
    didn't pass this flag" apart from "user passed a falsy value".
    """
    merged = _default_values()

    config_path = getattr(args, "config_file", None)
    if config_path:
        config_values = _load_config_file(config_path)
        unknown_keys = set(config_values) - set(_DEFAULTS)
        if unknown_keys:
            _logger.warning(
                "Ignoring unrecognised key(s) in config file %s: %s",
                config_path,
                ", ".join(sorted(unknown_keys)),
            )
        for key in _DEFAULTS:
            if key in config_values:
                merged[key] = config_values[key]

    for key, value in vars(args).items():
        if key in ("func", "config_file", "output_dir"):
            continue
        merged[key] = value

    return merged


def add_parser(parser):
    """add command line argument parser"""

    agegrid_cmd = parser.add_parser(
        "agegrid",
        aliases=("ag",),
        help=_help_str,
        add_help=True,
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # agegrid command arguments.
    #
    # Every optional argument below uses default=argparse.SUPPRESS instead of
    # a literal default value. This means a flag the user didn't type on the
    # command line simply won't appear in the parsed Namespace at all, which
    # is what lets _merge_config_and_args() tell "user didn't pass this" apart
    # from "user explicitly passed a falsy value" and correctly layer:
    #   built-in defaults (_DEFAULTS)  <  --config file  <  CLI flags
    # INPUT_FILES and OUTPUT_DIR are positional. OUTPUT_DIR stays required
    # (no default) rather than config-driven: argparse cannot unambiguously
    # split two consecutive variable-length positionals, so making OUTPUT_DIR
    # optional too would risk it silently swallowing INPUT_FILES tokens.
    agegrid_cmd.set_defaults(func=_run_seafloor_gridding)
    agegrid_cmd.add_argument(
        metavar="input_filenames",
        nargs="*",
        help="(optional) input rotation and topology files, or set 'input_filenames' "
        "in --config; prefer --rotations/--topologies instead, which take precedence "
        "over this positional argument if provided",
        dest="input_filenames",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        metavar="output_dir",
        help="(required) output directory (always required on the command line, even with --config)",
        dest="output_dir",
    )
    agegrid_cmd.add_argument(
        "-t",
        "--config",
        metavar="config_file",
        help="path to a TOML config file supplying any of the options below; "
        "command line flags override values from this file",
        default=None,
        dest="config_file",
    )
    agegrid_cmd.add_argument(
        "-m",
        "--model",
        metavar="model_name",
        help="reconstruction model name, this may suppress the INPUT_FILES, ROTATION_FILE, TOPOLOGY_FILE, and CONTINENTS_FILE arguments",
        dest="model_name",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "-c",
        "--continents",
        metavar="continents_filenames",
        nargs="*",
        help="input continental polygons files",
        dest="continents_filenames",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--rotations",
        metavar="rotation_filenames",
        nargs="*",
        help="alternative way to specify input rotation files",
        dest="rotation_filenames",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--topologies",
        metavar="topology_filenames",
        nargs="*",
        help="alternative way to specify input topology files",
        dest="topology_filenames",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--ridges",
        metavar="ridges_filenames",
        nargs="*",
        help="input ridges files",
        dest="ridges_filenames",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--isochrons",
        metavar="isochron_filenames",
        nargs="*",
        help="input isochron files",
        dest="isochron_filenames",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--isochrons-cob",
        metavar="isochron_cob_filenames",
        nargs="*",
        help="input isochron cob files",
        dest="isochron_cob_filenames",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "-r",
        "--resolution",
        metavar="grid_spacing",
        type=float,
        help=f"grid resolution (degrees); default: {_DEFAULTS['grid_spacing']}",
        dest="grid_spacing",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--interval-spacing",
        metavar="interval_spacing",
        type=float,
        help=f"The spacing between lines (in degrees) at which to generate interpolated isochrons.",
        dest="interval_spacing",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--refinement-levels",
        metavar="refinement_levels",
        type=int,
        help=f"mesh refinement levels; default: {_DEFAULTS['refinement_levels']}",
        dest="refinement_levels",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--ridge-sampling",
        metavar="ridge_sampling",
        type=float,
        help=f"MOR sampling resolution (degrees); default: {_DEFAULTS['ridge_sampling']}",
        dest="ridge_sampling",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--initial-spreadrate",
        metavar="initial_spreadrate",
        type=float,
        help=f"initial ocean spreading rate (km/Myr); default: {_DEFAULTS['initial_spreadrate']}",
        dest="initial_spreadrate",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--time-interval",
        metavar="time_interval",
        type=float,
        help=f"The time interval (Ma) to use for generating the age grid; default: {_DEFAULTS['time_interval']}",
        dest="time_interval",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--times",
        metavar="times",
        nargs="+",
        type=float,
        help=f"A list of times (Ma) to generate the age grid.",
        dest="times",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "-e",
        "--min-time",
        metavar="min_time",
        type=float,
        help=f"minimum time (Ma); default: {_DEFAULTS['min_time']}",
        dest="min_time",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "-s",
        "--max-time",
        metavar="max_time",
        type=float,
        help=f"maximum time (Ma); default: {_DEFAULTS['max_time']}",
        dest="max_time",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "-j",
        "--n-jobs",
        type=int,  # NOTE: previously missing, so a CLI-provided value stayed a
        # str and would misbehave when handed to multiprocessing/SeafloorGrid.
        help="number of processes to use; default: use all CPU available",
        metavar="n_jobs",
        dest="n_jobs",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "-f",
        "--file-collection",
        help="file collection name (optional)",
        metavar="file_collection",
        dest="file_collection",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "-u",
        "--include-unmasked",
        help="create unmasked grids in addition to masked ones",
        action="store_true",
        dest="include_unmasked",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--resume-from-checkpoints",
        help="flag to indicate whether to resume a previously interrupted gridding run "
        "(all parameters and input data should remain unchanged); "
        "default: re-run gridding from scratch",
        action="store_true",
        dest="resume_from_checkpoints",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--use-continent-contouring",
        help="flag to indicate using ptt's 'continent contouring' to generate continent masks",
        action="store_true",
        dest="use_continent_contouring",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "--use-isochron-interp",
        help="flag to indicate using 'isochron interpolation' to generate age grids",
        action="store_true",
        dest="use_isochron_interp",
        default=argparse.SUPPRESS,
    )
    agegrid_cmd.add_argument(
        "-a",
        "--anchor-plate-id",
        metavar="anchor_plate_id",
        type=int,
        help="anchor plate ID; default: 0",
        dest="anchor_plate_id",
        default=argparse.SUPPRESS,  # if unset, the default anchor plate ID will effectively be 0
    )


_help_str = "Generate seafloor grids, such as age grids and spreading rate grids."

__description__ = f"""{_help_str}

Example usage:
    - gplately ag output_dir -m merdith2021 -e 0 -s 10
    - gplately ag rotations.rot topologies.gpmlz output_dir -c continental_polygons.gpmlz -e 0 -s 10
    - gplately ag output_dir --config myconfig.toml
    - gplately ag output_dir --config myconfig.toml -e 0 -s 10   # CLI flags override the config file

Config file (--config) values are overridden by any matching command line
flag; anything not set in either place falls back to the built-in default.
OUTPUT_DIR is always required on the command line. Example myconfig.toml:

    [agegrid]
    model_name = "merdith2021"
    min_time = 0
    max_time = 10
    grid_spacing = 0.5
    unmasked = true
"""


def _run_topology_seafloor_gridding(
    model_name: str,
    input_filenames: List[str],
    continents_filenames: List[str],
    output_dir: str,
    min_time: float,
    max_time: float,
    ridge_time_step: float = 1,
    n_jobs: int = 1,
    refinement_levels: int = 6,
    grid_spacing: float = 0.1,
    ridge_sampling: float = 0.5,
    initial_spreadrate: float = 75,
    file_collection: str = "",
    include_unmasked: bool = False,
    plate_model_repo: str = "plate-model-repo",
    resume_from_checkpoints: bool = False,
    use_continent_contouring: bool = False,
    anchor_plate_id: Union[int, None] = None,
) -> None:
    """Create age grids for a plate model."""

    if model_name:
        if input_filenames or continents_filenames:
            raise Exception(
                f"The model name({model_name}) is provided with -m or --model. \n"
                + f"But the INPUT_FILES and/or CONTINENTS_FILE are also provided. \n"
                + f"The MODEL_NAME and (INPUT_FILES, CONTINENTS_FILE) are mutually exclusive. \n"
                + f"You either use (-m,--model) to specify a model name or use the first positional arguments and (-c,--continents) to specify the file paths. But not both.\n"
                + f"Fix your command line and try again."
            )
        try:
            # TODO: make "plate-model-repo" an optional parameter
            plate_model = PlateModelManager().get_model(
                model_name, data_dir=plate_model_repo
            )
        except:
            # try to use the plate model in the local folder if the network is down.
            plate_model = PlateModel(
                model_name, data_dir=plate_model_repo, readonly=True
            )

        if not plate_model:
            raise Exception(
                f"Unable to create PlateModel object for model {model_name}. "
                + f"Check your network connection and make sure the model name is correct. Run `pmm ls` to get a list of available model names."
            )

        rotation_files = plate_model.get_rotation_model()
        topology_files = plate_model.get_layer(
            "Topologies", return_none_if_not_exist=True
        )
        continent_files = plate_model.get_layer(
            "ContinentalPolygons", return_none_if_not_exist=True
        )
        if "Cratons" in plate_model.get_avail_layers():
            continent_files += plate_model.get_layer("Cratons")  # type: ignore

    else:
        rotation_files, topology_files = (
            _get_rotation_and_topology_files_from_input_filenames(input_filenames)
        )
        continent_files = continents_filenames

    if not rotation_files:
        raise Exception(
            "No rotation file(s) found. User must either provide rotation file(s) in the first positional arguments or specify a model name with (-m,--model)."
        )
    else:
        _logger.info(f"Using rotation files: {rotation_files}")

    if not topology_files:
        error_msg = "No topology file(s) found. "
        if model_name:
            error_msg += f"Make sure the model ({model_name}) has topology. Run `pmm ls {model_name}` to check the model."
        else:
            error_msg += f"You need to specify the topology files in the first positional arguments."
        raise Exception(error_msg)
    else:
        _logger.info(f"Using topology files: {topology_files}")

    if not continent_files:
        error_msg = "No continental polygon file(s) found. "
        if model_name:
            error_msg += f"Make sure the model ({model_name}) has continental polygons. Run `pmm ls {model_name}` to check the model."
        else:
            error_msg += f"You need to specify the continental polygon files with (-c,--continents). Run `gplately ag -h` to see help."
        raise Exception(error_msg)
    else:
        _logger.info(f"Using continent files: {continent_files}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ImportWarning)
        reconstruction = PlateReconstruction(
            rotation_model=rotation_files,
            topology_features=topology_files,
            anchor_plate_id=anchor_plate_id,
        )

        grid = TopologySeafloorGrid(
            reconstruction,
            min_time=min_time,
            max_time=max_time,
            save_directory=output_dir,
            ridge_time_step=ridge_time_step,
            refinement_levels=refinement_levels,
            grid_spacing=grid_spacing,
            ridge_sampling=ridge_sampling,
            initial_ocean_mean_spreading_rate=initial_spreadrate,
            file_collection=file_collection,
            resume_from_checkpoints=resume_from_checkpoints,
            continent_polygon_features=continent_files,
            use_continent_contouring=use_continent_contouring,
            nprocs=n_jobs,
        )

        grid.generate()
        for val in (
            SeafloorGrid.SEAFLOOR_AGE_KEY,
            SeafloorGrid.SPREADING_RATE_KEY,
        ):
            grid.lat_lon_z_to_netCDF(val, unmasked=include_unmasked)


def _run_seafloor_gridding(args):
    values = _merge_config_and_args(args)
    n_jobs = values["n_jobs"]
    if not n_jobs:
        try:
            n_jobs = multiprocessing.cpu_count()
        except NotImplementedError:
            n_jobs = 1
    start = time.time()

    file_collection = values["model_name"] or values["file_collection"]

    if values["use_isochron_interp"]:
        _logger.info("Using isochron interpolation to generate age grids.")
        _run_isochron_seafloor_gridding(args, values)
    else:
        _logger.info("Using topological reconstruction to generate age grids.")
        _run_topology_seafloor_gridding(
            model_name=values["model_name"],
            input_filenames=values["input_filenames"],
            continents_filenames=values["continents_filenames"],
            output_dir=args.output_dir,
            min_time=values["min_time"],
            max_time=values["max_time"],
            ridge_time_step=values["time_interval"],
            n_jobs=n_jobs,
            refinement_levels=values["refinement_levels"],
            grid_spacing=values["grid_spacing"],
            ridge_sampling=values["ridge_sampling"],
            initial_spreadrate=values["initial_spreadrate"],
            file_collection=file_collection,
            include_unmasked=values["include_unmasked"],
            resume_from_checkpoints=values["resume_from_checkpoints"],
            use_continent_contouring=values["use_continent_contouring"],
            anchor_plate_id=values["anchor_plate_id"],
        )

    end = time.time()
    hours_minutes_seconds = str(datetime.timedelta(seconds=end - start)).split(":")
    _logger.info(
        f"Completed creating age grids in {hours_minutes_seconds[0]} Hours, {hours_minutes_seconds[1]} Minutes, {hours_minutes_seconds[2].split('.')[0]} Seconds "
    )


def _run_isochron_seafloor_gridding(args, values):
    rotation_files, topology_files, continent_files = (
        _get_rotation_topology_and_continent_files(values)
    )
    if not topology_files:
        topology_files = None

    if not rotation_files:
        raise Exception(
            "No rotation file(s) found. User must either provide rotation file(s) with --rotation-filenames, or specify a model name with (-m,--model), or provide it in the first positional arguments."
        )
    else:
        _logger.info(f"Using rotation files: {rotation_files}")

    if not continent_files:
        continent_files = None
    else:
        _logger.info(f"Using continent files: {continent_files}")

    ridges = values["ridges_filenames"]
    isochrons = values["isochron_filenames"]
    isochron_cob = values["isochron_cob_filenames"]
    if not ridges:
        raise Exception(
            "No ridges file(s) found. User must either provide ridges file(s) with --ridges-filenames, or specify a model name with (-m,--model), or provide it in the first positional arguments."
        )
    else:
        _logger.info(f"Using ridges files: {ridges}")

    if not isochrons:
        raise Exception(
            "No isochron file(s) found. User must either provide isochron file(s) with --isochron-filenames, or specify a model name with (-m,--model), or provide it in the first positional arguments."
        )
    else:
        _logger.info(f"Using isochron files: {isochrons}")

    if not isochron_cob:
        raise Exception(
            "No isochron cob file(s) found. User must either provide isochron cob file(s) with --isochron-cob-filenames, or specify a model name with (-m,--model), or provide it in the first positional arguments."
        )
    else:
        _logger.info(f"Using isochron cob files: {isochron_cob}")

    plate_recon = PlateReconstruction(
        rotation_model=rotation_files,
        topology_features=topology_files,
        anchor_plate_id=values["anchor_plate_id"],
    )
    o = IsochronSeafloorGrid(
        plate_reconstruction=plate_recon,
        time_steps=_get_time_steps(values),
        ridges=ridges,
        isochrons=isochrons,
        iso_cob=isochron_cob,
        continental_polygons=continent_files,
        interval_spacing_degrees=values["interval_spacing"],
        grid_output_dir=args.output_dir,
    )
    o.generate(set([OutputScalarType.AGE, OutputScalarType.SPREADING_RATE]))


def _get_rotation_and_topology_files_from_input_filenames(input_filenames):
    rotation_files = []
    topology_files = []
    for fn in input_filenames:
        features = FeaturesFunctionArgument([fn]).get_features()
        if len(features) > 0:
            if (
                features[0].get_feature_type().to_qualified_string()
                == "gpml:TotalReconstructionSequence"
            ):
                rotation_files.append(fn)
            else:
                topology_files.append(fn)
    return rotation_files, topology_files


def _get_rotation_topology_and_continent_files(values):
    rotation_files = values["rotation_filenames"]
    topology_files = values["topology_filenames"]
    continent_files = values["continents_filenames"]
    if not isinstance(rotation_files, list):
        rotation_files = [rotation_files]
    if not isinstance(topology_files, list):
        topology_files = [topology_files]
    if not isinstance(continent_files, list):
        continent_files = [continent_files]

    r, t = _get_rotation_and_topology_files_from_input_filenames(
        values["input_filenames"]
    )
    rotation_files = list(set(rotation_files) | set(r))
    topology_files = list(set(topology_files) | set(t))

    model_name = values["model_name"]
    if model_name:
        plate_model = PlateModelManager().get_model(
            model_name, data_dir="plate-model-repo"
        )
        if not plate_model:
            _logger.error(
                f"Unable to create PlateModel object for model {model_name}. "
                + f"Check your network connection and make sure the model name is correct. Run `pmm ls` to get a list of available model names."
            )
        else:
            if not rotation_files:
                rotation_files = plate_model.get_rotation_model()

            all_layers = plate_model.get_avail_layers()
            if not topology_files and "Topologies" in all_layers:
                topology_files = plate_model.get_layer("Topologies")

            if not continent_files and "ContinentalPolygons" in all_layers:
                continent_files = plate_model.get_layer("ContinentalPolygons")
                if "Cratons" in all_layers:
                    continent_files += plate_model.get_layer("Cratons")  # type: ignore

    return rotation_files, topology_files, continent_files


def _get_time_steps(values):
    if values["times"]:
        return values["times"]
    else:
        return list(
            range(
                int(values["min_time"]),
                int(values["max_time"]) + 1,
                values["time_interval"],
            )
        )
