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
import datetime
import logging
import multiprocessing
import pathlib
import re
import time
from functools import partial
from typing import List, Optional, Union

from plate_model_manager import PlateModel, PlateModelManager

from ..grids import Raster

logger = logging.getLogger("gplately")

help_str = "Rotate a grid (or all grids in a folder) between plate-model reference frames."

__description__ = f"""{help_str}

The source ('from') and target ('to') rotation models can each be specified
either as a named plate model (downloaded via gplately's model registry) or
as one or more local rotation files.  The two options are mutually exclusive
within each side.

If the reconstruction time is not given explicitly with --time, the command
tries to deduce it from the filename by looking for a number immediately
followed by "Ma" (e.g. paleobathymetry_103Ma.nc → 103.0 Ma).

Example usage:
    - gplately rotate_grid input.nc output.nc --from-model Alfonso2024 --to-model Alfonso2024 --from-anchor 0 --to-anchor 701701 --time 100
    - gplately rotate_grid input_dir output_dir --from-model Alfonso2024 --to-model Alfonso2024 --from-anchor 0 --to-anchor 701701
    - gplately rotate_grid input.nc output.nc --from-rotation-files from.rot --to-rotation-files to.rot --time 100
"""


def add_parser(parser):
    """Add command line argument parser."""

    cmd = parser.add_parser(
        "rotate_grid",
        help=help_str,
        add_help=True,
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    cmd.set_defaults(func=_rotate_grid_cmd)

    cmd.add_argument(
        metavar="INPUT",
        help="input grid file (.nc) or directory containing .nc files",
        dest="input_path",
    )
    cmd.add_argument(
        metavar="OUTPUT",
        help="output grid file (.nc) or output directory",
        dest="output_path",
    )

    # --- source rotation model (mutually exclusive) ---
    from_group = cmd.add_mutually_exclusive_group(required=True)
    from_group.add_argument(
        "--from-model",
        metavar="MODEL_NAME",
        help="name of the source plate model (e.g. Alfonso2024, Muller2022)",
        dest="from_model_name",
        default=None,
    )
    from_group.add_argument(
        "--from-rotation-files",
        metavar="FILE",
        nargs="+",
        help="one or more source rotation files (.rot / .gpml)",
        dest="from_rotation_files",
        default=None,
    )

    # --- target rotation model (mutually exclusive) ---
    to_group = cmd.add_mutually_exclusive_group(required=True)
    to_group.add_argument(
        "--to-model",
        metavar="MODEL_NAME",
        help="name of the target plate model (e.g. Alfonso2024, Muller2022)",
        dest="to_model_name",
        default=None,
    )
    to_group.add_argument(
        "--to-rotation-files",
        metavar="FILE",
        nargs="+",
        help="one or more target rotation files (.rot / .gpml)",
        dest="to_rotation_files",
        default=None,
    )

    cmd.add_argument(
        "--from-anchor",
        metavar="PLATE_ID",
        type=int,
        help="anchor plate ID of the source reference frame (default: 0)",
        dest="from_anchor",
        default=0,
    )
    cmd.add_argument(
        "--to-anchor",
        metavar="PLATE_ID",
        type=int,
        help="anchor plate ID of the target reference frame (default: 0)",
        dest="to_anchor",
        default=0,
    )
    cmd.add_argument(
        "--non-reference-plate",
        metavar="PLATE_ID",
        type=int,
        help="arbitrary intermediate plate ID used during the rotation (default: 701)",
        dest="non_reference_plate",
        default=701,
    )
    cmd.add_argument(
        "-t",
        "--time",
        metavar="TIME_MA",
        type=float,
        help=(
            "reconstruction time in Ma; if omitted, the time is deduced from "
            "the filename (e.g. 103 from paleobathymetry_103Ma.nc)"
        ),
        dest="reconstruction_time",
        default=None,
    )
    cmd.add_argument(
        "-r",
        "--resolution",
        metavar="DEGREES",
        type=float,
        help="output grid spacing in degrees (default: 1.0)",
        default=1.0,
        dest="grid_spacing",
    )
    cmd.add_argument(
        "-j",
        "--n-jobs",
        metavar="N_JOBS",
        type=int,
        help="number of parallel worker processes (default: all available CPUs)",
        default=None,
        dest="n_jobs",
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_TIME_FROM_FILENAME_RE = re.compile(r"(\d+(?:\.\d+)?)Ma", re.IGNORECASE)


def _time_from_filename(path: pathlib.Path) -> Optional[float]:
    """Return the reconstruction time (Ma) encoded in *path*'s stem, or None."""
    match = _TIME_FROM_FILENAME_RE.search(path.stem)
    if match:
        return float(match.group(1))
    return None


def _get_rotation_files_for_model(
    model_name: str, data_dir: str = "plate-model-repo"
) -> List[str]:
    """Download (or locate from cache) and return rotation file paths for *model_name*."""
    try:
        plate_model = PlateModelManager().get_model(model_name, data_dir=data_dir)
    except Exception:
        plate_model = PlateModel(model_name, data_dir=data_dir, readonly=True)

    if not plate_model:
        raise ValueError(
            f"Could not find plate model '{model_name}'. "
            "Run `gplately list` to see available models."
        )
    return plate_model.get_rotation_model()


def rotate_single_grid(
    input_file: Union[str, pathlib.Path],
    output_file: Union[str, pathlib.Path],
    reconstruction_time: float,
    from_rotation_features_or_model,
    to_rotation_features_or_model,
    from_anchor: int = 0,
    to_anchor: int = 0,
    non_reference_plate: int = 701,
    grid_spacing: float = 1.0,
):
    """Rotate a single NetCDF grid file between two reference frames.

    Parameters
    ----------
    input_file : str or pathlib.Path
        Path to the input .nc grid file.
    output_file : str or pathlib.Path
        Path where the rotated .nc grid will be written.
    reconstruction_time : float
        Age (Ma) of the grid.
    from_rotation_features_or_model : str, list of str, or pygplates.RotationModel
        Rotation model / files for the source reference frame.
    to_rotation_features_or_model : str, list of str, or pygplates.RotationModel
        Rotation model / files for the target reference frame.
    from_anchor : int
        Anchor plate ID of the source reference frame.
    to_anchor : int
        Anchor plate ID of the target reference frame.
    non_reference_plate : int
        Intermediate plate used in the conversion (default 701).
    grid_spacing : float
        Output grid spacing in degrees.
    """
    raster = Raster(str(input_file))
    raster.rotate_reference_frames(
        grid_spacing_degrees=grid_spacing,
        reconstruction_time=reconstruction_time,
        from_rotation_features_or_model=from_rotation_features_or_model,
        to_rotation_features_or_model=to_rotation_features_or_model,
        from_rotation_reference_plate=from_anchor,
        to_rotation_reference_plate=to_anchor,
        non_reference_plate=non_reference_plate,
        output_name=str(output_file),
    )
    logger.info("  %s -> %s (%.1f Ma) done.", input_file, output_file, reconstruction_time)
    print(f"  {output_file} complete!")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def _rotate_grid_cmd(args):
    """Entry point called by argparse."""
    n_jobs = args.n_jobs
    if not n_jobs:
        try:
            n_jobs = multiprocessing.cpu_count()
        except NotImplementedError:
            n_jobs = 1

    start = time.time()

    # Resolve rotation file lists from model names (prime cache in main process
    # so parallel workers just do a fast local read).
    if args.from_model_name:
        from_rotation = _get_rotation_files_for_model(args.from_model_name)
    else:
        from_rotation = args.from_rotation_files

    if args.to_model_name:
        to_rotation = _get_rotation_files_for_model(args.to_model_name)
    else:
        to_rotation = args.to_rotation_files

    p_input = pathlib.Path(args.input_path)
    p_output = pathlib.Path(args.output_path)

    _worker = partial(
        _process_file,
        from_rotation=from_rotation,
        to_rotation=to_rotation,
        from_anchor=args.from_anchor,
        to_anchor=args.to_anchor,
        non_reference_plate=args.non_reference_plate,
        grid_spacing=args.grid_spacing,
        reconstruction_time=args.reconstruction_time,
    )

    if p_input.is_file():
        # --- single-file mode ---
        if p_output.suffix:
            output_file = p_output
        else:
            p_output.mkdir(parents=True, exist_ok=True)
            output_file = p_output / p_input.name
        _worker((p_input, output_file))

    elif p_input.is_dir():
        # --- directory mode ---
        if p_output.suffix:
            raise ValueError(
                "When INPUT is a directory, OUTPUT must also be a directory, not a file."
            )
        p_output.mkdir(parents=True, exist_ok=True)

        pairs = [
            (nc_file, p_output / nc_file.name)
            for nc_file in sorted(p_input.iterdir())
            if nc_file.suffix == ".nc" and nc_file.is_file() and not nc_file.name.startswith(".")
        ]

        if not pairs:
            raise ValueError(f"No .nc files found in '{p_input}'.")

        print(f"Found {len(pairs)} .nc file(s). Processing...")

        if n_jobs == 1:
            for pair in pairs:
                _worker(pair)
        else:
            with multiprocessing.Pool(n_jobs) as pool:
                pool.map(_worker, pairs)

    else:
        raise ValueError(f"Input path '{p_input}' does not exist.")

    end = time.time()
    h, m, s = str(datetime.timedelta(seconds=end - start)).split(":")
    logger.info(
        "Completed rotate_grid in %s Hours, %s Minutes, %s Seconds.",
        h,
        m,
        s.split(".")[0],
    )


def _process_file(
    pair,
    from_rotation,
    to_rotation,
    from_anchor,
    to_anchor,
    non_reference_plate,
    grid_spacing,
    reconstruction_time,
):
    """Worker function: rotate one (input_file, output_file) pair."""
    input_file, output_file = pair

    if reconstruction_time is None:
        t = _time_from_filename(pathlib.Path(input_file))
        if t is None:
            raise ValueError(
                f"Cannot deduce reconstruction time from filename '{input_file.name}'. "
                "Provide it explicitly with --time."
            )
    else:
        t = reconstruction_time

    rotate_single_grid(
        input_file=input_file,
        output_file=output_file,
        reconstruction_time=t,
        from_rotation_features_or_model=from_rotation,
        to_rotation_features_or_model=to_rotation,
        from_anchor=from_anchor,
        to_anchor=to_anchor,
        non_reference_plate=non_reference_plate,
        grid_spacing=grid_spacing,
    )
