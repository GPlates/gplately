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
import time
import warnings
from functools import partial
from typing import Optional, Sequence, Union

from ..grids import read_netcdf_grid, write_netcdf_grid

logger = logging.getLogger("gplately")


def add_parser(parser):
    """add command line argument parser"""

    grid_cmd = parser.add_parser(
        "regrid",
        aliases=("rg",),
        help=help_str,
        add_help=True,
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # regrid command arguments
    grid_cmd.set_defaults(func=_regrid_netcdf4)
    grid_cmd.add_argument(
        metavar="INPUT_FILE",
        help="input file or directory",
        dest="input_grid_filename",
    )
    grid_cmd.add_argument(
        metavar="OUTPUT_FILE",
        help="output file or directory",
        dest="output_grid_filename",
    )
    grid_cmd.add_argument(
        "-r",
        "--resolution",
        metavar="RESOLUTION",
        type=float,
        help="grid resolution (degrees)",
        default=None,
        dest="grid_spacing",
    )
    grid_cmd.add_argument(
        "-d",
        "--significant_digits",
        metavar="DIGITS",
        type=int,
        help="Round to specified number of significant digits",
        default=None,
        dest="significant_digits",
    )
    grid_cmd.add_argument(
        "-j",
        "--n_jobs",
        help="number of processes to use; default if specified: use all CPU available",
        metavar="N_JOBS",
        default=None,
        dest="n_jobs",
    )


help_str = "Compress and/or resample netCDF4 grids."

__description__ = f"""{help_str}

Example usage: 
    - gplately grid inputfile.nc -r 0.1 -d 2 outputfile.nc
    - gplately grid inputfile.nc --resolution 0.1 --significant-digits 2 outputfile.nc
    - gplately grid inputDirectory --significant-digits 2 outputDirectory
"""


def _batch_regrid_netcdf4(
    input_raster_filename,
    output_raster_filename,
    resample=None,
    significant_digits=None,
):

    grid = read_netcdf_grid(input_raster_filename, resample=resample)
    write_netcdf_grid(
        output_raster_filename, grid, significant_digits=significant_digits
    )

    print("  {} complete!".format(output_raster_filename))


def _regrid_netcdf4(args):
    n_jobs = args.n_jobs
    if not n_jobs:
        try:
            n_jobs = multiprocessing.cpu_count()
        except NotImplementedError:
            n_jobs = 1
    start = time.time()

    # determine if directory or filename
    p_input = pathlib.Path(args.input_grid_filename)
    p_output = pathlib.Path(args.output_grid_filename)

    if args.grid_spacing is None:
        grid_spacing = None
    else:
        grid_spacing = (args.grid_spacing, args.grid_spacing)

    if p_input.is_file():
        input_filename = p_input
        if p_output.suffix:
            output_filename = p_output
        else:
            p_output.mkdir(parents=False, exist_ok=True)
            output_filename = p_output.joinpath(p_input.name)

        _batch_regrid_netcdf4(
            input_filename,
            output_filename,
            resample=grid_spacing,
            significant_digits=args.significant_digits,
        )

    elif p_input.is_dir():
        if p_output.suffix:
            raise ValueError("Specify output directory, not a single output file")

        # find all .nc files in this directory
        input_raster_filenames = []
        output_raster_filenames = []
        for pathname in p_input.iterdir():
            if (
                pathname.suffix == ".nc"
                and pathname.is_file()
                and not pathname.name.startswith(".")
            ):
                input_raster_filenames.append(pathname)
                output_raster_filenames.append(p_output.joinpath(pathname.name))

        print(
            "Found {} netCDF files. Processing...".format(len(input_raster_filenames))
        )

        if n_jobs == 1:
            for in_name, out_name in zip(
                input_raster_filenames, output_raster_filenames
            ):
                _batch_regrid_netcdf4(
                    in_name,
                    out_name,
                    resample=grid_spacing,
                    significant_digits=args.significant_digits,
                )
        else:

            with multiprocessing.Pool(n_jobs) as pool:
                pool.starmap(
                    partial(
                        _batch_regrid_netcdf4,
                        resample=grid_spacing,
                        significant_digits=args.significant_digits,
                    ),
                    zip(input_raster_filenames, output_raster_filenames),
                )

    else:
        raise ValueError("Input filename or path does not exist.")

    end = time.time()
    hours_minutes_seconds = str(datetime.timedelta(seconds=end - start)).split(":")
    logger.info(
        f"Completed gridding in {hours_minutes_seconds[0]} Hours, {hours_minutes_seconds[1]} Minutes, {hours_minutes_seconds[2].split('.')[0]} Seconds "
    )
