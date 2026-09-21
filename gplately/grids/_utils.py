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
import numbers
import warnings

import numpy as np

# How many decimal places of reconstruction time appear in output filenames, when the
# caller does not say. These reproduce the naming of the workflows the code was ported
# from: the predicting-sediment-thickness distance grids use one decimal place, and the
# simple_paleobathymetry grids use none.
DEFAULT_DISTANCE_GRID_DECIMAL_PLACES_IN_TIME = 1
DEFAULT_DECIMAL_PLACES_IN_TIME = 0


def resolve_decimal_places_in_time(decimal_places_in_time, default):
    """Resolve and validate a ``decimal_places_in_time`` argument.

    ``None`` selects `default`, which is the value reproducing the original workflow's
    filenames for the output in question.
    """
    if decimal_places_in_time is None:
        decimal_places_in_time = default
    if (
        isinstance(decimal_places_in_time, bool)
        or not isinstance(decimal_places_in_time, numbers.Integral)
        or decimal_places_in_time < 0
    ):
        raise ValueError(
            "decimal_places_in_time must be a non-negative integer (or None), "
            f"but got {decimal_places_in_time!r}"
        )
    return decimal_places_in_time


def format_time_in_filename(time, decimal_places_in_time):
    """Format a reconstruction time for use in an output filename.

    Uses the same rule as pyBacktrack's ``output_file_decimal_places_in_time`` and
    ``merge_paleo_bathymetry_file_decimal_places_in_time``, which build their format as
    ``{time:.Nf}``, so grids written here and grids pyBacktrack goes looking for are named
    identically.
    """
    return "{:.{}f}".format(time, decimal_places_in_time)


def distance_grid_filename(grid_spacing, time, decimal_places_in_time):
    """Name of the mean-distance grid for one age grid.

    Written by :func:`gplately.generate_distance_grids` and read back by the
    ``generate-sediment-grids`` CLI subcommand, so the two must agree.
    """
    return "mean_distance_{:.1f}d_{}.nc".format(
        grid_spacing, format_time_in_filename(time, decimal_places_in_time)
    )


def check_times_are_distinct_in_filenames(
    times, decimal_places_in_time, filename_template
):
    """Raise `ValueError` if two different times would be written to the same file.

    Output filenames carry the time to a fixed number of decimal places, so a time step
    finer than that resolution makes consecutive times collide -- the grids are computed
    and then silently overwritten one another. Checking up front turns that into an error
    before any of the work is done.

    Parameters
    ----------
    times : sequence of float
        The times that will be written.
    decimal_places_in_time : int
        Decimal places of time in the filenames, as already resolved by
        :func:`resolve_decimal_places_in_time`.
    filename_template : str
        The filename with ``{}`` where the formatted time goes, e.g.
        ``"paleobathymetry_{}Ma.nc"``. Used only to make the error message concrete.

    Raises
    ------
    ValueError
        If two distinct times format to the same filename.
    """
    filenames = {}
    for time in times:
        filename = filename_template.format(
            format_time_in_filename(time, decimal_places_in_time)
        )
        if filename in filenames:
            clashing_time = filenames[filename]
            if clashing_time == time:
                continue  # the same time twice is the caller's business, not a collision
            # Find the coarsest resolution that would actually separate these times.
            sufficient = next(
                (
                    places
                    for places in range(
                        decimal_places_in_time + 1, decimal_places_in_time + 11
                    )
                    if format_time_in_filename(clashing_time, places)
                    != format_time_in_filename(time, places)
                ),
                None,
            )
            remedy = (
                f"pass decimal_places_in_time={sufficient} (or higher)"
                if sufficient is not None
                else "use times that differ by more than 1e-9"
            )
            raise ValueError(
                f"times {clashing_time} and {time} would both be written to "
                f"'{filename}', so one would silently overwrite the other -- {remedy}, "
                "or use a coarser time step"
            )
        filenames[filename] = time


def num_grid_points(
    spacing: float,
    start: float,
    stop: float,
) -> int:
    """
    Compute the number of grid points needed to cover [start, stop] at
    approximately the given spacing.

    If a requested spacing does not divide an extent evenly, the number
    of grid points must be snapped to an integer -- which in turn changes the
    *effective* spacing slightly. Warns if the
    effective spacing differs non-negligibly from what was requested, and
    returns the number of grid points (nodes) -- one more than the number of
    intervals.

    Parameters
    ----------
    spacing : float
        Desired spacing between adjacent grid points (must be > 0).
    start, stop : float
        Coordinate extent to cover (stop must be > start).

    Returns
    -------
    int
        Number of grid points (>= 2).
    """
    if spacing <= 0:
        raise ValueError(f"spacing must be positive, got {spacing}")
    extent = stop - start
    if extent <= 0:
        raise ValueError(f"stop ({stop}) must be greater than start ({start})")

    n_intervals = int(np.floor(extent / spacing))
    n_intervals = max(1, n_intervals)

    if not np.isclose(n_intervals * spacing, extent, atol=1e-6):
        warnings.warn(
            f"Spacing {spacing} does not divide extent {extent} into an integer "
            f"number of intervals; using effective spacing {extent / n_intervals} instead."
        )

    return n_intervals + 1
