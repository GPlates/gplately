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
import math
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


def time_index_at_or_after(time, time_increment):
    """Index of the first multiple of `time_increment` at or after `time`.

    A time that is a multiple of the increment must land on itself. Plain
    ``ceil(time / time_increment)`` does not guarantee that for floats: the times produced
    by a fractional step accumulate representation error, so 0.30000000000000004 / 0.1 is
    3.0000000000000004 and ceils to 4 -- snapping a 0.3 Ma grid to 0.4 Ma. Rounding first
    when the quotient is a whole number to within a tolerance keeps the intent.

    :func:`check_time_increment_covers_times` accepts exactly the times this function maps
    to themselves, so the check and the behaviour it guards cannot disagree.
    """
    multiples = time / time_increment
    nearest = round(multiples)
    if math.isclose(multiples, nearest, rel_tol=1e-9, abs_tol=1e-9):
        return int(nearest)
    return int(math.ceil(multiples))


def check_time_increment_covers_times(times, time_increment):
    """Validate `time_increment` for a backward reconstruction over `times`.

    Each age grid's walk starts at its own time snapped *up* to a multiple of
    `time_increment`, so that every age grid steps on one shared grid of times. A time that
    is not a multiple of the increment is therefore reconstructed from the wrong starting
    time -- silently, and with a plausible-looking result.

    Raises
    ------
    ValueError
        If `time_increment` is not positive, or if some time is not a multiple of it.
    """
    if time_increment <= 0:
        raise ValueError(f"time_increment must be positive, but got {time_increment!r}")

    for time in times:
        snapped_time = time_index_at_or_after(time, time_increment) * time_increment
        if not math.isclose(snapped_time, time, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(
                f"time {time} is not a multiple of time_increment {time_increment}, so it "
                f"would be reconstructed from {snapped_time} Ma instead -- use a "
                "time_increment that divides every time"
            )


def uniform_time_step(times, default=1.0):
    """The single time step separating `times`.

    Some APIs -- pyBacktrack's ``reconstruct_paleo_bathymetry_grids()`` among them -- do not
    take a list of times. They take a youngest time, an oldest time and one increment, and
    generate output at every step in between. Handing such an API a different increment
    makes it produce a different set of times than the caller asked for, so the increment
    has to be recovered from the times rather than borrowed from some other quantity.

    Parameters
    ----------
    times : sequence of float
        The times output is wanted at.
    default : float, default: 1.0
        Returned when there are fewer than two distinct times, where any positive increment
        produces the same single output.

    Raises
    ------
    ValueError
        If the times are not evenly spaced, and so cannot be described by one increment.
    """
    distinct_times = sorted(set(times))
    if len(distinct_times) < 2:
        return default

    steps = [
        later - earlier for earlier, later in zip(distinct_times, distinct_times[1:])
    ]
    if not all(
        math.isclose(step, steps[0], rel_tol=1e-9, abs_tol=1e-9) for step in steps
    ):
        raise ValueError(
            "times must be evenly spaced to be expressed as a single time increment, "
            f"but {distinct_times} step by {sorted(set(steps))}"
        )
    return steps[0]


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
