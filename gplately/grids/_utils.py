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
import warnings
import numpy as np


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
