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
import logging

import numpy as np  # pyright: ignore[reportMissingImports]

logger = logging.getLogger("gplately")


def _convert_longitude(grid_data, lons, map_fn, valid_max, seam_name, tol=1e-6):
    """
    Shared implementation for converting a grid's longitude convention.

    Parameters
    ----------
    grid_data : np.ndarray
        2D array of shape (len(lats), len(lons)).
    lons : array-like
        1D array of longitude values, monotonically increasing.
    map_fn : callable
        Function mapping raw longitude values into the target convention's
        base range, e.g. `lambda l: l % 360` or `lambda l: ((l + 180) % 360) - 180`.
    valid_max : float
        Upper bound of the target convention's range (360 or 180).
    seam_name : str
        Human-readable seam description used in error messages
        (e.g. "0/360" or "-180/180").

    Returns
    -------
    new_grid_data, new_lons

    Raises
    ------
    ValueError
        If grid_data's shape doesn't match lons, or if converting a *regional* grid would
        cross the seam and split the data into disconnected pieces.
    """
    lons = np.asarray(lons, dtype=float)
    grid_data = np.asarray(grid_data)

    if grid_data.shape[1] != lons.size:
        raise ValueError(
            f"grid_data has {grid_data.shape[1]} columns but lons has "
            f"{lons.size} values; they must match."
        )

    # Unwrap first: this accepts input that is already logically monotonic
    # but crosses the +/-180 (or 0/360) seam in mixed convention, e.g.
    # [160, 170, 180, -170, -160]. Any jump larger than half a period
    # (180 degrees) is corrected by adding/subtracting 360, turning the
    # array into a genuinely monotonically increasing sequence.
    lons = np.unwrap(lons, period=360)

    diffs = np.diff(lons)
    if not np.all(diffs > 0):
        bad_idx = np.where(diffs <= 0)[0]
        logger.error(
            "Non-monotonic longitude steps at indices %s: "
            "lons[i]=%s, lons[i+1]=%s, diff=%s",
            bad_idx.tolist(),
            lons[bad_idx],
            lons[bad_idx + 1],
            diffs[bad_idx],
        )
        raise ValueError(
            "The `lons` must be monotonically increasing (accounting for antimeridian wraparound)."
        )

    # Strip a redundant closing column (e.g. lons[0] and lons[-1] are the
    # same physical meridian, 360 degrees apart) -- handle separately at the end.
    has_duplicate_wrap = np.isclose(lons[-1] - lons[0], 360)
    if has_duplicate_wrap:
        core_lons, core_grid = lons[:-1], grid_data[:, :-1]
    else:
        core_lons, core_grid = lons, grid_data

    avg_spacing = np.mean(np.diff(core_lons))
    is_global = (core_lons[-1] - core_lons[0]) >= (360 - 1.5 * avg_spacing)

    # Step 1: simple relabeling (already monotonic after mapping)?
    simple = map_fn(core_lons)
    if np.all(np.diff(simple) > 0):
        new_core_lons, new_core_grid = simple, core_grid.copy()

    else:
        # Step 2: uniform shift -- handles boundary-only cases (e.g. a
        # regional grid starting/ending exactly at the seam).
        start = map_fn(core_lons[:1])[0]
        shifted = core_lons + (start - core_lons[0])

        if shifted[-1] <= valid_max + tol:
            new_core_lons, new_core_grid = shifted, core_grid.copy()

        elif is_global:
            # Step 3: genuine global wraparound -- roll columns.
            wrap_idx = np.argmax(np.diff(simple) < 0) + 1
            new_core_lons = np.concatenate([simple[wrap_idx:], simple[:wrap_idx]])
            new_core_grid = np.concatenate(
                [core_grid[:, wrap_idx:], core_grid[:, :wrap_idx]], axis=1
            )
            if not np.all(np.diff(new_core_lons) > 0):
                raise ValueError(
                    "Unexpected error: could not produce a monotonically "
                    "increasing longitude array after reordering."
                )
        else:
            raise ValueError(
                f"Converting this regional grid would cross the "
                f"{seam_name} seam, splitting the data into disconnected "
                f"pieces. This is only supported for global grids."
            )

    if has_duplicate_wrap:
        new_lons = np.append(new_core_lons, new_core_lons[0] + 360)
        new_grid = np.concatenate([new_core_grid, new_core_grid[:, :1]], axis=1)
    else:
        new_lons, new_grid = new_core_lons, new_core_grid

    return new_grid, new_lons


def to_longitude_positive_360(grid_data: np.ndarray, lons):
    """Convert a grid's longitude coordinates to the [0, 360] convention."""
    return _convert_longitude(
        grid_data,
        lons,
        map_fn=lambda l: l % 360,
        valid_max=360,
        seam_name="0/360",
    )


def to_longitude_signed_180(grid_data: np.ndarray, lons):
    """Convert a grid's longitude coordinates to the [-180, 180] convention."""
    return _convert_longitude(
        grid_data,
        lons,
        map_fn=lambda l: ((l + 180) % 360) - 180,
        valid_max=180,
        seam_name="-180/180",
    )


def upwrap_antimeridian_wraparound(lons):
    """
    Unwrap a longitude array which is wrapped around the antimeridian.

    Parameters
    ----------
    lons : array-like
        1D array of longitude values.

    Returns
    -------
    np.ndarray
        Unwrapped longitude array.
    """
    # Unwrap: this accepts input that is already logically monotonic
    # but crosses the +/-180 (or 0/360) seam in mixed convention, e.g.
    # [160, 170, 180, -170, -160], [340, 350, 360, 10, 20]. Any jump larger than half a period
    # (180 degrees) is corrected by adding/subtracting 360, turning the
    # array into a genuinely monotonically increasing sequence.
    lons = np.asarray(lons, dtype=float)
    if not np.all(np.diff(lons) > 0):
        lons = np.unwrap(lons, period=360)

    if np.any(lons < -180):
        lons += 360 * np.ceil((-180 - np.min(lons)) / 360)
    if np.any(lons > 360):
        lons -= 360 * np.ceil((np.max(lons) - 360) / 360)

    if not np.all(np.diff(lons) > 0):
        raise ValueError(
            "The input lons cannot be unwrapped to a monotonically increasing array. "
            "Please check the input values."
        )
    return lons
