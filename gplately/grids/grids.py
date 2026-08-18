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

"""This sub-module contains tools for working with MaskedArray, ndarray and netCDF4 rasters, as well as gridded-data."""

import logging
import math
import os
import warnings
from multiprocessing import cpu_count
from typing import Tuple, Union, cast, overload, Literal

import matplotlib.colors
import netCDF4
import numpy as np
import pygplates
from pygplates import (
    RotationModel as _RotationModel,
    FiniteRotation as _FiniteRotation,
    Feature as _Feature,
    FeaturesFunctionArgument as _FeaturesFunctionArgument,
    FeatureCollection as _FeatureCollection,
)
from rasterio.enums import MergeAlg
from rasterio.features import rasterize as _rasterize
from rasterio.transform import from_bounds as _from_bounds
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy.spatial import (
    cKDTree as _cKDTree,  # pyright: ignore[reportAttributeAccessIssue]
)
from scipy.spatial.transform import Rotation as _Rotation

from ..geometry import pygplates_to_shapely

# re-export, for backward compatibility, don't remove
from ..lib.regular_grid_interpolator import RegularGridInterpolator
from ..raster import Raster

logger = logging.getLogger("gplately")

__all__ = [
    "fill_raster",
    "read_netcdf_grid",
    "write_netcdf_grid",
    "default_netcdf_fill_value",
    "RegularGridInterpolator",
    "sample_grid",
    "reconstruct_grid",
    "rasterise",
    "rasterize",
    "Raster",
    # "TimeRaster",
]


def fill_raster(data, invalid=None):
    """Search a grid of ``data`` for invalid cells (i.e NaN-type entries) and fill each
    invalid cell with the value of its nearest valid neighbour.

    .. note::

        Uses scipy's ``distance_transform_edt`` function to perform an Exact Euclidean
        Distance Transform (EEDT). This locates the nearest valid neighbours of an invalid
        ``data`` cell.

        An optional parameter, ``invalid``, is a binary ndarray with the same dimensions
        as ``data`` and the following entries:

        * 1 if its corresponding entry in ``data`` is of NaN-type;
        * 0 if not NaN-type

        This will be used to locate nearest neighbour fill values during the Exact Euclidian
        Distance Transform. If ``invalid`` is not passed to ``fill_raster``, it will be created
        for the user.

    Parameters
    ----------
    data : MaskedArray
        A MaskedArray of data that may have invalid cells (i.e. entries of type NaN).

    invalid : ndarray, optional, default=None
        An ndarray with the same shape as ``data`` whose elements are 1 if its corresponding
        elements in ``data`` are of type ``NaN``, and 0 if its corresponding entries in ``data``
        are valid. An optional parameter - this will be created for the user if it isn't
        provided.

    Returns
    -------
    data : ndarray
        An updated ``data`` array where each invalid cell has been replaced with the value
        of its nearest valid neighbour.
    """
    masked_array = hasattr(data, "fill_value")
    mask_fill_value = None
    if masked_array:
        mask_fill_value = data.data == data.fill_value
        data = data.data.copy()
        data[mask_fill_value] = np.nan
    else:
        data = data.copy()

    if invalid is None:
        invalid = np.isnan(data)
        if masked_array:
            invalid += mask_fill_value
    ind = distance_transform_edt(invalid, return_distances=False, return_indices=True)
    assert ind is not None
    return data[tuple(ind)]


def _realign_grid(array, lons, lats):
    """realigns grid to -180/180 and flips the array if the latitudinal coordinates are decreasing."""
    lons = np.asarray(lons)
    lats = np.asarray(lats)

    # There must not be any duplicate longitudes or duplicate latitudes.
    lon_differences = np.diff(lons)
    if np.any(lon_differences == 0):
        raise ValueError("Longitudes contain duplicate values.")
    lat_differences = np.diff(lats)
    if np.any(lat_differences == 0):
        raise ValueError("Latitudes contain duplicate values.")

    # Check if longitudes and latitudes are in increasing order. If not then sort them.
    if not np.all(lon_differences > 0):
        sort_indices = np.argsort(lons)
        array = array[:, sort_indices]
        lons = lons[sort_indices]
    if not np.all(lat_differences > 0):
        sort_indices = np.argsort(lats)
        array = array[sort_indices, :]
        lats = lats[sort_indices]

    # If we need to wrap (180, 360) to (-180, 0).
    if lons[-1] > 180:
        mask_lon_gt_180 = lons > 180

        # If we have longitudes at 0 and 360 then we don't want to wrap the 360 column to 0 since we'd end up with two 0 columns.
        # Both columns should ideally be equal anyway (if input raster wraps 0->360 properly).
        if np.isclose(lons[0], 0.0) and np.isclose(lons[-1], 360.0):
            mask_lon_wrap = mask_lon_gt_180.copy()
            mask_lon_wrap[-1] = False  # drop the 360 column altogether
        else:
            mask_lon_wrap = mask_lon_gt_180

        # Wrap (180, 360) to (-180, 0).
        array = np.hstack([array[:, mask_lon_wrap], array[:, ~mask_lon_gt_180]])
        lons = np.hstack([lons[mask_lon_wrap] - 360.0, lons[~mask_lon_gt_180]])

        # If the input grid crossed the dateline (180) then create a matching column at -180 so that the output wraps -180->180 properly.
        #
        # The dateline was crossed if there are longitudes in (0, 180), noting that we already wrapped (180, 360) to (-180, 0).
        if lons[-1] > 0:
            if np.isclose(lons[-1], 180.0):
                # There's a longitude at 180, so duplicate it at -180.
                array = np.hstack([array[:, [-1]], array])
                lons = np.hstack([-180.0, lons])
            else:
                # There's no longitude at 180, so interpolate at 180 and insert at -180.
                #
                # The 360.0 accounts for the fact that lons[0] was wrapped above.
                interp_180_weight = (180.0 - lons[-1]) / (360.0 - (lons[-1] - lons[0]))
                interp_180_column = (
                    array[:, -1] * (1.0 - interp_180_weight)
                    + array[:, 0] * interp_180_weight
                )
                array = np.hstack([interp_180_column[:, np.newaxis], array])
                lons = np.hstack([-180.0, lons])

    return array, lons, lats


def _guess_data_variable_name(cdf: netCDF4.Dataset, x_name: str, y_name: str) -> Union[str, None]:  # type: ignore
    """best effort to find out the data variable name"""
    vars = cdf.variables.keys()
    for var in vars:
        dimensions = cdf.variables[var].dimensions
        if len(dimensions) != 2:  # only consider two-dimensional data
            continue
        else:
            if dimensions[0] == y_name and dimensions[1] == x_name:
                return var
    return None


def _is_a_common_name_for_longitude(name: str) -> bool:
    """Return True if the `name` parameter is a possible common name for longitude."""
    return name in ["lon", "lons", "longitude", "x", "east", "easting", "eastings"]


def _is_a_common_name_for_latitude(name: str) -> bool:
    """Return True if the `name` parameter is a possible common name for latitude."""
    return name in ["lat", "lats", "latitude", "y", "north", "northing", "northings"]


def _spaced_axis(start, stop, step):
    """Build an inclusive coordinate axis from `start` to `stop`, sampled every `step`.

    Equivalent to ``np.arange(start, stop + step, step)``, but with exact endpoints.
    Accumulated floating-point error in ``np.arange`` can push the final sample past
    `stop` — a 0.2-degree global latitude axis ends at 90.00000000000256, which trips
    the pole-clipping guard in `sample_grid`. Deriving the sample count and handing it
    to `np.linspace` keeps both endpoints exact.

    `step` must share the sign of ``stop - start``, so descending axes (as produced by
    an ``upper`` origin) are handled the same way as ascending ones. As with
    ``np.arange``, a mis-signed `step` yields an empty axis.
    """
    n = int(round((stop - start) / step)) + 1
    return np.linspace(start, stop, max(n, 0))


def _find_extent_from_data(
    data, origin
) -> Union[Tuple[float, float, float, float], None]:
    """Try to find the extent from data. Return None if data doesn't contain coordinates.
    As of 2025-12-10, only support xarray.DataArray."""
    extent = None
    lons = None
    lats = None
    try:
        for name in data.coords:
            if not lats and _is_a_common_name_for_latitude(name):
                lats = data.coords[name]
            elif not lons and _is_a_common_name_for_longitude(name):
                lons = data.coords[name]
        if lons is not None and lats is not None:
            extent = (
                float(lons.min()),
                float(lons.max()),
                float(lats.min()),
                float(lats.max()),
            )
    except Exception as ex:
        logger.debug(ex)
        return None

    return _adjust_extent_for_origin(extent, origin)


def read_netcdf_grid(
    filename,
    return_grids: bool = False,
    realign: bool = False,
    resample=None,
    resize=None,
    x_dimension_name: str = "",
    y_dimension_name: str = "",
    data_variable_name: str = "",
) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray], np.ndarray]:
    """Read grid data from a NetCDF (.nc) file.

    Parameters
    ----------
    filename : str
        Full path to the ``netCDF`` raster file.
    return_grids : bool, optional, default=False
        If set to ``True``, returns lon, lat arrays associated with the grid data.
    realign : bool, optional, default=False
        if set to ``True``, realigns grid to -180/180 and flips the array if the latitudinal coordinates are decreasing.
    resample : tuple, optional, default=None
        If provided as ``resample = (spacingX, spacingY)``, the grid data will be resampled with these x and y resolutions.
    resize : tuple, optional, default=None
        If provided as ``resample = (resX, resY)``, the grid data will be resized to the number of columns (resX) and rows (resY).
    x_dimension_name : str, optional, default=""
        If the grid file uses the comman names, such as ``x``, ``lon``, ``lons`` or ``longitude``,
        you need not to provide this parameter. Otherwise, you need to tell us what the x dimension name is.
    y_dimension_name : str, optional, default=""
        If the grid file uses the comman names, such as ``y``, ``lat``, ``lats`` or ``latitude``,
        you need not to provide this parameter. Otherwise, you need to tell us what the y dimension name is.
    data_variable_name : str, optional, default=""
        GPlately will try its best to guess the data variable name.
        However, it would be much better if you tell us what the data variable name is.
        Otherwise, GPlately's guess may/may not be correct.

    Returns
    -------
    grid_z : `MaskedArray`_
        A `MaskedArray`_ object containing the grid data. The longitudes are re-aligned between -180 and 180 degrees.
    lon, lat : `MaskedArray`_
        When ``return_grids`` is ``True``, return two additional `MaskedArray`_ objects containing the longitudes and latitudes of the grid data.


    .. _MaskedArray: https://numpy.org/doc/stable/reference/maskedarray.generic.html
    """

    def find_label(keys, labels):
        for label in labels:
            if label in keys:
                return label
        return None

    # possible permutations of lon/lat/z
    label_lon = ["lon", "lons", "longitude", "x", "east", "easting", "eastings"]
    label_lat = ["lat", "lats", "latitude", "y", "north", "northing", "northings"]
    label_z = ["z", "data", "values", "Band1", "__xarray_dataarray_variable__"]

    # add capitalise and upper case permutations
    label_lon = (
        label_lon
        + [label.capitalize() for label in label_lon]
        + [label.upper() for label in label_lon]
    )
    label_lat = (
        label_lat
        + [label.capitalize() for label in label_lat]
        + [label.upper() for label in label_lat]
    )
    label_z = (
        label_z
        + [label.capitalize() for label in label_z]
        + [label.upper() for label in label_z]
    )

    # open netCDF file and re-align from -180, 180 degrees
    with netCDF4.Dataset(filename, "r") as cdf:
        keys = cdf.variables.keys()

        # find the names of variables
        if data_variable_name:
            key_z = data_variable_name
        else:
            key_z = find_label(keys, label_z)
        if x_dimension_name:
            key_lon = x_dimension_name
        else:
            key_lon = find_label(keys, label_lon)
        if y_dimension_name:
            key_lat = y_dimension_name
        else:
            key_lat = find_label(keys, label_lat)

        if key_lon is None or key_lat is None:
            raise ValueError(
                f"Cannot find x,y or lon/lat coordinates in netcdf. The dimensions in the file are {cdf.dimensions.keys()}"
            )

        if key_z is None:
            key_z = _guess_data_variable_name(cdf, key_lon, key_lat)

        if key_z is None:
            raise ValueError(
                f"Cannot find z data in netcdf. The variables in the file are {cdf.variables.keys()}"
            )

        # extract data from cdf variables
        # TODO: the dimensions of data may not be (lat, lon). It is possible(but unlikely?) that the dimensions are(lon, lat).
        # just note you may need numpy.swapaxes() here.
        if len(cdf[key_z].dimensions) != 2:
            raise Exception(
                f"The data in the netcdf file is not two-dimensional. This function can only handle two-dimensional data."
                + f"The dimensions in the file are {cdf[key_z].dimensions.keys()}"
            )
        cdf_grid = cdf[key_z][:]
        cdf_lon = cdf[key_lon][:]
        cdf_lat = cdf[key_lat][:]

        # fill missing values
        if np.issubdtype(cdf_grid.dtype, np.floating):
            if hasattr(cdf[key_z], "missing_value"):
                fill_value = cdf[key_z].missing_value
                cdf_grid[np.isclose(cdf_grid, fill_value, rtol=0.1)] = np.nan
            elif hasattr(cdf[key_z], "_FillValue"):
                fill_value = cdf[key_z]._FillValue
                cdf_grid[np.isclose(cdf_grid, fill_value, rtol=0.1)] = np.nan

        # convert to boolean array
        if np.issubdtype(cdf_grid.dtype, np.integer):
            unique_grid = np.unique(cdf_grid)
            if len(unique_grid) == 2:
                if (unique_grid == [0, 1]).all():
                    cdf_grid = cdf_grid.astype(bool)

    # we realign the grid to -180/180 when the longitudes are from 0 to 360
    # this is a temporary fix. we need a more sophisticated solution.
    if np.max(cdf_lon) > 180:
        # realign longitudes to -180/180 dateline
        cdf_grid_z, cdf_lon, cdf_lat = _realign_grid(cdf_grid, cdf_lon, cdf_lat)
    else:
        cdf_grid_z = cdf_grid

    # resample
    if resample is not None:
        spacingX, spacingY = resample

        # don't resample if already the same resolution
        dX = np.diff(cdf_lon).mean()
        dY = np.diff(cdf_lat).mean()

        if not np.isclose(dX, spacingX) or not np.isclose(dY, spacingY):
            lon_grid = _spaced_axis(cdf_lon.min(), cdf_lon.max(), spacingX)
            lat_grid = _spaced_axis(cdf_lat.min(), cdf_lat.max(), spacingY)
            lonq, latq = np.meshgrid(lon_grid, lat_grid)
            original_extent = (
                cdf_lon[0],
                cdf_lon[-1],
                cdf_lat[0],
                cdf_lat[-1],
            )
            cdf_grid_z = sample_grid(
                lonq,
                latq,
                cdf_grid_z,
                method="nearest",
                extent=original_extent,
                return_indices=False,
            )
            cdf_lon = lon_grid
            cdf_lat = lat_grid

    # resize
    if resize is not None:
        resX, resY = resize

        # don't resize if already the same shape
        if resX != cdf_grid_z.shape[1] or resY != cdf_grid_z.shape[0]:  # type: ignore
            original_extent = (
                cdf_lon[0],
                cdf_lon[-1],
                cdf_lat[0],
                cdf_lat[-1],
            )
            lon_grid = np.linspace(original_extent[0], original_extent[1], resX)
            lat_grid = np.linspace(original_extent[2], original_extent[3], resY)
            lonq, latq = np.meshgrid(lon_grid, lat_grid)

            cdf_grid_z = sample_grid(
                lonq,
                latq,
                cdf_grid_z,
                method="nearest",
                extent=original_extent,
                return_indices=False,
            )
            cdf_lon = lon_grid
            cdf_lat = lat_grid

    # Fix grids with 9e36 as the fill value for nan.
    # cdf_grid_z.fill_value = float('nan')
    # cdf_grid_z.data[cdf_grid_z.data > 1e36] = cdf_grid_z.fill_value

    if return_grids:
        return cdf_grid_z, cdf_lon, cdf_lat
    else:
        return cdf_grid_z


def write_netcdf_grid(
    filename,
    grid,
    extent: Union[tuple, str] = "global",
    significant_digits=None,
    fill_value: Union[str, float, bool, None] = None,
    metadata: Union[dict, None] = None,
    title: Union[str, None] = None,
):
    """Write geological data contained in a ``grid`` to a netCDF4 grid with a specified ``filename``.

    Notes
    -----
    The written netCDF4 grid has the same latitudinal and longitudinal (row and column) dimensions as ``grid``.
    It has three variables:

    * Latitudes of ``grid`` data
    * Longitudes of ``grid`` data
    * The data stored in ``grid``

    However, the latitudes and longitudes of the grid returned to the user are constrained to those
    specified in ``extent``.
    By default, ``extent`` assumes a global latitudinal and longitudinal span: `extent=[-180,180,-90,90]`.

    Parameters
    ----------
    filename : str
        The full path (including a filename and the ".nc" extension) to save the created netCDF4 ``grid`` to.

    grid : array-like
        An ndarray grid containing data to be written into a `netCDF` (.nc) file. Note: Rows correspond to
        the data's latitudes, while the columns correspond to the data's longitudes.

    extent : list, default=[-180,180,-90,90]
        Four elements that specify the [min lon, max lon, min lat, max lat] to constrain the lat and lon
        variables of the netCDF grid to. If no extents are supplied, full global extent `[-180, 180, -90, 90]`
        is assumed.

    significant_digits : int, optional
        Optionally applies lossy data compression up to a specified number of significant digits.
        This significantly reduces file size, but make sure the required precision is preserved in the
        saved netcdf file.

    fill_value : scalar or False or None, default=None
        Value used to fill in missing data.

        If ``False`` is specified then no fill value is used, and you must ensure that data was written to *all* elements of ``grid``.
        And any NaN elements will be written as the raw bit pattern for NaN (rather than automatically converted to a fill value).

        If ``None`` is specified then a default fill value is used, as follows:

        * If ``significant_digits`` is NOT specified and the grid data type is floating-point then the default is `np.nan`.
          This essentially means that `np.nan` is only used (as the default) when *losslessly* compressing *floating-point* data.
          This is because *lossy* compression with a NaN fill value appears to *not* always mask out NaN regions.
        * In all other cases the default is determined by `netCDF` based on the grid type.
          For example, the default for floating-point types is 9.969209968386869e+36 (see `netCDF4.default_fillvals`) and
          the default for *signed* integers is the largest negative value supported by the integer type (for *unsigned* its largest value).
          In this case, please ensure the default is outside the range of your valid grid data, otherwise specify a custom fill value.

        ...and to query the default value associated with ``None`` you can call :func:`default_netcdf_fill_value`.

    metadata : dict, default=None
        Optional metadata to store as global netCDF attributes.

    title : str, default=None
        Title to store as the global ``title`` netCDF attribute. If ``None``,
        defaults to ``"Grid produced by gplately <version>"``.

    Returns
    -------
    A netCDF grid will be saved to the path specified in ``filename``.
    """
    from gplately import __version__ as _version

    if extent == "global":
        extent = (-180, 180, -90, 90)
    else:
        extent = tuple(extent)
        assert len(extent) == 4, "specify the [min lon, max lon, min lat, max lat]"

    nrows, ncols = np.shape(grid)

    assert isinstance(extent, tuple)
    lon_grid = np.linspace(extent[0], extent[1], ncols)
    lat_grid = np.linspace(extent[2], extent[3], nrows)

    data_kwds = {"compression": "zlib", "complevel": 6}

    def _compute_coordinate_bounds(coords):
        coords = np.asarray(coords, dtype=float)
        if coords.size == 1:
            edges = np.array([coords[0] - 0.5, coords[0] + 0.5], dtype=float)
        else:
            edges = np.empty(coords.size + 1, dtype=float)
            mids = 0.5 * (coords[:-1] + coords[1:])
            edges[1:-1] = mids
            edges[0] = coords[0] - (mids[0] - coords[0])
            edges[-1] = coords[-1] + (coords[-1] - mids[-1])
        return np.column_stack((edges[:-1], edges[1:]))

    with netCDF4.Dataset(filename, "w", driver=None) as cdf:
        if title is None:
            cdf.title = "Grid produced by gplately " + str(_version)
        else:
            cdf.title = str(title)
        if metadata:
            for key, value in metadata.items():
                if value is None:
                    continue
                attr_name = str(key).strip()
                if not attr_name:
                    continue
                if isinstance(value, np.generic):
                    value = value.item()
                elif isinstance(value, tuple):
                    value = list(value)
                elif isinstance(value, (list, dict, str, int, float, bool)):
                    pass
                else:
                    value = str(value)
                cdf.setncattr(attr_name, value)

        # ACDD-style geospatial discovery metadata derived from extent.
        lon_min = float(min(extent[0], extent[1]))
        lon_max = float(max(extent[0], extent[1]))
        lat_min = float(min(extent[2], extent[3]))
        lat_max = float(max(extent[2], extent[3]))
        cdf.geospatial_lon_min = lon_min
        cdf.geospatial_lon_max = lon_max
        cdf.geospatial_lat_min = lat_min
        cdf.geospatial_lat_max = lat_max
        cdf.geospatial_bounds = (
            f"POLYGON (({lon_min} {lat_min}, {lon_max} {lat_min}, "
            f"{lon_max} {lat_max}, {lon_min} {lat_max}, {lon_min} {lat_min}))"
        )

        cdf.createDimension("lon", lon_grid.size)
        cdf.createDimension("lat", lat_grid.size)
        cdf.createDimension("bnds", 2)
        cdf_lon = cdf.createVariable("lon", lon_grid.dtype, ("lon",), **data_kwds)
        cdf_lat = cdf.createVariable("lat", lat_grid.dtype, ("lat",), **data_kwds)
        cdf_lon[:] = lon_grid
        cdf_lat[:] = lat_grid

        cdf_lon.units = "degrees_east"
        cdf_lon.standard_name = "lon"
        cdf_lon.bounds = "lon_bnds"
        cdf_lon.actual_range = [lon_grid[0], lon_grid[-1]]

        cdf_lat.units = "degrees_north"
        cdf_lat.standard_name = "lat"
        cdf_lat.bounds = "lat_bnds"
        cdf_lat.actual_range = [lat_grid[0], lat_grid[-1]]

        # create container variable for CRS: lon/lat WGS84 datum
        crso = cdf.createVariable("crs", "i4")
        crso.long_name = "Lon/Lat Coords in WGS84"
        crso.grid_mapping_name = "latitude_longitude"
        crso.longitude_of_prime_meridian = 0.0
        crso.semi_major_axis = 6378137.0
        crso.inverse_flattening = 298.257223563
        crso.spatial_ref = """GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]"""

        # add more keyword arguments for quantizing data
        if significant_digits:
            data_kwds["significant_digits"] = int(significant_digits)
            data_kwds["quantize_mode"] = "GranularBitRound"

        # The fill value can be False, but not True.
        if isinstance(fill_value, bool):
            if fill_value:
                raise ValueError(
                    "'fill_value' cannot be True; it should be False, None or a number"
                )

        if grid.dtype is np.dtype(bool):
            # Boolean arrays need to be converted to integers since
            # there's no such thing as a mask on a boolean array.
            grid = grid.astype("i1")
            fill_value = False  # no pre-filling

        if fill_value is None:
            # The fill value was not specified, so we'll set it to the default.
            fill_value = default_netcdf_fill_value(grid, significant_digits)
            if fill_value is None:
                raise ValueError(
                    "Grid type does not have a default fill value (according to netCDF4)"
                )

        # Set the fill value keyword argument.
        #
        # If this is None then netCDF4 will pre-fill using its default fill value (for the grid type).
        # If this is False then netCDF4 will not pre-fill. Note that True is not really a valid value (it behaves as the value 1).
        # Otherwise netCDF4 will pre-fill using the specified value.
        #
        # Note: It seems better to set the 'fill_value' keyword argument (when creating z variable)
        #       rather than set the 'missing_value' attribute (on the z variable after creating it).
        #       This translates to setting the '_FillValue' attribute instead of the 'missing_value' attribute.
        #       And this appears to work better when compressing/quantizing (eg, with 'significant_digits').
        data_kwds["fill_value"] = fill_value

        cdf_data = cdf.createVariable("z", grid.dtype, ("lat", "lon"), **data_kwds)

        # Ensure min and max z values are properly registered.
        if isinstance(fill_value, bool):
            # Fill value is False (note that True is not a valid value/type).
            # So all values are expected to be valid (note that NaN is a valid floating-point type)
            cdf_data.actual_range = [np.nanmin(grid), np.nanmax(grid)]
        elif np.isnan(fill_value):
            # Fill value is NaN.
            cdf_data.actual_range = [np.nanmin(grid), np.nanmax(grid)]
        else:
            # Fill value is a non-NaN number, so create a grid mask using it.
            grid_mask = grid != fill_value
            cdf_data.actual_range = [
                # Note: grid elements could still contain NaN values (though unlikely)...
                np.nanmin(grid[grid_mask]),
                np.nanmax(grid[grid_mask]),
            ]

        cdf_data.standard_name = "z"

        # cdf_data.add_offset = 0.0
        cdf_data.grid_mapping = "crs"
        # cdf_data.set_auto_maskandscale(False)

        #
        # NOTE: Create the lon/lat bounds variables AFTER creating the z variable so that 'gmt grdinfo' reports correctly.
        #
        #       Otherwise it reports, for example, (x_min, x_max) = (0, 1) and (y_min, y_max) = (179.9, -179.9) instead of
        #       (x_min, x_max) = (-180, 180) and (y_min, y_max) = (-90, 90) for a grid with 0.2 degree resolution.
        #       Apparently GMT scans variables based on the order they were written to the binary file structure, so if
        #       'z' is created before 'lon_bnds' then GMT processes 'z' first. And since 'lon_bnds' is a 2D variable
        #       with a column dimension of 2, we don't want GMT to get confused and try to read it as a grid plane.
        #
        cdf_lon_bnds = cdf.createVariable("lon_bnds", "f8", ("lon", "bnds"))
        cdf_lon_bnds[:, :] = _compute_coordinate_bounds(lon_grid)
        #
        cdf_lat_bnds = cdf.createVariable("lat_bnds", "f8", ("lat", "bnds"))
        cdf_lat_bnds[:, :] = _compute_coordinate_bounds(lat_grid)

        # write data
        cdf_data[:, :] = grid


def default_netcdf_fill_value(grid, significant_digits=None):
    """Return the default fill value that would be used when calling ``write_netcdf_grid`` with ``fill_value=None``.

    Notes
    -----
    This is useful when you need to set some values in ``grid`` to a fill value (eg, masking out continents from a seafloor age grid)
    before writing out the grid with ``write_netcdf_grid``.

    If ``significant_digits`` is NOT specified (ie, None) and the grid data type is floating-point then the default is `np.nan`.
    This essentially means that `np.nan` is only used (as the default) when *losslessly* compressing *floating-point* data.
    This is because *lossy* compression with a NaN fill value appears to *not* always mask out NaN regions.

    In all other cases the default is determined by `netCDF` based on the grid type.
    For example, the default for floating-point types is 9.969209968386869e+36 (see `netCDF4.default_fillvals`).
    In this case, please ensure the default is outside the range of your valid grid data, otherwise you should specify
    a custom fill value when calling ``write_netcdf_grid``.

    Parameters
    ----------
    grid : ndarray
        The grid that the default fill value will be used for (when writing the grid data to a `netCDF` file).

    significant_digits : int, optional
        Whether lossy data compression will be applied when writing the grid data to a `netCDF` file.
        This should be the same value that you will pass to ``write_netcdf_grid``.

    Returns
    -------
    The default fill value, or None if the grid type does not have a default value (according to netCDF4).
    """
    # If we're NOT using *lossy* compression and the grid type is floating-point
    # then set the fill value to np.nan.
    if significant_digits is None and np.issubdtype(grid.dtype, np.floating):
        return np.nan

    # When using *lossy* compression, we can't seem to use NaN as a fill value
    # because reading the resultant grid does not seem to mask out the NaN regions.
    # It was reported in https://github.com/GPlates/gplately/pull/125
    # that 2 significant digits was enough to preserve NaN masks, but
    # it doesn't seem to work for me (using netCDF4 1.7.2 on Windows).
    # Even 7 significant digits doesn't work.
    #
    # Instead we just use the default '_FillValue' provided by netCDF4.
    # These default values differ depending on the grid type, and can be accessed with 'netCDF4.default_fillvals'.
    # For example, the default for floating-point types is 9.969209968386869e+36.
    #
    # Note: This will return None if netCDF4 does not have a default fill value for the grid type.
    return netCDF4.default_fillvals.get(
        grid.dtype.str.lstrip("<>=|")  # Remove endianness/alignment prefixes
    )


@overload
def sample_grid(
    lon,
    lat,
    grid,
    method: str = "linear",
    extent: Union[tuple, str] = "global",
    origin=None,
    *,
    return_indices: Literal[False] = False,
) -> np.ndarray: ...
@overload
def sample_grid(
    lon,
    lat,
    grid,
    method: str = "linear",
    extent: Union[tuple, str] = "global",
    origin=None,
    *,
    return_indices: Literal[True],
) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray]]: ...
def sample_grid(
    lon,
    lat,
    grid,
    method="linear",
    extent: Union[tuple, str] = "global",
    origin=None,
    *,
    return_indices=False,
):
    """Sample point data with given `lon` and `lat` coordinates onto a `grid`
    using spline interpolation.

    Parameters
    ----------
    lon, lat : array_like
        The longitudes and latitudes of the points to interpolate onto the
        gridded data. Must be broadcastable to a common shape.
    grid : Raster or array_like
        An array whose elements define a grid. The number of rows corresponds
        to the number of point latitudes, while the number of columns
        corresponds to the number of point longitudes.
    method : str or int; default: 'linear'
        The order of spline interpolation. Must be an integer in the range
        0-5. 'nearest', 'linear', and 'cubic' are aliases for 0, 1, and 3,
        respectively.
    extent : str or 4-tuple, default: 'global'
        4-tuple to specify (min_lon, max_lon, min_lat, max_lat) extents
        of the raster. If no extents are supplied, full global extent
        [-180,180,-90,90] is assumed (equivalent to `extent='global'`).
        For array data with an upper-left origin, make sure `min_lat` is
        greater than `max_lat`, or specify `origin` parameter.
    origin : {'lower', 'upper'}, optional
        When `data` is an array, use this parameter to specify the origin
        (upper left or lower left) of the data (overriding `extent`).
    return_indices : bool, default=False
        Whether to return the row and column indices of the nearest grid
        points.

    Returns
    -------
    numpy.ndarray
        The values interpolated at the input points.
    indices : 2-tuple of numpy.ndarray
        The i- and j-indices of the nearest grid points to the input
        points, only present if `return_indices=True`.

    Raises
    ------
    ValueError
        If an invalid `method` is provided.
    RuntimeWarning
        If `lat` contains any invalid values outside of the interval
        [-90, 90]. Invalid values will be clipped to this interval.

    Notes
    -----
    If `return_indices` is set to `True`, the nearest array indices
    are returned as a tuple of arrays, in (i, j) or (lat, lon) format.

    An example output:

        # The first array holds the rows of the raster where point data spatially falls near.
        # The second array holds the columns of the raster where point data spatially falls near.
        sampled_indices = (array([1019, 1019, 1019, ..., 1086, 1086, 1087]), array([2237, 2237, 2237, ...,  983,  983,  983]))
    """
    order = {
        "nearest": 0,
        "linear": 1,
        "cubic": 3,
    }.get(method, method)
    if order not in {0, 1, 2, 3, 4, 5}:
        raise ValueError(f"Invalid `method` parameter: {method}")

    if isinstance(grid, Raster):
        extent = grid.extent
        if np.ma.isMaskedArray(grid.data):
            grid = np.ma.asarray(grid.data, dtype=float).filled(np.nan)
        else:
            grid = np.array(grid.data)
    else:
        extent = _parse_extent(extent, origin)
        if np.ma.isMaskedArray(grid):
            grid = np.ma.asarray(grid, dtype=float).filled(np.nan)
        grid = _check_grid(grid)

    # Do not wrap from North to South Pole (or vice versa)
    if np.any(np.abs(lat) > 90.0):
        warnings.warn(
            "Invalid values encountered in lat; clipping to [-90, 90]",
            RuntimeWarning,
        )
        lat = np.clip(lat, -90.0, 90.0)

    dx = (extent[1] - extent[0]) / (np.shape(grid)[1] - 1)
    dy = (extent[3] - extent[2]) / (np.shape(grid)[0] - 1)
    point_i = (lat - extent[2]) / dy
    point_j = (lon - extent[0]) / dx

    point_coords = np.vstack(
        (
            np.ravel(point_i),
            np.ravel(point_j),
        )
    )
    if np.ndim(grid) == 2:
        interpolated = map_coordinates(
            np.array(grid, dtype="float"),
            point_coords,
            order=order,
            mode="grid-wrap",
            prefilter=order > 1,
        )
        interpolated = np.reshape(interpolated, np.shape(lon))
    else:  # ndim(grid) == 3
        depth = np.shape(grid)[2]
        interpolated = []
        interpolated_k = np.array([])
        for k in range(depth):
            interpolated_k = map_coordinates(
                grid[..., k],
                point_coords,
                order=order,
                mode="grid-wrap",
                prefilter=order > 1,
            )
            interpolated_k = np.reshape(
                interpolated_k,
                np.shape(lon),
            )
            interpolated.append(interpolated_k)
        del interpolated_k
        interpolated = np.stack(interpolated, axis=-1)

    interpolated = interpolated.astype(grid.dtype)
    if return_indices:
        indices = (
            np.rint(np.ravel(point_i)).astype(np.int_),
            np.rint(np.ravel(point_j)).astype(np.int_),
        )
        return interpolated, indices
    return interpolated


def reconstruct_grid(
    grid,
    partitioning_features,
    rotation_model,
    to_time,
    from_time=0.0,
    extent: Union[tuple, str] = "global",
    origin=None,
    fill_value=None,
    threads=1,
    anchor_plate_id=None,
    x_dimension_name: str = "",
    y_dimension_name: str = "",
    data_variable_name: str = "",
):
    """Reconstruct a gridded dataset to a given reconstruction time.

    .. note::

        Use :meth:`Raster.reconstruct` whenever is possible. This :func:`reconstruct_grid` is better to be private.

    Parameters
    ----------
    grid : array_like, or str
        The grid to be reconstructed. If ``grid`` is a filename, it will be loaded using :meth:`read_netcdf_grid`.
    partitioning_features : valid argument to pygplates.FeaturesFunctionArgument
        Features used to partition the ``grid`` by plate ID, usually a static
        polygons file. The ``partitioning_features`` may be a single
        ``pygplates.Feature`` object, a ``pygplates.FeatureCollection``, a filename (:class:`str`), or a (potentially
        nested) sequence of any combination of the above types.
    rotation_model : valid argument to pygplates.RotationModel
        The rotation model used to reconstruct the ``grid``.
        The ``rotation_model`` may be a ``pygplates.RotationModel`` object, a rotation ``pygplates.FeatureCollection``, a rotation filename
        (:class:`str`), a rotation ``pygplates.Feature``, a sequence of
        rotation features, or a (potentially nested) sequence of any combination of the above types.
    to_time : float
        Time to which ``grid`` will be reconstructed.
    from_time : float, default=0.0
        Time from which to reconstruct the ``grid``.
    extent : tuple or str, default="global"
        Extent of the ``grid``. Valid arguments are a tuple of the form (xmin, xmax, ymin, ymax), or the string "global",
        equivalent to (-180.0, 180.0, -90.0, 90.0).
    origin : {"upper", "lower"}, optional
        Origin of the ``grid`` - either lower-left or upper-left. By default, determined from `extent`.
    fill_value : float, int, or tuple, optional, default=None
        The value to be used for regions outside of ``partitioning_features``
        at ``to_time``. If not provided, this value will be determined based on the input.
    threads : int, default=1
        Number of threads to use for certain computationally heavy routines.
    anchor_plate_id : int, optional, default=None
        ID of the anchored plate. By default, use the default anchor plate ID of ``rotation_model``
        if it's a ``pygplates.RotationModel`` (otherwise zero).
    x_dimension_name : str, optional, default=""
        If the grid file uses comman names, such as "x", "lon", "lons" or "longitude", you need not set this parameter.
        Otherwise, you need to tell us what the x dimension name is.
    y_dimension_name : str, optional, default=""
        If the grid file uses comman names, such as "y", "lat", "lats" or "latitude", you need not set this parameter.
        Otherwise, you need to tell us what the y dimension name is.
    data_variable_name : str, optional, default=""
        The program will try its best to determine the data variable name.
        However, it would be better if you could tell us what the data variable name is.
        Otherwise, the program will guess. The result may/may not be correct.

    Returns
    -------
    numpy.ndarray
        The reconstructed grid. Areas for which no plate ID could be
        determined from ``partitioning_features`` will be filled with ``fill_value``.


    .. note::

        For two-dimensional grids, ``fill_value`` should be a single
        number. The default value will be ``np.nan`` for float or
        complex types, the minimum value for integer types, and the
        maximum value for unsigned types.
        For RGB image grids, ``fill_value`` should be a 3-tuple RGB
        colour code or a matplotlib colour name. The default value
        will be black (0.0, 0.0, 0.0).
        For RGBA image grids, ``fill_value`` should be a 4-tuple RGBA
        colour code or a matplotlib colour name. The default fill
        value will be transparent black (0.0, 0.0, 0.0, 0.0).
    """
    if math.isclose(to_time, from_time):
        warnings.warn(
            "Reconstruction time is the same as the original time; returning input grid unchanged",
            UserWarning,
        )
        return grid

    assert rotation_model is not None, "`rotation_model` cannot be None."

    # first, try and see if the `grid` is a file path.
    if isinstance(grid, (str, bytes, os.PathLike)) and os.path.isfile(grid):
        grid = np.array(
            read_netcdf_grid(
                grid,
                x_dimension_name=x_dimension_name,
                y_dimension_name=y_dimension_name,
                data_variable_name=data_variable_name,
            )
        )
    else:
        # If the grid is not a file, we assume it is already an array-like object and proceed without loading.
        # convert grid data to numpy array. This will make a copy of the input grid.
        grid = np.array(grid)

    extent = _parse_extent(extent, origin)
    dtype = grid.dtype

    # Determine number of threads to use
    if isinstance(threads, str):
        if threads.lower() in {"all", "max"}:
            threads = cpu_count()
        else:
            raise ValueError(f"Invalid `threads` value: {threads}")
    threads = min([int(threads), cpu_count()])
    threads = max([threads, 1])

    grid = grid.squeeze()
    grid = _check_grid(grid)

    # Determine fill_value
    if fill_value is None:
        if grid.ndim == 2:
            if dtype.kind == "i":
                fill_value = np.iinfo(dtype).min
            elif dtype.kind == "u":
                fill_value = np.iinfo(dtype).max
            else:  # dtype.kind in ("f", "c")
                fill_value = np.nan
        else:  # grid.ndim == 3
            if dtype.kind in ("i", "u"):
                fill_value = tuple([0] * grid.shape[2])
            else:  # dtype.kind == "f"
                fill_value = tuple([0.0] * grid.shape[2])
    if isinstance(fill_value, str):
        if grid.ndim == 2:
            raise TypeError(f"Invalid fill_value for 2D grid: {fill_value}")
        fill_value = np.array(matplotlib.colors.to_rgba(fill_value))
        if dtype.kind == "u":
            fill_value = (fill_value * 255.0).astype("u1")
            fill_value = np.clip(fill_value, 0, 255)
        fill_value = tuple(fill_value)[: grid.shape[2]]

    if (
        grid.ndim == 3
        and grid.shape[2] == 4
        and hasattr(fill_value, "__len__")
        and len(fill_value) == 3  # type: ignore
    ):  # give fill colour maximum alpha value if not specified
        fill_alpha = 255 if dtype.kind in ("i", "u") else 1.0
        fill_value = (*fill_value, fill_alpha)  # type: ignore
    if np.size(fill_value) != np.atleast_3d(grid).shape[-1]:
        raise ValueError(
            f"Shape mismatch: fill_value size: {np.size(fill_value)}, grid shape: {np.shape(grid)}"
        )

    xmin, xmax, ymin, ymax = extent
    ny, nx = grid.shape[:2]

    if isinstance(partitioning_features, _FeaturesFunctionArgument):
        partitioning_features = _FeatureCollection(partitioning_features.get_features())
    elif not isinstance(partitioning_features, _FeatureCollection):
        partitioning_features = _FeatureCollection(
            _FeaturesFunctionArgument(partitioning_features).get_features()
        )

    if not isinstance(rotation_model, _RotationModel):
        rotation_model = _RotationModel(rotation_model)

    lons = np.linspace(xmin, xmax, nx)
    lats = np.linspace(ymin, ymax, ny)
    m_lons, m_lats = np.meshgrid(lons, lats)

    valid_partitioning_features = [
        i
        for i in partitioning_features
        if i.is_valid_at_time(from_time) and i.is_valid_at_time(to_time)
    ]
    plate_ids = rasterise(
        features=valid_partitioning_features,
        rotation_model=rotation_model,
        key="plate_id",
        time=from_time,
        extent=extent,
        shape=grid.shape[:2],
        origin=origin,
        anchor_plate_id=anchor_plate_id,
    )
    valid_output_mask = (
        rasterise(
            features=valid_partitioning_features,
            rotation_model=rotation_model,
            key="plate_id",
            time=to_time,
            extent=extent,
            shape=grid.shape[:2],
            origin=origin,
            anchor_plate_id=anchor_plate_id,
        )
        != -1
    )

    valid_mask = plate_ids != -1
    valid_m_lons = m_lons[valid_mask]
    valid_m_lats = m_lats[valid_mask]
    assert plate_ids is not None
    valid_plate_ids = plate_ids[valid_mask]
    if grid.ndim == 2:
        valid_data = grid[valid_mask]
    else:
        valid_data = np.empty(
            (grid.shape[2], np.sum(valid_mask)),
            dtype=dtype,
        )
        for k in range(grid.shape[2]):
            valid_data[k, :] = grid[..., k][valid_mask]

    if grid.ndim == 2:
        output_grid = np.full(grid.shape, fill_value)
    else:
        output_grid = np.empty(grid.shape, dtype=dtype)
        for k in range(grid.shape[2]):
            output_grid[..., k] = fill_value[k]  # type: ignore
    output_lons = m_lons[valid_output_mask]
    output_lats = m_lats[valid_output_mask]

    unique_plate_ids, inv = np.unique(valid_plate_ids, return_inverse=True)
    rotations_dict = {}
    for plate in unique_plate_ids:
        rot = rotation_model.get_rotation(
            to_time=float(to_time),
            from_time=float(from_time),
            moving_plate_id=int(plate),
            anchor_plate_id=anchor_plate_id,  # if None then uses default anchor plate of 'rotation_model'
        )
        if not isinstance(rot, _FiniteRotation):
            raise ValueError(f"No rotation found for plate ID: {plate}")
        lat, lon, angle = rot.get_lat_lon_euler_pole_and_angle_degrees()
        angle = np.deg2rad(angle)
        vec = _lat_lon_to_vector(lat, lon, degrees=True)
        rotations_dict[plate] = vec * angle
    rotations_array = np.array([rotations_dict[x] for x in unique_plate_ids])[inv]
    combined_rotations = _Rotation.from_rotvec(rotations_array)

    point_vecs = _lat_lon_to_vector(
        np.ravel(valid_m_lats),
        np.ravel(valid_m_lons),
        degrees=True,
    )
    rotated_vecs = combined_rotations.apply(point_vecs)

    tree = _cKDTree(rotated_vecs)
    output_vecs = _lat_lon_to_vector(
        output_lats,
        output_lons,
        degrees=True,
    )
    # Compatibility with older versions of SciPy:
    # 'n_jobs' argument was replaced with 'workers'
    try:
        _, indices = tree.query(output_vecs, k=1, workers=threads)
    except TypeError:
        _, indices = tree.query(output_vecs, k=1, n_jobs=threads)

    if grid.ndim == 2:
        output_data = valid_data[indices]
        output_grid[valid_output_mask] = output_data
    else:
        for k in range(grid.shape[2]):
            output_data = valid_data[k, indices]
            output_grid[..., k][valid_output_mask] = output_data

    return output_grid


def rasterise(
    features,
    rotation_model=None,
    key: Union[str, float, int, list] = "plate_id",
    time=None,
    resx=1.0,
    resy=1.0,
    shape=None,
    extent: Union[tuple, str] = "global",
    origin=None,
    tessellate_degrees=0.1,
    anchor_plate_id=None,
):
    """Rasterise geometries or GPlates features at a given reconstruction time.

    This function is particularly useful for rasterising static polygons
    to extract a grid of plate IDs.

    Parameters
    ----------
    features : geometries or features
        `features` may be a single `pygplates.Feature`, a
        `pygplates.FeatureCollection`, a `str` filename,
        or a (potentially nested) sequence of any combination of the
        above types.
        Alternatively, `features` may also be a sequence of geometry types
        (`pygplates.GeometryOnSphere` or `pygplates.ReconstructionGeometry`).
        In this case, `rotation_model` and `time` will be ignored, and
        `key` must be an array_like of the same length as `features`.
    rotation_model : valid argument for pygplates.RotationModel, optional
        `rotation_model` may be a `pygplates.RotationModel`, a rotation
        feature collection (pygplates.FeatureCollection), a rotation filename
        (`str`), a rotation feature (`pygplates.Feature`), a sequence of
        rotation features, or a (potentially nested) sequence of any
        combination of the above types.
        Alternatively, if time not given, a rotation model is
        not usually required.
    key : str or array_like, default "plate_id"
        The value used to create the rasterised grid. May be any of
        the following values:
        - "plate_id"
        - "conjugate_plate_id"
        - "from_age"
        - "to_age"
        - "left_plate"
        - "right_plate"
        Alternatively, `key` may be a sequence of the same length as
        `features`.
    time : float, optional
        Reconstruction time at which to perform rasterisation. If given,
        `rotation_model` must also be specified.
    resx, resy : float, default 1.0
        Resolution (in degrees) of the rasterised grid.
    shape : tuple, optional
        If given, the output grid will have the specified shape,
        overriding `resx` and `resy`.
    extent : tuple or "global", default "global"
        Extent of the rasterised grid. Valid arguments are a tuple of
        the form (xmin, xmax, ymin, ymax), or the string "global",
        equivalent to (-180.0, 180.0, -90.0, 90.0).
    origin : {"upper", "lower"}, optional
        Origin (upper-left or lower-left) of the output array. By default,
        determined from `extent`.
    tessellate_degrees : float, default 0.1
        Densify pyGPlates geometries to this resolution before conversion.
        Can be disabled by specifying `tessellate_degrees=None`, but this
        may provide inaccurate results for low-resolution input geometries.

    Returns
    -------
    grid : numpy.ndarray
        The output array will have the shape specified in `shape`, if given.
        The origin of the array will be in the lower-left corner of
        the area specified in `extent`, unless `resx` or `resy` is negative.

    Raises
    ------
    ValueError
        If an invalid `key` value is passed.
    TypeError
        If `rotation_model` is not supplied and `time` is not `None`.

    Notes
    -----
    This function is used by gplately.grids.reconstruct_grids to rasterise
    static polygons in order to extract their plate IDs.
    """
    valid_keys = {
        "plate_id",
        "conjugate_plate_id",
        "from_age",
        "to_age",
        "left_plate",
        "right_plate",
    }
    if isinstance(key, str):
        key = key.lower()
        if key not in valid_keys:
            raise ValueError(
                "Invalid key: {}".format(key)
                + "\nkey must be one of {}".format(valid_keys)
            )

    extent = _parse_extent(extent, origin)
    minx, maxx, miny, maxy = extent

    if minx > maxx:
        resx = -1.0 * np.abs(resx)
    if miny > maxy:
        resy = -1.0 * np.abs(resy)

    if shape is not None:
        lons = np.linspace(minx, maxx, shape[1], endpoint=True)
        lats = np.linspace(miny, maxy, shape[0], endpoint=True)
    else:
        lons = _spaced_axis(minx, maxx, resx)
        lats = _spaced_axis(miny, maxy, resy)
    nx = lons.size
    ny = lats.size

    try:
        features = _FeaturesFunctionArgument(features).get_features()
        geometries = None
    except Exception as err:
        if not str(err).startswith("Python argument types in"):
            # Not a Boost.Python.ArgumentError
            raise err
        geometries = pygplates_to_shapely(
            features,
            tessellate_degrees=tessellate_degrees,
        )

    reconstructed = []
    if geometries is None:
        if rotation_model is None:
            if time is not None:
                raise TypeError(
                    "Rotation model must be provided if `time` is not `None`"
                )
            rotation_model = _RotationModel(_Feature())
            time = 0.0
        features = _FeaturesFunctionArgument(features).get_features()
        if time is None:
            time = 0.0
        time = float(time)

        pygplates.reconstruct(  # type: ignore
            features,
            rotation_model,
            reconstructed,
            time,
            anchor_plate_id=anchor_plate_id,
        )
        geometries = pygplates_to_shapely(
            reconstructed,
            tessellate_degrees=tessellate_degrees,
        )
    if not isinstance(geometries, list):
        geometries = [geometries]

    if isinstance(key, str):
        values, fill_value, dtype = _get_rasterise_values(key, reconstructed)
    else:
        if isinstance(key, (int, float)):
            key = [key] * len(geometries)
        if len(key) != len(geometries):
            raise ValueError(
                f"Shape mismatch: len(key) = {len(key)}, len(geometries) = {len(geometries)}"
            )
        values = np.array(key)
        dtype = values.dtype
        if dtype.kind == "u":
            fill_value = np.iinfo(dtype).max
        elif dtype.kind == "i":
            fill_value = -1
        elif dtype.kind == "f":
            fill_value = np.nan
        else:
            raise TypeError("Unrecognised dtype for `key`: {}".format(dtype))

    return _rasterise_geometries(
        geometries=geometries,
        values=values,
        out_shape=(ny, nx),
        fill_value=fill_value,
        dtype=dtype,
        merge_alg=MergeAlg.replace,
        transform=_from_bounds(minx, miny, maxx, maxy, nx, ny),
    )


def _get_rasterise_values(
    key,
    reconstructed,
):
    valid_keys = {
        "plate_id",
        "conjugate_plate_id",
        "from_age",
        "to_age",
        "left_plate",
        "right_plate",
    }
    if key == "plate_id":
        values = [i.get_feature().get_reconstruction_plate_id() for i in reconstructed]
        fill_value = -1
        dtype = np.int32
    elif key == "conjugate_plate_id":
        values = [i.get_feature().get_conjugate_plate_id() for i in reconstructed]
        fill_value = -1
        dtype = np.int32
    elif key == "from_age":
        values = [i.get_feature().get_valid_time()[0] for i in reconstructed]
        fill_value = np.nan
        dtype = np.float32
    elif key == "to_age":
        values = [i.get_feature().get_valid_time()[1] for i in reconstructed]
        fill_value = np.nan
        dtype = np.float32
    elif key == "left_plate":
        values = [i.get_feature().get_left_plate() for i in reconstructed]
        fill_value = -1
        dtype = np.int32
    elif key == "right_plate":
        values = [i.get_feature().get_right_plate() for i in reconstructed]
        fill_value = -1
        dtype = np.int32
    else:
        raise ValueError(
            "Invalid key: {}".format(key) + "\nkey must be one of {}".format(valid_keys)
        )
    return values, fill_value, dtype


def _rasterise_geometries(
    geometries,
    values,
    out_shape,
    fill_value,
    dtype,
    transform,
    merge_alg=MergeAlg.replace,
):
    shapes = zip(geometries, values)
    out = _rasterize(
        shapes=shapes,
        out_shape=out_shape,
        fill=fill_value,
        dtype=dtype,
        merge_alg=merge_alg,
        transform=transform,
    )
    return np.flipud(out)


rasterize = rasterise  # alias for American English spelling


def _lat_lon_to_vector(lat, lon, degrees=False):
    """Convert (lat, lon) coordinates (degrees or radians) to vectors on
    the unit sphere. Returns a vector of shape (3,) if `lat` and `lon` are
    single values, else an array of shape (N, 3) containing N (x, y, z)
    row vectors, where N is the size of `lat` and `lon`.
    """
    lon = np.atleast_1d(lon).flatten()
    lat = np.atleast_1d(lat).flatten()
    if degrees:
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)

    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)

    size = x.size
    if size == 1:
        x = np.atleast_1d(np.squeeze(x))[0]
        y = np.atleast_1d(np.squeeze(y))[0]
        z = np.atleast_1d(np.squeeze(z))[0]
        return np.array((x, y, z))

    x = x.reshape((-1, 1))
    y = y.reshape((-1, 1))
    z = z.reshape((-1, 1))
    return np.hstack((x, y, z))


def _vector_to_lat_lon(
    x,
    y,
    z,
    degrees=False,
    return_array=False,
):
    """Convert one or more (x, y, z) vectors (on the unit sphere) to
    (lat, lon) coordinate pairs, in degrees or radians.
    """
    x = np.atleast_1d(x).flatten()
    y = np.atleast_1d(y).flatten()
    z = np.atleast_1d(z).flatten()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lat = np.arcsin(z)
        lon = np.arctan2(y, x)
        if degrees:
            lat = np.rad2deg(lat)
            lon = np.rad2deg(lon)

    if lat.size == 1 and not return_array:
        lat = np.atleast_1d(np.squeeze(lat))[0]
        lon = np.atleast_1d(np.squeeze(lon))[0]
        return (lat, lon)

    lat = lat.reshape((-1, 1))
    lon = lon.reshape((-1, 1))
    return lat, lon


def _check_grid_shape(data):
    """Check data is a 2D grid or a 3D RGB(A) image."""
    ndim = np.ndim(data)
    shape = np.shape(data)
    valid = True
    if ndim not in (2, 3):
        # ndim == 2: greyscale image/grid
        # ndim == 3: colour RGB(A) image
        valid = False
    if ndim == 3 and shape[2] not in (3, 4):
        # shape[2] == 3: colour image (RGB)
        # shape[2] == 4: colour image w/ transparency (RGBA)
        valid = False

    if not valid:
        raise ValueError("Invalid grid shape: {}".format(shape))


def _check_image_values(data):
    """Check values are within correct range for an RGB(A) image."""
    dtype = data.dtype
    if dtype.kind == "i":
        data = data.astype("u1")
        dtype = data.dtype
    min_value = np.nanmin(data)
    max_value = np.nanmax(data)
    if min_value < 0:
        raise ValueError("Invalid value for RGB(A) image: {}".format(min_value))
    if (dtype.kind == "f" and max_value > 1.0) or (
        dtype.kind == "u" and max_value > 255
    ):
        raise ValueError("Invalid value for RGB(A) image: {}".format(max_value))
    return data


def _check_grid(data):
    """Check grid shape and values make sense."""
    if not isinstance(data, np.ndarray):
        data = np.array(data)

    _check_grid_shape(data)
    if data.ndim == 3:
        # data is an RGB(A) image
        data = _check_image_values(data)
    return data


def _parse_extent(extent, origin) -> Tuple[float, float, float, float]:
    """Default values: extent='global', origin=None"""
    if hasattr(extent, "lower"):  # i.e. a string
        extent = extent.lower()

    if extent is None or extent == "global":
        extent = (-180.0, 180.0, -90.0, 90.0)
    elif len(extent) != 4:
        raise TypeError("`extent` must be a four-element tuple, 'global', or None")

    extent = tuple(float(i) for i in extent)
    return cast(
        Tuple[float, float, float, float], _adjust_extent_for_origin(extent, origin)
    )


def _adjust_extent_for_origin(
    extent, origin
) -> Union[Tuple[float, float, float, float], None]:
    """Adjust upper/lower bounds of extent according to origin."""
    if extent is None:
        return None

    if origin is None:
        return extent

    origin = str(origin).lower()
    if origin == "lower" and extent[2] > extent[3]:
        extent = (
            extent[0],
            extent[1],
            extent[3],
            extent[2],
        )
    if origin == "upper" and extent[2] < extent[3]:
        extent = (
            extent[0],
            extent[1],
            extent[3],
            extent[2],
        )

    return extent


# class TimeRaster(Raster):
#     """A class for the temporal manipulation of raster data. To be added soon!"""
#     def __init__(self, plate_reconstruction=None, filename=None, array=None, extent=None, resample=None):
#         raise NotImplementedError(
#             "This class has not been implemented; use `Raster` instead"
#         )
#         super(TimeRaster, self).__init__(plate_reconstruction)
