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

"""
This sub-module contains tools for working with MaskedArray, ndarray and netCDF4 rasters, as well as
gridded-data.

Some methods available in `grids`:

* Point data can be interpolated onto a raster or grid with Scipy using linear or
nearest-neighbour interpolation.
* Rasters can be resampled with a set of X and Y-direction spacings, and can be resized
using given X and Y resolutions.
* Grids with invalid (NaN-type) data cells can have their NaN entries replaced
with the values of their nearest valid neighbours.

Classes
-------
* RegularGridInterpolator
* Raster
"""

import copy
import logging
import math
import warnings
from multiprocessing import cpu_count
from typing import Tuple, Union

import matplotlib.colors
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pygplates
from cartopy.crs import PlateCarree as _PlateCarree
from cartopy.mpl.geoaxes import GeoAxes as _GeoAxes
from rasterio.enums import MergeAlg
from rasterio.features import rasterize as _rasterize
from rasterio.transform import from_bounds as _from_bounds
from scipy.interpolate import RegularGridInterpolator as _RGI
from scipy.interpolate import griddata
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy.spatial import cKDTree as _cKDTree  # type: ignore
from scipy.spatial.transform import Rotation as _Rotation

from .geometry import pygplates_to_shapely
from .reconstruction import PlateReconstruction as _PlateReconstruction
from .tools import _deg2pixels, griddata_sphere

logger = logging.getLogger("gplately")

__all__ = [
    "fill_raster",
    "read_netcdf_grid",
    "write_netcdf_grid",
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
    return data[tuple(ind)]


def _realign_grid(array, lons, lats):
    """realigns grid to -180/180 and flips the array if the latitudinal coordinates are decreasing."""
    mask_lons = lons > 180

    # realign to -180/180
    if mask_lons.any():
        dlon = np.diff(lons).mean()
        array = np.hstack([array[:, mask_lons], array[:, ~mask_lons]])
        lons = np.hstack([lons[mask_lons] - 360 - dlon, lons[~mask_lons]])

    if lats[0] > lats[-1]:
        array = np.flipud(array)
        lats = lats[::-1]

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


def read_netcdf_grid(
    filename,
    return_grids: bool = False,
    realign: bool = False,
    resample=None,
    resize=None,
    x_dimension_name: str = "",
    y_dimension_name: str = "",
    data_variable_name: str = "",
) -> Union[
    Tuple[np.ma.MaskedArray, np.ma.MaskedArray, np.ma.MaskedArray], np.ma.MaskedArray
]:
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
        if hasattr(cdf[key_z], "missing_value") and np.issubdtype(
            cdf_grid.dtype, np.floating
        ):
            fill_value = cdf[key_z].missing_value
            cdf_grid[np.isclose(cdf_grid, fill_value, rtol=0.1)] = np.nan

        # convert to boolean array
        if np.issubdtype(cdf_grid.dtype, np.integer):
            unique_grid = np.unique(cdf_grid)
            if len(unique_grid) == 2:
                if (unique_grid == [0, 1]).all():
                    cdf_grid = cdf_grid.astype(bool)

    if realign:
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

        if spacingX != dX or spacingY != dY:
            lon_grid = np.arange(cdf_lon.min(), cdf_lon.max() + spacingX, spacingX)
            lat_grid = np.arange(cdf_lat.min(), cdf_lat.max() + spacingY, spacingY)
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
        if resX != cdf_grid_z.shape[1] or resY != cdf_grid_z.shape[0]:
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
    fill_value: Union[float, None] = np.nan,
):
    """Write geological data contained in a `grid` to a netCDF4 grid with a specified `filename`.

    Notes
    -----
    The written netCDF4 grid has the same latitudinal and longitudinal (row and column) dimensions as `grid`.
    It has three variables:

    * Latitudes of `grid` data
    * Longitudes of `grid` data
    * The data stored in `grid`

    However, the latitudes and longitudes of the grid returned to the user are constrained to those
    specified in `extent`.
    By default, `extent` assumes a global latitudinal and longitudinal span: `extent=[-180,180,-90,90]`.

    Parameters
    ----------
    filename : str
        The full path (including a filename and the ".nc" extension) to save the created netCDF4 `grid` to.

    grid : array-like
        An ndarray grid containing data to be written into a `netCDF` (.nc) file. Note: Rows correspond to
        the data's latitudes, while the columns correspond to the data's longitudes.

    extent : list, default=[-180,180,-90,90]
        Four elements that specify the [min lon, max lon, min lat, max lat] to constrain the lat and lon
        variables of the netCDF grid to. If no extents are supplied, full global extent `[-180, 180, -90, 90]`
        is assumed.

    significant_digits : int
        Applies lossy data compression up to a specified number of significant digits.
        This significantly reduces file size, but make sure the required precision is preserved in the
        saved netcdf file.

    fill_value : scalar, NoneType, default: np.nan
        Value used to fill in missing data. By default this is np.nan.

    Returns
    -------
    A netCDF grid will be saved to the path specified in `filename`.
    """
    import netCDF4

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

    with netCDF4.Dataset(filename, "w", driver=None) as cdf:
        cdf.title = "Grid produced by gplately " + str(_version)
        cdf.createDimension("lon", lon_grid.size)
        cdf.createDimension("lat", lat_grid.size)
        cdf_lon = cdf.createVariable("lon", lon_grid.dtype, ("lon",), **data_kwds)
        cdf_lat = cdf.createVariable("lat", lat_grid.dtype, ("lat",), **data_kwds)
        cdf_lon[:] = lon_grid
        cdf_lat[:] = lat_grid

        # Units for Geographic Grid type
        cdf_lon.units = "degrees_east"
        cdf_lon.standard_name = "lon"
        cdf_lon.actual_range = [lon_grid[0], lon_grid[-1]]
        cdf_lat.units = "degrees_north"
        cdf_lat.standard_name = "lat"
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
            # significant_digits needs to be >= 2 so that NaNs are preserved
            data_kwds["significant_digits"] = max(2, int(significant_digits))
            data_kwds["quantize_mode"] = "GranularBitRound"

        # boolean arrays need to be converted to integers
        # no such thing as a mask on a boolean array
        if grid.dtype is np.dtype(bool):
            grid = grid.astype("i1")
            fill_value = None

        cdf_data = cdf.createVariable("z", grid.dtype, ("lat", "lon"), **data_kwds)

        # netCDF4 uses the missing_value attribute as the default _FillValue
        # without this, _FillValue defaults to 9.969209968386869e+36
        if fill_value is not None:
            cdf_data.missing_value = fill_value
            grid_mask = grid != fill_value

            cdf_data.actual_range = [
                np.nanmin(grid[grid_mask]),
                np.nanmax(grid[grid_mask]),
            ]

        else:
            # ensure min and max z values are properly registered
            cdf_data.actual_range = [np.nanmin(grid), np.nanmax(grid)]

        cdf_data.standard_name = "z"

        # cdf_data.add_offset = 0.0
        cdf_data.grid_mapping = "crs"
        # cdf_data.set_auto_maskandscale(False)

        # write data
        cdf_data[:, :] = grid


class RegularGridInterpolator(_RGI):
    """A class to sample gridded data at a set of point coordinates using either linear or nearest-neighbour
    interpolation methods. It is a child class of `scipy 1.10`'s [`RegularGridInterpolator`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html) class.

    This will only work for scipy version 1.10 onwards.

    Attributes
    ----------
    points : tuple of ndarrays of float with shapes (m1, ), …, (mn, )
        Each array contains point coordinates that define the regular grid in n dimensions.
    values : ndarray
        The data on a regular grid. Note: the number of rows corresponds to the number of point latitudes, while the number
        of columns corresponds to the number of point longitudes.
    method : str, default=’linear’
        The method of interpolation to perform. Supported are "linear" and "nearest". Assumes “linear” by default.
    bounds_error : bool, default=false
        Choose whether to return a ValueError and terminate the interpolation if any provided sample points are out
        of grid bounds. By default, it is set to `False`. In this case, all out-of-bound point values are replaced
        with the `fill_value` (defined below) if supplied.
    fill_value : float, default=np.nan
        Used to replace point values that are out of grid bounds, provided that ‘bounds_error’ is false.

    """

    def __init__(
        self, points, values, method="linear", bounds_error=False, fill_value=np.nan
    ):
        super(RegularGridInterpolator, self).__init__(
            points, values, method, bounds_error, fill_value
        )

    def __call__(self, xi, method=None, return_indices=False, return_distances=False):
        """Samples gridded data at a set of point coordinates. Uses either a linear or nearest-neighbour interpolation `method`.

        Uses the gridded data specified in the sample_grid method parameter. Note: if any provided sample points are out of
        grid bounds and a corresponding error message was suppressed (by specifying bounds_error=False), all out-of-bound
        point values are replaced with the self.fill_value attribute ascribed to the RegularGridInterpolator object (if it
        exists). Terminates otherwise.

        This is identical to scipy 1.10's RGI object.

        Parameters
        ----------
        xi : ndarray of shape (..., ndim)
            The coordinates of points to sample the gridded data at.

        method : str, default=None
            The method of interpolation to perform. Supported are "linear" and "Nearest". Assumes “linear” interpolation
            if None provided.

        return_indices : bool, default=False
            Choose whether to return indices of neighbouring sampling points.

        return_distances : bool, default=False
            Choose whether to return normal distances between interpolated points and neighbouring sampling points.

        Returns
        -------
        output_tuple : tuple of ndarrays
            The first ndarray in the output tuple holds the interpolated grid data. If sample point distances and indices are
            required, these are returned as subsequent tuple elements.

        Raises
        ------
        ValueError
            * Raised if the string method supplied is not “linear” or “nearest”.
            * Raised if the provided sample points for interpolation (xi) do not have the same dimensions as the supplied grid.
            * Raised if the provided sample points for interpolation include any point out of grid bounds. Alerts user which
            dimension (index) the point is located. Only raised if the RegularGridInterpolator attribute bounds_error is set
            to True. If suppressed, out-of-bound points are replaced with a set fill_value.
        """
        method = self.method if method is None else method
        if method not in ["linear", "nearest"]:
            raise ValueError("Method '%s' is not defined" % method)

        xi, xi_shape, ndim, nans, out_of_bounds = self._prepare_xi(xi)

        indices, norm_distances = self._find_indices(xi.T)

        if method == "linear":
            result = self._evaluate_linear(indices, norm_distances)
        elif method == "nearest":
            result = self._evaluate_nearest(indices, norm_distances)
        if not self.bounds_error and self.fill_value is not None:
            result[out_of_bounds] = self.fill_value

        interp_output = result.reshape(xi_shape[:-1] + self.values.shape[ndim:])
        output_tuple = [interp_output]

        if return_indices:
            output_tuple.append(indices)
        if return_distances:
            output_tuple.append(norm_distances)

        if return_distances or return_indices:
            return tuple(output_tuple)
        else:
            return output_tuple[0]

    def _prepare_xi(self, xi):
        try:
            from scipy.interpolate.interpnd import _ndim_coords_from_arrays
        except ImportError:
            # SciPy 1.15 renamed interpnd to _interpnd (see https://github.com/scipy/scipy/pull/21754).
            from scipy.interpolate._interpnd import _ndim_coords_from_arrays

        ndim = len(self.grid)
        xi = _ndim_coords_from_arrays(xi, ndim=ndim)
        if xi.shape[-1] != len(self.grid):
            raise ValueError(
                "The requested sample points xi have dimension "
                f"{xi.shape[-1]} but this "
                f"RegularGridInterpolator has dimension {ndim}"
            )

        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])

        # find nans in input
        nans = np.any(np.isnan(xi), axis=-1)

        if self.bounds_error:
            for i, p in enumerate(xi.T):
                if not np.logical_and(
                    np.all(self.grid[i][0] <= p), np.all(p <= self.grid[i][-1])
                ):
                    raise ValueError(
                        "One of the requested xi is out of bounds "
                        "in dimension %d" % i
                    )
            out_of_bounds = None
        else:
            out_of_bounds = self._find_out_of_bounds(xi.T)

        return xi, xi_shape, ndim, nans, out_of_bounds

    def _find_out_of_bounds(self, xi):
        # check for out of bounds xi
        out_of_bounds = np.zeros((xi.shape[1]), dtype=bool)
        # iterate through dimensions
        for x, grid in zip(xi, self.grid):
            out_of_bounds += x < grid[0]
            out_of_bounds += x > grid[-1]
        return out_of_bounds

    def _find_indices(self, xi):
        """Index identifier outsourced from scipy 1.9's
        RegularGridInterpolator to ensure stable
        operations with all versions of scipy >1.0.
        """
        # find relevant edges between which xi are situated
        indices = []
        # compute distance to lower edge in unity units
        norm_distances = []
        # iterate through dimensions
        for x, grid in zip(xi, self.grid):
            i = np.searchsorted(grid, x) - 1
            i[i < 0] = 0
            i[i > grid.size - 2] = grid.size - 2
            indices.append(i)

            # compute norm_distances, incl length-1 grids,
            # where `grid[i+1] == grid[i]`
            denom = grid[i + 1] - grid[i]
            with np.errstate(divide="ignore", invalid="ignore"):
                norm_dist = np.where(denom != 0, (x - grid[i]) / denom, 0)
            norm_distances.append(norm_dist)

        return indices, norm_distances

    def _evaluate_linear(self, indices, norm_distances):
        """Linear interpolator outsourced from scipy 1.9's
        RegularGridInterpolator to ensure stable
        operations with all versions of scipy >1.0.
        """
        import itertools

        # slice for broadcasting over trailing dimensions in self.values
        vslice = (slice(None),) + (None,) * (self.values.ndim - len(indices))

        # Compute shifting up front before zipping everything together
        shift_norm_distances = [1 - yi for yi in norm_distances]
        shift_indices = [i + 1 for i in indices]

        # The formula for linear interpolation in 2d takes the form:
        # values = self.values[(i0, i1)] * (1 - y0) * (1 - y1) + \
        #          self.values[(i0, i1 + 1)] * (1 - y0) * y1 + \
        #          self.values[(i0 + 1, i1)] * y0 * (1 - y1) + \
        #          self.values[(i0 + 1, i1 + 1)] * y0 * y1
        # We pair i with 1 - yi (zipped1) and i + 1 with yi (zipped2)
        zipped1 = zip(indices, shift_norm_distances)
        zipped2 = zip(shift_indices, norm_distances)

        # Take all products of zipped1 and zipped2 and iterate over them
        # to get the terms in the above formula. This corresponds to iterating
        # over the vertices of a hypercube.
        hypercube = itertools.product(*zip(zipped1, zipped2))
        values = 0.0
        for h in hypercube:
            edge_indices, weights = zip(*h)
            weight = 1.0
            for w in weights:
                weight *= w
            values += np.asarray(self.values[edge_indices]) * weight[vslice]
        return values

    def _evaluate_nearest(self, indices, norm_distances):
        """Nearest neighbour interpolator outsourced from scipy 1.9's
        RegularGridInterpolator to ensure stable
        operations with all versions of scipy >1.0.
        """
        idx_res = [
            np.where(yi <= 0.5, i, i + 1) for i, yi in zip(indices, norm_distances)
        ]
        return self.values[tuple(idx_res)]


def sample_grid(
    lon,
    lat,
    grid,
    method="linear",
    extent: Union[tuple, str] = "global",
    origin=None,
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
        raise ValueError("Invalid `method` parameter: {}".format(method))

    if isinstance(grid, Raster):
        extent = grid.extent
        grid = np.array(grid.data)
    else:
        extent = _parse_extent_origin(extent, origin)
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

    point_coords = np.row_stack(
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
    try:
        grid = np.array(
            read_netcdf_grid(
                grid,
                x_dimension_name=x_dimension_name,
                y_dimension_name=y_dimension_name,
                data_variable_name=data_variable_name,
            )
        )  # load grid data from file
    except Exception:
        grid = np.array(grid)  # copy grid data to array
    if to_time == from_time:
        return grid
    elif rotation_model is None:
        raise TypeError("`rotation_model` must be provided if `to_time` != `from_time`")

    extent = _parse_extent_origin(extent, origin)
    dtype = grid.dtype

    if isinstance(threads, str):
        if threads.lower() in {"all", "max"}:
            threads = cpu_count()
        else:
            raise ValueError("Invalid `threads` value: {}".format(threads))
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
            raise TypeError("Invalid fill_value for 2D grid: {}".format(fill_value))
        fill_value = np.array(matplotlib.colors.to_rgba(fill_value))
        if dtype.kind == "u":
            fill_value = (fill_value * 255.0).astype("u1")
            fill_value = np.clip(fill_value, 0, 255)
        fill_value = tuple(fill_value)[: grid.shape[2]]

    if (
        grid.ndim == 3
        and grid.shape[2] == 4
        and hasattr(fill_value, "__len__")
        and len(fill_value) == 3
    ):  # give fill colour maximum alpha value if not specified
        fill_alpha = 255 if dtype.kind in ("i", "u") else 1.0
        fill_value = (*fill_value, fill_alpha)
    if np.size(fill_value) != np.atleast_3d(grid).shape[-1]:
        raise ValueError(
            "Shape mismatch: "
            + "fill_value size: {}".format(np.size(fill_value))
            + ", grid shape: {}".format(np.shape(grid))
        )

    xmin, xmax, ymin, ymax = extent
    ny, nx = grid.shape[:2]

    if isinstance(partitioning_features, pygplates.FeaturesFunctionArgument):
        partitioning_features = pygplates.FeatureCollection(
            partitioning_features.get_features()
        )
    elif not isinstance(partitioning_features, pygplates.FeatureCollection):
        partitioning_features = pygplates.FeatureCollection(
            pygplates.FeaturesFunctionArgument(partitioning_features).get_features()
        )

    if not isinstance(rotation_model, pygplates.RotationModel):
        rotation_model = pygplates.RotationModel(rotation_model)

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
            output_grid[..., k] = fill_value[k]
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
        if not isinstance(rot, pygplates.FiniteRotation):
            raise ValueError("No rotation found for plate ID: {}".format(plate))
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
        _, indices = tree.query(
            output_vecs,
            k=1,
            workers=threads,
        )
    except TypeError as err:
        if "Unexpected keyword argument" in err.args[0] and "workers" in err.args[0]:
            _, indices = tree.query(
                output_vecs,
                k=1,
                n_jobs=threads,
            )
        else:
            raise err

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
    key="plate_id",
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

    extent = _parse_extent_origin(extent, origin)
    minx, maxx, miny, maxy = extent

    if minx > maxx:
        resx = -1.0 * np.abs(resx)
    if miny > maxy:
        resy = -1.0 * np.abs(resy)

    if shape is not None:
        lons = np.linspace(minx, maxx, shape[1], endpoint=True)
        lats = np.linspace(miny, maxy, shape[0], endpoint=True)
    else:
        lons = np.arange(minx, maxx + resx, resx)
        lats = np.arange(miny, maxy + resy, resy)
    nx = lons.size
    ny = lats.size

    try:
        features = pygplates.FeaturesFunctionArgument(features).get_features()
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
            rotation_model = pygplates.RotationModel(pygplates.Feature())
            time = 0.0
        features = pygplates.FeaturesFunctionArgument(features).get_features()
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
        if not hasattr(key, "__len__"):
            key = [key] * len(geometries)
        if len(key) != len(geometries):
            raise ValueError(
                "Shape mismatch: len(key) = {}, ".format(len(key))
                + "len(geometries) = {}".format(len(geometries))
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


rasterize = rasterise


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
    ndim = data.ndim
    dtype = data.dtype
    _check_grid_shape(data)

    if ndim == 3:
        # data is an RGB(A) image
        data = _check_image_values(data)

    return data


def _parse_extent_origin(extent, origin):
    """Default values: extent='global', origin=None"""
    if hasattr(extent, "lower"):  # i.e. a string
        extent = extent.lower()

    if extent is None or extent == "global":
        extent = (-180.0, 180.0, -90.0, 90.0)
    elif len(extent) != 4:
        raise TypeError("`extent` must be a four-element tuple, 'global', or None")
    extent = tuple(float(i) for i in extent)

    if origin is not None:
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


class Raster(object):
    """The functionalities include sampling data at points using spline
    interpolation, resampling rasters with new X and Y-direction spacings and
    resizing rasters using new X and Y grid pixel resolutions. NaN-type data
    in rasters can be replaced with the values of their nearest valid neighbours.
    """

    def __init__(
        self,
        data=None,
        plate_reconstruction=None,
        extent: Union[str, tuple] = "global",
        realign=False,
        resample=None,
        resize=None,
        time=0.0,
        origin=None,
        x_dimension_name: str = "",
        y_dimension_name: str = "",
        data_variable_name: str = "",
        **kwargs,
    ):
        """Constructor. Create a :class:`Raster` object.

        Parameters
        ----------
        data : str or array-like
            The raster data, either as a file path (:class:`str`) or array data.

        plate_reconstruction : PlateReconstruction
            A :class:`PlateReconstruction` object to provide the following essential components for reconstructing points.

            * :py:attr:`PlateReconstruction.rotation_model`
            * :py:attr:`PlateReconstruction.topology_featues`
            * :py:attr:`PlateReconstruction.static_polygons`

        extent : str or 4-tuple, default: 'global'
            4-tuple to specify (min_lon, max_lon, min_lat, max_lat) extents
            of the raster. If no extents are supplied, full global extent
            (-180, 180, -90, 90) is assumed (equivalent to ``extent='global'``).
            For array data with an upper-left origin, make sure ``min_lat`` is
            greater than ``max_lat``, or specify ``origin`` parameter.

        resample : 2-tuple, optional
            Optionally resample grid, pass spacing in X and Y direction as a
            2-tuple e.g. resample=(spacingX, spacingY).

        resize : 2-tuple, optional
            Optionally resample grid to X-columns, Y-rows as a
            2-tuple e.g. resample=(resX, resY).

        time : float, default: 0.0
            The geological time the time-dependant raster data.

        origin : {'lower', 'upper'}, optional
            When ``data`` is an array, use this parameter to specify the origin
            (upper left or lower left) of the data (overriding ``extent``).

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

        **kwargs
            Handle deprecated arguments such as ``PlateReconstruction_object``, ``filename``, and ``array``.
        """
        if isinstance(data, self.__class__):
            self._data = data._data.copy()
            self.plate_reconstruction = data.plate_reconstruction
            self._lons = data._lons
            self._lats = data._lats
            self._time = data._time
            return

        if "PlateReconstruction_object" in kwargs.keys():
            warnings.warn(
                "`PlateReconstruction_object` keyword argument has been "
                + "deprecated, use `plate_reconstruction` instead",
                DeprecationWarning,
            )
            if plate_reconstruction is None:
                plate_reconstruction = kwargs.pop("PlateReconstruction_object")
        if "filename" in kwargs.keys() and "array" in kwargs.keys():
            raise TypeError(
                "Both `filename` and `array` were provided; use "
                + "one or the other, or use the `data` argument"
            )
        if "filename" in kwargs.keys():
            warnings.warn(
                "`filename` keyword argument has been deprecated, "
                + "use `data` instead",
                DeprecationWarning,
            )
            if data is None:
                data = kwargs.pop("filename")
        if "array" in kwargs.keys():
            warnings.warn(
                "`array` keyword argument has been deprecated, " + "use `data` instead",
                DeprecationWarning,
            )
            if data is None:
                data = kwargs.pop("array")
        for key in kwargs.keys():
            raise TypeError(
                "Raster.__init__() got an unexpected keyword argument "
                + "'{}'".format(key)
            )
        self.plate_reconstruction = plate_reconstruction

        if time < 0.0:
            raise ValueError("Invalid time: {}".format(time))
        time = float(time)
        self._time = time

        if data is None:
            raise TypeError("`data` argument (or `filename` or `array`) is required")
        if isinstance(data, str):
            # Filename
            self._filename = data
            self._data, lons, lats = read_netcdf_grid(
                data,
                return_grids=True,
                realign=realign,
                resample=resample,
                resize=resize,
                x_dimension_name=x_dimension_name,
                y_dimension_name=y_dimension_name,
                data_variable_name=data_variable_name,
            )
            self._lons = lons
            self._lats = lats

        else:
            # numpy array
            self._filename = None
            extent = _parse_extent_origin(extent, origin)
            data = _check_grid(data)
            self._data = np.array(data)
            self._lons = np.linspace(extent[0], extent[1], self.data.shape[1])
            self._lats = np.linspace(extent[2], extent[3], self.data.shape[0])
            if realign:
                # realign to -180,180 and flip grid
                self._data, self._lons, self._lats = _realign_grid(
                    self._data, self._lons, self._lats
                )

        if (not isinstance(data, str)) and (resample is not None):
            self.resample(*resample, inplace=True)

        if (not isinstance(data, str)) and (resize is not None):
            self.resize(*resize, inplace=True)

    @property
    def time(self):
        """The geological time of the time-dependant raster data.

        :type: float
        """
        return self._time

    @time.setter
    def time(self, new_time: float):
        """Set a new reconstruction time."""
        try:
            new_time_f = float(new_time)
        except ValueError:
            raise ValueError(f"Invalid new reconstruction time: {new_time}")
        if new_time_f < 0.0:
            raise ValueError(
                f"The reconstruction time ({new_time_f}) must be greater than 0."
            )
        if not math.isclose(self._time, new_time_f):
            self._time = new_time_f
            self.reconstruct(new_time_f, inplace=True)

    @property
    def data(self):
        """Array containing the raster data. This attribute can be modified after creating the :class:`Raster` object.

        :type: ndarray, shape (ny, nx)
        """
        return self._data

    @data.setter
    def data(self, z):
        z = np.array(z)
        if z.shape != np.shape(self.data):
            raise ValueError(
                "Shape mismatch: old dimensions are {}, new are {}".format(
                    np.shape(self.data),
                    z.shape,
                )
            )
        self._data = z

    @property
    def lons(self):
        """The x-coordinates of the raster data. This attribute can be modified after creating the :class:`Raster` object.

        :type: ndarray, shape (nx,)
        """
        return self._lons

    @lons.setter
    def lons(self, x):
        x = np.array(x).ravel()
        if x.size != np.shape(self.data)[1]:
            raise ValueError(
                "Shape mismatch: data x-dimension is {}, new value is {}".format(
                    np.shape(self.data)[1],
                    x.size,
                )
            )
        self._lons = x

    @property
    def lats(self):
        """The y-coordinates of the raster data. This attribute can be modified after creating the :class:`Raster` object.

        :type: ndarray, shape (ny,)
        """
        return self._lats

    @lats.setter
    def lats(self, y):
        y = np.array(y).ravel()
        if y.size != np.shape(self.data)[0]:
            raise ValueError(
                "Shape mismatch: data y-dimension is {}, new value is {}".format(
                    np.shape(self.data)[0],
                    y.size,
                )
            )
        self._lats = y

    @property
    def extent(self):
        """The spatial extent ``(x0, x1, y0, y1)`` of the data. If not supplied, global extent ``(-180, 180, -90, 90)`` is assumed.

        If y0 < y1, the origin is the lower-left corner; else the upper-left.

        :type:  tuple of 4 floats
        """
        return (
            float(self.lons[0]),
            float(self.lons[-1]),
            float(self.lats[0]),
            float(self.lats[-1]),
        )

    @property
    def origin(self):
        """The origin (``lower`` or ``upper``) of the data array.

        :type: str
        """
        if self.lats[0] < self.lats[-1]:
            return "lower"
        else:
            return "upper"

    @property
    def shape(self):
        """The shape of the data array."""
        return np.shape(self.data)

    @property
    def size(self):
        """The size of the data array."""
        return np.size(self.data)

    @property
    def dtype(self):
        """The data type of the array."""
        return self.data.dtype

    @property
    def ndim(self):
        """The number of dimensions in the array."""
        return np.ndim(self.data)

    @property
    def filename(self):
        """The filename used to create the :class:`Raster` object.
        If the object was created directly from an array, this attribute is ``None``.

        :type: str or None
        """
        return self._filename

    @property
    def plate_reconstruction(self):
        """A :class:`PlateReconstruction` object to provide the following essential components for reconstructing points.

            * :py:attr:`PlateReconstruction.rotation_model`
            * :py:attr:`PlateReconstruction.topology_featues`
            * :py:attr:`PlateReconstruction.static_polygons`

        :type: PlateReconstruction
        """
        return self._plate_reconstruction

    @plate_reconstruction.setter
    def plate_reconstruction(self, reconstruction):
        if reconstruction is None:
            # Remove `plate_reconstruction` attribute
            pass
        elif not isinstance(reconstruction, _PlateReconstruction):
            # Convert to a `PlateReconstruction` if possible
            try:
                reconstruction = _PlateReconstruction(*reconstruction)
            except Exception:
                reconstruction = _PlateReconstruction(reconstruction)
        self._plate_reconstruction = reconstruction

    def copy(self):
        """Return a copy of the :class:`Raster` object.

        Returns
        -------
        Raster
            A copy of the current :class:`Raster` object.
        """
        return Raster(
            self.data.copy(), self.plate_reconstruction, self.extent, time=self.time
        )

    def interpolate(
        self,
        lons,
        lats,
        method="linear",
        return_indices=False,
    ):
        """Sample grid data at a set of points using spline interpolation.

        Parameters
        ----------
        lons, lats : array_like
            The longitudes and latitudes of the points to interpolate onto the
            gridded data. Must be broadcastable to a common shape.
        method : str or int; default: 'linear'
            The order of spline interpolation. Must be an integer in the range
            0-5. ``nearest``, ``linear``, and ``cubic`` are aliases for 0, 1, and 3,
            respectively.
        return_indices : bool, default=False
            Whether to return the row and column indices of the nearest grid
            points.

        Returns
        -------
        numpy.ndarray
            The values interpolated at the input points.
        indices : 2-tuple of numpy.ndarray
            The i- and j-indices of the nearest grid points to the input
            points, only present if ``return_indices=True``.

        Raises
        ------
        ValueError
            If an invalid ``method`` is provided.
        RuntimeWarning
            If ``lats`` contains any invalid values outside of the interval
            [-90, 90]. Invalid values will be clipped to this interval.


        .. note::

            If ``return_indices`` is set to ``True``, the nearest array indices
            are returned as a tuple of arrays, in ``(i, j)`` or ``(lat, lon)`` format.


        An example output:

        .. code:: console

            # The first array holds the rows of the raster where point data spatially falls near.
            # The second array holds the columns of the raster where point data spatially falls near.
            sampled_indices = (array([1019, 1019, 1019, ..., 1086, 1086, 1087]), array([2237, 2237, 2237, ...,  983,  983,  983]))
        """
        return sample_grid(
            lon=lons,
            lat=lats,
            grid=self,
            method=method,
            return_indices=return_indices,
        )

    def resample(self, spacingX, spacingY, method="linear", inplace=False):
        """Resamples the grid with a new ``spacingX`` and ``spacingY``, meshed with linear interpolation.

        .. note::

            Ultimately, the :meth:`resample` changes the lat-lon resolution of the gridded data. The
            larger the x and y spacings given are, the larger the pixellation of raster data.

            The :meth:`resample` creates new latitude and longitude arrays with specified spacings in the
            X and Y directions (``spacingX`` and ``spacingY``). These arrays are linearly interpolated
            into a new raster. If ``inplace`` is set to ``True``, the respaced latitude array, longitude
            array and raster will inplace the ones currently attributed to the :class:`Raster` object.

        Parameters
        ----------
        spacingX, spacingY : ndarray
            Specify the spacing in the X and Y directions with which to resample. The larger
            ``spacingX`` and ``spacingY`` are, the larger the raster pixels become (less resolved).
            Note: to keep the size of the raster consistent, set ``spacingX = spacingY``;
            otherwise, if for example ``spacingX > spacingY``, the raster will appear stretched
            longitudinally.

        method : str or int; default: 'linear'
            The order of spline interpolation. Must be an integer in the range
            0-5. 'nearest', 'linear', and 'cubic' are aliases for 0, 1, and 3,
            respectively.

        inplace : bool, default=False
            Choose to overwrite the data (the ``self.data`` attribute), latitude array
            (``self.lats``) and longitude array (``self.lons``) currently attributed to the
            :class:`Raster` object.

        Returns
        -------
        Raster
            The resampled grid. If ``inplace`` is set to ``True``, this raster overwrites the
            one attributed to ``data``.
        """
        spacingX = np.abs(spacingX)
        spacingY = np.abs(spacingY)
        if self.origin == "upper":
            spacingY *= -1.0

        lons = np.arange(self.extent[0], self.extent[1] + spacingX, spacingX)
        lats = np.arange(self.extent[2], self.extent[3] + spacingY, spacingY)
        lonq, latq = np.meshgrid(lons, lats)

        data = self.interpolate(lonq, latq, method=method)
        if inplace:
            self._data = data
            self._lons = lons
            self._lats = lats
        else:
            return Raster(data, self.plate_reconstruction, self.extent, self.time)

    def resize(self, resX, resY, inplace=False, method="linear", return_array=False):
        """Resize the grid with a new resolution (``resX`` and ``resY``) using linear interpolation.

        .. note::

            Ultimately, The :meth:`resize` "stretches" a raster in the x and y directions. The larger
            the resolutions in x and y, the more stretched the raster appears in x and y.

            It creates new latitude and longitude arrays with specific resolutions in
            the X and Y directions (``resX`` and ``resY``). These arrays are linearly interpolated
            into a new raster. If ``inplace`` is set to ``True``, the resized latitude, longitude
            arrays and raster will inplace the ones currently attributed to the :class:`Raster` object.

        Parameters
        ----------
        resX, resY : ndarray
            Specify the resolutions with which to resize the raster. The larger ``resX`` is,
            the more longitudinally-stretched the raster becomes. The larger ``resY`` is, the
            more latitudinally-stretched the raster becomes.

        method : str or int; default: 'linear'
            The order of spline interpolation. Must be an integer in the range
            0-5. 'nearest', 'linear', and 'cubic' are aliases for 0, 1, and 3,
            respectively.

        inplace : bool, default=False
            Choose to overwrite the data (the ``self.data`` attribute), latitude array
            (``self.lats``) and longitude array (``self.lons``) currently attributed to the
            :class:`Raster` object.

        return_array : bool, default False
            Return a ``numpy.ndarray``, rather than a :class:`Raster` object.

        Returns
        -------
        Raster
            The resized grid. If ``inplace`` is set to ``True``, the data in :attr:`Raster.data` will be overwritten.
        """
        # construct grid
        lons = np.linspace(self.extent[0], self.extent[1], resX)
        lats = np.linspace(self.extent[2], self.extent[3], resY)
        lonq, latq = np.meshgrid(lons, lats)

        data = self.interpolate(lonq, latq, method=method)
        if inplace:
            self._data = data
            self._lons = lons
            self._lats = lats
        if return_array:
            return data
        else:
            return Raster(data, self.plate_reconstruction, self.extent, time=self.time)

    def fill_NaNs(self, inplace=False, return_array=False):
        """Search for the invalid ``data`` cells containing NaN-type entries and
        replaces NaNs with the value of the nearest valid data cell.

        Parameters
        ---------
        inplace : bool, default=False
            Choose whether to overwrite the grid currently held in the ``data`` attribute with the filled grid.

        return_array : bool, default False
            Return a ``numpy.ndarray``, rather than a :class:`Raster`.

        Returns
        --------
        Raster
            The resized grid. If ``inplace`` is set to ``True``, the data in :attr:`Raster.data` will be overwritten.
        """
        data = fill_raster(self.data)
        if inplace:
            self._data = data
        if return_array:
            return data
        else:
            return Raster(data, self.plate_reconstruction, self.extent, time=self.time)

    def save_to_netcdf4(self, filename, significant_digits=None, fill_value=np.nan):
        """Saves the grid attributed to the :class:`Raster` object to the given ``filename`` (including
        the ".nc" extension) in netCDF4 format."""
        write_netcdf_grid(
            str(filename), self.data, self.extent, significant_digits, fill_value
        )

    def reconstruct(
        self,
        time,
        fill_value=None,
        partitioning_features=None,
        threads=1,
        anchor_plate_id=None,
        inplace=False,
        return_array=False,
    ):
        """Reconstruct the raster from its initial time (``self.time``) to a new time.

        Parameters
        ----------
        time : float
            Time to which the data will be reconstructed.
        fill_value : float, int, str, or tuple, optional
            The value to be used for regions outside of the static polygons
            at ``time``. By default (``fill_value=None``), this value will be
            determined based on the input.
        partitioning_features : sequence of Feature or str, optional
            The features used to partition the raster grid and assign plate
            IDs. By default, ``self.plate_reconstruction.static_polygons``
            will be used, but alternatively any valid argument to
            ``pygplates.FeaturesFunctionArgument`` can be specified here.
        threads : int, default 1
            Number of threads to use for certain computationally heavy routines.
        anchor_plate_id : int, optional
            ID of the anchored plate. By default, reconstructions are made with respect to
            the anchor plate ID specified in the :class:`PlateReconstruction` object.
        inplace : bool, default False
            Perform the reconstruction in-place (replace the raster's data with the reconstructed data).
        return_array : bool, default False
            Return a ``numpy.ndarray``, rather than a :class:`Raster`.

        Returns
        -------
        Raster or np.ndarray
            The reconstructed grid. Areas for which no plate ID could be determined will be filled with ``fill_value``.

        .. note::

            For two-dimensional grids, ``fill_value`` should be a single
            number. The default value will be ``np.nan`` for float or
            complex types, the minimum value for integer types, and the
            maximum value for unsigned types.
            For RGB image grids, ``fill_value`` should be a 3-tuple RGB
            colour code or a matplotlib colour string. The default value
            will be black (0.0, 0.0, 0.0) or (0, 0, 0).
            For RGBA image grids, ``fill_value`` should be a 4-tuple RGBA
            colour code or a matplotlib colour string. The default fill
            value will be transparent black (0.0, 0.0, 0.0, 0.0) or
            (0, 0, 0, 0).
        """
        try:
            to_time_f = float(time)
        except ValueError:
            raise ValueError(f"Invalid reconstruction time: {time}")
        if to_time_f < 0.0:
            raise ValueError(
                f"The reconstruction time ({to_time_f}) must be greater than 0."
            )

        # A valid PlateReconstruction object is required!
        assert self.plate_reconstruction is not None

        if partitioning_features is None:
            partitioning_features = self.plate_reconstruction.static_polygons

        result = reconstruct_grid(
            grid=self.data,
            partitioning_features=partitioning_features,
            rotation_model=self.plate_reconstruction.rotation_model,
            from_time=self.time,
            to_time=to_time_f,
            extent=self.extent,
            origin=self.origin,
            fill_value=fill_value,
            threads=threads,
            anchor_plate_id=anchor_plate_id,
        )

        raster_rotation_model = self.plate_reconstruction.rotation_model
        # use the new reconstructed raster data to replace the current Raster obj
        # TODO: maybe need to put anchor_plate_id into rotation_model if it is not None
        if inplace:
            self.data = result
            self._time = to_time_f
            if (
                anchor_plate_id is not None
                and raster_rotation_model
                and raster_rotation_model.get_default_anchor_plate_id()
                != anchor_plate_id
            ):
                self.plate_reconstruction.rotation_model = pygplates.RotationModel(
                    raster_rotation_model, default_anchor_plate_id=anchor_plate_id
                )
            if return_array:
                return result
            return self

        # create a new Raster obj to return
        if not return_array:
            result = Raster(
                data=result,
                plate_reconstruction=copy.deepcopy(self.plate_reconstruction),
                extent=self.extent,
                time=to_time_f,
                origin=self.origin,
            )
            if (
                anchor_plate_id is not None
                and raster_rotation_model
                and raster_rotation_model.get_default_anchor_plate_id()
                != anchor_plate_id
            ):
                result.plate_reconstruction.rotation_model = pygplates.RotationModel(
                    raster_rotation_model, default_anchor_plate_id=anchor_plate_id
                )
        return result

    def imshow(self, ax=None, projection=None, **kwargs):
        """Display raster data.

        A pre-existing matplotlib ``Axes`` instance is used if available,
        else a new one is created. The ``origin`` and ``extent`` of the image
        are determined automatically and should not be specified.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            If specified, the image will be drawn within these axes.
        projection : cartopy.crs.Projection, optional
            The map projection to be used. If both ``ax`` and ``projection``
            are specified, this will be checked against the ``projection``
            attribute of ``ax``, if it exists.
        **kwargs : dict, optional
            Any further keyword arguments are passed to
            ``matplotlib.pyplot.imshow`` or ``matplotlib.axes.Axes.imshow``,
            where appropriate.

        Returns
        -------
        matplotlib.image.AxesImage

        Raises
        ------
        ValueError
            If ``ax`` and ``projection`` are both specified, but do not match (i.e. ``ax.projection != projection``).
        """
        for kw in ("origin", "extent"):
            if kw in kwargs.keys():
                raise TypeError(
                    "imshow got an unexpected keyword argument: {}".format(kw)
                )
        if ax is None:
            existing_figure = len(plt.get_fignums()) > 0
            current_axes = plt.gca()
            if projection is None:
                ax = current_axes
            elif (
                isinstance(current_axes, _GeoAxes)
                and current_axes.projection == projection
            ):
                ax = current_axes
            else:
                if not existing_figure:
                    current_axes.remove()
                ax = plt.axes(projection=projection)
        elif projection is not None:
            # projection and ax both specified
            if isinstance(ax, _GeoAxes) and ax.projection == projection:
                pass  # projections match
            else:
                raise ValueError(
                    "Both `projection` and `ax` were specified, but"
                    + " `projection` does not match `ax.projection`"
                )

        if isinstance(ax, _GeoAxes) and "transform" not in kwargs.keys():
            kwargs["transform"] = _PlateCarree()
        extent = self.extent
        if self.origin == "upper":
            extent = (
                extent[0],
                extent[1],
                extent[3],
                extent[2],
            )
        im = ax.imshow(self.data, origin=self.origin, extent=extent, **kwargs)
        return im

    plot = imshow

    def rotate_reference_frames(
        self,
        grid_spacing_degrees,
        reconstruction_time,
        from_rotation_features_or_model=None,  # filename(s), or pyGPlates feature(s)/collection(s) or a RotationModel
        to_rotation_features_or_model=None,  # filename(s), or pyGPlates feature(s)/collection(s) or a RotationModel
        from_rotation_reference_plate=0,
        to_rotation_reference_plate=0,
        non_reference_plate=701,
        output_name=None,
    ):
        """Rotate a grid defined in one plate model reference frame
        within a :class:`Raster` object to another plate reconstruction model reference frame.

        Parameters
        ----------
        grid_spacing_degrees : float
            The spacing (in degrees) for the output rotated grid.
        reconstruction_time : float
            The time at which to rotate the input grid.
        from_rotation_features_or_model : str, list of str, or instance of pygplates.RotationModel
            A filename, or a list of filenames, or a pyGPlates
            RotationModel object that defines the rotation model
            that the input grid is currently associated with.
        to_rotation_features_or_model : str, list of str, or instance of pygplates.RotationModel
            A filename, or a list of filenames, or a pyGPlates
            RotationModel object that defines the rotation model
            that the input grid shall be rotated with.
        from_rotation_reference_plate : int, default = 0
            The current reference plate for the plate model the grid
            is defined in. Defaults to the anchor plate 0.
        to_rotation_reference_plate : int, default = 0
            The desired reference plate for the plate model the grid
            is being rotated to. Defaults to the anchor plate 0.
        non_reference_plate : int, default = 701
            An arbitrary placeholder reference frame with which
            to define the "from" and "to" reference frames.
        output_name : str, default None
            If passed, the rotated grid is saved as a netCDF grid to this filename.

        Returns
        -------
        Raster
            An instance of the :class:`Raster` object containing the rotated grid.
        """

        if from_rotation_features_or_model is None:
            if self.plate_reconstruction is None:
                raise ValueError("Set a plate reconstruction model")
            from_rotation_features_or_model = self.plate_reconstruction.rotation_model
        if to_rotation_features_or_model is None:
            if self.plate_reconstruction is None:
                raise ValueError("Set a plate reconstruction model")
            to_rotation_features_or_model = self.plate_reconstruction.rotation_model

        # Create the pygplates.FiniteRotation that rotates
        # between the two reference frames.
        from_rotation_model = pygplates.RotationModel(from_rotation_features_or_model)
        to_rotation_model = pygplates.RotationModel(to_rotation_features_or_model)
        from_rotation = from_rotation_model.get_rotation(
            reconstruction_time,
            non_reference_plate,
            anchor_plate_id=from_rotation_reference_plate,
        )
        to_rotation = to_rotation_model.get_rotation(
            reconstruction_time,
            non_reference_plate,
            anchor_plate_id=to_rotation_reference_plate,
        )
        reference_frame_conversion_rotation = to_rotation * from_rotation.get_inverse()

        # Resize the input grid to the specified output resolution before rotating
        resX = _deg2pixels(grid_spacing_degrees, self.extent[0], self.extent[1])
        resY = _deg2pixels(grid_spacing_degrees, self.extent[2], self.extent[3])
        resized_input_grid = self.resize(resX, resY, inplace=False)

        # Get the flattened lons, lats
        llons, llats = np.meshgrid(resized_input_grid.lons, resized_input_grid.lats)
        llons = llons.ravel()
        llats = llats.ravel()

        # Convert lon-lat points of Raster grid to pyGPlates points
        input_points = pygplates.MultiPointOnSphere(
            (lat, lon) for lon, lat in zip(llons, llats)
        )
        # Get grid values of the resized Raster object
        values = np.array(resized_input_grid.data).ravel()

        # Rotate grid nodes to the other reference frame
        output_points = reference_frame_conversion_rotation * input_points

        # Assemble rotated points with grid values.
        out_lon = np.empty_like(llons)
        out_lat = np.empty_like(llats)
        zdata = np.empty_like(values)
        for i, point in enumerate(output_points):
            out_lat[i], out_lon[i] = point.to_lat_lon()
            zdata[i] = values[i]

        # Create a regular grid on which to interpolate lats, lons and zdata
        # Use the extent of the original Raster object
        extent_globe = self.extent

        resX = (
            int(np.floor((extent_globe[1] - extent_globe[0]) / grid_spacing_degrees))
            + 1
        )
        resY = (
            int(np.floor((extent_globe[3] - extent_globe[2]) / grid_spacing_degrees))
            + 1
        )

        grid_lon = np.linspace(extent_globe[0], extent_globe[1], resX)
        grid_lat = np.linspace(extent_globe[2], extent_globe[3], resY)

        X, Y = np.meshgrid(grid_lon, grid_lat)

        # Interpolate lons, lats and zvals over a regular grid using nearest
        # neighbour interpolation
        Z = griddata_sphere((out_lon, out_lat), zdata, (X, Y), method="nearest")

        # Write output grid to netCDF if requested.
        if output_name:
            write_netcdf_grid(output_name, Z, extent=extent_globe)

        return Raster(data=Z)

    def query(self, lons, lats, region_of_interest=None):
        """Given a set of location coordinates, return the grid values at these locations.

        Parameters
        ----------
        lons: list
            a list of longitudes of the location coordinates
        lats: list
            a list of latitude of the location coordinates
        region_of_interest: float
            the radius of the region of interest in km
            this is the arch length. we need to calculate the straight distance between the two points in 3D space from this arch length.

        Returns
        -------
        list
            a list of grid values for the given locations.
        """

        if not hasattr(self, "spatial_cKDTree"):
            # build the spatial tree if the tree has not been built yet
            x0 = self.extent[0]
            x1 = self.extent[1]
            y0 = self.extent[2]
            y1 = self.extent[3]
            yn = self.data.shape[0]
            xn = self.data.shape[1]
            # we assume the grid is Grid-line Registration, not Pixel Registration
            # http://www.soest.hawaii.edu/pwessel/courses/gg710-01/GMT_grid.pdf
            # TODO: support both Grid-line and Pixel Registration
            grid_x, grid_y = np.meshgrid(
                np.linspace(x0, x1, xn), np.linspace(y0, y1, yn)
            )
            # in degrees
            self.grid_cell_radius = (
                math.sqrt(math.pow(((y0 - y1) / yn), 2) + math.pow(((x0 - x1) / xn), 2))
                / 2
            )
            self.data_mask = ~np.isnan(self.data)
            grid_points = [
                pygplates.PointOnSphere((float(p[1]), float(p[0]))).to_xyz()
                for p in np.dstack((grid_x, grid_y))[self.data_mask]
            ]
            logger.debug("building the spatial tree...")
            self.spatial_cKDTree = _cKDTree(grid_points)

        query_points = [
            pygplates.PointOnSphere((float(p[1]), float(p[0]))).to_xyz()
            for p in zip(lons, lats)
        ]

        if region_of_interest is None:
            # convert the arch length(in degrees) to direct length in 3D space
            roi = 2 * math.sin(math.radians(self.grid_cell_radius / 2.0))
        else:
            roi = 2 * math.sin(
                region_of_interest / pygplates.Earth.mean_radius_in_kms / 2.0
            )

        dists, indices = self.spatial_cKDTree.query(
            query_points, k=1, distance_upper_bound=roi
        )
        # print(dists, indices)
        return np.concatenate((self.data[self.data_mask], [math.nan]))[indices]

    def clip_by_extent(self, extent):
        """Clip the raster according to a given extent ``(x_min, x_max, y_min, y_max)``.
        The extent of the returned raster may be slightly bigger than the given extent.
        This happens when the border of the given extent fall between two gird lines.

        Parameters
        ----------
        extent: tuple
            A tuple of 4 (min_lon, max_lon, min_lat, max_lat) extent.

        Returns
        --------
        Raster
            The clipped grid.
        """
        if (
            extent[0] >= extent[1]
            or extent[2] >= extent[3]
            or extent[0] < -180
            or extent[1] > 180
            or extent[2] < -90
            or extent[3] > 90
        ):
            raise Exception(f"Invalid extent: {extent}")
        if (
            extent[0] < self.extent[0]
            or extent[1] > self.extent[1]
            or extent[2] < self.extent[2]
            or extent[3] > self.extent[3]
        ):
            raise Exception(
                f"The given extent is out of scope. {extent} -- {self.extent}"
            )
        y_len, x_len = self.data.shape
        logger.debug(f"the shape of raster data x:{x_len} y:{y_len}")

        x0 = math.floor(
            (extent[0] - self.extent[0])
            / (self.extent[1] - self.extent[0])
            * (x_len - 1)
        )
        x1 = math.ceil(
            (extent[1] - self.extent[0])
            / (self.extent[1] - self.extent[0])
            * (x_len - 1)
        )
        # print(x0, x1)
        y0 = math.floor(
            (extent[2] - self.extent[2])
            / (self.extent[3] - self.extent[2])
            * (y_len - 1)
        )
        y1 = math.ceil(
            (extent[3] - self.extent[2])
            / (self.extent[3] - self.extent[2])
            * (y_len - 1)
        )
        # print(y0, y1)
        new_extent = (
            x0 / (x_len - 1) * (self.extent[1] - self.extent[0]) - 180,
            x1 / (x_len - 1) * (self.extent[1] - self.extent[0]) - 180,
            y0 / (y_len - 1) * (self.extent[3] - self.extent[2]) - 90,
            y1 / (y_len - 1) * (self.extent[3] - self.extent[2]) - 90,
        )
        # print(new_extent)
        # print(self.data[y0 : y1 + 1, x0 : x1 + 1].shape)
        return Raster(
            data=self.data[y0 : y1 + 1, x0 : x1 + 1],
            extent=new_extent,
        )

    def _clip_by_polygon(self, polygon):
        """TODO:"""
        pass

    def __array__(self):
        return np.array(self.data)

    def __add__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return self.data + other.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = self.data + other
        new_raster.data = new_data
        return new_raster

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return self.data - other.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = self.data - other
        new_raster.data = new_data
        return new_raster

    def __rsub__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return other.data - self.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = other - self.data
        new_raster.data = new_data
        return new_raster

    def __mul__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return self.data * other.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = self.data * other
        new_raster.data = new_data
        return new_raster

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return self.data / other.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = self.data / other
        new_raster.data = new_data
        return new_raster

    def __rtruediv__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return other.data / self.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = other / self.data
        new_raster.data = new_data
        return new_raster

    def __floordiv__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return self.data // other.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = self.data // other
        new_raster.data = new_data
        return new_raster

    def __rfloordiv__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return other.data // self.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = other // self.data
        new_raster.data = new_data
        return new_raster

    def __mod__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return self.data % other.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = self.data % other
        new_raster.data = new_data
        return new_raster

    def __rmod__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return other.data % self.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = other % self.data
        new_raster.data = new_data
        return new_raster

    def __pow__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return self.data**other.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = self.data**other
        new_raster.data = new_data
        return new_raster

    def __rpow__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster
            # to take properties from
            return other.data**self.data

        # Return Raster with new data
        new_raster = self.copy()
        new_data = other**self.data
        new_raster.data = new_data
        return new_raster


# class TimeRaster(Raster):
#     """A class for the temporal manipulation of raster data. To be added soon!"""
#     def __init__(self, PlateReconstruction_object=None, filename=None, array=None, extent=None, resample=None):
#         raise NotImplementedError(
#             "This class has not been implemented; use `Raster` instead"
#         )
#         super(TimeRaster, self).__init__(PlateReconstruction_object)
