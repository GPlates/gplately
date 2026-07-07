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

import copy
import logging
import math
from os import name
from pathlib import Path
import warnings
from multiprocessing import cpu_count
from typing import List, Tuple, Union, overload, Literal
from enum import Enum

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pygplates
from cartopy.crs import PlateCarree as _PlateCarree
from cartopy.mpl.geoaxes import GeoAxes as _GeoAxes
from rasterio.enums import MergeAlg
from rasterio.features import rasterize as _rasterize
from rasterio.transform import from_bounds as _from_bounds
from scipy.spatial import (
    cKDTree as _cKDTree,  # pyright: ignore[reportAttributeAccessIssue]
)
from scipy.spatial import KDTree

from .lib.exceptions import NegativeReconstructionTime

from .geometry import pygplates_to_shapely
from .reconstruction import PlateReconstruction
from .tools import _deg2pixels, griddata_sphere
from .utils.io_utils import load_feature_collection

logger = logging.getLogger("gplately")


__all__ = [
    "GridRegistration",
    "Raster",
]


# The numbers representing the registration types are the same with PyGMT's grid registration enum values
class GridRegistration(Enum):
    # good for geoscience data because sampled at locations(points) on the grid
    Gridline = 0
    # good for image data, the data represents the average value of an area
    Pixel = 1


class Raster(object):
    """A class to represent a raster grid with time-dependent reconstruction capabilities."""

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
        grid_registration=GridRegistration.Gridline,
        x_dimension_name: str = "",
        y_dimension_name: str = "",
        data_variable_name: str = "",
        **kwargs,
    ):
        """Constructor. Create a :class:`Raster` object.

        Parameters
        ----------
        data : str or array-like or :class:`Raster`
            The raster data, either as a file path (:class:`str`) or array data or a ``Raster`` object.
            If a ``Raster`` object is specified then all other arguments are ignored except ``plate_reconstruction``
            which, if it is not ``None``, will override the plate reconstruction of the ``Raster`` object.
            The data parameter accepts `numpy.ndarray`, `xarray.DataArray` or or any object that can be converted to a `numpy.ndarray`.
            Use `xarray.DataArray` if you want to specify the longitudes and latitudes of the raster data. If you use `numpy.ndarray`, then you must specify the `extent` parameter to tell us the longitudes and latitudes of the raster data.
            The default value is ``None``, which is for backwards compatibility only. In the future, this parameter will be required and the default value will be removed.
        plate_reconstruction : PlateReconstruction
            A :class:`PlateReconstruction` object for raster reconstruction.
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
            The geological time of the time-dependant raster data.
        origin : {'lower', 'upper'}, optional
            When ``data`` is an array, use this parameter to specify the origin
            (upper left or lower left) of the data (overriding ``extent``).
        cell_registration : {'gridline', 'pixel'}, optional, default: 'gridline'
            Specify whether the raster data is gridline-registered or pixel-registered.
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
        # initialise the attribute to None to avoid potential issues with the setter method
        self._plate_reconstruction = None
        # set the initial reconstruction time of the Raster. This will not reconstruct the raster data.
        # the initial reconstruction time is used to indicate the geological time of the provided raster data.
        # later, when the user calls the setter method of the `time` property, the raster data will be reconstructed to the new time.
        self._time = self._get_valid_reconstruction_time(time)
        self._grid_registration = grid_registration
        self._data_var_name = data_variable_name
        self._lons = None
        self._lats = None
        self._default_value = None

        # deal with deprecated arguments, such as ``PlateReconstruction_object``, ``filename``, and ``array``
        # if, in some exceptional cases, the user has to use the deprecated arguments,
        # we will still allow them to do so only when the new `data` and `plate_reconstruction` parameters are not provided,
        # but we will raise warnings and errors as needed.
        # and this function will check for unexpected keyword arguments.
        data, plate_reconstruction = self._handle_deprecated_args(
            data, plate_reconstruction, kwargs
        )

        # if the "data" parameter is a "Raster" object, we do a copy from the other Raster object
        # we also allow the user to override the plate reconstruction of the other Raster object by
        # providing a new plate reconstruction object.
        if isinstance(data, self.__class__):
            self._copy_constructor(data, plate_reconstruction)
            return

        self.plate_reconstruction = plate_reconstruction
        assert (
            data is not None
        ), "`data` argument (or `filename` or `array`) is required."

        # handle the data parameter is a path to a NetCDF file
        if isinstance(data, str):
            if not Path(data).is_file():
                raise FileNotFoundError(f"File not found: {data}")
            self._filename = data
            self._data, lons, lats = self.read_netcdf_grid(
                data,
                return_grids=True,
                realign=realign,
                resample=resample,
                resize=resize,
                x_dimension_name=x_dimension_name,
                y_dimension_name=y_dimension_name,
                data_variable_name=data_variable_name,
            )
            if np.ma.isMaskedArray(self._data):
                self._data = np.ma.asarray(self._data, dtype=float).filled(np.nan)
            self._lons = lons
            self._lats = lats
        else:
            # if the "data" parameter is a numpy array or xarray.DataArray object
            self._filename = None
            extent = self._parse_extent(extent, origin)

            # if the "data" parameter is an xarray.Dataset object, we will try to get a DataArray by the data_variable_name
            # or the first data variable in the Dataset will be used.
            if isinstance(data, xr.Dataset):
                if data_variable_name and data_variable_name in data.data_vars:
                    data = data[data_variable_name]
                else:
                    first_var = next(iter(data.data_vars))
                    data = data[first_var]

            # try to extract the extent from input data if it is an xarray.DataArray object
            if isinstance(data, xr.DataArray):
                extent_from_data = self._find_extent_from_data(data, origin)
                if extent_from_data is not None and extent != extent_from_data:
                    extent = extent_from_data
                    logger.info(
                        f"Raster.__init__(): Use the extent extracted from xarray.DataArray: {extent}."
                    )

            if np.ma.isMaskedArray(data):
                data = np.ma.asarray(data, dtype=float).filled(np.nan)
            data = self._check_grid(data)
            self._data = np.array(data)  # copy to avoid modifying original data

            # get lons and lats from the input data if it is an xarray.DataArray object
            if isinstance(data, xr.DataArray):
                for name in data.coords:
                    if name == x_dimension_name or self._is_a_common_name_for_latitude(
                        str(name)
                    ):
                        self._lats = data.coords[name]
                    if (
                        name == y_dimension_name
                        or self._is_a_common_name_for_longitude(str(name))
                    ):
                        self._lons = data.coords[name]

            # if we cannot find reliable lons and lats from the input data, we will generate them based on the extent and the shape of the data
            if not self._lons:
                self._lons = np.linspace(extent[0], extent[1], self.data.shape[1])
            if not self._lats:
                self._lats = np.linspace(extent[2], extent[3], self.data.shape[0])

            # we realign the grid to -180/180 when the longitudes are from 0 to 360
            # this is a temporary fix. we need a more sophisticated solution.
            # for example, some people may use (-360-0) or some other ranges for longitudes. It is unlikely, but possible.
            if np.max(self._lons) > 180:
                # realign to -180,180 and flip grid if needed
                self._data, self._lons, self._lats = self._realign_grid(
                    self._data, self._lons, self._lats
                )

        if (not isinstance(data, str)) and (resample is not None):
            self.resample(*resample, inplace=True)

        if (not isinstance(data, str)) and (resize is not None):
            self.resize(*resize, inplace=True)

    def to_data_array(self, name=""):
        """Convert the raster to an xarray DataArray with latitude and longitude coordinates."""
        if not name:
            name = self._data_var_name if self._data_var_name else "z"
        da = xr.DataArray(
            self.data,
            coords={
                "lat": (
                    "lat",
                    self.lats,
                    {
                        "standard_name": "latitude",
                        "long_name": "latitude",
                        "units": "degrees_north",
                    },
                ),
                "lon": (
                    "lon",
                    self.lons,
                    {
                        "standard_name": "longitude",
                        "long_name": "longitude",
                        "units": "degrees_east",
                    },
                ),
            },
            dims=["lat", "lon"],
            name=name,
        )
        return da

    @property
    def time(self) -> float:
        """The geological time of the time-dependant raster data."""
        return self._time

    @property
    def default_value(self):
        """The default value used for the raster data."""
        if self._default_value is None:
            dtype = self.data.dtype
            if np.issubdtype(dtype, np.floating):
                self._default_value = np.nan
            elif np.issubdtype(dtype, np.signedinteger):
                self._default_value = np.iinfo(dtype).min
            elif np.issubdtype(dtype, np.unsignedinteger):
                self._default_value = np.iinfo(dtype).max
            else:
                self._default_value = 0
        return self._default_value

    @time.setter
    def time(self, new_time: float):
        """Set a new reconstruction time and reconstruct the raster if necessary."""
        new_time_f = self._get_valid_reconstruction_time(new_time)
        if not math.isclose(self._time, new_time_f):
            self._time = new_time_f
            logger.info(
                f"Raster.time: Reconstructing raster data to new time {new_time_f} Ma."
            )
            self.reconstruct(new_time_f, inplace=True)

    @property
    def data(self) -> np.ndarray:
        """Numpy array containing the raster data."""
        return self._data

    @data.setter
    def data(self, z):
        z = np.array(z)
        if z.shape != np.shape(self.data):
            raise ValueError(
                f"Shape mismatch error: old dimensions are {np.shape(self.data)}, new dimensions are {z.shape}"
            )
        self._data = z

    @property
    def lons(self) -> np.ndarray:
        """The longitude coordinates of the raster data."""
        assert self._lons is not None, "Longitude coordinates are not set."
        return self._lons

    @lons.setter
    def lons(self, x):
        x = np.array(x).ravel()
        if x.size != np.shape(self.data)[1]:
            raise ValueError(
                f"The length of the new longitude array ({x.size}) does not match the number of columns in the data array ({np.shape(self.data)[1]})"
            )
        self._lons = x

    @property
    def lats(self) -> np.ndarray:
        """The latitude coordinates of the raster data."""
        assert self._lats is not None, "Latitude coordinates are not set."
        return self._lats

    @lats.setter
    def lats(self, y):
        y = np.array(y).ravel()
        if y.size != np.shape(self.data)[0]:
            raise ValueError(
                f"The length of the new latitude array ({y.size}) does not match the number of rows in the data array ({np.shape(self.data)[0]})"
            )
        self._lats = y

    @property
    def extent(self) -> Tuple[float, float, float, float]:
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
    def origin(self) -> Literal["lower", "upper"]:
        """The origin (``lower`` or ``upper``) of the data array."""
        if self.lats[0] < self.lats[-1]:
            return "lower"
        else:
            return "upper"

    @property
    def shape(self) -> Tuple[int, int]:
        """The shape of the data array."""
        return np.shape(self.data)

    @property
    def size(self) -> int:
        """The size of the data array."""
        return np.size(self.data)

    @property
    def dtype(self) -> np.dtype:
        """The data type of the array."""
        return self.data.dtype

    @property
    def ndim(self) -> int:
        """The number of dimensions in the array."""
        return np.ndim(self.data)

    @property
    def filename(self) -> Union[str, None]:
        """The filename used to create the :class:`Raster` object.
        If the object was created directly from an array, this attribute is ``None``."""
        return self._filename

    @property
    def plate_reconstruction(self) -> Union[PlateReconstruction, None]:
        """A :class:`PlateReconstruction` object for raster reconstruction."""
        return self._plate_reconstruction

    @plate_reconstruction.setter
    def plate_reconstruction(self, new_plate_recon_obj):
        if new_plate_recon_obj is not None and not isinstance(
            new_plate_recon_obj, PlateReconstruction
        ):
            # Convert to a `PlateReconstruction` if possible
            try:
                new_plate_recon_obj = PlateReconstruction(*new_plate_recon_obj)
            except Exception:
                new_plate_recon_obj = PlateReconstruction(new_plate_recon_obj)
        self._plate_reconstruction = new_plate_recon_obj

    def copy(self) -> "Raster":
        """Return a copy of the :class:`Raster` object."""
        return Raster(
            self.data.copy(),
            copy.deepcopy(self.plate_reconstruction),
            self.extent,
            time=self.time,
        )

    def gap_filling(self, method="nearest", inplace=False) -> "Raster":
        """Fill gaps in the raster data using the specified method.

        Parameters
        ----------
        method : str, default: 'nearest'
            The method to use for gap filling. Options include 'nearest', 'linear', and 'cubic'.
        inplace : bool, default: False
            If True, modify the current Raster object. If False, return a new Raster object.

        Returns
        -------
        Raster
            A new Raster object with gaps filled, or the current object if inplace is True.
        """
        # TODO: generated code, not tested yet, do not use
        # Create a mask of the valid (non-NaN) data points
        valid_mask = ~np.isnan(self.data)
        valid_points = np.column_stack(
            (self.lons[valid_mask.any(axis=0)], self.lats[valid_mask.any(axis=1)])
        )
        valid_values = self.data[valid_mask]

        # Create a grid of points for interpolation
        lon_grid, lat_grid = np.meshgrid(self.lons, self.lats)
        grid_points = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))

        # Perform interpolation to fill gaps
        from scipy.interpolate import griddata

        filled_data = griddata(
            valid_points, valid_values, grid_points, method=method
        ).reshape(self.data.shape)

        if inplace:
            self.data = filled_data
            return self
        else:
            return Raster(
                filled_data,
                copy.deepcopy(self.plate_reconstruction),
                self.extent,
                time=self.time,
            )

    @overload
    def interpolate(
        self,
        lons,
        lats,
        method="linear",
        *,
        return_indices: Literal[False] = False,
    ) -> np.ndarray: ...

    @overload
    def interpolate(
        self,
        lons,
        lats,
        method="linear",
        *,
        return_indices: Literal[True],
    ) -> tuple[np.ndarray, tuple[np.ndarray, np.ndarray]]: ...

    def interpolate(
        self,
        lons,
        lats,
        method="linear",
        *,
        return_indices=False,
    ) -> Union[np.ndarray, tuple[np.ndarray, tuple[np.ndarray, np.ndarray]]]:
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
        return self.sample_grid(
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
            The resampled grid. If ``inplace`` is set to ``True`` the returned raster is ``self``,
            otherwise a new :class:`Raster` object is returned.
        """
        spacingX = np.abs(spacingX)
        spacingY = np.abs(spacingY)
        if self.origin == "upper":
            spacingY *= -1.0

        # Don't need to resample if the spacings are the same.
        dX = np.diff(self.lons).mean()
        dY = np.diff(self.lats).mean()
        if np.isclose(dX, spacingX) and np.isclose(dY, spacingY):
            if inplace:
                # Return current Raster.
                return self
            else:
                # Return a COPY of the current Raster.
                return Raster(
                    self.data,  # note: this gets copied in Raster constructor
                    self.plate_reconstruction,
                    self.extent,
                    time=self.time,
                )

        lons = np.arange(self.extent[0], self.extent[1] + spacingX, spacingX)
        lats = np.arange(self.extent[2], self.extent[3] + spacingY, spacingY)
        lonq, latq = np.meshgrid(lons, lats)

        data = self.interpolate(lonq, latq, method=method)
        if inplace:
            self._data = data
            self._lons = lons
            self._lats = lats
            # Return current Raster.
            return self
        else:
            # Return a new Raster.
            return Raster(data, self.plate_reconstruction, self.extent, time=self.time)

    @overload
    def resize(
        self, resX, resY, inplace=False, method="linear", *, return_array: Literal[True]
    ) -> np.ndarray: ...
    @overload
    def resize(
        self,
        resX,
        resY,
        inplace=False,
        method="linear",
        *,
        return_array: Literal[False] = False,
    ) -> "Raster": ...
    def resize(
        self, resX, resY, inplace=False, method="linear", *, return_array=False
    ) -> Union[np.ndarray, "Raster"]:
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
        Raster or numpy.ndarray
            The resized grid. If ``inplace`` is set to ``True`` the returned raster is ``self``
            (or returned array is ``self.data`` if ``return_array`` is ``True``), otherwise a new :class:`Raster`
            object is returned (or a new data array if ``return_array`` is ``True``).
        """
        # Don't need to resize if the shape is the same.
        if resX == len(self.lons) and resY == len(self.lats):
            if inplace:
                # Return the current data or Raster.
                return self.data if return_array else self
            else:
                # Return a COPY of the current data or Raster.
                return (
                    self.data.copy()
                    if return_array
                    else Raster(
                        self.data,  # note: this gets copied in Raster constructor
                        self.plate_reconstruction,
                        self.extent,
                        time=self.time,
                    )
                )

        # construct grid
        lons = np.linspace(self.extent[0], self.extent[1], resX)
        lats = np.linspace(self.extent[2], self.extent[3], resY)
        lonq, latq = np.meshgrid(lons, lats)

        data: np.ndarray = self.interpolate(lonq, latq, method=method)

        if inplace:
            self._data = data
            self._lons = lons
            self._lats = lats

            # Return the current data or Raster.
            if return_array:
                return self.data
            else:
                return self

        else:
            # Return the new data or Raster.
            if return_array:
                return data
            else:
                return Raster(
                    data,
                    self.plate_reconstruction,
                    self.extent,
                    time=self.time,
                )

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
        data = self.fill_raster(self.data)
        if inplace:
            self._data = data
        if return_array:
            return data
        else:
            return Raster(data, self.plate_reconstruction, self.extent, time=self.time)

    def save_to_netcdf4(self, filename, significant_digits=None, fill_value=None):
        """Saves the grid attributed to the :class:`Raster` object to the given ``filename`` (including
        the ".nc" extension) in netCDF4 format."""
        self.write_netcdf_grid(
            str(filename), self.data, self.extent, significant_digits, fill_value
        )

    @overload
    def reconstruct(
        self,
        time,
        *,
        fill_value=None,
        partitioning_features=None,
        threads=1,
        anchor_plate_id=None,
        inplace=False,
        return_array: Literal[False] = False,
    ) -> "Raster": ...

    @overload
    def reconstruct(
        self,
        time,
        *,
        fill_value=None,
        partitioning_features=None,
        threads=1,
        anchor_plate_id=None,
        inplace=False,
        return_array: Literal[True],
    ) -> np.ndarray: ...

    def reconstruct(
        self,
        time,
        *,
        fill_value=None,
        partitioning_features=None,
        threads=1,
        anchor_plate_id=None,
        inplace=False,
        return_array=False,
    ) -> Union["Raster", np.ndarray]:
        """Reconstruct the raster from its current time to a new time.

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
        to_time_f = self._get_valid_reconstruction_time(time)

        assert (
            self.plate_reconstruction is not None
        ), "A valid PlateReconstruction object is required!"

        assert (
            self.plate_reconstruction.rotation_model is not None
        ), "A valid RotationModel object is required!"

        if partitioning_features is None:
            partitioning_features = self.plate_reconstruction.static_polygons

        assert (
            partitioning_features is not None
        ), "No partitioning features, such as static polygons, provided!"

        use_old_implementation = False
        if use_old_implementation:
            result = self.reconstruct_grid(
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
        else:
            # prepare rotation model for reconstruction
            rotation_model = None
            if anchor_plate_id is not None:
                # the self._get_rotation_model_with_a_different_default_anchor_plate_id() return None
                # if the anchor_plate_id is the same as the default anchor plate ID in the current rotation model
                # but it doesn't matter for the function call self._recconstruct_raster() because it will use the default rotation model in self.plate_reconstruction.rotation_model if rotation_model is None
                rotation_model = (
                    self._get_rotation_model_with_a_different_default_anchor_plate_id(
                        anchor_plate_id
                    )
                )

            result = self._recconstruct_raster(
                to_time=to_time_f,
                rotation_model=rotation_model,
                partitioning_features=partitioning_features,
                threads=threads,
            )

        # use the new reconstructed raster data to replace the current Raster obj
        # put anchor_plate_id into rotation_model if it is not None
        if inplace:
            self.data = result
            self._time = to_time_f
            if (
                anchor_plate_id is not None
                and (
                    rot_model := self._get_rotation_model_with_a_different_default_anchor_plate_id(
                        anchor_plate_id
                    )
                )
                is not None
            ):
                self.plate_reconstruction.rotation_model = rot_model
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
                result.plate_reconstruction is not None
                and anchor_plate_id is not None
                and (
                    rot_model := self._get_rotation_model_with_a_different_default_anchor_plate_id(
                        anchor_plate_id
                    )
                )
                is not None
            ):
                result.plate_reconstruction.rotation_model = rot_model
        return result

    def _recconstruct_raster(
        self,
        to_time: float,
        rotation_model: Union[pygplates.RotationModel, None] = None,
        partitioning_features: Union[pygplates.FeatureCollection, None] = None,
        threads: Union[int, None] = None,
    ) -> np.ndarray:
        """Reconstruct the raster from its current time to a new time.

        .. note::
            The `RotationModel` object associated with the orginal raster data in this object is self.plate_reconstruction.rotation_model.
            User may provide a different `RotationModel` object to reconstruct the raster data to a new time, such as useing a different anchor plate ID.

        This is a private method and not intended to be called directly by users.
        Instead, users should call the public method :meth:`reconstruct` to reconstruct the raster data.
        The doc is provided here for developers/AI assistants who want to understand the implementation of the raster reconstruction.

        Parameters
        ----------
        to_time : float
            Time to which the data will be reconstructed.
        rotation_model : Union[pygplates.RotationModel, None], default None
            The rotation model used to reconstruct the raster data to the `to_time`.
            If None is provided, the default rotation model in self.plate_reconstruction.rotation_model will be used.
        partitioning_features : Union[pygplates.FeatureCollection, None], default None
            The features used to partition the raster grid and assign plate IDs, such as static polygons.
            If None is provided, the `self.plate_reconstruction.static_polygons` will be used.
        threads : int, default 1
            Number of threads to use for certain computationally heavy routines.

        Returns
        -------
        np.ndarray
            The reconstructed raster data. Areas outside of the partitioning polygons will be filled with invalid/default values.
        """
        # The original time of the raster data, which is the time when the raster data was created or last reconstructed.
        from_time = self.time

        # If the intended reconstruction time is the same as the original time,
        # and no different rotation model is provided, we will not perform any reconstruction and return the raster data unchanged.
        # When user provides a different rotation model, but the same `to_time` and `from_time`,
        # it means user wants to reconstruct the raster data using a different anchor plate ID or
        # even a different compatible rotation model, we go ahead and reconstruct it.
        if rotation_model is None and math.isclose(to_time, from_time):
            warnings.warn(
                "Reconstruction time is the same as the original time and no other rotation model is provided; returning input grid unchanged",
                UserWarning,
            )
            return self.data

        # Get the partitioning polygons at the `from_time` and their corresponding features(we need the features to query plate IDs, or other attributes).
        # Only get those polygons that are valid at both the `from_time` and `to_time`,
        # and if they intersect with the raster data extent(not global raster), they will be cut to fit into the extent.
        # This will make sure the reconstructed raster doesn't contain interpolated data which shouldn't be there.
        # For example, if the orignal raster covers only half of Australia, you may not want the interpolated data to show up in the other half of Australia after reconstruction.
        if partitioning_features is None:
            assert (
                self.plate_reconstruction is not None
            ), "A valid PlateReconstruction object is required here!"
            partitioning_features = load_feature_collection(
                self.plate_reconstruction.static_polygons
            )
        (
            partitioning_polygons_at_from_time,
            features_of_partitioning_polygons_at_from_time,
        ) = self._get_partitioning_polygons_at_from_time(
            partitioning_features=partitioning_features,
            from_time=from_time,
            to_time=to_time,
        )
        logger.debug(
            f"Number of partitioning polygons at from_time {from_time}: {len(partitioning_polygons_at_from_time)}"
        )
        assert len(partitioning_polygons_at_from_time) == len(
            features_of_partitioning_polygons_at_from_time
        ), "The number of partitioning polygons and their corresponding features must be the same."

        # Partition the original raster data points using the partitioning polygons
        # and get the plate IDs from the associated features for each data point.
        m_lons, m_lats = np.meshgrid(self.lons, self.lats)
        plate_id_grid = self._get_plate_id_grid(
            m_lons,
            m_lats,
            partitioning_polygons_at_from_time,
            features_of_partitioning_polygons_at_from_time,
        )
        assert (
            plate_id_grid.shape == m_lons.shape
        ), "The shape of the plate ID grid should match the shape of the raster data grid."

        logger.debug(
            f"Number of valid raster data points inside partitioning polygons at from_time {from_time}: {np.count_nonzero(~np.isnan(plate_id_grid))}"
        )

        # Reconstruct the raster data points to the `to_time` using the given rotation model and plate IDs.
        # The returned `reconstructed_original_sample_point_lat_lon_array` and `original_sample_point_row_col_index_array` are 1D numpy arrays of the same length,
        # containing the reconstructed latitudes and longitudes of the raster data points and
        # the corresponding row and column indices of the raster data points in the original raster.
        # Note: only the raster data points that are inside the partitioning polygons at `from_time` will be reconstructed.
        (
            reconstructed_original_sample_point_lat_lon_array,
            original_sample_point_row_col_index_array,
        ) = self._reconstruct_raster_data_points(
            from_time=from_time,
            to_time=to_time,
            m_lons=m_lons,
            m_lats=m_lats,
            plate_id_grid=plate_id_grid,
            rotation_model=rotation_model,
        )

        # We also need to reconstruct the partitioning polygons to the `to_time`.
        partitioning_polygons_at_to_time: List[pygplates.PolygonOnSphere] = (
            self._get_partitioning_polygons_at_to_time(
                partitioning_polygons_at_from_time,
                features_of_partitioning_polygons_at_from_time,
                from_time,
                to_time,
                rotation_model,
            )
        )

        # Initialize the output raster and get the sample points that are inside the partitioning polygons at `to_time`.
        (
            output_raster,
            output_sample_points_lons,
            output_sample_points_lats,
            output_sample_points_row_col_indices,
        ) = self._get_output_raster_and_sample_points(partitioning_polygons_at_to_time)

        # build a KDTree from the reconstructed original raster data points
        # and query the nearest neighbor for the sample points of output raster
        reconstructed_original_sample_points_lats = (
            reconstructed_original_sample_point_lat_lon_array[:, 0]
        )
        reconstructed_original_sample_points_lons = (
            reconstructed_original_sample_point_lat_lon_array[:, 1]
        )

        reconstructed_original_sample_points_vecs = self._lat_lon_to_vector(
            reconstructed_original_sample_points_lats,
            reconstructed_original_sample_points_lons,
            degrees=True,
        )
        # Build a KDTree from the reconstructed original raster data points
        tree = KDTree(reconstructed_original_sample_points_vecs)

        output_vecs = self._lat_lon_to_vector(
            output_sample_points_lats,
            output_sample_points_lons,
            degrees=True,
        )
        assert (
            len(output_vecs) > 1
        ), "No enough output sample points to proceed. Please check your partitioning features to make sure there are valid polygons at `to_time`."

        # Compatibility with older versions of SciPy:
        # 'n_jobs' argument was replaced with 'workers'
        if threads is None:
            cpu_no = cpu_count()
            if cpu_no > 2:
                threads = cpu_count() - 2  # leave 2 cores for other tasks
            else:
                threads = 1
        try:
            _, indices = tree.query(output_vecs, k=1, workers=threads)
        except TypeError:
            _, indices = tree.query(output_vecs, k=1, n_jobs=threads)  # type: ignore

        assert isinstance(indices, np.ndarray)

        # Fill the output grid with the values from the original raster data.
        for out_idx, src_idx in enumerate(indices):
            row, col = output_sample_points_row_col_indices[out_idx]
            output_raster[row, col] = self.data[
                original_sample_point_row_col_index_array[src_idx][0],
                original_sample_point_row_col_index_array[src_idx][1],
            ]

        return output_raster

    def _get_output_raster_and_sample_points(
        self, partitioning_polygons_at_to_time
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get the output raster and the sample points that are inside the partitioning polygons at `to_time`.

        Parameters
        ----------
        partitioning_polygons_at_to_time : list[pygplates.PolygonOnSphere]
            A list of partitioning polygons at the `to_time`.

        Returns
        -------
        output_raster : np.ndarray
            A 2D array of the output raster, initialized with None values.
        output_sample_points_lons : np.ndarray
            A 1D array of longitudes of the sample points that are inside the partitioning polygons at `to_time`.
        output_sample_points_lats : np.ndarray
            A 1D array of latitudes of the sample points that are inside the partitioning polygons at `to_time`.
        output_sample_points_row_col_indices : np.ndarray
            A 1D array of tuples, where each tuple contains the row and column indices of the sample points in the output raster that are inside the partitioning polygons at `to_time`.
        """
        # If the original raster extent is global, the output raster extent will also be global.
        left, right, bottom, top = self.extent
        if math.isclose(abs(left - right), 360) or math.isclose(abs(bottom - top), 180):
            output_raster_extent = self.extent
        else:
            # Get output raster extent from the partitioning polygons at the `to_time`.
            coords = np.vstack(
                [
                    polygon.to_lat_lon_array()
                    for polygon in partitioning_polygons_at_to_time
                ]
            )
            output_raster_extent = (
                coords[:, 1].min(),  # min lon
                coords[:, 1].max(),  # max lon
                coords[:, 0].min(),  # min lat
                coords[:, 0].max(),  # max lat
            )

        # Build the initial output mesh grid from the output raster extent
        output_lon_mesh, output_lat_mesh = np.meshgrid(
            np.linspace(
                output_raster_extent[0], output_raster_extent[1], len(self.lons)
            ),
            np.linspace(
                output_raster_extent[2], output_raster_extent[3], len(self.lats)
            ),
        )

        # Preserve trailing dimensions (eg RGB/RGBA channels) so per-cell assignment
        # accepts vector pixel values from multiband source rasters.
        output_shape = output_lon_mesh.shape
        if self.data.ndim > 2:
            output_shape = output_shape + self.data.shape[2:]
        output_raster = np.full(output_shape, self.default_value, dtype=self.data.dtype)

        # Rasterise to get a boolean mask of which output grid cells fall inside any partitioning polygon.
        assert (
            len(partitioning_polygons_at_to_time) > 0
        ), "No partitioning polygons at `to_time` to rasterize."
        shapely_geoms = pygplates_to_shapely(
            partitioning_polygons_at_to_time,
            tessellate_degrees=0.1,
        )
        if not isinstance(shapely_geoms, list):
            shapely_geoms = [shapely_geoms]
        ny, nx = output_lon_mesh.shape
        minx, maxx, miny, maxy = output_raster_extent
        polygon_mask = _rasterize(
            shapes=zip(shapely_geoms, [1] * len(shapely_geoms)),
            out_shape=(ny, nx),
            fill=0,
            dtype=np.uint8,
            merge_alg=MergeAlg.replace,
            transform=_from_bounds(minx, miny, maxx, maxy, nx, ny),
        )
        assert (
            polygon_mask is not None
        ), "Rasterization of partitioning polygons failed. This should never happen."
        polygon_mask = np.flipud(polygon_mask).astype(bool)

        # Get a list of lons and lats of the points that are inside the reconstructed partitioning polygons.
        # Only query data for the points that are inside the reconstructed partitioning polygons.
        output_sample_points_lons = []
        output_sample_points_lats = []
        output_sample_points_row_col_indices = []
        rows, cols = np.where(polygon_mask)
        for row, col in zip(rows, cols):
            output_sample_points_lons.append(output_lon_mesh[row, col])
            output_sample_points_lats.append(output_lat_mesh[row, col])
            output_sample_points_row_col_indices.append((row, col))

        return (
            output_raster,
            np.array(output_sample_points_lons),
            np.array(output_sample_points_lats),
            np.array(output_sample_points_row_col_indices),
        )

    def _get_partitioning_polygons_at_to_time(
        self,
        partitioning_polygons_at_from_time: list[pygplates.PolygonOnSphere],
        features_of_partitioning_polygons_at_from_time: list[pygplates.Feature],
        from_time: float,
        to_time: float,
        rotation_model: pygplates.RotationModel,
    ) -> List[pygplates.PolygonOnSphere]:
        """Get the partitioning polygons at the `to_time` by reconstructing the partitioning polygons at the `from_time`."""
        polygon_feature_collection = pygplates.FeatureCollection()
        for p, f in zip(
            partitioning_polygons_at_from_time,
            features_of_partitioning_polygons_at_from_time,
        ):
            feature = pygplates.Feature()  # type: ignore
            feature.set_reconstruction_plate_id(f.get_reconstruction_plate_id())  # type: ignore
            feature.set_geometry(p)
            polygon_feature_collection.add(feature)

        # If the from_time is not 0, we need to reverse reconstruct the partitioning polygons
        # to the present-day coordinates before reconstructing them to the to_time.
        # Note: this is done by using the original rotation model in self.plate_reconstruction.rotation_model.
        assert self.plate_reconstruction is not None
        if not math.isclose(from_time, 0.0):
            pygplates.reverse_reconstruct(  # type: ignore
                polygon_feature_collection,
                self.plate_reconstruction.rotation_model,
                from_time,
            )
        if rotation_model is None:
            rotation_model = self.plate_reconstruction.rotation_model
        assert rotation_model is not None, "A valid RotationModel object is required!"
        rfgs = []
        pygplates.reconstruct(polygon_feature_collection, rotation_model, rfgs, to_time)  # type: ignore
        polygons_at_to_time = []
        for rfg in rfgs:
            geom = rfg.get_reconstructed_geometry()
            if geom is not None and isinstance(geom, pygplates.PolygonOnSphere):
                polygons_at_to_time.append(geom)
        return polygons_at_to_time

    def _get_plate_id_grid(
        self,
        lon_mesh: np.ndarray,
        lat_mesh: np.ndarray,
        partitioning_polygons: list[pygplates.PolygonOnSphere],
        features_of_partitioning_polygons: list[pygplates.Feature],
    ) -> np.ndarray:
        """Get the plate ID grid for the raster data points based on the partitioning polygons and their corresponding features.

        Parameters
        ----------
        lon_mesh : np.ndarray
            A 2D array of longitudes of the raster data points.
        lat_mesh : np.ndarray
            A 2D array of latitudes of the raster data points.
        partitioning_polygons : list[pygplates.PolygonOnSphere]
            A list of partitioning polygons.
        features_of_partitioning_polygons : list[pygplates.Feature]
            A list of features corresponding to the partitioning polygons.

        Returns
        -------
        np.ndarray
            A 2D array of plate IDs corresponding to the raster data points. Each number in the returned 2D array is the plate ID of the polygon that contains that lat_lon point,
            or -1 if no polygon contains the point.
        """
        assert (
            len(partitioning_polygons) > 0
        ), "At least one partitioning polygon is required! Do not call this function if there are no partitioning polygons at all."
        plate_ids = [
            f.get_reconstruction_plate_id()  # type: ignore
            for f in features_of_partitioning_polygons
        ]

        # Convert pygplates PolygonOnSphere objects to shapely geometries for rasterio
        shapely_geoms = pygplates_to_shapely(
            partitioning_polygons,
            tessellate_degrees=0.1,
        )

        ny, nx = lon_mesh.shape
        minx = float(lon_mesh[0, 0])
        maxx = float(lon_mesh[0, -1])
        miny = float(lat_mesh[0, 0])
        maxy = float(lat_mesh[-1, 0])

        if not isinstance(shapely_geoms, list):
            shapely_geoms = [shapely_geoms]

        return np.flipud(
            _rasterize(
                shapes=zip(shapely_geoms, plate_ids),
                out_shape=(ny, nx),
                fill=-1,
                dtype=np.int32,
                merge_alg=MergeAlg.replace,
                transform=_from_bounds(minx, miny, maxx, maxy, nx, ny),
            )
        )

    def _reconstruct_raster_data_points(
        self,
        from_time: float,
        to_time: float,
        m_lons: np.ndarray,
        m_lats: np.ndarray,
        plate_id_grid: np.ndarray,
        rotation_model: pygplates.RotationModel,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Reconstruct the raster data points from `from_time` to `to_time` using the given rotation model and plate IDs.

        .. note::
            Only the raster data points have valid plate IDs will be reconstructed.
            The raster data points with invalid plate IDs (i.e., -1) will be ignored.

        Parameters
        ----------
        from_time : float
            The time from which to reconstruct the raster data points.
        to_time : float
            The time to which to reconstruct the raster data points.
        m_lons : np.ndarray
            A 2D array of longitudes of the raster data points.
        m_lats : np.ndarray
            A 2D array of latitudes of the raster data points.
        plate_id_grid : np.ndarray
            A 2D array of plate IDs corresponding to the raster data points.
        rotation_model : pygplates.RotationModel
            The rotation model to use for the reconstruction.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            A tuple containing two 2D arrays:
            - The first array contains the reconstructed latitudes and longitudes of the raster data points.
            - The second array contains the row and column indices of the raster data points in the original raster.
        """
        # m_lons, m_lats and plate_id_grid are 2D arrays of the same shape,
        # representing the longitudes, latitudes and plate IDs of the raster data points.
        assert (
            m_lons.shape == m_lats.shape == plate_id_grid.shape
        ), "Input arrays must have the same shape."

        if rotation_model is None:
            rotation_model = self.plate_reconstruction.rotation_model
        assert rotation_model is not None, "A valid RotationModel object is required!"

        reconstructed_lat_lon_array = []
        row_col_index_array = []
        # Group the raster data points by their plate IDs and reconstruct them to the new time
        unique_plate_ids = np.unique(plate_id_grid)
        for plate_id in unique_plate_ids:
            if plate_id == -1:
                # Ignore the raster data points that are outside of the partitioning polygons at `from_time`
                continue
            rows, cols = np.where(plate_id_grid == plate_id)
            for i, j in zip(rows, cols):
                if (
                    not math.isclose(from_time, 0.0)
                    and self.plate_reconstruction is not None
                    and self.plate_reconstruction.rotation_model is not None
                    and rotation_model is not self.plate_reconstruction.rotation_model
                ):
                    # If the from_time is not 0, and the rotation model used for reconstruction is different from the original rotation model in self.plate_reconstruction.rotation_model,
                    # we need to first reconstruct the raster data points to the present-day coordinates using the original rotation model in self.plate_reconstruction.rotation_model, and then reconstruct them to the `to_time` using the given rotation model.
                    assert self.plate_reconstruction.rotation_model is not None
                    rotation_to_present_day = (
                        self.plate_reconstruction.rotation_model.get_rotation(
                            to_time=float(0),
                            from_time=float(from_time),
                            moving_plate_id=int(plate_id),
                        )
                    )
                    rotation_to_time = rotation_model.get_rotation(
                        to_time=float(to_time),
                        from_time=0.0,
                        moving_plate_id=int(plate_id),
                    )
                    rotation = rotation_to_time * rotation_to_present_day
                else:
                    rotation = rotation_model.get_rotation(
                        to_time=float(to_time),
                        from_time=float(from_time),
                        moving_plate_id=int(plate_id),
                    )
                if not isinstance(rotation, pygplates.FiniteRotation):
                    raise ValueError(f"No rotation found for plate ID: {plate_id}")

                rotated_point = rotation * pygplates.PointOnSphere(
                    m_lats[i, j], m_lons[i, j]
                )
                reconstructed_lat_lon_array.append(rotated_point.to_lat_lon())
                row_col_index_array.append((i, j))
        assert len(reconstructed_lat_lon_array) == len(row_col_index_array)
        return np.array(reconstructed_lat_lon_array), np.array(row_col_index_array)

    def _partition_raster_data_points(
        self,
        lon_mesh: np.ndarray,
        lat_mesh: np.ndarray,
        polygons: list[pygplates.PolygonOnSphere],
    ) -> np.ndarray:
        """Find which partitioning polygon each point in the lat_lon mesh belongs to.
        Return a 2D array of the same shape as the input lat_lon mesh.
        Each number in the returned 2D array is the index of the polygon that contains that lat_lon point,
        or -1 if no polygon contains the point.
        """
        polygon_idx_grid = np.full(lon_mesh.shape, -1, dtype=int)
        # Check each grid point to see if it is inside the polygon.
        for (row, col), _ in np.ndenumerate(polygon_idx_grid):
            lat = lat_mesh[row, col]
            lon = lon_mesh[row, col]
            for idx, polygon in enumerate(polygons):
                # If a point is found inside a polygon, assign the polygon index to that grid point
                # and break the loop for that point.
                if polygon.is_point_in_polygon((lat, lon)):
                    polygon_idx_grid[row, col] = idx
                    break  # No need to check other polygons for this point.
        return polygon_idx_grid

    def _get_partitioning_polygons_at_from_time(
        self,
        partitioning_features: Union[
            pygplates.FeatureCollection, List[pygplates.Feature]
        ],
        *,
        from_time: float,
        to_time: float,
    ) -> Tuple[List[pygplates.PolygonOnSphere], List[pygplates.Feature]]:
        """Return the partitioning polygons at `from_time` and their corresponding features.

        - The polygons must be valid at both `from_time` and `to_time`.
        - The polygons will be reconstructed to the `from_time` if `from_time` is greater than 0(not present-day).
        - The polygons will be cut to fit into the extent of the raster. If the raster is global, the cutting step is unnecessary.
        """
        # The polygons must be valid at both from_time and to_time.
        valid_partitioning_features: list[pygplates.Feature] = [
            f
            for f in partitioning_features
            if f.is_valid_at_time(from_time) and f.is_valid_at_time(to_time)
        ]

        partitioning_polygons: list[pygplates.PolygonOnSphere] = []
        associated_features: list[pygplates.Feature] = []
        if from_time > 0.0:
            # Reconstruct the valid partitioning polygons to the `from_time` using the rotation model in this object.
            rfgs = []
            assert self.plate_reconstruction is not None
            pygplates.reconstruct(  # type: ignore
                valid_partitioning_features,
                self.plate_reconstruction.rotation_model,
                rfgs,
                from_time,
            )
            for rfg in rfgs:
                geom = rfg.get_reconstructed_geometry()
                if geom is not None and isinstance(geom, pygplates.PolygonOnSphere):
                    partitioning_polygons.append(geom)
                    associated_features.append(rfg.get_feature())
        else:
            # If from_time is 0, we don't need to reconstruct the partitioning polygons, just use the polygons from the valid partitioning features directly.
            for f in valid_partitioning_features:
                for geom in f.get_all_geometries():
                    if geom is not None and isinstance(geom, pygplates.PolygonOnSphere):
                        partitioning_polygons.append(geom)
                        associated_features.append(f)

        left, right, bottom, top = self.extent
        # If the extent is global, no need to cut the partitioning polygons.
        if math.isclose(abs(left - right), 360) or math.isclose(abs(bottom - top), 180):
            return partitioning_polygons, associated_features

        # Otherwise, define the rectangular extent polygon to cut the partitioning polygons
        extent_polygon = pygplates.PolygonOnSphere((lat, lon) for lon, lat in [(left, bottom), (left, top), (right, top), (right, bottom)])  # type: ignore
        extent_feature = pygplates.Feature()  # type: ignore
        extent_feature.set_geometry(extent_polygon)  # type: ignore
        partitioner = pygplates.PlatePartitioner(  # type: ignore
            pygplates.FeatureCollection([extent_feature]),
            pygplates.RotationModel([]),
        )

        # cut the partitioning polygons with the extent polygon and return the new polygons and their corresponding features
        inside_geometries = []
        outside_geometries = []
        polygons_within_extent: list[pygplates.PolygonOnSphere] = []
        polygons_within_extent_features: list[pygplates.Feature] = []
        for polygon, feature in zip(partitioning_polygons, associated_features):
            partitioner.partition_geometry(polygon, inside_geometries, outside_geometries)  # type: ignore
            for inside_geom in inside_geometries:
                if isinstance(inside_geom, pygplates.PolygonOnSphere):
                    polygons_within_extent.append(inside_geom)
                else:
                    # The PlatePartitioner may cut a polygon into polylines. Convert the polylines into polygons by connecting the endpoints of the polylines to form a closed polygon.
                    polygons_within_extent.append(
                        pygplates.PolygonOnSphere(inside_geom)
                    )
                polygons_within_extent_features.append(feature)
        assert len(polygons_within_extent) == len(polygons_within_extent_features)
        return polygons_within_extent, polygons_within_extent_features

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
        self.data[self.data is None] = 0
        im = ax.imshow(self.data, origin=self.origin, extent=extent, **kwargs)
        return im

    plot = imshow

    def rotate_reference_frames(
        self,
        grid_spacing_degrees,
        reconstruction_time,
        from_rotation_features_or_model=None,
        to_rotation_features_or_model=None,
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
        from_rotation_features_or_model : str, list of str, instance of pygplates.RotationModel, filename(s), or pyGPlates feature(s)/collection(s)
            A filename, or a list of filenames, or a pyGPlates
            RotationModel object that defines the rotation model
            that the input grid is currently associated with.
        to_rotation_features_or_model : str, list of str, instance of pygplates.RotationModel, filename(s), or pyGPlates feature(s)/collection(s)
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
            self.write_netcdf_grid(output_name, Z, extent=extent_globe)

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

    def _get_rotation_model_with_a_different_default_anchor_plate_id(
        self, anchor_plate_id: int
    ) -> Union[pygplates.RotationModel, None]:
        """Check if the specified anchor plate id is different from the default anchor plate id in the current rotation model.
        If they are different, create and return a new rotation model with the specified default anchor plate id.
        Otherwise, return None."""
        assert (
            self.plate_reconstruction and self.plate_reconstruction.rotation_model
        ), "The self.plate_reconstruction.rotation_model must be valid here!"
        if (
            anchor_plate_id
            != self.plate_reconstruction.rotation_model.get_default_anchor_plate_id()
        ):
            return pygplates.RotationModel(
                self.plate_reconstruction.rotation_model,
                default_anchor_plate_id=anchor_plate_id,
            )
        return None

    def _get_valid_reconstruction_time(self, new_time) -> float:
        """Validate the new reconstruction time and return it as a float.
        Raise a ValueError if the new reconstruction time is invalid (e.g., negative or not a number).
        """
        try:
            new_time_f = float(new_time)
            if new_time_f < 0.0:
                raise NegativeReconstructionTime(new_time_f)
            return new_time_f
        except ValueError:
            raise ValueError(
                f"Invalid reconstruction time: {new_time}. Must be a float number greater than 0."
            )

    def _copy_constructor(self, other, plate_reconstruction):
        self._data = other._data.copy()  # type: ignore
        # Use specified plate reconstruction (if specified),
        # otherwise use the plate reconstruction from 'data'.
        if plate_reconstruction is not None:
            self.plate_reconstruction = plate_reconstruction
        else:
            self.plate_reconstruction = copy.deepcopy(other.plate_reconstruction)
        self._lons = other._lons.copy()
        self._lats = other._lats.copy()
        self._time = other._time
        self._filename = other._filename

    def _handle_deprecated_args(self, data, plate_reconstruction, kwargs):
        _data = data
        _plate_reconstruction = plate_reconstruction

        if not plate_reconstruction and "PlateReconstruction_object" in kwargs.keys():
            warnings.warn(
                "`PlateReconstruction_object` keyword argument is deprecated, use `plate_reconstruction` instead",
                DeprecationWarning,
            )
            _plate_reconstruction = kwargs.pop("PlateReconstruction_object")

        if "filename" in kwargs.keys() and "array" in kwargs.keys():
            raise TypeError(
                "The `filename` and `array` arguments are mutually exclusive and both are deprecated. Use `data` instead."
            )

        if data is None and "filename" in kwargs.keys():
            warnings.warn(
                "The `filename` keyword argument is deprecated, use `data` instead",
                DeprecationWarning,
            )
            _data = kwargs.pop("filename")

        if data is None and "array" in kwargs.keys():
            warnings.warn(
                "The `array` keyword argument is deprecated, use `data` instead",
                DeprecationWarning,
            )
            _data = kwargs.pop("array")

        for key in kwargs.keys():
            raise TypeError(
                f"Raster.__init__() got an unexpected keyword argument '{key}'."
            )

        return _data, _plate_reconstruction

    def _lat_lon_to_vector(self, lat, lon, degrees=False):
        from .grids import _lat_lon_to_vector as _impl

        return _impl(lat, lon, degrees=degrees)

    def _is_a_common_name_for_longitude(self, name: str) -> bool:
        """Return True if the `name` parameter is a possible common name for longitude."""
        return name in ["lon", "lons", "longitude", "x", "east", "easting", "eastings"]

    def _is_a_common_name_for_latitude(self, name: str) -> bool:
        """Return True if the `name` parameter is a possible common name for latitude."""
        return name in [
            "lat",
            "lats",
            "latitude",
            "y",
            "north",
            "northing",
            "northings",
        ]

    def __array__(self):
        return np.array(self.data)

    def __add__(self, other):
        if isinstance(other, Raster):
            # Return array, since we don't know which Raster to take properties from
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

    # NOTE: Raster/GridRegistration were moved from grids.py into this module.
    # A number of low-level helper functions are still implemented in grids.py.
    # We provide lazy wrappers here to avoid import-time circular dependencies.
    # TODO: sort out this issue later
    def fill_raster(self, data, invalid=None):
        from .grids import fill_raster as _impl

        return _impl(data, invalid=invalid)

    def _realign_grid(self, array, lons, lats):
        from .grids import _realign_grid as _impl

        return _impl(array, lons, lats)

    def _find_extent_from_data(self, data, origin):
        from .grids import _find_extent_from_data as _impl

        return _impl(data, origin)

    def read_netcdf_grid(self, *args, **kwargs):
        from .grids import read_netcdf_grid as _impl

        return _impl(*args, **kwargs)

    def write_netcdf_grid(self, *args, **kwargs):
        from .grids import write_netcdf_grid as _impl

        return _impl(*args, **kwargs)

    def sample_grid(self, *args, **kwargs):
        from .grids import sample_grid as _impl

        return _impl(*args, **kwargs)

    def reconstruct_grid(self, *args, **kwargs):
        from .grids import reconstruct_grid as _impl

        return _impl(*args, **kwargs)

    def _check_grid(self, data) -> np.ndarray:
        from .grids import _check_grid as _impl

        return _impl(data)

    def _parse_extent(self, extent, origin) -> Tuple[float, float, float, float]:
        from .grids import _parse_extent as _impl

        return _impl(extent, origin)
