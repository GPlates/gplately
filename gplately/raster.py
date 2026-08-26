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
import numbers
from pathlib import Path
import warnings
from multiprocessing import cpu_count
from typing import List, Tuple, Union, cast, overload, Literal
import pygmt
from xarray.core.types import InterpOptions

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr
import pygplates
from pygplates import (
    RotationModel as _RotationModel,
    PolygonOnSphere as _PolygonOnSphere,
    FeatureCollection as _FeatureCollection,
    Feature as _Feature,
    FiniteRotation as _FiniteRotation,
    MultiPointOnSphere as _MultiPointOnSphere,
)
from cartopy.crs import PlateCarree as _PlateCarree
from cartopy.mpl.geoaxes import GeoAxes as _GeoAxes
from rasterio.enums import MergeAlg
from rasterio.features import rasterize as _rasterize
from rasterio.transform import from_bounds as _from_bounds
from scipy.spatial import (
    cKDTree as _cKDTree,  # pyright: ignore[reportAttributeAccessIssue]
)
from scipy.spatial import KDTree
from scipy.spatial.transform import Rotation as _Rotation
from scipy.ndimage import map_coordinates

from .utils.longitude_convert import _convert_longitude, upwrap_antimeridian_wraparound
from .lib.exceptions import NegativeReconstructionTime
from .lib.enums import GridRegistration, LongitudeConvention, InterpMethod

from .geometry import pygplates_to_shapely
from .reconstruction import PlateReconstruction
from .tools import griddata_sphere
from .utils.io_utils import load_feature_collection, load_data_array_from_netcdf

logger = logging.getLogger("gplately")


__all__ = [
    "GridRegistration",
    "Raster",
]


class Raster(object):
    """A class to represent a raster grid with time-dependent reconstruction capabilities."""

    def __init__(
        self,
        data=None,
        plate_reconstruction=None,
        extent: Union[str, tuple] = "global",
        resample=None,
        resize=None,
        time=0.0,
        origin=None,
        *,
        lons=None,
        lats=None,
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
            Warning: The coordinates embeded in the data or the ``lons`` and ``lats`` parameters will override ``extent``.
        resample : 2-tuple, optional
            Optionally resample grid, pass spacing in X and Y direction as a
            2-tuple e.g. resample=(spacingX, spacingY).
        resize : 2-tuple, optional
            Optionally resample grid to X-columns, Y-rows as a
            2-tuple e.g. resample=(resX, resY).
        time : float, default: 0.0
            The geological time of the time-dependant raster data.
        origin : {'lower', 'upper'}, optional
            When ``data`` is a plain numpy array, use this parameter to specify the origin
            (upper left or lower left) of the data.
            Warning: The coordinates embeded in the data or the ``lons`` and ``lats`` parameters will override ``origin``.
        lons : array-like, optional
            1D array of longitude values. If not provided, will be inferred from the extent, origin and the shape of the data.
        lats : array-like, optional
            1D array of latitude values. If not provided, will be inferred from the extent, origin and the shape of the data.
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
        self._lon_name = x_dimension_name
        self._lat_name = y_dimension_name
        self._data_var_name = data_variable_name
        self._lons = None
        self._lats = None
        self._fill_value = None
        self._filename = None
        self._data_mask = None
        self._spatial_cKDTree = None
        self._spatial_cKDTree_values = None

        # deal with deprecated arguments, such as ``PlateReconstruction_object``, ``filename``, and ``array``
        # if, in some exceptional cases, the user has to use the deprecated arguments,
        # we will still allow them to do so only when the new `data` and `plate_reconstruction` parameters are not provided,
        # but we will raise warnings and errors as needed.
        # and this function will check for unexpected keyword arguments.
        data, plate_reconstruction = self._handle_deprecated_args(
            data, plate_reconstruction, kwargs
        )

        # If the "data" parameter is a "Raster" object, we do a copy from the other Raster object.
        # We also allow the user to override the plate reconstruction of the other Raster object by
        # providing a new plate reconstruction object.
        if isinstance(data, self.__class__):
            self._copy_from_other(data)
            if plate_reconstruction is not None:
                self.plate_reconstruction = plate_reconstruction
            return

        self.plate_reconstruction = plate_reconstruction
        assert (
            data is not None
        ), "`data` argument (or `filename` or `array`) is required."

        # Now we load the data. The data can be from the several sources:
        # - a path to a NetCDF file
        # - a numpy array object
        # - an xarray.DataArray object
        # - an xarray.Dataset object
        load_data_success = False
        try:
            # first, try and see if we can load the data from a file, such as a NetCDF file or a GeoTIFF file.
            self._load_data_from_file(filename=data)
            load_data_success = True
        except Exception as e:
            logger.debug(
                f"Raster.__init__(): Failed to load data from file. Error: {e}. Will try to load data from array."
            )
            # if we cannot load the data from a file, we will try to load the data in other ways, such as from a numpy array or an xarray.DataArray object.

        if not load_data_success and isinstance(data, xr.Dataset):
            # if the "data" parameter is an xarray.Dataset object, we will try to get a DataArray by the data_variable_name
            # or the first data variable in the Dataset will be used.
            if isinstance(data, xr.Dataset):
                if data_variable_name and data_variable_name in data.data_vars:
                    data = data[data_variable_name]
                else:
                    first_var = next(iter(data.data_vars))
                    data = data[first_var]
                    logger.debug(
                        f"Raster.__init__(): Using the first data variable({first_var}) in the xarray.Dataset object as the raster data."
                    )
                # at this point, the data variable is now an xarray.DataArray object, and we will continue to load the data from the DataArray object.

        if not load_data_success and isinstance(data, xr.DataArray):
            logger.debug(
                "Raster.__init__(): Loading data from xarray.DataArray object."
            )
            # load self._data, self._lons, self._lats, self._extent from the input data if it is an xarray.DataArray object
            self._load_data_from_data_array(data)
            load_data_success = True

        if np.ma.isMaskedArray(data):
            logger.debug("Raster.__init__(): Loading data from masked array.")
            # If the input data is a masked array, we will convert it to a regular numpy array and fill the masked values with np.nan
            # TODO: Check the input data type and fill the masked array properly
            self._data_masked = data
            self._data = np.ma.asarray(data, dtype=float).filled(np.nan)
            load_data_success = True

        if not load_data_success:
            logger.debug(
                "Raster.__init__(): Loading data from numpy array or other array-like object."
            )
            self._data = np.array(data)  # copy to avoid modifying original data

        self._validate_data()

        # if we cannot find reliable lons and lats from the input data, try the user provided `lons` and `lats`.
        if self._lons is None or self._lats is None:
            self._lons = lons
            self._lats = lats
        # if user did not provide lons and lats, we will generate them based on the extent and the shape of the data
        if self._lons is None or self._lats is None:
            self._lats, self._lons = self._get_lats_lons_from_extent_origin(
                extent, origin
            )

        assert (
            self._lons is not None and self._lats is not None
        ), "Failed to determine the longitudes and latitudes of the raster data."
        assert (
            len(self._lons) == self._data.shape[1]
        ), "The length of the longitude array does not match the number of columns in the data array."
        assert (
            len(self._lats) == self._data.shape[0]
        ), "The length of the latitude array does not match the number of rows in the data array."

        if (not isinstance(data, str)) and (resample is not None):
            self.resample(*resample, inplace=True)

        if (not isinstance(data, str)) and (resize is not None):
            self.resize(*resize, inplace=True)

    @property
    def time(self) -> float:
        """The geological time of the time-dependant raster data."""
        return self._time

    @time.setter
    def time(self, new_time: float):
        """Set a new reconstruction time and reconstruct the raster data."""
        new_time_f = self._get_valid_reconstruction_time(new_time)
        if not math.isclose(self._time, new_time_f):
            self._time = new_time_f
            logger.info(
                f"Reconstructing raster data inplace from current time {self._time} to a new time {new_time_f} Ma."
            )
            self.reconstruct(new_time_f, inplace=True)
        else:
            logger.warning(
                f"The given reconstruction time {new_time_f} Ma is the same as the current time {self._time} Ma. Nothing will be done."
            )

    @property
    def fill_value(self):
        """The fill value used for the raster data.

        This property is being set when this Raster object is being created with Raster reconstruction.
        The `fill_value` means there is no valid data at the corresponding location.

        The value of this property could be:

            - None, which means this raster data was never created by reconstruction.
            - a single number for 2D scalar rasters, such as np.nan, minimum value for signed integer, and the
                maximum value for unsigned interger.
            - a 3-tuple RGB colour code, such as black (0.0, 0.0, 0.0) or (0, 0, 0).
            - a 4-tuple RGBA colour code, such as transparent black (0.0, 0.0, 0.0, 0.0) or (0, 0, 0, 0).
        """
        return self._fill_value

    @fill_value.setter
    def fill_value(self, value):
        assert value is None or isinstance(
            value, (numbers.Number, tuple)
        ), f"fill_value must be None, a number, or a tuple. Got {type(value)}: {value}."
        if isinstance(value, tuple):
            assert len(value) in (3, 4), "fill_value tuple must be of length 3 or 4."
            for v in value:
                assert isinstance(
                    v, numbers.Number
                ), f"fill_value tuple must contain only numbers. Got {type(value)}: {value}."
        self._fill_value = value

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
        self._invalidate_spatial_cache()

    @property
    def masked_data(self) -> np.ma.masked_array:
        """Masked Numpy array for the raster data."""
        if self._data_mask is None:
            self._data_mask = np.isnan(self._data)
        return np.ma.masked_array(self._data, mask=self._data_mask)

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
        self._invalidate_spatial_cache()

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
        self._invalidate_spatial_cache()

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
    def longitude_convention(self) -> LongitudeConvention:
        """The longitude convention of the raster data.

        :type:  LongitudeConvention
        """
        if np.all((self.lons >= 0) & (self.lons <= 180)):
            return LongitudeConvention.POSITIVE_180
        elif np.all((self.lons >= 0) & (self.lons <= 360)) and np.any(
            (self.lons > 180)
        ):
            return LongitudeConvention.POSITIVE_360
        elif np.all((self.lons >= -180) & (self.lons <= 180)) and np.any(
            (self.lons < 0)
        ):
            return LongitudeConvention.SIGNED_180
        else:
            return LongitudeConvention.OTHER

    @property
    def normalized_extent(self) -> Tuple[float, float, float, float]:
        """The conventional spatial extent of the data.
        The "extent" property above may not be in the conventional order,
        especially when the origin is the upper-left corner.
        The "normalized_extent" property always returns the extent in the conventional order.

        Regardless of origin, extent is always:
            extent = [left, right, bottom, top]
            or
            extent = [xmin, xmax, ymin, ymax]

        The format never changes — bottom always means the smaller y value, top always means the larger y value.

        :type:  tuple of 4 floats
        """
        _lon_diffs = np.diff(self.lons)
        increasing_lon = np.all(_lon_diffs > 0)
        decreasing_lon = np.all(_lon_diffs < 0)
        assert (
            increasing_lon or decreasing_lon
        ), "Longitudes are neither strictly increasing nor decreasing"

        _lat_diffs = np.diff(self.lats)
        increasing_lat = np.all(_lat_diffs > 0)
        decreasing_lat = np.all(_lat_diffs < 0)
        assert (
            increasing_lat or decreasing_lat
        ), "Latitudes are neither strictly increasing nor decreasing"

        return (
            np.min(self.lons),
            np.max(self.lons),
            np.min(self.lats),
            np.max(self.lats),
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

    @classmethod
    def from_points(
        cls,
        lon,
        lat,
        values,
        region=None,
        spacing="0.1d",
        tension=0.35,
        preprocess="blockmean",
        **surface_kwargs,
    ) -> "Raster":
        """
        Class method to create a :class:`Raster` object from scattered geographic points.
        Interpolate scattered geographic points onto a regular grid using
        PyGMT's `surface` (Green's-function-based minimum-curvature gridding),
        with optional block-averaging preprocessing to avoid duplicate/
        near-duplicate point errors.

        Parameters
        ----------
        lon, lat, values : array-like
            1D arrays (or lists) of equal length giving the longitude,
            latitude, and data value of each scattered point.
        region : str or list, optional
            GMT-style region specification [xmin, xmax, ymin, ymax].
            If None, it is inferred from the data extent (with no padding).
        spacing : str or float, optional
            Grid spacing passed to `surface`/`blockmean` (e.g. "0.1d" for
            0.1 degree, or "10k" for 10 km). Default "0.1d".
        tension : float, optional
            Tension factor for `surface`, between 0 (minimum curvature,
            can overshoot) and 1 (harmonic, no overshoot). Default 0.35.
        preprocess : {"blockmean", "blockmedian", None}, optional
            Whether to pre-bin the scattered points onto the target grid
            spacing before running `surface`. `surface` requires at most
            one point per grid cell, so this is standard practice for
            real-world (noisy / clustered / duplicate) data. Set to None
            to skip and pass points to `surface` directly. Default
            "blockmean".
        **surface_kwargs :
            Any additional keyword arguments forwarded to `pygmt.surface`
            (e.g. `maxradius`, `convergence`, etc.).

        Returns
        -------
        Raster
            A Raster object containing the gridded surface.

        Example
        -------
        >>> import numpy as np
        >>> lon = np.random.uniform(120, 130, 500)
        >>> lat = np.random.uniform(30, 40, 500)
        >>> val = np.random.uniform(0, 100, 500)
        >>> raster = Raster.from_points(lon, lat, val, spacing="0.05d")
        """
        lon = np.asarray(lon, dtype=float)
        lat = np.asarray(lat, dtype=float)
        values = np.asarray(values, dtype=float)

        if not (len(lon) == len(lat) == len(values)):
            raise ValueError("lon, lat, and values must have the same length")

        df = pd.DataFrame({"x": lon, "y": lat, "z": values})

        if region is None:
            # Snap bounding box outward to an integer multiple of `spacing`
            region = pygmt.info(data=df[["x", "y"]], spacing=spacing)

        region = tuple(region)  # Ensure region is a tuple for PyGMT

        if preprocess == "blockmean":
            df = pygmt.blockmean(data=df, region=region, spacing=spacing)
        elif preprocess == "blockmedian":
            df = pygmt.blockmedian(data=df, region=region, spacing=spacing)
        elif preprocess is not None:
            raise ValueError(
                "preprocess must be one of 'blockmean', 'blockmedian', or None"
            )
        with pygmt.config(GMT_VERBOSE="q"):
            grid = pygmt.surface(
                data=df,
                region=region,
                spacing=spacing,
                tension=tension,
                **surface_kwargs,
            )

        assert (
            grid is not None
        ), "PyGMT surface gridding returned None. Maybe `outgrid` is set (grid output will be stored in file set by `outgrid`)"
        return cls(grid)

    def to_data_array(self, name=""):
        """Convert the raster to an xarray DataArray with spatial coordinates.

        Supports both:

        - 2D scalar rasters with dimensions ``(lat, lon)``
        - 3D RGB/RGBA rasters with dimensions ``(lat, lon, band)``
        """
        if not name:
            name = self._data_var_name if self._data_var_name else "z"

        coords = {
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
        }

        if self.data.ndim == 2:
            dims = ["lat", "lon"]
        elif self.data.ndim == 3:
            if self.data.shape[2] == 3:
                band_labels = np.array(["r", "g", "b"], dtype=object)
            elif self.data.shape[2] == 4:
                band_labels = np.array(["r", "g", "b", "a"], dtype=object)
            else:
                raise ValueError(
                    "3D raster data must have 3 (RGB) or 4 (RGBA) channels."
                )

            coords["band"] = (
                "band",
                band_labels,
                {
                    "long_name": "color channel",
                },
            )
            dims = ["lat", "lon", "band"]
        else:
            raise ValueError(
                f"Unsupported raster dimensionality {self.data.ndim}; expected 2D or 3D data."
            )

        return xr.DataArray(
            self.data,
            coords=coords,
            dims=dims,
            name=name,
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
        use_spatial_tree: bool = False,
        use_old_implementation: bool = False,
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
        use_spatial_tree: bool = False,
        use_old_implementation: bool = False,
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
        use_spatial_tree: bool = False,
        use_old_implementation: bool = False,
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
        use_spatial_tree : bool, default False
            Whether to use a spatial tree for faster feature lookup.

        Returns
        -------
        Raster or np.ndarray
            The reconstructed Raster as a Raster object or a numpy array.
            Areas with no valid data will be filled with ``fill_value``.


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
        ), "Partitioning features, such as static polygons, are required! They can be provided either in the PlateReconstruction object or as an argument to this function."

        # Flag to roll back to the old implementation quickly if needed.
        # Keep this flag here for a while until the dust settles. -- MC 2026-07-20
        if not use_old_implementation:
            rotation_model = None  # if this is None, the default rotation model in self.plate_reconstruction.rotation_model will be used
            if anchor_plate_id is not None:
                # the self._get_rotation_model_with_a_different_default_anchor_plate_id() returns None if the anchor_plate_id is the same as the default anchor plate ID in the current rotation model
                # but it doesn't matter for the function call self._reconstruct_raster() because it will use the default rotation model in self.plate_reconstruction.rotation_model if rotation_model is None
                rotation_model = (
                    self._get_rotation_model_with_a_different_default_anchor_plate_id(
                        anchor_plate_id
                    )
                )
            _parsed_fill_value = self._parse_fill_value(fill_value)
            # Raster normalization is expensive. Normalize your raster beforehand if you need to call reconstruct() repeatedly.
            if not self.is_normalized():
                # If the raster is in 0-360 convention, we need to convert it to -180/180 convention for reconstruction,
                # and then convert it back to 0-360 convention after reconstruction.
                # The reason is that the static polygons are usually defined in -180/180 convention.
                # The mismatch of longitude conventions may cause problems when determining the plate IDs for the raster cells.

                # NOTE: See the warning below
                # TODO: The code below assumes the partitioning features (static polygons) are in -180/180 convention.
                # If they are in 0-360 convention, the code may not work correctly. We need to check the longitude convention of the partitioning features and convert them if necessary.
                _normalized_raster = self.normalized()
                _reconstructed_raster = _normalized_raster._reconstruct_impl(
                    to_time=to_time_f,
                    rotation_model=rotation_model,
                    partitioning_features=partitioning_features,
                    threads=threads,
                    fill_value=_parsed_fill_value,
                    use_spatial_tree=use_spatial_tree,
                )

                if self.longitude_convention == LongitudeConvention.POSITIVE_360:
                    _reconstructed_raster = (
                        _reconstructed_raster.to_longitude_positive_360()
                    )

                if self.origin == "upper":
                    _reconstructed_raster.data = np.flipud(_reconstructed_raster.data)
                    _reconstructed_raster.lats = np.flip(_reconstructed_raster.lats)
            else:
                _reconstructed_raster = self._reconstruct_impl(
                    to_time=to_time_f,
                    rotation_model=rotation_model,
                    partitioning_features=partitioning_features,
                    threads=threads,
                    fill_value=_parsed_fill_value,
                    use_spatial_tree=use_spatial_tree,
                )

            if inplace:
                self.data = _reconstructed_raster.data
                self._time = to_time_f
                self.fill_value = _parsed_fill_value
                if rotation_model is not None:
                    self.plate_reconstruction.rotation_model = rotation_model
                if return_array:
                    return self.data
                return self
            else:
                if return_array:
                    return _reconstructed_raster.data
                return _reconstructed_raster
        else:
            # Code below uses the old implementation for reconstruction
            _result = self._reconstruct_grid(
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
            # Use the new reconstructed raster data to replace the data in the current Raster object.
            # Put anchor_plate_id into rotation_model if it is not None.
            if inplace:
                self.data = _result
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
                    return _result
                return self

            # If `inplace` is False, create a new Raster object to return
            if not return_array:
                _result = Raster(
                    data=_result,
                    plate_reconstruction=copy.copy(self.plate_reconstruction),
                    extent=self.extent,
                    time=to_time_f,
                    origin=self.origin,
                )
                if (
                    _result.plate_reconstruction is not None
                    and anchor_plate_id is not None
                    and (
                        rot_model := self._get_rotation_model_with_a_different_default_anchor_plate_id(
                            anchor_plate_id
                        )
                    )
                    is not None
                ):
                    _result.plate_reconstruction.rotation_model = rot_model
            return _result

    def plot(self, ax_or_fig=None, projection=None, use_gmt=False, **kwargs):
        """Plot the raster data using either matplotlib or pygmt.

        Parameters
        ----------
        ax_or_fig : matplotlib.axes.Axes or matplotlib.figure.Figure, optional
            If specified, the image will be drawn within these axes or figure.
        projection : cartopy.crs.Projection, optional
            The map projection to be used. If both ``ax_or_fig`` and ``projection`
            are specified, this will be checked against the ``projection``
            attribute of ``ax_or_fig``, if it exists.
        use_gmt : bool, default False
            If True, use pygmt to plot the raster data. If False, use matplotlib.
        **kwargs : dict, optional
            Any further keyword arguments are passed to the plotting function.
        """
        if not use_gmt:
            ax = kwargs.pop("ax", None)
            if ax_or_fig and ax:
                raise ValueError("Only one of `ax_or_fig` or `ax` can be specified.")
            if ax and ax_or_fig is None:
                ax_or_fig = ax
            return self.imshow(ax=ax_or_fig, projection=projection, **kwargs)
        else:
            from .plot.pygmt_plot import PygmtPlotEngine

            fig = kwargs.pop("fig", None)
            if ax_or_fig and fig:
                raise ValueError("Only one of `ax_or_fig` or `fig` can be specified.")
            if fig and ax_or_fig is None:
                ax_or_fig = fig
            if not ax_or_fig:
                raise ValueError(
                    "When `use_gmt` is True, either `ax_or_fig` or `fig` must be specified."
                )
            return PygmtPlotEngine().plot_grid(
                ax_or_fig=ax_or_fig, grid=self, projection=projection, **kwargs
            )

    def copy(self) -> "Raster":
        """Return a copy of the :class:`Raster` object."""
        new_obj = self.__class__.__new__(self.__class__)
        new_obj._copy_from_other(self)
        return new_obj

    def is_global(self) -> bool:
        if not math.isclose(abs(self.lats[-1] - self.lats[0]), 180):
            return False

        if np.isclose(self.lons[-1] - self.lons[0], 360):
            # something like [0, 1, 2, ..., 359, 360] or [-180, -179, ..., 179, 180]
            return True
        avg_spacing = np.mean(np.diff(self.lons))
        # something like [0, 1, 2, ..., 358, 359] or [-180, -179, ..., 178, 179], without the duplicate 360 or 180 at the end, but still covering the whole globe
        return (self.lons[-1] - self.lons[0]) >= (360 - 1.5 * avg_spacing)

    def is_normalized(self) -> bool:
        """Check if the raster is normalized.

        For now, the normalized raster means longitude: [-180, 180] and latitude: [-90, 90].
        """
        return (
            self.longitude_convention == LongitudeConvention.SIGNED_180
            and self.origin == "lower"
        )

    def is_longitude_wrapped(self) -> bool:
        """Check if the longitude coordinates have wrapped around the antimeridian."""
        if np.any(np.abs(np.diff(self.lons)) > 180):
            return True
        return False

    def unwrap_longitude(self):
        """Unwrap the longitude coordinates around the antimeridian.
        For example, if the longitude coordinates are [170, 175, 180, -175, -170] in range[-180, 180], they will be unwrapped to [170, 175, 180, 185, 190] in range[0, 360].

        .. note:: There is no need to rearrange the raster data because unwrapping the longitude coordinates does not change the order of the data.
            The data will still be in the same order as before, but the longitude coordinates will be adjusted to avoid wrapping around the antimeridian.
        """
        if self.is_longitude_wrapped():
            self.lons = upwrap_antimeridian_wraparound(self.lons)

    def normalized(self) -> "Raster":
        """Return a normalized Raster object.

        For now, the normalized raster means longitude: [-180, 180] and latitude: [-90, 90].

        Returns
        -------
        Raster
            A new Raster object with normalized longitude and latitude.
            If the current Raster object is already normalized, it will return itself.
        """
        if (
            self.longitude_convention == LongitudeConvention.SIGNED_180
            and self.origin == "lower"
        ):
            return self

        if self.longitude_convention != LongitudeConvention.SIGNED_180:
            _normalized_raster = self.to_longitude_signed_180()
        else:
            _normalized_raster = self.copy()

        if _normalized_raster.origin == "upper":
            _normalized_raster.data = np.flipud(_normalized_raster.data)
            _normalized_raster.lats = np.flip(_normalized_raster.lats)
        return _normalized_raster

    def fill_gaps(
        self,
        *,
        method="nearest",
        use_gmt=False,
        use_spatial_tree=False,
        inplace=False,
        valid_mask=None,
        invalid_value=None,
    ) -> "Raster":
        """Fill invalid cells in a raster using interpolation.

        Supports both scalar rasters (``ny x nx``) and image rasters
        (``ny x nx x channels``) where ``channels`` is 3 (RGB) or 4 (RGBA).

        Parameters
        ----------
        method : str, default: 'nearest'
            Interpolation method passed to :func:`scipy.interpolate.griddata`.
            Typical options are ``'nearest'``, ``'linear'`` and ``'cubic'``.
        use_gmt : bool, default: False
            If ``True``, use PyGMT ``grdfill``/``fillgrd`` for filling gaps.
            This option currently supports only 2D scalar rasters.
        use_spatial_tree : bool, default: False
            If ``True``, use method :meth:`Raster._query_by_KDTree` to fill the gaps.
            Mutually exclusive with ``use_gmt``. This option support both 2D scalar rasters and 3D RGB/RGBA rasters.
        inplace : bool, default: False
            If ``True``, modify and return the current object. Otherwise,
            return a new :class:`Raster`.
        valid_mask : array-like of bool, optional
            A 2D mask with shape ``(len(lats), len(lons))`` indicating valid
            source pixels (``True`` = valid). If omitted, validity is inferred
            from finite values and, if provided, ``invalid_value``.
        invalid_value : scalar or sequence, optional
            Additional marker for invalid pixels.

            - For scalar rasters, cells equal to this value are treated as invalid.
            - For RGB/RGBA rasters, this can be a scalar (applied to all channels)
              or a sequence with length equal to the channel count.

        Returns
        -------
        Raster
            A raster with gaps filled.

        Raises
        ------
        ValueError
            If raster dimensionality/channel count is unsupported, mask shape is
            invalid, or no valid points are available for interpolation.
        """
        data = np.asarray(self.data)
        is_image = self._validate_fill_gaps_data(data)

        if use_gmt and use_spatial_tree:
            raise ValueError("`use_gmt` and `use_spatial_tree` are mutually exclusive.")
        if use_spatial_tree and str(method).lower() != "nearest":
            raise ValueError(
                "use_spatial_tree=True currently supports only method='nearest'."
            )

        effective_valid_mask = self._build_fill_gaps_valid_mask(
            data,
            is_image=is_image,
            valid_mask=valid_mask,
            invalid_value=invalid_value,
        )

        if np.all(effective_valid_mask):
            return self if inplace else self.copy()

        if not np.any(effective_valid_mask):
            raise ValueError("No valid points available to interpolate gaps.")

        gap_mask = ~effective_valid_mask

        if use_gmt:
            filled_data = self._fill_gaps_with_gmt(data, gap_mask, is_image, method)
        elif use_spatial_tree:
            filled_data = self._fill_gaps_with_spatial_tree(
                data,
                gap_mask,
                effective_valid_mask,
            )
        else:
            filled_data = self._fill_gaps_with_griddata(
                data,
                gap_mask,
                effective_valid_mask,
                is_image,
                method,
            )

        if np.issubdtype(data.dtype, np.integer):
            info = np.iinfo(data.dtype)
            filled_data = np.clip(np.rint(filled_data), info.min, info.max).astype(
                data.dtype
            )
        else:
            filled_data = filled_data.astype(data.dtype, copy=False)

        if inplace:
            self.data = filled_data
            return self
        ret = self.copy()
        ret.data = filled_data
        return ret

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

    def interp(
        self,
        *,
        lons,
        lats,
        data: Union[xr.DataArray, None] = None,
        method: InterpMethod = InterpMethod.LINEAR,
        fill_value=np.nan,
        pointwise=True,
    ) -> np.ndarray:
        """
        Interpolate `data` at given longitude/latitude locations.

        `data` may be:
        - 2D, with dims ('lat', 'lon')
        - 3D, with dims ('lat', 'lon', <channel>) for RGB/RGBA images

        Parameters
        ----------
        pointwise : bool
            If True (default), `lons`/`lats` are treated as paired query
            points (lons[i], lats[i]) — e.g. sampling at scattered station
            locations. Output shape: (N,) or (N, C).

            If False, `lons`/`lats` define a new rectangular grid (outer
            product) — e.g. resampling the whole image onto a new lat/lon
            mesh. Output shape: (len(lats), len(lons)) or
            (len(lats), len(lons), C).

        Returns
        -------
        np.ndarray
        """
        lons = np.asarray(lons)
        lats = np.asarray(lats)

        if data is None:
            data = self.to_data_array()

        assert data is not None and isinstance(
            data, xr.DataArray
        ), "The `data` variable must be an xarray DataArray"

        orig_dtype = data.dtype

        if "lat" not in data.dims or "lon" not in data.dims:
            raise ValueError(
                f"The `data` variable must have 'lat' and 'lon' dimensions, got {data.dims}"
            )

        if pointwise:
            if lons.shape != lats.shape:
                raise ValueError(
                    "lons and lats must have the same shape for pointwise interp"
                )

            # wrap query points as DataArrays sharing a common dim so xarray
            # does pointwise (vectorized) interpolation instead of building
            # an outer-product grid
            lon_da = xr.DataArray(lons, dims="points")
            lat_da = xr.DataArray(lats, dims="points")

            result = data.interp(
                lon=lon_da,
                lat=lat_da,
                method=cast(InterpOptions, method.value),
                kwargs={"fill_value": fill_value, "bounds_error": False},
            )

            # move 'points' to the front, keep any remaining dims (e.g. channel) after it
            extra_dims = [d for d in result.dims if d != "points"]
            if extra_dims:
                result = result.transpose("points", *extra_dims)

        else:
            # grid mode: lons/lats are new 1D coordinate axes, xarray
            # broadcasts them as an outer product, replacing 'lon' and 'lat'
            # in place — original dim order (lat, lon, [channel]) is preserved
            result = data.interp(
                lon=lons,
                lat=lats,
                method=cast(InterpOptions, method.value),
                kwargs={"fill_value": fill_value, "bounds_error": False},
            )

        values = result.values

        if np.issubdtype(orig_dtype, np.integer):
            info = np.iinfo(orig_dtype)
            # NaNs can appear from fill_value / out-of-bounds points; interp
            # output can't be safely rounded to int while NaNs are present,
            # so fill them with 0 first (adjust default if 0 isn't a sane nodata value)
            values = np.nan_to_num(values, nan=0)
            values = np.clip(np.round(values), info.min, info.max).astype(orig_dtype)
        else:
            values = values.astype(orig_dtype)

        return values

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
        values, indices = self.sample_values(
            lons=lons,
            lats=lats,
            method=method,
        )
        if return_indices:
            return values, indices
        return values

    def resample(
        self,
        spacingX: float,
        spacingY: float,
        method="linear",
        inplace=False,
        strict_spacing: bool = True,
    ):
        """Resample raster data onto a new lon-lat grid.

        .. note::

            This method changes the lat-lon resolution of the gridded data. Larger
            spacing values produce a coarser grid.

            New latitude and longitude arrays are created from ``spacingX`` and
            ``spacingY``, and data values are interpolated onto that target grid.
            If ``inplace`` is ``True``, the current :class:`Raster` object is updated.

        Parameters
        ----------
        spacingX, spacingY : float
            Target spacing in degrees in the X (longitude) and Y (latitude)
            directions. Both must be positive.

        method : str or int; default: 'linear'
            The order of spline interpolation. Must be an integer in the range
            0-5. 'nearest', 'linear', and 'cubic' are aliases for 0, 1, and 3,
            respectively.

        inplace : bool, default=False
            Choose to overwrite the data (the ``self.data`` attribute), latitude array
            (``self.lats``) and longitude array (``self.lons``) currently attributed to the
            :class:`Raster` object.

        strict_spacing : bool, default=True
            Controls whether spacing or extent is preserved exactly.

            - If ``True``, output spacing is exactly ``spacingX``/``spacingY`` and
                the output extent may differ slightly from the input extent.
            - If ``False``, output extent is preserved exactly and the effective
                spacing may differ slightly from ``spacingX``/``spacingY``.


        Returns
        -------
        Raster
            The resampled grid. If ``inplace`` is set to ``True`` the returned raster is ``self``,
            otherwise a new :class:`Raster` object is returned.
        """
        if not (np.isfinite(spacingX) and np.isfinite(spacingY)):
            raise ValueError("Spacing values must be finite numbers.")
        if spacingX <= 0.0 or spacingY <= 0.0:
            raise ValueError("Spacing must be positive.")

        # Don't need to resample if the spacings are the same.
        dX = np.abs(np.diff(self.lons)).mean()
        dY = np.abs(np.diff(self.lats)).mean()
        if np.isclose(dX, spacingX) and np.isclose(dY, spacingY):
            if inplace:
                return self
            else:
                return Raster(
                    self.data,  # note: this gets copied in Raster constructor
                    self.plate_reconstruction,
                    self.extent,
                    time=self.time,
                )

        x0, x1, y0, y1 = self.extent
        nx = round(np.abs(x1 - x0) / spacingX) + 1
        ny = round(np.abs(y1 - y0) / spacingY) + 1
        if nx < 2 or ny < 2:
            raise ValueError(
                f"Resampling spacings are too large for the given extent. Resulting grid would have shape ({ny}, {nx})."
            )
        if nx > 3601 or ny > 1801:
            logger.warning(
                f"The resulting grid after resampling would have shape ({ny}, {nx}). It may take too long to reconstruct. Consider using smaller spacings."
            )
        if strict_spacing:
            x_direction = 1.0 if x1 >= x0 else -1.0
            y_direction = 1.0 if y1 >= y0 else -1.0
            lons = x0 + x_direction * spacingX * np.arange(nx, dtype=float)
            lats = y0 + y_direction * spacingY * np.arange(ny, dtype=float)
        else:
            lons = np.linspace(x0, x1, nx)
            lats = np.linspace(y0, y1, ny)
        lonq, latq = np.meshgrid(lons, lats)

        data = self.interpolate(lonq, latq, method=method)
        resampled_extent = (
            float(lons[0]),
            float(lons[-1]),
            float(lats[0]),
            float(lats[-1]),
        )
        if inplace:
            self._data = data
            self._lons = lons
            self._lats = lats
            self._invalidate_spatial_cache()
            return self
        else:
            return Raster(
                data,
                self.plate_reconstruction,
                resampled_extent,
                time=self.time,
                origin=self.origin,
            )

    @overload
    def resize(
        self,
        resX: int,
        resY: int,
        inplace=False,
        method="linear",
        *,
        return_array: Literal[True],
    ) -> np.ndarray: ...
    @overload
    def resize(
        self,
        resX: int,
        resY: int,
        inplace=False,
        method="linear",
        *,
        return_array: Literal[False] = False,
    ) -> "Raster": ...
    def resize(
        self,
        resX: int,
        resY: int,
        inplace=False,
        method="linear",
        *,
        return_array=False,
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
        resX, resY : int
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
            self._invalidate_spatial_cache()
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
        """Deprecated. Use :meth:`fill_gaps` instead. This method will be removed in a future version of the library."""
        r = self.fill_gaps(inplace=inplace, method="nearest")
        if return_array:
            return r.data
        else:
            return r

    def save_to_netcdf4(self, filename, significant_digits=None, fill_value=None):
        """Saves the :class:`Raster` object as a netCDF4 file."""
        self._write_netcdf_grid(
            str(filename), self.data, self.extent, significant_digits, fill_value
        )

    def save_as_image(
        self,
        filename,
        projection=None,
        cmap="viridis",
        vmin=None,
        vmax=None,
        dpi=300,
        colorbar_label="",
        data_only=False,
    ):
        """Saves the :class:`Raster` object as a image."""
        if not data_only:
            if projection is None:
                projection = _PlateCarree()

            fig = plt.figure(figsize=(12, 6), dpi=dpi)
            ax = fig.add_subplot(111, projection=projection)
            im = self.plot(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax)
            if im is not None:
                fig.colorbar(
                    im,
                    orientation="horizontal",
                    shrink=0.4,
                    pad=0.05,
                    label=colorbar_label,
                )
            plt.savefig(filename, dpi=dpi, bbox_inches="tight")
        else:
            plt.imsave(filename, np.flipud(self.data), cmap=cmap, vmin=vmin, vmax=vmax)

        plt.close()

    def rotate_reference_frames(
        self,
        reconstruction_time=None,
        from_rotation_features_or_model=None,
        to_rotation_features_or_model=None,
        grid_spacing_degrees=None,
        from_rotation_reference_plate=0,
        to_rotation_reference_plate=0,
        non_reference_plate=701,
        output_name=None,
    ):
        """Rotate a grid defined in one plate model reference frame
        within a :class:`Raster` object to another plate reconstruction model reference frame.

        Parameters
        ----------
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
        grid_spacing_degrees : float, optional
            The spacing (in degrees) for the output rotated grid.
            If not specified, the output grid keeps the input raster shape.
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

        if reconstruction_time is None:
            raise ValueError("`reconstruction_time` is required.")

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
        from_rotation_model = _RotationModel(from_rotation_features_or_model)
        to_rotation_model = _RotationModel(to_rotation_features_or_model)
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

        # Resize the input grid to the specified output resolution before rotating.
        if grid_spacing_degrees is None:
            resX = len(self.lons)
            resY = len(self.lats)
        else:
            if not np.isfinite(grid_spacing_degrees) or grid_spacing_degrees <= 0:
                raise ValueError(
                    "`grid_spacing_degrees` must be a positive finite number."
                )
            from .grids import num_grid_points as _num_grid_points

            resX = _num_grid_points(
                grid_spacing_degrees, self.extent[0], self.extent[1]
            )
            resY = _num_grid_points(
                grid_spacing_degrees, self.extent[2], self.extent[3]
            )
        resized_input_grid = self.resize(resX, resY, inplace=False)

        # Get the flattened lons, lats
        llons, llats = np.meshgrid(resized_input_grid.lons, resized_input_grid.lats)
        llons = llons.ravel()
        llats = llats.ravel()

        # Convert lon-lat points of Raster grid to pyGPlates points
        input_points = _MultiPointOnSphere((lat, lon) for lon, lat in zip(llons, llats))
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

        grid_lon = np.linspace(extent_globe[0], extent_globe[1], resX)
        grid_lat = np.linspace(extent_globe[2], extent_globe[3], resY)

        X, Y = np.meshgrid(grid_lon, grid_lat)

        # Interpolate lons, lats and zvals over a regular grid using nearest
        # neighbour interpolation
        Z = griddata_sphere((out_lon, out_lat), zdata, (X, Y), method="nearest")

        # Write output grid to netCDF if requested.
        if output_name:
            self._write_netcdf_grid(output_name, Z, extent=extent_globe)

        return Raster(data=Z, extent=extent_globe)

    def sample_values(self, *, lons, lats, method="linear"):
        order = {
            "nearest": 0,
            "linear": 1,
            "cubic": 3,
        }.get(method, method)
        if order not in {0, 1, 2, 3, 4, 5}:
            raise ValueError(f"Invalid `method` parameter: {method}")

        extent = self.extent
        grid = self.data

        # Do not wrap from North to South Pole (or vice versa)
        if np.any(np.abs(lats) > 90.0):
            if np.any(np.abs(lats) > 90.0 + 1e-8):
                # Only raise a warning when the values are really invalid, not just slightly out of bounds due to floating point errors.
                warnings.warn(
                    f"Invalid values({lats[np.abs(lats) > 90.0]}) encountered in latitudes; clipping to [-90, 90]",
                    RuntimeWarning,
                )
            lats = np.clip(lats, -90.0, 90.0)

        dx = (extent[1] - extent[0]) / (np.shape(grid)[1] - 1)
        dy = (extent[3] - extent[2]) / (np.shape(grid)[0] - 1)
        point_i = (lats - extent[2]) / dy
        point_j = (lons - extent[0]) / dx

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
            interpolated = np.reshape(interpolated, np.shape(lons))
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
                    np.shape(lons),
                )
                interpolated.append(interpolated_k)
            del interpolated_k
            interpolated = np.stack(interpolated, axis=-1)

        interpolated = interpolated.astype(grid.dtype)

        indices = (
            np.rint(np.ravel(point_i)).astype(np.int_),
            np.rint(np.ravel(point_j)).astype(np.int_),
        )
        return interpolated, indices

    def query(
        self,
        lons: np.ndarray,
        lats: np.ndarray,
        *,
        interpolation_method: str = "nearest",
        region_of_interest: Union[None, float] = None,
        pointwise: bool = True,
    ):
        """Query raster values at given longitude/latitude coordinates.

        Parameters
        ----------
        lons, lats : np.ndarray
            Longitude and latitude coordinates.
        interpolation_method : str, default "nearest"
            Interpolation method, such as "linear" or "nearest". See `Raster.InterpMethod` for details.
        pointwise : bool, default True
            If True, sample paired points `(lons[i], lats[i])`. If False,
            treat `lons` and `lats` as 1D axes of an output grid.


        Returns
        -------
        tuple
            Sampled values at the given coordinates.
        """
        if region_of_interest is not None:
            try:
                roi_float = float(region_of_interest)
            except (ValueError, TypeError):
                raise ValueError(
                    f"Invalid value for region_of_interest: {region_of_interest}"
                )
            if roi_float < 0:
                raise ValueError(
                    f"Invalid value for region_of_interest: {region_of_interest}. It must be non-negative."
                )
            return self._query_by_KDTree(
                lons,
                lats,
                roi_float,
                pointwise=pointwise,
            )

        return self.interp(
            lons=lons,
            lats=lats,
            method=InterpMethod(interpolation_method),
            pointwise=pointwise,
        )

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

    def clip_by_polygons(self, polygons: list[_PolygonOnSphere]):
        """TODO:"""
        pass

    def to_longitude_positive_360(self, inplace=False):
        """Convert a grid's longitude coordinates to the [0, 360] convention."""
        _grid_data, _lons = _convert_longitude(
            self.data,
            self.lons,
            map_fn=lambda l: l % 360,
            valid_max=360,
            seam_name="0/360",
        )
        if inplace:
            self._data = _grid_data
            self._lons = _lons
            self._invalidate_spatial_cache()
            return self
        else:
            new_raster = self.copy()
            new_raster._data = _grid_data
            new_raster._lons = _lons
            return new_raster

    def to_longitude_signed_180(self, inplace=False):
        """Convert a grid's longitude coordinates to the [-180, 180] convention."""
        if self.longitude_convention == LongitudeConvention.SIGNED_180:
            logger.warning(
                "The raster is already in the [-180, 180] longitude convention. No conversion is performed. Returning the original raster object."
            )
            return self
        _grid_data, _lons = _convert_longitude(
            self.data,
            self.lons,
            map_fn=lambda l: ((l + 180) % 360) - 180,
            valid_max=180,
            seam_name="-180/180",
        )
        if inplace:
            self._data = _grid_data
            self._lons = _lons
            self._invalidate_spatial_cache()
            return self
        else:
            new_raster = self.copy()
            new_raster._data = _grid_data
            new_raster._lons = _lons
            return new_raster

    def imshow(self, ax=None, projection=None, **kwargs):
        """Deprecated. Use :meth:`plot` instead. Plot the raster data using matplotlib."""
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

    def _reconstruct_impl(
        self,
        to_time: float,
        rotation_model: Union[_RotationModel, None] = None,
        partitioning_features: Union[_FeatureCollection, None] = None,
        *,
        fill_value: Union[float, int, str, tuple],
        threads: Union[int, None] = None,
        use_spatial_tree: bool = False,
    ) -> "Raster":
        """Reconstruct the raster from its current time to a new time.

        .. note::
            The `RotationModel` object associated with the orginal raster data in this object is self.plate_reconstruction.rotation_model.
            User may provide a different `RotationModel` object to reconstruct the raster data to a new time, such as useing a different anchor plate ID.

        This is a private method and not intended to be called directly by users.
        Instead, users should call the public method :meth:`reconstruct` to reconstruct the raster data.
        The doc is provided here for developers/AI assistants who want to understand the implementation of the raster reconstruction.

        .. note::

            There are two options to reconstruct the raster data.

            Option 1: Use spatial tree (use_spatial_tree=True)
                1. Get the partitioning polygons at the `from_time`.
                2. Reconstruct the partitioning polygons to the `to_time`.
                3. Initialize the output raster and get the output sample points that are inside the partitioning polygons at `to_time`.
                4. Reconstruct the input sample points to the `to_time`.
                5. Use spatial tree to find the nearest input sample point for each output sample point and get the value.

            Option 2: Use raster interpolation (use_spatial_tree=False)
                1. Get the partitioning polygons at the `from_time`.
                2. Reconstruct the partitioning polygons to the `to_time`.
                3. Initialize the output raster and get the output sample points that are inside the partitioning polygons at `to_time`.
                4. Reverse reconstruct the output sample points to the `from_time` and query the raster values at those points.

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
        fill_value : Union[float, int, str, tuple]
            The value to be used for regions outside of the static polygons at ``to_time``.
            The fill value must have been parsed, see :meth:`_parse_fill_value`.
        threads : int, default 1
            Number of threads to use for certain computationally heavy routines.
        use_spatial_tree : bool, default False
            Whether to use a spatial tree for faster feature lookup.

        Returns
        -------
        Raster
            The reconstructed raster. Areas outside of the partitioning polygons will be filled with invalid/default values.
        """
        assert (
            self.is_normalized()
        ), "Raster data must be normalized before calling Raster._reconstruct_impl()."

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
            return self

        # Get the partitioning polygons at the `from_time` and their corresponding features(we need the features to query plate IDs, or other attributes).
        # Only get those polygons that are valid at both the `from_time` and `to_time`,
        # and if they intersect with the raster data extent(not global raster), they will be cut to fit into the extent.
        # This will make sure the reconstructed raster doesn't contain interpolated data which shouldn't be there.
        # For example, if the orignal raster covers only half of Australia, you may not want the interpolated data to show up in the other half of Australia after reconstruction.
        if partitioning_features is None:
            assert (
                self.plate_reconstruction is not None
            ), "A valid PlateReconstruction object is required here!"
            partitioning_features = self.plate_reconstruction.static_polygons

        # Accept all supported forms (filenames, FeatureCollection, Feature, lists)
        # and normalize them before feature-level operations.
        partitioning_features = load_feature_collection(partitioning_features)
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

        # We also need to reconstruct the partitioning polygons to the `to_time`.
        (
            partitioning_polygons_at_to_time,
            features_of_partitioning_polygons_at_to_time,
        ) = self._get_partitioning_polygons_at_to_time(
            partitioning_polygons_at_from_time,
            features_of_partitioning_polygons_at_from_time,
            from_time,
            to_time,
            rotation_model,
        )

        # Initialize the output raster and get the output sample points that are inside the partitioning polygons at `to_time`.
        (
            output_raster,
            output_raster_extent,
            output_sample_points_lons,
            output_sample_points_lats,
            output_sample_points_pids,
            output_sample_points_row_col_indices,
        ) = self._get_output_raster_and_sample_points(
            partitioning_polygons_at_to_time,
            features_of_partitioning_polygons_at_to_time,
            fill_value=fill_value,
        )

        if not use_spatial_tree:
            logger.debug(
                f"Using raster interpolation to reconstruct the raster data from time {from_time} to {to_time}."
            )
            (
                reverse_reconstructed_lons,
                reverse_reconstructed_lats,
            ) = self._reverse_reconstruct_data_points(
                from_time=from_time,
                to_time=to_time,
                lons=output_sample_points_lons,
                lats=output_sample_points_lats,
                plate_ids=output_sample_points_pids,
                rotation_model=rotation_model,
            )

            values = self.interp(
                lons=reverse_reconstructed_lons,
                lats=reverse_reconstructed_lats,
                method=InterpMethod.NEAREST,
            )

            # The interpolated values may contain a few old fill value, which is a tiny problem by itself(just a small number near the polgyon boundaries).
            # However, the issue will be amplified when use the method `Raster.fill_gaps()` to fill the gaps in the raster, which will propagate the old fill value to a much larger area.
            # We need to replace the old fill value with the new fill value.
            if self.fill_value is not None and self.fill_value != fill_value:
                logger.debug(
                    f"Replacing fill value {self.fill_value} with new fill value {fill_value}."
                )
                Raster._replace_values(
                    values,
                    self.fill_value,
                    fill_value,
                )

            output_raster[
                output_sample_points_row_col_indices[:, 0],
                output_sample_points_row_col_indices[:, 1],
            ] = values
        else:
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

            # build a KDTree from the reconstructed original raster data points
            # and query the nearest neighbor for the sample points of output raster
            reconstructed_original_sample_points_vecs = self._lat_lon_to_vector(
                reconstructed_original_sample_point_lat_lon_array[:, 0],
                reconstructed_original_sample_point_lat_lon_array[:, 1],
                degrees=True,
            )
            # Build a KDTree from the reconstructed original raster data points
            tree = _cKDTree(reconstructed_original_sample_points_vecs)

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
            src_rows = original_sample_point_row_col_index_array[indices, 0]
            src_cols = original_sample_point_row_col_index_array[indices, 1]
            output_raster[
                output_sample_points_row_col_indices[:, 0],
                output_sample_points_row_col_indices[:, 1],
            ] = self.data[src_rows, src_cols]

        if self.is_global():
            ret_raster_obj = Raster(
                data=output_raster,
                extent=(-180, 180, -90, 90),
                time=to_time,
            )
        else:
            ret_raster_obj = Raster(
                data=output_raster,
                extent=output_raster_extent,
                time=to_time,
            )
        assert (
            self.plate_reconstruction is not None
        ), "The PlateReconstruction object in this Raster object should not be None."
        ret_raster_obj.plate_reconstruction = self.plate_reconstruction.copy()
        if rotation_model is not None:
            ret_raster_obj.plate_reconstruction.rotation_model = rotation_model
        ret_raster_obj.fill_value = fill_value
        return ret_raster_obj

    def _get_output_raster_and_sample_points(
        self,
        partitioning_polygons_at_to_time,
        features_of_partitioning_polygons_at_to_time,
        *,
        fill_value,
    ) -> Tuple[np.ndarray, tuple, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get the output raster and the sample points that are inside the partitioning polygons at `to_time`.

        Parameters
        ----------
        partitioning_polygons_at_to_time : list[pygplates.PolygonOnSphere]
            A list of partitioning polygons at the `to_time`.
        features_of_partitioning_polygons_at_to_time : list[pygplates.Feature]
            A list of features corresponding to the partitioning polygons at the `to_time`.
        fill_value : float, int, str, or tuple
            The value to be used for regions outside of the partitioning polygons at ``to_time``.
            The fill value must have been parsed, see :meth:`_parse_fill_value`.

        Returns
        -------
        output_raster : np.ndarray
            A 2D array of the output raster, initialized with the specified fill value.
        output_raster_extent : tuple
            The extent of the output raster.
        output_sample_points_lons : np.ndarray
            A 1D array of longitudes of the sample points that are inside the partitioning polygons at `to_time`.
        output_sample_points_lats : np.ndarray
            A 1D array of latitudes of the sample points that are inside the partitioning polygons at `to_time`.
        output_sample_points_pids : np.ndarray
            A 1D array of plate IDs corresponding to the sample points that are inside the partitioning polygons at `to_time`.
        output_sample_points_row_col_indices : np.ndarray
            A 1D array of tuples, where each tuple contains the row and column indices of the sample points in the output raster that are inside the partitioning polygons at `to_time`.

        """
        # If the original raster extent is global, the output raster extent will also be global.
        if self.is_global():
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
        output_raster = np.full(
            output_shape,
            fill_value,
            dtype=self.data.dtype,
        )

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
        pids = [
            f.get_reconstruction_plate_id()  # type: ignore
            for f in features_of_partitioning_polygons_at_to_time
        ]
        assert len(shapely_geoms) == len(
            pids
        ), "The number of partitioning polygons and their corresponding plate IDs must be the same."
        ny, nx = output_lon_mesh.shape
        minx, maxx, miny, maxy = output_raster_extent
        polygon_mask = _rasterize(
            shapes=zip(shapely_geoms, pids),
            out_shape=(ny, nx),
            fill=0,
            dtype=np.uint32,
            merge_alg=MergeAlg.replace,
            transform=_from_bounds(minx, miny, maxx, maxy, nx, ny),
        )
        assert (
            polygon_mask is not None
        ), "Rasterization of partitioning polygons failed. This should never happen."
        polygon_mask = np.flipud(polygon_mask)

        # Extract output sample points in a vectorized way for better performance.
        rows, cols = np.where(polygon_mask > 0)
        output_sample_points_lons = output_lon_mesh[rows, cols]
        output_sample_points_lats = output_lat_mesh[rows, cols]
        output_sample_points_pids = polygon_mask[rows, cols]
        output_sample_points_row_col_indices = np.column_stack((rows, cols))

        return (
            output_raster,
            output_raster_extent,
            output_sample_points_lons,
            output_sample_points_lats,
            output_sample_points_pids,
            output_sample_points_row_col_indices,
        )

    def _get_partitioning_polygons_at_to_time(
        self,
        partitioning_polygons_at_from_time: list[_PolygonOnSphere],
        features_of_partitioning_polygons_at_from_time: list[_Feature],
        from_time: float,
        to_time: float,
        rotation_model: _RotationModel,
    ) -> Tuple[list[_PolygonOnSphere], list[_Feature]]:
        """Get the partitioning polygons at the `to_time` by reconstructing the partitioning polygons at the `from_time`."""
        polygon_feature_collection = _FeatureCollection()
        for p, f in zip(
            partitioning_polygons_at_from_time,
            features_of_partitioning_polygons_at_from_time,
        ):
            feature = _Feature()
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
        features_of_partitioning_polygons = []
        for rfg in rfgs:
            geom = rfg.get_reconstructed_geometry()
            if geom is not None and isinstance(geom, _PolygonOnSphere):
                polygons_at_to_time.append(geom)
                features_of_partitioning_polygons.append(rfg.get_feature())
        return polygons_at_to_time, features_of_partitioning_polygons

    def _get_plate_id_grid(
        self,
        lon_mesh: np.ndarray,
        lat_mesh: np.ndarray,
        partitioning_polygons: list[_PolygonOnSphere],
        features_of_partitioning_polygons: list[_Feature],
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

    def _reverse_reconstruct_data_points(
        self,
        *,
        from_time: float,
        to_time: float,
        lons: np.ndarray,
        lats: np.ndarray,
        plate_ids: np.ndarray,
        rotation_model: _RotationModel,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Reverse reconstruct the raster data points from `to_time` to `from_time` using the given rotation model and plate IDs."""
        assert (
            lons.shape == lats.shape == plate_ids.shape
        ), "Input arrays must have the same shape."

        assert np.all(plate_ids > 0), "All plate IDs must be positive integers."
        assert (
            self.plate_reconstruction is not None
        ), "A valid PlateReconstruction object is required!"

        if rotation_model is None:
            rotation_model = self.plate_reconstruction.rotation_model
        assert rotation_model is not None, "A valid RotationModel object is required!"

        unique_plate_ids, inv = np.unique(plate_ids, return_inverse=True)

        # Compute one finite rotation per unique plate, then broadcast by inverse index.
        rotations_dict: dict[int, np.ndarray] = {}
        for plate_id in unique_plate_ids:
            if (
                not math.isclose(to_time, 0.0)
                and self.plate_reconstruction is not None
                and self.plate_reconstruction.rotation_model is not None
                and rotation_model is not self.plate_reconstruction.rotation_model
            ):
                # If to_time is not 0 and a different rotation model is used,
                # rotate from to_time to present-day with the provided model,
                # then from present-day to from_time with the original model.
                assert self.plate_reconstruction.rotation_model is not None
                rotation_to_present_day = rotation_model.get_rotation(
                    to_time=0.0,
                    from_time=float(to_time),
                    moving_plate_id=int(plate_id),
                )
                rotation_to_from_time = (
                    self.plate_reconstruction.rotation_model.get_rotation(
                        to_time=float(from_time),
                        from_time=0.0,
                        moving_plate_id=int(plate_id),
                    )
                )
                rotation = rotation_to_from_time * rotation_to_present_day
            else:
                rotation = rotation_model.get_rotation(
                    to_time=float(from_time),
                    from_time=float(to_time),
                    moving_plate_id=int(plate_id),
                )

            if not isinstance(rotation, _FiniteRotation):
                raise ValueError(f"No rotation found for plate ID: {plate_id}")

            lat, lon, angle = rotation.get_lat_lon_euler_pole_and_angle_degrees()
            angle = np.deg2rad(angle)
            axis_vec = self._lat_lon_to_vector(lat, lon, degrees=True)
            rotations_dict[int(plate_id)] = np.asarray(axis_vec, dtype=float) * float(
                angle
            )

        rotation_vectors = np.array(
            [rotations_dict[int(pid)] for pid in unique_plate_ids]
        )
        combined_rotations = _Rotation.from_rotvec(rotation_vectors[inv])

        point_vecs = self._lat_lon_to_vector(
            lats,
            lons,
            degrees=True,
        )
        rotated_vecs = combined_rotations.apply(point_vecs)

        # Convert vectors back to (lat, lon) in degrees.
        rotated_lats = np.rad2deg(np.arcsin(np.clip(rotated_vecs[:, 2], -1.0, 1.0)))
        rotated_lons = np.rad2deg(np.arctan2(rotated_vecs[:, 1], rotated_vecs[:, 0]))

        return rotated_lons, rotated_lats

    def _reconstruct_raster_data_points(
        self,
        from_time: float,
        to_time: float,
        m_lons: np.ndarray,
        m_lats: np.ndarray,
        plate_id_grid: np.ndarray,
        rotation_model: _RotationModel,
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
            assert (
                self.plate_reconstruction is not None
            ), "A valid PlateReconstruction object is required!"
            rotation_model = self.plate_reconstruction.rotation_model
        assert rotation_model is not None, "A valid RotationModel object is required!"

        # Reconstruct only the points with valid plate IDs.
        valid_mask = plate_id_grid != -1
        rows, cols = np.where(valid_mask)
        valid_plate_ids = plate_id_grid[valid_mask]

        # Keep output shape compatible with downstream code.
        if valid_plate_ids.size == 0:
            return np.empty((0, 2), dtype=float), np.empty((0, 2), dtype=int)

        unique_plate_ids, inv = np.unique(valid_plate_ids, return_inverse=True)

        # Compute one finite rotation per unique plate, then broadcast by inverse index.
        rotations_dict: dict[int, np.ndarray] = {}
        for plate_id in unique_plate_ids:
            if (
                not math.isclose(from_time, 0.0)
                and self.plate_reconstruction is not None
                and self.plate_reconstruction.rotation_model is not None
                and rotation_model is not self.plate_reconstruction.rotation_model
            ):
                # If from_time is not 0 and a different rotation model is used,
                # rotate from from_time to present-day with the original model,
                # then from present-day to to_time with the provided model.
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

            if not isinstance(rotation, _FiniteRotation):
                raise ValueError(f"No rotation found for plate ID: {plate_id}")

            lat, lon, angle = rotation.get_lat_lon_euler_pole_and_angle_degrees()
            angle = np.deg2rad(angle)
            axis_vec = self._lat_lon_to_vector(lat, lon, degrees=True)
            rotations_dict[int(plate_id)] = np.asarray(axis_vec, dtype=float) * float(
                angle
            )

        rotation_vectors = np.array(
            [rotations_dict[int(pid)] for pid in unique_plate_ids]
        )
        combined_rotations = _Rotation.from_rotvec(rotation_vectors[inv])

        point_vecs = self._lat_lon_to_vector(
            m_lats[valid_mask],
            m_lons[valid_mask],
            degrees=True,
        )
        rotated_vecs = combined_rotations.apply(point_vecs)

        # Convert vectors back to (lat, lon) in degrees.
        rotated_lats = np.rad2deg(np.arcsin(np.clip(rotated_vecs[:, 2], -1.0, 1.0)))
        rotated_lons = np.rad2deg(np.arctan2(rotated_vecs[:, 1], rotated_vecs[:, 0]))

        reconstructed_lat_lon_array = np.column_stack((rotated_lats, rotated_lons))
        row_col_index_array = np.column_stack((rows, cols))

        assert len(reconstructed_lat_lon_array) == len(row_col_index_array)
        return reconstructed_lat_lon_array, row_col_index_array

    def _partition_raster_data_points(
        self,
        lon_mesh: np.ndarray,
        lat_mesh: np.ndarray,
        polygons: list[_PolygonOnSphere],
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
        partitioning_features: Union[_FeatureCollection, List[_Feature]],
        *,
        from_time: float,
        to_time: float,
    ) -> Tuple[List[_PolygonOnSphere], List[_Feature]]:
        """Return the partitioning polygons at `from_time` and their corresponding features.

        - The polygons must be valid at both `from_time` and `to_time`.
        - The polygons will be reconstructed to the `from_time` if `from_time` is greater than 0(not present-day).
        - The polygons will be cut to fit into the extent of the raster. If the raster is global, the cutting step is unnecessary.
        """
        # The polygons must be valid at both from_time and to_time.
        valid_partitioning_features = [
            f
            for f in partitioning_features
            if f.is_valid_at_time(from_time) and f.is_valid_at_time(to_time)
        ]

        partitioning_polygons: list[_PolygonOnSphere] = []
        associated_features: list[_Feature] = []
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
                if geom is not None and isinstance(geom, _PolygonOnSphere):
                    partitioning_polygons.append(geom)
                    associated_features.append(rfg.get_feature())
        else:
            # If from_time is 0, we don't need to reconstruct the partitioning polygons, just use the polygons from the valid partitioning features directly.
            for f in valid_partitioning_features:
                for geom in f.get_all_geometries():
                    if geom is not None and isinstance(geom, _PolygonOnSphere):
                        partitioning_polygons.append(geom)
                        associated_features.append(f)

        # If the raster is global, no need to cut the partitioning polygons.
        if self.is_global():
            return partitioning_polygons, associated_features

        # Otherwise, define the rectangular extent polygon to cut the partitioning polygons
        left, right, bottom, top = self.extent
        extent_polygon = _PolygonOnSphere((lat, lon) for lon, lat in [(left, bottom), (left, top), (right, top), (right, bottom)])  # type: ignore
        extent_feature = _Feature()
        extent_feature.set_geometry(extent_polygon)
        partitioner = (
            pygplates.PlatePartitioner(  # pyright: ignore[reportAttributeAccessIssue]
                _FeatureCollection([extent_feature]),
                _RotationModel([]),
            )
        )

        # cut the partitioning polygons with the extent polygon and return the new polygons and their corresponding features
        inside_geometries = []
        outside_geometries = []
        polygons_within_extent: list[_PolygonOnSphere] = []
        polygons_within_extent_features: list[_Feature] = []
        for polygon, feature in zip(partitioning_polygons, associated_features):
            partitioner.partition_geometry(polygon, inside_geometries, outside_geometries)  # type: ignore
            for inside_geom in inside_geometries:
                if isinstance(inside_geom, _PolygonOnSphere):
                    polygons_within_extent.append(inside_geom)
                else:
                    # The PlatePartitioner may cut a polygon into polylines. Convert the polylines into polygons by connecting the endpoints of the polylines to form a closed polygon.
                    polygons_within_extent.append(_PolygonOnSphere(inside_geom))
                polygons_within_extent_features.append(feature)
        assert len(polygons_within_extent) == len(polygons_within_extent_features)
        return polygons_within_extent, polygons_within_extent_features

    def _query_by_KDTree(
        self,
        lons: np.ndarray,
        lats: np.ndarray,
        region_of_interest: float = np.inf,
        pointwise: bool = True,
    ):
        """Given a list of location coordinates, return the nearest valid raster values at these locations.
        If `region_of_interest` is given, don't consider the cells out of the region of interest.

        Parameters
        ----------
        lons: list
            a list of longitudes of the location coordinates
        lats: list
            a list of latitude of the location coordinates
        region_of_interest: float
            the radius of the region of interest in km. This is the great arch length.
        pointwise : bool
            If True (default), `lons`/`lats` are treated as paired query
            points (lons[i], lats[i]) — e.g. sampling at scattered station
            locations. Output shape: (N,) or (N, C).

            If False, `lons`/`lats` define a new rectangular grid (outer
            product) — e.g. resampling the whole image onto a new lat/lon
            mesh. Output shape: (len(lats), len(lons)) or
            (len(lats), len(lons), C).

        Returns
        -------
        list
            a list of grid values for the given locations.
        """

        lons_arr = np.asarray(lons, dtype=float)
        lats_arr = np.asarray(lats, dtype=float)

        if pointwise:
            if lons_arr.shape != lats_arr.shape:
                raise ValueError(
                    "lons and lats must have the same shape for pointwise query"
                )
            output_shape = (lons_arr.size,)
            query_lons = lons_arr.ravel()
            query_lats = lats_arr.ravel()
        else:
            if lons_arr.ndim != 1 or lats_arr.ndim != 1:
                raise ValueError(
                    "When pointwise is False, lons and lats must be 1D arrays."
                )
            output_shape = (lats_arr.size, lons_arr.size)
            lon_array, lat_array = np.meshgrid(lons_arr, lats_arr)
            query_lons = lon_array.ravel()
            query_lats = lat_array.ravel()

        result_dtype = np.result_type(self.data.dtype, float)

        if self._spatial_cKDTree is None or self._spatial_cKDTree_values is None:
            # The tree has not been built yet. Build the spatial tree.
            masked_data = np.ma.asarray(self.masked_data)
            mask = np.ma.getmaskarray(masked_data)

            if masked_data.ndim == 2:
                valid_mask = ~mask
            else:
                valid_mask = ~np.any(mask, axis=2)

            if not np.any(valid_mask):
                empty_shape = output_shape
                if self.data.ndim == 3:
                    empty_shape = empty_shape + (self.data.shape[2],)
                return np.full(empty_shape, np.nan, dtype=result_dtype)

            lon_grid, lat_grid = np.meshgrid(self.lons, self.lats)
            valid_grid_points = self._lat_lon_to_vector(
                lat_grid[valid_mask], lon_grid[valid_mask], degrees=True
            )
            valid_grid_points = np.atleast_2d(valid_grid_points)
            self._spatial_cKDTree = _cKDTree(valid_grid_points)
            self._spatial_cKDTree_values = np.asarray(self.data)[valid_mask]

        query_points = self._lat_lon_to_vector(query_lats, query_lons, degrees=True)

        if region_of_interest < 0:
            raise ValueError("region_of_interest must be non-negative.")

        distance_upper_bound = np.inf
        if np.isfinite(region_of_interest):
            distance_upper_bound = 2 * math.sin(
                region_of_interest
                / pygplates.Earth.mean_radius_in_kms  # pyright: ignore[reportAttributeAccessIssue]
                / 2.0
            )

        _, indices = self._spatial_cKDTree.query(
            query_points, k=1, distance_upper_bound=distance_upper_bound
        )
        indices = np.asarray(indices, dtype=np.int_)
        if indices.ndim == 0:
            indices = indices.reshape(1)
        values = self._spatial_cKDTree_values
        assert values is not None

        valid_hits = indices < len(values)

        if values.ndim == 1:
            result = np.full(indices.shape, np.nan, dtype=result_dtype)
            result[valid_hits] = values[indices[valid_hits]]
            return result.reshape(output_shape)

        channel_count = values.shape[1]
        result = np.full((indices.size, channel_count), np.nan, dtype=result_dtype)
        result[valid_hits] = values[indices[valid_hits]]
        return result.reshape(output_shape + (channel_count,))

    def _get_rotation_model_with_a_different_default_anchor_plate_id(
        self, anchor_plate_id: int
    ) -> Union[_RotationModel, None]:
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
            return _RotationModel(
                self.plate_reconstruction.rotation_model,
                default_anchor_plate_id=anchor_plate_id,
            )
        return None

    def _validate_fill_gaps_data(self, data: np.ndarray) -> bool:
        if data.ndim == 2:
            return False
        if data.ndim == 3 and data.shape[2] in (3, 4):
            return True
        raise ValueError(
            "fill_gaps supports either 2D scalar rasters or 3D RGB/RGBA rasters."
        )

    def _build_fill_gaps_valid_mask(
        self,
        data: np.ndarray,
        *,
        is_image: bool,
        valid_mask,
        invalid_value,
    ) -> np.ndarray:
        ny, nx = data.shape[:2]
        finite_mask = np.isfinite(data).all(axis=2) if is_image else np.isfinite(data)

        invalid_mask = np.zeros((ny, nx), dtype=bool)
        if invalid_value is not None:
            if is_image:
                invalid_arr = np.asarray(invalid_value)
                if invalid_arr.ndim == 0:
                    invalid_mask = np.all(data == invalid_arr.item(), axis=2)
                else:
                    invalid_arr = invalid_arr.ravel()
                    if invalid_arr.size != data.shape[2]:
                        raise ValueError(
                            "For RGB/RGBA rasters, invalid_value must be a scalar "
                            "or have one value per channel."
                        )
                    invalid_mask = np.all(
                        data == invalid_arr.reshape((1, 1, -1)), axis=2
                    )
            else:
                invalid_mask = data == invalid_value

        inferred_valid_mask = finite_mask & (~invalid_mask)
        if valid_mask is None:
            return inferred_valid_mask

        user_valid_mask = np.asarray(valid_mask, dtype=bool)
        if user_valid_mask.shape != (ny, nx):
            raise ValueError(
                f"valid_mask shape mismatch: expected {(ny, nx)}, got {user_valid_mask.shape}."
            )
        return user_valid_mask & inferred_valid_mask

    def _fill_gaps_with_gmt(
        self,
        data: np.ndarray,
        gap_mask: np.ndarray,
        is_image: bool,
        method,
    ) -> np.ndarray:
        if is_image:
            raise ValueError("use_gmt=True is only supported for 2D scalar rasters.")

        try:
            import pygmt
        except Exception as exc:
            raise ImportError(
                "PyGMT is required when use_gmt=True, but it could not be imported."
            ) from exc

        gmt_fill = getattr(pygmt, "grdfill", None)
        if gmt_fill is None:
            gmt_fill = getattr(pygmt, "fillgrd", None)
        if gmt_fill is None:
            raise AttributeError(
                "PyGMT does not provide grdfill/fillgrd in this environment."
            )

        method_l = str(method).lower()
        from packaging.version import parse as _parse_version

        if _parse_version(pygmt.__version__) < _parse_version("0.16.0"):
            fill_option_map = {
                "nearest": {"neighborfill": True},
                "linear": {"splinefill": True},
                "cubic": {"splinefill": True},
                "spline": {"splinefill": True},
            }
        else:
            fill_option_map = {
                "nearest": {"neighbor_fill": True},
                "linear": {"spline_fill": True},
                "cubic": {"spline_fill": True},
                "spline": {"spline_fill": True},
            }
        fill_kwargs = fill_option_map.get(method_l)
        if fill_kwargs is None:
            raise ValueError(
                "Unsupported method for use_gmt=True. "
                "Use one of: nearest, linear, cubic, spline."
            )

        gmt_input = np.asarray(data, dtype=np.float64).copy()
        gmt_input[gap_mask] = np.nan
        gmt_grid = xr.DataArray(
            gmt_input,
            dims=("lat", "lon"),
            coords={"lat": self.lats, "lon": self.lons},
            name=self._data_var_name if self._data_var_name else "z",
        )

        try:
            gmt_filled = gmt_fill(grid=gmt_grid, **fill_kwargs)
        except TypeError:
            legacy_mode_map = {
                "nearest": "n",
                "linear": "s",
                "cubic": "s",
                "spline": "s",
            }
            gmt_filled = gmt_fill(grid=gmt_grid, mode=legacy_mode_map[method_l])

        filled_data = np.asarray(data, dtype=np.float64).copy()
        gmt_filled = np.asarray(gmt_filled)
        filled_data[gap_mask] = gmt_filled[gap_mask]
        return filled_data

    def _fill_gaps_with_griddata(
        self,
        data: np.ndarray,
        gap_mask: np.ndarray,
        effective_valid_mask: np.ndarray,
        is_image: bool,
        method,
    ) -> np.ndarray:
        lon_grid, lat_grid = np.meshgrid(self.lons, self.lats)
        valid_points = np.column_stack(
            (lon_grid[effective_valid_mask], lat_grid[effective_valid_mask])
        )

        from scipy.interpolate import griddata

        filled_data = np.asarray(data, dtype=np.float64).copy()
        if is_image:
            for channel_idx in range(data.shape[2]):
                channel = data[:, :, channel_idx]
                valid_values = channel[effective_valid_mask]
                interpolated = griddata(
                    valid_points,
                    valid_values,
                    (lon_grid, lat_grid),
                    method=method,
                )
                if np.any(np.isnan(interpolated[gap_mask])):
                    nearest = griddata(
                        valid_points,
                        valid_values,
                        (lon_grid, lat_grid),
                        method="nearest",
                    )
                    interpolated = np.where(
                        np.isnan(interpolated), nearest, interpolated
                    )
                filled_data[:, :, channel_idx][gap_mask] = interpolated[gap_mask]
        else:
            valid_values = data[effective_valid_mask]
            interpolated = griddata(
                valid_points,
                valid_values,
                (lon_grid, lat_grid),
                method=method,
            )
            if np.any(np.isnan(interpolated[gap_mask])):
                nearest = griddata(
                    valid_points,
                    valid_values,
                    (lon_grid, lat_grid),
                    method="nearest",
                )
                interpolated = np.where(np.isnan(interpolated), nearest, interpolated)
            filled_data[gap_mask] = interpolated[gap_mask]

        return filled_data

    def _fill_gaps_with_spatial_tree(
        self,
        data: np.ndarray,
        gap_mask: np.ndarray,
        effective_valid_mask: np.ndarray,
    ) -> np.ndarray:
        if not np.any(gap_mask):
            return np.asarray(data, dtype=np.float64).copy()

        lon_grid, lat_grid = np.meshgrid(self.lons, self.lats)
        query_lons = lon_grid[gap_mask]
        query_lats = lat_grid[gap_mask]

        old_data_mask = self._data_mask
        old_tree = self._spatial_cKDTree
        old_tree_values = self._spatial_cKDTree_values
        try:
            invalid_mask = ~effective_valid_mask
            if data.ndim == 3:
                invalid_mask = np.repeat(
                    invalid_mask[:, :, np.newaxis], data.shape[2], axis=2
                )

            self._data_mask = invalid_mask
            self._spatial_cKDTree = None
            self._spatial_cKDTree_values = None

            sampled = self._query_by_KDTree(
                query_lons,
                query_lats,
                region_of_interest=np.inf,
                pointwise=True,
            )
        finally:
            self._data_mask = old_data_mask
            self._spatial_cKDTree = old_tree
            self._spatial_cKDTree_values = old_tree_values

        filled_data = np.asarray(data, dtype=np.float64).copy()
        filled_data[gap_mask] = sampled
        return filled_data

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

    def _copy_from_other(self, other: "Raster"):
        _other_dict = other.__dict__.copy()
        _other_dict.pop("_plate_reconstruction", None)

        self.__dict__.update(_other_dict)

        self.data = other.data.copy()
        self.lons = other.lons.copy()
        self.lats = other.lats.copy()
        self._plate_reconstruction = copy.copy(other.plate_reconstruction)

    def _load_data_from_file(self, filename):
        """Load raster data from a file, such as a NetCDF file or an image file.

        Parameters
        ----------
        filename : str
            The path to the file containing the raster data.

        Raises
        ------
        FileNotFoundError
            If the specified file does not exist.
        ValueError
            If the file format is not supported or if there is an error reading the file.
        """
        if not Path(filename).is_file():
            raise FileNotFoundError(
                f"Raster._load_data_from_file(): File not found: {filename}."
            )

        x_data_array = load_data_array_from_netcdf(filename, self._data_var_name)
        self._load_data_from_data_array(x_data_array)
        self._filename = filename

        if self._lons is None or self._lats is None:
            logger.warning(
                f"Raster._load_data_from_file(): Could not find longitude and latitude coordinates in the file: {filename}."
                + f"The coordinate names in the file are {list(x_data_array.coords.keys())}."
            )

    def _load_data_from_data_array(self, x_data_array: xr.DataArray):
        """Load raster data from an xarray DataArray."""
        # get the data
        self._data = x_data_array.to_numpy()
        # get the lons and lats
        if self._lon_name in x_data_array.coords:
            self._lons = x_data_array.coords[self._lon_name].to_numpy()
        else:
            for name in x_data_array.coords:
                if self._is_a_common_name_for_longitude(str(name)):
                    self._lons = x_data_array.coords[name].to_numpy()
                    break
        if self._lat_name in x_data_array.coords:
            self._lats = x_data_array.coords[self._lat_name].to_numpy()
        else:
            for name in x_data_array.coords:
                if self._is_a_common_name_for_latitude(str(name)):
                    self._lats = x_data_array.coords[name].to_numpy()
                    break

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

    def _invalidate_spatial_cache(self):
        self._data_mask = None
        self._spatial_cKDTree = None
        self._spatial_cKDTree_values = None

    def _is_a_common_name_for_longitude(self, name: str) -> bool:
        """Return True if the `name` parameter is a possible common name for longitude."""
        return name.lower() in [
            "lon",
            "lons",
            "longitude",
            "x",
            "east",
            "easting",
            "eastings",
        ]

    def _is_a_common_name_for_latitude(self, name: str) -> bool:
        """Return True if the `name` parameter is a possible common name for latitude."""
        return name.lower() in [
            "lat",
            "lats",
            "latitude",
            "y",
            "north",
            "northing",
            "northings",
        ]

    def _parse_fill_value(self, fill_value):
        """Normalize ``fill_value`` to match raster dimensionality and dtype.

        Rules
        -----
        - 2D rasters return a scalar.
        - 3D rasters return a tuple with one value per channel.
        - ``None`` chooses a dtype-aware default:
          - 2D integer: dtype minimum (signed) or maximum (unsigned)
          - 2D float/complex: ``np.nan``
          - 3D integer: all zeros
          - 3D float/complex: all zeros as float
        - Color strings (for example ``"black"`` or ``"#ff0000"``) are
          supported for 3D rasters only and converted through RGBA.
        - For RGBA rasters, RGB input is accepted and an opaque alpha is appended.
        - For 3D rasters, scalar numeric input is broadcast to all channels.

        Parameters
        ----------
        fill_value : scalar, str, sequence, or None
            Fill value provided by caller.

        Returns
        -------
        scalar or tuple
            Scalar for 2D rasters, tuple for 3D rasters.

        Raises
        ------
        ValueError
            If dimensionality is unsupported, the value shape is invalid,
            or conversion to target dtype is impossible.
        TypeError
            If a color string is supplied for a 2D raster.
        """
        grid = self.data
        dtype = grid.dtype

        if grid.ndim not in (2, 3):
            raise ValueError(
                f"Unsupported raster dimensionality {grid.ndim}; expected 2D or 3D data."
            )

        expected_size = 1 if grid.ndim == 2 else int(grid.shape[2])
        dtype_kind = dtype.kind

        def _coerce_scalar(value):
            if dtype_kind == "b":
                return bool(value)
            if dtype_kind in ("i", "u"):
                if isinstance(value, numbers.Real) and not math.isfinite(value):
                    raise ValueError(
                        f"Non-finite fill_value {value} is not allowed for integer dtype {dtype}."
                    )
                info = np.iinfo(dtype)
                ivalue = int(np.rint(float(value)))
                return int(np.clip(ivalue, info.min, info.max))
            if dtype_kind == "f":
                return float(value)
            if dtype_kind == "c":
                raise ValueError(
                    "Complex raster dtypes are not supported for fill_value normalization."
                )
            raise ValueError(f"Unsupported dtype for fill_value normalization: {dtype}")

        if fill_value is None:
            if grid.ndim == 2:
                if dtype_kind == "i":
                    fill_value = np.iinfo(dtype).min
                elif dtype_kind == "u":
                    fill_value = np.iinfo(dtype).max
                else:
                    fill_value = np.nan
            else:
                base = 0 if dtype_kind in ("i", "u", "b") else 0.0
                fill_value = tuple([base] * expected_size)
            return fill_value

        if isinstance(fill_value, str):
            if grid.ndim == 2:
                raise TypeError(f"Invalid fill_value for 2D grid: {fill_value}")

            import matplotlib.colors

            try:
                rgba = np.asarray(matplotlib.colors.to_rgba(fill_value), dtype=float)
            except ValueError as exc:
                raise ValueError(f"Invalid color specification: {fill_value}") from exc

            if dtype_kind in ("i", "u"):
                max_value = np.iinfo(dtype).max
                rgba = np.rint(rgba * max_value)
            elif dtype_kind == "b":
                rgba = rgba >= 0.5

            return tuple(rgba[:expected_size])

        if grid.ndim == 2:
            if hasattr(fill_value, "__len__") and not isinstance(
                fill_value, (str, bytes)
            ):
                fill_array = np.asarray(fill_value)
                if fill_array.size != 1:
                    raise ValueError(
                        f"Shape mismatch: fill_value size: {fill_array.size}, expected: 1, grid shape: {np.shape(grid)}"
                    )
                fill_value = fill_array.reshape(-1)[0]
            return _coerce_scalar(fill_value)

        # From here on: 3D raster handling.
        if np.isscalar(fill_value):
            fill_items = [fill_value] * expected_size
        else:
            fill_array = np.asarray(fill_value).ravel()
            if expected_size == 4 and fill_array.size == 3:
                alpha_default = np.iinfo(dtype).max if dtype_kind in ("i", "u") else 1.0
                fill_array = np.append(fill_array, alpha_default)
            elif fill_array.size == 1:
                fill_array = np.repeat(fill_array, expected_size)

            if fill_array.size != expected_size:
                raise ValueError(
                    f"Shape mismatch: fill_value size: {fill_array.size}, expected: {expected_size}, grid shape: {np.shape(grid)}"
                )
            fill_items = fill_array.tolist()

        return tuple(_coerce_scalar(v) for v in fill_items)

    def _validate_data(self):
        _data = self.data
        assert isinstance(_data, np.ndarray), "Raster data must be a numpy ndarray"

        _ndim = np.ndim(_data)
        _shape = np.shape(_data)
        if _ndim not in (2, 3):
            # ndim == 2: greyscale image/grid
            # ndim == 3: colour RGB(A) image
            raise ValueError(
                f"Raster data must be 2D or 3D, but got {_ndim}D with shape {_shape}"
            )
        if _ndim == 3 and _shape[2] not in (3, 4):
            # shape[2] == 3: colour image (RGB)
            # shape[2] == 4: colour image w/ transparency (RGBA)
            raise ValueError(
                f"Image data must have 3 (RGB) or 4 (RGBA) channels, but got {_shape[2]} channels"
            )
        if _data.ndim == 3:
            # data is an RGB(A) image
            _dtype = _data.dtype
            if _dtype.kind == "i":
                _data = _data.astype("u1")
                _dtype = _data.dtype
            _min_value = np.nanmin(_data)
            _max_value = np.nanmax(_data)
            if _min_value < 0:
                raise ValueError(f"Invalid value for RGB(A) image: {_min_value}")
            if (_dtype.kind == "f" and _max_value > 1.0) or (
                _dtype.kind == "u" and _max_value > 255
            ):
                raise ValueError(f"Invalid value for RGB(A) image: {_max_value}")

    def _get_lats_lons_from_extent_origin(
        self,
        extent: Union[str, Tuple[float, float, float, float], None],
        origin: Union[Literal["lower", "upper"], str, None],
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Build latitude and longitude coordinate vectors from raster extent.

        Parameters
        ----------
        extent : {'global', None} or 4-tuple
            Spatial bounds as ``(min_lon, max_lon, min_lat, max_lat)``.
            ``None`` or ``'global'`` maps to ``(-180, 180, -90, 90)``.
        origin : {'lower', 'upper'} or None
            Desired vertical axis direction for latitude values.
            If provided, latitude bounds are reordered to match:
            - ``'lower'`` -> ascending latitudes
            - ``'upper'`` -> descending latitudes

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            ``(lats, lons)`` matching ``self.data.shape`` as
            ``(n_rows, n_cols)``.

        Raises
        ------
        TypeError
            If ``extent`` is not ``None``, ``'global'``, or a 4-element sequence.
        ValueError
            If ``origin`` is not one of ``'lower'`` or ``'upper'``.
        """
        # Normalise extent input.
        if extent is None:
            min_lon, max_lon, min_lat, max_lat = (
                -180.0,
                180.0,
                -90.0,
                90.0,
            )
        elif isinstance(extent, str):
            if extent.lower() != "global":
                raise ValueError(
                    "`extent` string must be 'global' (or use a 4-element tuple)."
                )
            min_lon, max_lon, min_lat, max_lat = (-180.0, 180.0, -90.0, 90.0)
        else:
            if len(extent) != 4:
                raise TypeError(
                    "`extent` must be a four-element tuple, 'global', or None"
                )
            min_lon, max_lon, min_lat, max_lat = (
                float(extent[0]),
                float(extent[1]),
                float(extent[2]),
                float(extent[3]),
            )

        if origin is not None:
            origin_str = str(origin).lower()
            if origin_str not in ("lower", "upper"):
                raise ValueError("`origin` must be one of: 'lower', 'upper', or None")

            if origin_str == "lower" and min_lat > max_lat:
                # increasing latitudes for 'lower' origin
                min_lat, max_lat = max_lat, min_lat
            elif origin_str == "upper" and min_lat < max_lat:
                # decreasing latitudes for 'upper' origin
                min_lat, max_lat = max_lat, min_lat

        lons = np.linspace(min_lon, max_lon, self.data.shape[1])
        lats = np.linspace(min_lat, max_lat, self.data.shape[0])

        return lats, lons

    @staticmethod
    def _replace_values(values: np.ndarray, search_value, replace_value):
        """Replace all occurrences of `search_value` in `values` with `replace_value`.

        Parameters
        ----------
        values : numpy.ndarray
            The array in which to search and replace values.
        search_value : scalar or array-like
            The value(s) to search for in `values`. Can be a scalar or an array of the same shape as `values`.
        replace_value : scalar or array-like
            The value(s) to replace `search_value` with. Can be a scalar or an array of the same shape as `values`.
        """
        if values.ndim == 1:
            # for a np array of scalar values.
            if np.isscalar(search_value) and np.isnan(search_value):
                search_mask = np.isnan(values)
            else:
                search_mask = values == search_value
            values[search_mask] = replace_value
        else:
            # for a np array of multiple values.
            search_arr = np.asarray(search_value)
            if search_arr.ndim == 0:
                search_mask = np.all(values == search_arr.item(), axis=1)
            else:
                search_mask = np.all(values == search_arr.reshape((1, -1)), axis=1)

            target_fill_arr = np.asarray(replace_value)
            if target_fill_arr.ndim == 0:
                values[search_mask, :] = target_fill_arr.item()
            else:
                values[search_mask, :] = target_fill_arr.reshape((1, -1))
        return

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
    def _write_netcdf_grid(self, *args, **kwargs):
        from .grids import write_netcdf_grid as _impl

        return _impl(*args, **kwargs)

    def _reconstruct_grid(self, *args, **kwargs):
        from .grids import reconstruct_grid as _impl

        return _impl(*args, **kwargs)
