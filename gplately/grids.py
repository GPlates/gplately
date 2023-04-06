"""Tools for working with MaskedArray, ndarray and netCDF4 rasters, as well as
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
import warnings
from multiprocessing import cpu_count

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pygplates
from cartopy.crs import PlateCarree as _PlateCarree
from cartopy.mpl.geoaxes import GeoAxes as _GeoAxes
from rasterio.enums import MergeAlg
from rasterio.features import rasterize as _rasterize
from rasterio.transform import from_bounds as _from_bounds
from scipy.interpolate import RegularGridInterpolator as _RGI
from scipy.ndimage import (
    distance_transform_edt,
    map_coordinates,
)
from scipy.spatial import cKDTree as _cKDTree
from scipy.spatial.transform import Rotation as _Rotation

from .geometry import pygplates_to_shapely
from .reconstruction import PlateReconstruction as _PlateReconstruction

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


def fill_raster(data,invalid=None):
    """Search a grid of `data` for invalid cells (i.e NaN-type entries) and fill each
    invalid cell with the value of its nearest valid neighbour. 

    Notes
    -----
    Uses `scipy`'s `distance_transform_edt` function to perform an Exact Euclidean 
    Distance Transform (EEDT). This locates the nearest valid neighbours of an invalid 
    `data` cell. 

    An optional parameter, `invalid`, is a binary ndarray with the same dimensions 
    as `data` and the following entries:

    * 1 if its corresponding entry in `data` is of NaN-type;
    * 0 if not NaN-type

    This will be used to locate nearest neighbour fill values during the Exact Euclidian 
    Distance Transform. If `invalid` is not passed to `fill_raster`, it will be created 
    for the user.

    Parameters
    ----------
    data : MaskedArray
        A MaskedArray of data that may have invalid cells (i.e. entries of type NaN).

    invalid : ndarray, optional, default=None
        An ndarray with the same shape as `data` whose elements are 1 if its corresponding 
        elements in `data` are of type `NaN`, and 0 if its corresponding entries in `data` 
        are valid. An optional parameter - this will be created for the user if it isn’t 
        provided.

    Returns
    -------
    data : ndarray
        An updated `data` array where each invalid cell has been replaced with the value 
        of its nearest valid neighbour. 
    """
    masked_array = hasattr(data, "fill_value")
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


def realign_grid(array, lons, lats):
    mask_lons = lons > 180

    # realign to -180/180
    if mask_lons.any():
        dlon = np.diff(lons).mean()
        array = np.hstack([array[:,mask_lons], array[:,~mask_lons]])
        lons  = np.hstack([lons[mask_lons] - 360 - dlon, lons[~mask_lons]])

    if lats[0] > lats[-1]:
        array = np.flipud(array)
        lats = lats[::-1]

    return array, lons, lats

def read_netcdf_grid(filename, return_grids=False, realign=False, resample=None):
    """Read a `netCDF` (.nc) grid from a given `filename` and return its data as a
    `MaskedArray`.

    Notes
    -----
    If a `resample` tuple is passed with X and Y spacings (`spacingX`, `spacingY`), 
    the gridded data in `filename` will be resampled with these resolutions. 

    By default, only the `MaskedArray` is returned to the user. However, if `return_grids` is 
    set to `True`, the `MaskedArray` will be returned along with two additional arrays 
    in a `tuple`:

    * A 1d `MaskedArray` containing the longitudes of the `netCDF` gridded data
    * A 1d `MaskedArray` containing the latitudes of the `netCDF` gridded data 
    
    Parameters
    ----------
    filename : str
        Full path to the `netCDF` raster file.
        
    return_grids : bool, optional, default=False
        If set to `True`, returns lon, lat arrays associated with the grid data.

    realign : bool, optional, default=False
        if set to `True`, realigns grid to -180/180 and flips the array if the
        latitudinal coordinates are decreasing.
        
    resample : tuple, optional, default=None
        If passed as `resample = (spacingX, spacingY)`, the given `netCDF` grid is resampled 
        with these x and y resolutions.

    Returns
    -------
    grid_z : MaskedArray
        A `MaskedArray` containing the gridded data from the supplied netCDF4 `filename`. 
        Entries' longitudes are re-aligned between -180 and 180 degrees.

    lon, lat : 1d MaskedArrays
        `MaskedArrays` encasing longitude and latitude variables belonging to the 
        supplied netCDF4 file. Longitudes are rescaled between -180 and 180 degrees. 
        An example output of `cdf_lat` is:

            masked_array(data=[-90. , -89.9, -89.8, ...,  89.8,  89.9,  90. ], mask=False, fill_value=1e+20)
    """

    def find_label(keys, labels):
        for label in labels:
            if label in keys:
                return label
        return None


    import netCDF4

    # possible permutations of lon/lat/z
    label_lon = ['lon', 'lons', 'longitude', 'x', 'east', 'easting', 'eastings']
    label_lat = ['lat', 'lats', 'latitude', 'y', 'north', 'northing', 'northings']
    label_z   = ['z', 'data', 'values']

    # add capitalise and upper case permutations
    label_lon = label_lon + [label.capitalize() for label in label_lon] + [label.upper() for label in label_lon]
    label_lat = label_lat + [label.capitalize() for label in label_lat] + [label.upper() for label in label_lat]
    label_z = label_z + [label.capitalize() for label in label_z] + [label.upper() for label in label_z]

    # open netCDF file and re-align from -180, 180 degrees
    with netCDF4.Dataset(filename, 'r') as cdf:
        keys = cdf.variables.keys()
        
        # find the names of variables
        key_z   = find_label(keys, label_z)
        key_lon = find_label(keys, label_lon)
        key_lat = find_label(keys, label_lat)

        if key_lon is None or key_lat is None:
            raise ValueError("Cannot find x,y or lon/lat coordinates in netcdf")
        if key_z is None:
            raise ValueError("Cannot find z data in netcdf")

        # extract data from cdf variables
        cdf_grid = cdf[key_z][:]
        cdf_lon  = cdf[key_lon][:]
        cdf_lat  = cdf[key_lat][:]

    if realign:
        # realign longitudes to -180/180 dateline
        cdf_grid_z, cdf_lon, cdf_lat = realign_grid(cdf_grid, cdf_lon, cdf_lat)
    else:
        cdf_grid_z = cdf_grid

    # resample
    if resample is not None:
        spacingX, spacingY = resample
        lon_grid = np.arange(cdf_lon.min(), cdf_lon.max()+spacingX, spacingX)
        lat_grid = np.arange(cdf_lat.min(), cdf_lat.max()+spacingY, spacingY)
        lonq, latq = np.meshgrid(lon_grid, lat_grid)
        original_extent = (
            cdf_lon[0],
            cdf_lon[-1],
            cdf_lat[0],
            cdf_lat[-1],
        )
        cdf_grid_z = sample_grid(
            lonq, latq,
            cdf_grid_z,
            method="nearest",
            extent=original_extent,
            return_indices=False,
        )
        cdf_lon = lon_grid
        cdf_lat = lat_grid
            
    # Fix grids with 9e36 as the fill value for nan. 
    #cdf_grid_z.fill_value = float('nan')
    #cdf_grid_z.data[cdf_grid_z.data > 1e36] = cdf_grid_z.fill_value
    
    if return_grids:
        return cdf_grid_z, cdf_lon, cdf_lat
    else:
        return cdf_grid_z
    
def write_netcdf_grid(filename, grid, extent=[-180,180,-90,90]):
    """ Write geological data contained in a `grid` to a netCDF4 grid with a specified `filename`.

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

    extent : 1D numpy array, default=[-180,180,-90,90]
        Four elements that specify the [min lon, max lon, min lat, max lat] to constrain the lat and lon 
        variables of the netCDF grid to. If no extents are supplied, full global extent `[-180, 180, -90, 90]` 
        is assumed. 

    Returns
    -------
    A netCDF grid will be saved to the path specified in `filename`.
    """
    import netCDF4
    
    nrows, ncols = np.shape(grid)
    
    lon_grid = np.linspace(extent[0], extent[1], ncols)
    lat_grid = np.linspace(extent[2], extent[3], nrows)
    
    with netCDF4.Dataset(filename, 'w', driver=None) as cdf:
        cdf.title = "Grid produced by gplately"
        cdf.createDimension('lon', lon_grid.size)
        cdf.createDimension('lat', lat_grid.size)
        cdf_lon = cdf.createVariable('lon', lon_grid.dtype, ('lon',), zlib=True)
        cdf_lat = cdf.createVariable('lat', lat_grid.dtype, ('lat',), zlib=True)
        cdf_lon[:] = lon_grid
        cdf_lat[:] = lat_grid

        # Units for Geographic Grid type
        cdf_lon.units = "degrees_east"
        cdf_lon.standard_name = 'lon'
        cdf_lon.actual_range = [lon_grid[0], lon_grid[-1]]
        cdf_lat.units = "degrees_north"
        cdf_lat.standard_name = 'lat'
        cdf_lat.actual_range = [lat_grid[0], lat_grid[-1]]

        cdf_data = cdf.createVariable('z', grid.dtype, ('lat','lon'), zlib=True)
        # netCDF4 uses the missing_value attribute as the default _FillValue
        # without this, _FillValue defaults to 9.969209968386869e+36
        cdf_data.missing_value = np.nan
        cdf_data.standard_name = 'z'
        #Ensure pygmt registers min and max z values properly
        cdf_data.actual_range = [np.nanmin(grid), np.nanmax(grid)]

        cdf_data[:,:] = grid


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
    def __init__(self, points, values, method="linear", bounds_error=False, fill_value=np.nan):
        super(RegularGridInterpolator, self).__init__(points, values, method, bounds_error, fill_value)

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
            result = self._evaluate_linear(indices,
                                           norm_distances)
        elif method == "nearest":
            result = self._evaluate_nearest(indices,
                                            norm_distances)
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
        from scipy.interpolate.interpnd import _ndim_coords_from_arrays
        ndim = len(self.grid)
        xi = _ndim_coords_from_arrays(xi, ndim=ndim)
        if xi.shape[-1] != len(self.grid):
            raise ValueError("The requested sample points xi have dimension "
                             f"{xi.shape[-1]} but this "
                             f"RegularGridInterpolator has dimension {ndim}")

        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])

        # find nans in input
        nans = np.any(np.isnan(xi), axis=-1)

        if self.bounds_error:
            for i, p in enumerate(xi.T):
                if not np.logical_and(np.all(self.grid[i][0] <= p),
                                      np.all(p <= self.grid[i][-1])):
                    raise ValueError("One of the requested xi is out of bounds "
                                     "in dimension %d" % i)
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
            with np.errstate(divide='ignore', invalid='ignore'):
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
        vslice = (slice(None),) + (None,)*(self.values.ndim - len(indices))

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
        values = 0.
        for h in hypercube:
            edge_indices, weights = zip(*h)
            weight = 1.
            for w in weights:
                weight *= w
            values += np.asarray(self.values[edge_indices]) * weight[vslice]
        return values


    def _evaluate_nearest(self, indices, norm_distances):
        """Nearest neighbour interpolator outsourced from scipy 1.9's 
        RegularGridInterpolator to ensure stable 
        operations with all versions of scipy >1.0. 
        """
        idx_res = [np.where(yi <= .5, i, i + 1)
                   for i, yi in zip(indices, norm_distances)]
        return self.values[tuple(idx_res)]


def sample_grid(
    lon,
    lat,
    grid,
    method="linear",
    extent="global",
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
    extent="global",
    origin=None,
    fill_value=None,
    threads=1,
    anchor_plate_id=0,
):
    """Reconstruct a gridded dataset to a given reconstruction time.

    Parameters
    ----------
    grid : array_like, or str
        The grid to be reconstructed. If `grid` is a filename, it will be
        loaded using `gplately.grids.read_netcdf_grid`.
    partitioning_features : valid argument to pygplates.FeaturesFunctionArgument
        Features used to partition `grid` by plate ID, usually a static
        polygons file. `partitioning_features` may be a single
        feature (`pygplates.Feature`), a feature collection
        (`pygplates.FeatureCollection`), a filename (`str`), or a (potentially
        nested) sequence of any combination of the above types.
    rotation_model : valid argument to pygplates.RotationModel
        The rotation model used to reconstruct `grid`.
        `rotation_model` may be a rotation model object
        (`pygplates.RotationModel`), a rotation feature collection
        (`pygplates.FeatureCollection`), a rotation filename
        (`str`), a rotation feature (`pygplates.Feature`), a sequence of
        rotation features, or a (potentially nested) sequence of any
        combination of the above types.
    to_time : float
        Time to which `grid` will be reconstructed.
    from_time : float, default 0.0
        Time from which to reconstruct `grid`.
    extent : tuple or "global", default "global"
        Extent of `grid`. Valid arguments are a tuple of
        the form (xmin, xmax, ymin, ymax), or the string "global",
        equivalent to (-180.0, 180.0, -90.0, 90.0).
    origin : {"upper", "lower"}, optional
        Origin of `grid` - either lower-left or upper-left. By default,
        determined from `extent`.
    fill_value : float, int, or tuple, optional
        The value to be used for regions outside of `partitioning_features`
        at `to_time`. By default (`fill_value=None`), this value will be
        determined based on the input.
    threads : int, default 1
        Number of threads to use for certain computationally heavy routines.
    anchor_plate_id : int, default 0
        ID of the anchored plate.

    Returns
    -------
    numpy.ndarray
        The reconstructed grid. Areas for which no plate ID could be
        determined from `partitioning_features` will be filled with
        `fill_value`.

    Notes
    -----
    For two-dimensional grids, `fill_value` should be a single
    number. The default value will be `np.nan` for float or
    complex types, the minimum value for integer types, and the
    maximum value for unsigned types.
    For RGB image grids, `fill_value` should be a 3-tuple RGB
    colour code or a matplotlib colour string. The default value
    will be black (0.0, 0.0, 0.0) or (0, 0, 0).
    For RGBA image grids, `fill_value` should be a 4-tuple RGBA
    colour code or a matplotlib colour string. The default fill
    value will be transparent black (0.0, 0.0, 0.0, 0.0) or
    (0, 0, 0, 0).
    """
    try:
        grid = np.array(read_netcdf_grid(grid))  # load grid data from file
    except Exception:
        grid = np.array(grid)  # copy grid data to array
    if to_time == from_time:
        return grid
    elif rotation_model is None:
        raise TypeError(
            "`rotation_model` must be provided if `to_time` != `from_time`"
        )

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
            raise TypeError(
                "Invalid fill_value for 2D grid: {}".format(fill_value)
            )
        fill_value = np.array(matplotlib.colors.to_rgba(fill_value))
        if dtype.kind == "u":
            fill_value = (fill_value * 255.0).astype("u1")
            fill_value = np.clip(fill_value, 0, 255)
        fill_value = tuple(fill_value)[:grid.shape[2]]

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
            pygplates.FeaturesFunctionArgument(
                partitioning_features
            ).get_features()
        )

    if not isinstance(rotation_model, pygplates.RotationModel):
        rotation_model = pygplates.RotationModel(rotation_model)

    lons = np.linspace(xmin, xmax, nx)
    lats = np.linspace(ymin, ymax, ny)
    m_lons, m_lats = np.meshgrid(lons, lats)

    valid_partitioning_features = [
        i for i in partitioning_features
        if i.is_valid_at_time(from_time)
        and i.is_valid_at_time(to_time)
    ]
    plate_ids = rasterise(
        features=valid_partitioning_features,
        rotation_model=rotation_model,
        key="plate_id",
        time=from_time,
        extent=extent,
        shape=grid.shape[:2],
        origin=origin,
    )
    valid_output_mask = rasterise(
        features=valid_partitioning_features,
        rotation_model=rotation_model,
        key="plate_id",
        time=to_time,
        extent=extent,
        shape=grid.shape[:2],
        origin=origin,
    ) != -1

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
            anchor_plate_id=int(anchor_plate_id),
        )
        if not isinstance(rot, pygplates.FiniteRotation):
            raise ValueError("No rotation found for plate ID: {}".format(plate))
        lat, lon, angle = rot.get_lat_lon_euler_pole_and_angle_degrees()
        angle = np.deg2rad(angle)
        vec = _lat_lon_to_vector(lat, lon, degrees=True)
        rotations_dict[plate] = vec * angle
    rotations_array = np.array(
        [rotations_dict[x] for x in unique_plate_ids]
    )[inv]
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
        if (
            "Unexpected keyword argument" in err.args[0]
            and "workers" in err.args[0]
        ):
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
    extent="global",
    origin=None,
    tessellate_degrees=0.1,
):
    """Rasterise GPlates objects at a given reconstruction time.

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
        reconstructed = None

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

        reconstructed = []
        pygplates.reconstruct(
            features,
            rotation_model,
            reconstructed,
            time,
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
            "Invalid key: {}".format(key)
            + "\nkey must be one of {}".format(valid_keys)
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
        raise ValueError(
            "Invalid value for RGB(A) image: {}".format(min_value)
        )
    if (
        (dtype.kind == "f" and max_value > 1.0)
        or (dtype.kind == "u" and max_value > 255)
    ):
        raise ValueError(
            "Invalid value for RGB(A) image: {}".format(max_value)
        )
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
        raise TypeError(
            "`extent` must be a four-element tuple, 'global', or None"
        )
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
    """A class for working with raster data.

    `Raster`'s functionalities inclue sampling data at points using spline
    interpolation, resampling rasters with new X and Y-direction spacings and
    resizing rasters using new X and Y grid pixel resolutions. NaN-type data
    in rasters can be replaced with the values of their nearest valid
    neighbours.

    Parameters
    ----------
    data : str or array-like
        The raster data, either as a filename (`str`) or array.

    plate_reconstruction : PlateReconstruction
        Allows for the accessibility of PlateReconstruction object attributes.
        Namely, PlateReconstruction object attributes rotation_model,
        topology_features and static_polygons can be used in the `Raster`
        object if called using “self.plate_reconstruction.X”, where X is the
        attribute.

    extent : str or 4-tuple, default: 'global'
        4-tuple to specify (min_lon, max_lon, min_lat, max_lat) extents
        of the raster. If no extents are supplied, full global extent
        [-180,180,-90,90] is assumed (equivalent to `extent='global'`).
        For array data with an upper-left origin, make sure `min_lat` is
        greater than `max_lat`, or specify `origin` parameter.

    resample : 2-tuple, optional
        Optionally resample grid, pass spacing in X and Y direction as a
        2-tuple e.g. resample=(spacingX, spacingY).

    time : float, default: 0.0
        The time step represented by the raster data. Used for raster
        reconstruction.

    origin : {'lower', 'upper'}, optional
        When `data` is an array, use this parameter to specify the origin
        (upper left or lower left) of the data (overriding `extent`).

    **kwargs
        Handle deprecated arguments such as `PlateReconstruction_object`,
        `filename`, and `array`.

    Attributes
    ----------
    data : ndarray, shape (ny, nx)
        Array containing the underlying raster data. This attribute can be
        modified after creation of the `Raster`.
    plate_reconstruction : PlateReconstruction
        An object of GPlately's `PlateReconstruction` class, like the
        `rotation_model`, a set of reconstructable `topology_features` and
        `static_polygons` that belong to a particular plate model. These
        attributes can be used in the `Raster` object if called using
        “self.plate_reconstruction.X”, where X is the attribute. This
        attribute can be modified after creation of the `Raster`.
    extent : tuple of floats
        Four-element array to specify [min lon, max lon, min lat, max lat]
        extents of any sampling points. If no extents are supplied, full
        global extent [-180,180,-90,90] is assumed.
    lons : ndarray, shape (nx,)
        The x-coordinates of the raster data. This attribute can be modified
        after creation of the `Raster`.
    lats : ndarray, shape (ny,)
        The y-coordinates of the raster data. This attribute can be modified
        after creation of the `Raster`.
    origin : {'lower', 'upper'}
        The origin (lower or upper left) or the data array.
    filename : str or None
        The filename used to create the `Raster` object. If the object was
        created directly from an array, this attribute is `None`.

    Methods
    -------
    interpolate(lons, lats, method='linear', return_indices=False)
        Sample gridded data at a set of points using spline interpolation.

    resample(spacingX, spacingY, overwrite=False)
        Resamples the grid using X & Y-spaced lat-lon arrays, meshed with
        linear interpolation.

    resize(resX, resY, overwrite=False)
        Resizes the grid with a specific resolution and samples points
        using linear interpolation.

    fill_NaNs(overwrite=False)
        Searches for invalid 'data' cells containing NaN-type entries and
        replaces NaNs with the value of the nearest valid data cell.

    reconstruct(time, fill_value=None, partitioning_features=None,
                threads=1, anchor_plate_id=0, inplace=False)
        Reconstruct the raster from its initial time (`self.time`) to a new
        time.
    """
    def __init__(
        self,
        data=None,
        plate_reconstruction=None,
        extent="global",
        realign=False,
        resample=None,
        time=0.0,
        origin=None,
        **kwargs
    ):
        """Constructs all necessary attributes for the raster object.

        Note: either a str path to a netCDF file OR an ndarray representing a grid must be specified.

        Parameters
        ----------
        data : str or array-like
            The raster data, either as a filename (`str`) or array.

        plate_reconstruction : PlateReconstruction
            Allows for the accessibility of PlateReconstruction object attributes. Namely, PlateReconstruction object
            attributes rotation_model, topology_featues and static_polygons can be used in the points object if called using
            “self.plate_reconstruction.X”, where X is the attribute.

        extent : str or 4-tuple, default: 'global'
            4-tuple to specify (min_lon, max_lon, min_lat, max_lat) extents
            of the raster. If no extents are supplied, full global extent
            [-180,180,-90,90] is assumed (equivalent to `extent='global'`).
            For array data with an upper-left origin, make sure `min_lat` is
            greater than `max_lat`, or specify `origin` parameter.

        resample : 2-tuple, optional
            Optionally resample grid, pass spacing in X and Y direction as a
            2-tuple e.g. resample=(spacingX, spacingY).

        time : float, default: 0.0
            The time step represented by the raster data. Used for raster
            reconstruction.

        origin : {'lower', 'upper'}, optional
            When `data` is an array, use this parameter to specify the origin
            (upper left or lower left) of the data (overriding `extent`).

        **kwargs
            Handle deprecated arguments such as `PlateReconstruction_object`,
            `filename`, and `array`.
        """
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
                "`array` keyword argument has been deprecated, "
                + "use `data` instead",
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
            raise TypeError(
                "`data` argument (or `filename` or `array`) is required"
            )
        if isinstance(data, str):
            # Filename
            self._filename = data
            self._data, lons, lats = read_netcdf_grid(
                data,
                return_grids=True,
                realign=realign,
                resample=resample,
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
                self._data, self._lons, self._lats = realign_grid(self._data, self._lons, self._lats)

        if (not isinstance(data, str)) and (resample is not None):
            self.resample(*resample, overwrite=True)

    @property
    def time(self):
        """The time step represented by the raster data."""
        return self._time

    @property
    def data(self):
        """The object's raster data.

        Can be modified.
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
        """The x-coordinates of the raster data.

        Can be modified.
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
        """The y-coordinates of the raster data.

        Can be modified.
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
        """The spatial extent (x0, x1, y0, y1) of the data.

        If y0 < y1, the origin is the lower-left corner; else the upper-left.
        """
        return (
            float(self.lons[0]),
            float(self.lons[-1]),
            float(self.lats[0]),
            float(self.lats[-1]),
        )

    @property
    def origin(self):
        """The origin of the data array, used for e.g. plotting."""
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
        """The filename of the raster file used to create the object.

        If a NumPy array was used instead, this attribute is `None`.
        """
        return self._filename

    @property
    def plate_reconstruction(self):
        """The `PlateReconstruction` object to be used for raster
        reconstruction.
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
        """ Returns a copy of the Raster
        
        Returns
        -------
        Raster
            A copy of the current Raster object
        """
        return Raster(self.data.copy(), self.plate_reconstruction, self.extent, self.time)

    def interpolate(
        self,
        lons,
        lats,
        method="linear",
        return_indices=False,
    ):
        """Interpolate a set of point data onto the gridded data provided
        to the `Raster` object.

        Parameters
        ----------
        lons, lats : array_like
            The longitudes and latitudes of the points to interpolate onto the
            gridded data. Must be broadcastable to a common shape.
        method : str or int; default: 'linear'
            The order of spline interpolation. Must be an integer in the range
            0-5. 'nearest', 'linear', and 'cubic' are aliases for 0, 1, and 3,
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
            points, only present if `return_indices=True`.

        Raises
        ------
        ValueError
            If an invalid `method` is provided.
        RuntimeWarning
            If `lats` contains any invalid values outside of the interval
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
        return sample_grid(
            lon=lons,
            lat=lats,
            grid=self,
            method=method,
            return_indices=return_indices,
        )

    def resample(self, spacingX, spacingY, method="linear", inplace=False):
        """Resample the `grid` passed to the `Raster` object with a new `spacingX` and 
        `spacingY` using linear interpolation.

        Notes
        -----
        Ultimately, `resample` changes the lat-lon resolution of the gridded data. The
        larger the x and y spacings given are, the larger the pixellation of raster data. 

        `resample` creates new latitude and longitude arrays with specified spacings in the
        X and Y directions (`spacingX` and `spacingY`). These arrays are linearly interpolated 
        into a new raster. If `inplace` is set to `True`, the respaced latitude array, longitude 
        array and raster will inplace the ones currently attributed to the `Raster` object.

        Parameters
        ----------
        spacingX, spacingY : ndarray
            Specify the spacing in the X and Y directions with which to resample. The larger 
            `spacingX` and `spacingY` are, the larger the raster pixels become (less resolved).
            Note: to keep the size of the raster consistent, set `spacingX = spacingY`; 
            otherwise, if for example `spacingX > spacingY`, the raster will appear stretched 
            longitudinally. 

        method : str or int; default: 'linear'
            The order of spline interpolation. Must be an integer in the range
            0-5. 'nearest', 'linear', and 'cubic' are aliases for 0, 1, and 3,
            respectively.

        inplace : bool, default=False
            Choose to overwrite the data (the `self.data` attribute), latitude array 
            (`self.lats`) and longitude array (`self.lons`) currently attributed to the 
            `Raster` object. 

        Returns
        -------
        Raster
            The resampled grid. If `inplace` is set to `True`, this raster overwrites the
            one attributed to `data`.
        """
        spacingX = np.abs(spacingX)
        spacingY = np.abs(spacingY)
        if self.origin == "upper":
            spacingY *= -1.0

        lons = np.arange(self.extent[0], self.extent[1]+spacingX, spacingX)
        lats = np.arange(self.extent[2], self.extent[3]+spacingY, spacingY)
        lonq, latq = np.meshgrid(lons, lats)

        data = self.interpolate(lonq, latq, method=method)
        if inplace:
            self._data = data
            self._lons = lons
            self._lats = lats
        else:
            return Raster(data, self.plate_reconstruction, self.extent, self.time)


    def resize(self, resX, resY, inplace=False, method="linear", return_array=False):
        """Resize the grid passed to the `Raster` object with a new x and y resolution 
        (`resX` and `resY`) using linear interpolation. 

        Notes
        -----
        Ultimately, `resize` "stretches" a raster in the x and y directions. The larger
        the resolutions in x and y, the more stretched the raster appears in x and y.

        It creates new latitude and longitude arrays with specific resolutions in 
        the X and Y directions (`resX` and `resY`). These arrays are linearly interpolated
        into a new raster. If `inplace` is set to `True`, the resized latitude, longitude 
        arrays and raster will inplace the ones currently attributed to the `Raster` object.

        Parameters
        ----------
        resX, resY : ndarray
            Specify the resolutions with which to resize the raster. The larger `resX` is,
            the more longitudinally-stretched the raster becomes. The larger `resY` is, the
            more latitudinally-stretched the raster becomes.

        method : str or int; default: 'linear'
            The order of spline interpolation. Must be an integer in the range
            0-5. 'nearest', 'linear', and 'cubic' are aliases for 0, 1, and 3,
            respectively.

        inplace : bool, default=False
            Choose to overwrite the data (the `self.data` attribute), latitude array 
            (`self.lats`) and longitude array (`self.lons`) currently attributed to the 
            `Raster` object. 

        return_array : bool, default False
            Return a `numpy.ndarray`, rather than a `Raster`.

        Returns
        -------
        Raster
            The resized grid. If `inplace` is set to `True`, this raster overwrites the
            one attributed to `data`.
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
            return Raster(data, self.plate_reconstruction, self.extent, self.time)


    def fill_NaNs(self, inplace=False, return_array=False):
        """Search raster for invalid ‘data’ cells containing NaN-type entries replaces them 
        with the value of their nearest valid data cells.

        Parameters
        ---------
        inplace : bool, default=False
            Choose whether to overwrite the grid currently held in the `data` attribute with
            the filled grid.

        return_array : bool, default False
            Return a `numpy.ndarray`, rather than a `Raster`.

        Returns
        --------
        Raster
            The resized grid. If `inplace` is set to `True`, this raster overwrites the
            one attributed to `data`.
        """
        data = fill_raster(self.data)
        if inplace:
            self._data = data
        if return_array:
            return data
        else:
            return Raster(data, self.plate_reconstruction, self.extent, self.time)


    def save_to_netcdf4(self, filename):
        """ Saves the grid attributed to the `Raster` object to the given `filename` (including
        the ".nc" extension) in netCDF4 format."""
        write_netcdf_grid(str(filename), self.data, self.extent)


    def reconstruct(
        self,
        time,
        fill_value=None,
        partitioning_features=None,
        threads=1,
        anchor_plate_id=0,
        inplace=False,
        return_array=False,
    ):
        """Reconstruct raster data to a given time.

        Parameters
        ----------
        time : float
            Time to which the data will be reconstructed.
        fill_value : float, int, str, or tuple, optional
            The value to be used for regions outside of the static polygons
            at `time`. By default (`fill_value=None`), this value will be
            determined based on the input.
        partitioning_features : sequence of Feature or str, optional
            The features used to partition the raster grid and assign plate
            IDs. By default, `self.plate_reconstruction.static_polygons`
            will be used, but alternatively any valid argument to
            `pygplates.FeaturesFunctionArgument` can be specified here.
        threads : int, default 1
            Number of threads to use for certain computationally heavy
            routines.
        anchor_plate_id : int, default 0
            ID of the anchored plate.
        inplace : bool, default False
            Perform the reconstruction in-place (replace the raster's data
            with the reconstructed data).
        return_array : bool, default False
            Return a `numpy.ndarray`, rather than a `Raster`.

        Returns
        -------
        Raster or np.ndarray
            The reconstructed grid. Areas for which no plate ID could be
            determined will be filled with `fill_value`.

        Raises
        ------
        TypeError
            If this `Raster` has no `plate_reconstruction` set.

        Notes
        -----
        For two-dimensional grids, `fill_value` should be a single
        number. The default value will be `np.nan` for float or
        complex types, the minimum value for integer types, and the
        maximum value for unsigned types.
        For RGB image grids, `fill_value` should be a 3-tuple RGB
        colour code or a matplotlib colour string. The default value
        will be black (0.0, 0.0, 0.0) or (0, 0, 0).
        For RGBA image grids, `fill_value` should be a 4-tuple RGBA
        colour code or a matplotlib colour string. The default fill
        value will be transparent black (0.0, 0.0, 0.0, 0.0) or
        (0, 0, 0, 0).
        """
        if time < 0.0:
            raise ValueError("Invalid time: {}".format(time))
        time = float(time)
        if self.plate_reconstruction is None:
            raise TypeError(
                "Cannot perform reconstruction - "
                + "`plate_reconstruction` has not been set"
            )
        if partitioning_features is None:
            partitioning_features = self.plate_reconstruction.static_polygons
        result =  reconstruct_grid(
            grid=self.data,
            partitioning_features=partitioning_features,
            rotation_model=self.plate_reconstruction.rotation_model,
            from_time=self.time,
            to_time=time,
            extent=self.extent,
            origin=self.origin,
            fill_value=fill_value,
            threads=threads,
            anchor_plate_id=anchor_plate_id,
        )

        if inplace:
            self.data = result
            self._time = time
            if return_array:
                return result
            return self

        if not return_array:
            result = type(self)(
                data=result,
                plate_reconstruction=self.plate_reconstruction,
                extent=self.extent,
                time=time,
                origin=self.origin,
            )
        return result


    def imshow(self, ax=None, projection=None, **kwargs):
        """Display raster data.

        A pre-existing matplotlib `Axes` instance is used if available,
        else a new one is created. The `origin` and `extent` of the image
        are determined automatically and should not be specified.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            If specified, the image will be drawn within these axes.
        projection : cartopy.crs.Projection, optional
            The map projection to be used. If both `ax` and `projection`
            are specified, this will be checked against the `projection`
            attribute of `ax`, if it exists.
        **kwargs : dict, optional
            Any further keyword arguments are passed to
            `matplotlib.pyplot.imshow` or `matplotlib.axes.Axes.imshow`,
            where appropriate.

        Returns
        -------
        matplotlib.image.AxesImage

        Raises
        ------
        ValueError
            If `ax` and `projection` are both specified, but do not match
            (i.e. `ax.projection != projection`).
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


# class TimeRaster(Raster):
#     """A class for the temporal manipulation of raster data. To be added soon!"""
#     def __init__(self, PlateReconstruction_object=None, filename=None, array=None, extent=None, resample=None):
#         raise NotImplementedError(
#             "This class has not been implemented; use `Raster` instead"
#         )
#         super(TimeRaster, self).__init__(PlateReconstruction_object)
