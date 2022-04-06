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
* TimeRaster
"""
import concurrent.futures
from multiprocessing import cpu_count
import warnings

import pygplates
import numpy as np
from rasterio.enums import MergeAlg
from rasterio.features import rasterize as _rasterize
from rasterio.transform import from_bounds as _from_bounds
from scipy.interpolate import RegularGridInterpolator as _RGI
from scipy.ndimage import distance_transform_edt, map_coordinates

from .geometry import pygplates_to_shapely
from .reconstruction import Points as _Points

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
    "TimeRaster",
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

def read_netcdf_grid(filename, return_grids=False, resample=None):
    """Read a `netCDF` (.nc) grid from a given `filename` and return its data as a
    `MaskedArray`. Re-align longitudes of raster data from -180 to 180 degrees.

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
        
    resample : tuple, optional, default=None
        If passed as `resample = (spacingX, spacingY)`, the given `netCDF` grid is resampled 
        with these x and y resolutions.

    Returns
    -------
    cdf_grid_z : MaskedArray
        A `MaskedArray` containing the gridded data from the supplied netCDF4 `filename`. 
        Entries' longitudes are re-aligned between -180 and 180 degrees.

    cdf_lon, cdf_lat : 1d MaskedArrays
        `MaskedArrays` encasing longitude and latitude variables belonging to the 
        supplied netCDF4 file. Longitudes are rescaled between -180 and 180 degrees. 
        An example output of `cdf_lat` is:

            masked_array(data=[-90. , -89.9, -89.8, ...,  89.8,  89.9,  90. ], mask=False, fill_value=1e+20)
    """
    import netCDF4
    
    # open netCDF file and re-align from -180, 180 degrees
    with netCDF4.Dataset(filename, 'r') as cdf:
        cdf_grid = cdf["z"]
        try:
            cdf_lon = cdf['lon'][:]
            cdf_lat = cdf['lat'][:]
        except:
            cdf_lon = cdf['x'][:]
            cdf_lat = cdf['y'][:]
            
        cdf_lon_mask = cdf_lon[:] > 180
        dlon = np.diff(cdf_lon[:]).mean()
        
        if cdf_lon_mask.any():
            cdf_grid_z = np.hstack([cdf_grid[:,cdf_lon_mask], cdf_grid[:,~cdf_lon_mask]])
            cdf_lon = np.hstack([cdf_lon[cdf_lon_mask]-360-dlon, cdf_lon[~cdf_lon_mask]])
        else:
            cdf_grid_z = cdf_grid[:]

    # resample
    if resample is not None:
        spacingX, spacingY = resample
        lon_grid = np.arange(cdf_lon.min(), cdf_lon.max()+spacingX, spacingX)
        lat_grid = np.arange(cdf_lat.min(), cdf_lat.max()+spacingY, spacingY)
        lonq, latq = np.meshgrid(lon_grid, lat_grid)
        interp = RegularGridInterpolator((cdf_lat, cdf_lon), cdf_grid_z, method='nearest', bounds_error=False)
        cdf_grid_z = interp((latq, lonq))
        cdf_lon = lon_grid
        cdf_lat = lat_grid
            
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
    
    with netCDF4.Dataset(filename, 'w') as cdf:
        cdf.createDimension('x', lon_grid.size)
        cdf.createDimension('y', lat_grid.size)
        cdf_lon = cdf.createVariable('x', lon_grid.dtype, ('x',), zlib=True)
        cdf_lat = cdf.createVariable('y', lat_grid.dtype, ('y',), zlib=True)
        cdf_lon[:] = lon_grid
        cdf_lat[:] = lat_grid
        cdf_lon.units = "degrees"
        cdf_lat.units = "degrees"

        cdf_data = cdf.createVariable('z', grid.dtype, ('y','x'), zlib=True)
        cdf_data[:,:] = grid


class RegularGridInterpolator(_RGI):
    """A class to sample gridded data at a set of point coordinates using either linear or nearest-neighbour 
    interpolation methods. It is a child class of `scipy`'s [`RegularGridInterpolator`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html) class. 

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
        from scipy.interpolate.interpnd import _ndim_coords_from_arrays
        method = self.method if method is None else method
        if method not in ["linear", "nearest"]:
            raise ValueError("Method '%s' is not defined" % method)

        ndim = len(self.grid)
        xi = _ndim_coords_from_arrays(xi, ndim=ndim)
        if xi.shape[-1] != len(self.grid):
            raise ValueError("The requested sample points xi have dimension "
                             "%d, but this RegularGridInterpolator has "
                             "dimension %d" % (xi.shape[1], ndim))

        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])

        if self.bounds_error:
            for i, p in enumerate(xi.T):
                if not np.logical_and(np.all(self.grid[i][0] <= p),
                                      np.all(p <= self.grid[i][-1])):
                    raise ValueError("One of the requested xi is out of bounds "
                                     "in dimension %d" % i)

        indices, norm_distances, out_of_bounds = self._find_indices(xi.T)
        if method == "linear":
            result = self._evaluate_linear(indices,
                                           norm_distances,
                                           out_of_bounds)
        elif method == "nearest":
            result = self._evaluate_nearest(indices,
                                            norm_distances,
                                            out_of_bounds)
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


def sample_grid(lon, lat, grid, extent=[-180,180,-90,90], return_indices=False, return_distances=False, method='linear'):
    """Sample point data with given `lon` and `lat` coordinates onto a `grid` using either a linear or nearest-neighbour 
    interpolation `method`.

    Notes
    -----
    If `return_indices` is set to `True`, the indices of raster points that were used as neighbouring sampling points
    are returned as an array containing two arrays:

    * array [0] is for the raster row coordinate (lat), and 
    * array [1] is for the raster column (lon) coordinate.

    An example output:

        # The first array holds the rows of the raster where point data spatially falls near.
        # The second array holds the columns of the raster where point data spatially falls near.
        sampled_indices = [array([1019, 1019, 1019, ..., 1086, 1086, 1087]), array([2237, 2237, 2237, ...,  983,  983,  983])]

    If `return_distances` is set to `True`, the distances between the raster sampling points and interpolated points 
    are returned as an array containing two arrays:

    * array [0] is for the latitudinal component of distance between the raster sampling point and the interpolated point.
    * array [1] is for the longitudinal component of distance between the raster sampling point and the interpolated point.

    An example output:

        # The first array holds the lat-component of the normal dist, while the second array holds the lon-component.
        sampled_dist = [array([5.30689060e-05, 3.47557804e-02, 1.03967049e-01, ..., 3.46526690e-02, 5.77772021e-01, 1.20890767e-01]), 
        array([4.41756600e-04, 2.89440621e-01, 8.66576791e-01, ..., 4.08341107e-01, 3.74526858e-01, 3.40690957e-01])]


    Parameters
    ----------
    lon, lat : 1d arrays
        Two arrays each specifying the longitudes and latitudes of the points to interpolate on the grid.

    grid : ndarray or MaskedArray
        An array whose elements define a grid. The number of rows corresponds to the number of point latitudes, while
        the number of columns corresponds to the number of point longitudes.

    extent : 1D numpy array, default=[-180,180,-90,90]
        A four-element array to specify the [min lon, max lon, min lat, max lat] with which to constrain lat and lon sampling
        points with respect to the given grid. If no extents are supplied, full global extent is assumed. 

    return_indices : bool, default=False
        Choose whether to return the row and column indices of points on the `grid` used to interpolate the point data. 

    return_distances : bool, default=False
        Choose whether to return the row and column normal distances between interpolated points and neighbouring 
        sampling points.

    method : str, default=’linear’
        The method of interpolation to perform. Supported are "linear" and "nearest". Assumes “linear” by default.

    Returns
    -----
    output_tuple : tuple of ndarrays
        By default, `output_tuple` has one ndarray - this holds the values of the grid data where interpolated points lie. 
        If sample point indices and/or distances have been requested (by setting `return_indices` and/or `return_distances`
        to `True`), these are returned as subsequent tuple elements. 

    Raises
    ------
    ValueError
        * Raised if the string method supplied is not “linear” or “nearest”.
        * Raised if the provided sample points for interpolation (xi) do not have the same dimensions as the supplied grid. 
        * Raised if the provided sample points for interpolation include any point out of grid bounds. Alerts user which 
        dimension (index) the point is located. Only raised if the RegularGridInterpolator attribute bounds_error is set 
        to True. If suppressed, out-of-bound points are replaced with a set fill_value. 

    """
    interpolator = RegularGridInterpolator((np.linspace(extent[2], extent[3], grid.shape[0]),
                                            np.linspace(extent[0], extent[1], grid.shape[1])),
                                            grid, method=method)

    return interpolator(np.c_[lat, lon], return_indices=return_indices, return_distances=return_distances)


def reconstruct_grid(
    grid,
    partitioning_features=None,
    rotation_model=None,
    to_time=0.0,
    from_time=0.0,
    extent="global",
    origin="upper",
    plate_ids=None,
    threads=1,
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
    origin : {"upper", "lower"}
        Origin of `grid` - either lower-left or upper-left.
    plate_ids : array_like, optional
        If a rasterised grid of plate IDs has already been obtained
        (e.g. using `gplately.grids.rasterise`), it can be provided here
        in order to avoid repeatedly rasterising `partitioning_features`.
        `plate_ids` must be of the same shape as `grid`.
    threads : int, default 1
        Number of threads to use for certain computationally heavy subroutines.

    Returns
    -------
    numpy.ndarray
        The reconstructed grid. Areas for which no plate ID could be obtained
        from `partitioning_features` will be filled with either `-1` or
        `np.nan`, depending on the dtype of `grid`.
    """
    try:
        grid = np.array(read_netcdf_grid(grid))
    except Exception:
        pass
    if to_time == from_time:
        return grid
    elif rotation_model is None:
        raise TypeError(
            "`rotation_model` must be provided if `to_time` != `from_time`"
        )

    if plate_ids is not None:
        plate_ids = np.array(plate_ids)
    if plate_ids is not None and plate_ids.shape != grid.shape:
        raise ValueError(
            "Shape mismatch: "
            + "`grid.shape` == {}, ".format(grid.shape)
            + "`plate_ids.shape` == {}".format(plate_ids.shape)
        )
    if origin.lower() not in {"lower", "upper"}:
        raise ValueError("Invalid `origin` value: {}".format(origin))
    origin = origin.lower()
    dtype = grid.dtype
    if dtype.kind in ("b", "u"):
        grid = grid.astype(int)
        dtype = grid.dtype

    if isinstance(threads, str):
        if threads.lower() in {"all", "max"}:
            threads = cpu_count()
        else:
            raise ValueError("Invalid `threads` value: {}".format(threads))
    threads = min([int(threads), cpu_count()])
    threads = max([threads, 1])

    grid = grid.squeeze()
    if grid.ndim != 2:
        raise ValueError("`grid` has invalid shape {}".format(grid.shape))
    if extent == "global":
        extent = (-180, 180, -90, 90)
    xmin, xmax, ymin, ymax = extent
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    if ymin > ymax:
        ymin, ymax = ymax, ymin
    ny, nx = grid.shape
    resx = (xmax - xmin) / nx
    resy = (ymax - ymin) / ny

    if not isinstance(partitioning_features, pygplates.FeatureCollection):
        partitioning_features = pygplates.FeatureCollection(
            pygplates.FeaturesFunctionArgument(
                partitioning_features
            ).get_features()
        )
    if not isinstance(rotation_model, pygplates.RotationModel):
        rotation_model = pygplates.RotationModel(rotation_model)

    lats = np.arange(ymin + resy * 0.5, ymax, resy)
    lons = np.arange(xmin + resx * 0.5, xmax, resx)
    lons, lats = np.meshgrid(lons, lats)
    if plate_ids is None:
        plate_ids = rasterise(
            features=partitioning_features,
            rotation_model=rotation_model,
            key="plate_id",
            time=None if to_time == 0.0 and rotation_model is None else to_time,
            extent=extent,
            shape=grid.shape,
            origin=origin,
        )
    plate_ids = plate_ids.flatten()

    unique_plate_ids = np.unique(plate_ids)
    rotations_dict = {}
    for plate in unique_plate_ids:
        if plate == -1:
            continue
        rot = rotation_model.get_rotation(
            float(from_time),
            int(plate),
            float(to_time),
        )
        if not isinstance(rot, pygplates.FiniteRotation):
            continue
        lat, lon, angle = rot.get_lat_lon_euler_pole_and_angle_degrees()
        angle = np.deg2rad(angle)
        vec = _lat_lon_to_vector(lat, lon, degrees=True)
        rotations_dict[plate] = (vec, angle)

    point_vecs = _lat_lon_to_vector(
        lats,
        lons,
        degrees=True,
        threads=threads,
    )

    rotated_vecs = np.full_like(point_vecs, np.nan)
    if threads > 1:
        executor = concurrent.futures.ThreadPoolExecutor(threads)
        plate_ids_divided = np.array_split(unique_plate_ids, threads)

        def _fill(ids, out):
            for id in ids:
                if id == -1:
                    continue
                index = plate_ids == id
                vec_subset = point_vecs[index, :]
                rotation, angle = rotations_dict[id]
                rotated = _rotate(vec_subset, rotation, angle)
                out[index] = rotated

        futures = {}
        for i in range(threads):
            args = (
                _fill,
                plate_ids_divided[i],
                rotated_vecs,
            )
            futures[executor.submit(*args)] = i
        concurrent.futures.wait(futures)
        executor.shutdown(False)
    else:
        for plate_id in unique_plate_ids:
            if plate_id == -1:
                continue
            index = plate_ids == plate_id
            vec_subset = point_vecs[index, :]
            rotation, angle = rotations_dict[plate_id]
            rotated = _rotate(vec_subset, rotation, angle)
            rotated_vecs[index] = rotated

    x = rotated_vecs[:, 0]
    y = rotated_vecs[:, 1]
    z = rotated_vecs[:, 2]
    rotated_lats, rotated_lons = _vector_to_lat_lon(
        x,
        y,
        z,
        degrees=True,
        return_array=True,
        threads=threads,
    )
    if origin == "upper":
        rotated_y = (ymax - rotated_lats) / resy
    else:
        rotated_y = (rotated_lats - ymin) / resy
    rotated_x = np.abs((rotated_lons - xmin) / resx)

    mask = plate_ids != -1
    interp_coords = np.vstack(
        (
            rotated_y.reshape((1, -1)),
            rotated_x.reshape((1, -1)),
        )
    )
    # data = np.full(rotated_lats.size, np.nan)
    if dtype.kind == "i":
        fill_value = -1
    elif dtype.kind in ("f", "c"):
        fill_value = np.nan
    else:
        fill_value = np.nan
    data = np.full(rotated_lats.size, fill_value, dtype=dtype)
    tmp = map_coordinates(
        grid,
        interp_coords[:, mask],
        mode="grid-wrap",
        order=0,
    ).squeeze()
    data[mask] = tmp
    data = data.reshape(grid.shape)
    if origin == "upper":
        data = np.flipud(data)
    return data


def rasterise(
    features,
    rotation_model=None,
    key="plate_id",
    time=None,
    resx=1.0,
    resy=1.0,
    shape=None,
    extent="global",
    origin="upper",
):
    """Rasterise GPlates objects at a given reconstruction time.

    This function is particularly useful for rasterising static polygons
    to extract a grid of plate IDs.

    Parameters
    ----------
    features : valid argument for pygplates.FeaturesFunctionArgument
        `features` may be a single `pygplates.Feature`, a
        `pygplates.FeatureCollection`, a `str` filename,
        or a (potentially nested) sequence of any combination of the
        above types.
    rotation_model : valid argument for pygplates.RotationModel, optional
        `rotation_model` may be a `pygplates.RotationModel`, a rotation
        feature collection (pygplates.FeatureCollection), a rotation filename
        (`str`), a rotation feature (`pygplates.Feature`), a sequence of
        rotation features, or a (potentially nested) sequence of any
        combination of the above types.
        Alternatively, if time not given, a rotation model is
        not usually required.
    key : str, default "plate_id"
        The value used to create the rasterised grid. May be any of
        the following values:
        - "plate_id"
        - "conjugate_plate_id"
        - "from_age"
        - "to_age"
        - "left_plate"
        - "right_plate"
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
    origin : {"upper", "lower"}
        Origin (upper-left or lower-left) of the output array.

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
    if origin.lower() not in {"upper", "lower"}:
        raise ValueError("Invalid `origin`: {}".format(origin))
    origin = origin.lower()
    valid_keys = {
        "plate_id",
        "conjugate_plate_id",
        "from_age",
        "to_age",
        "left_plate",
        "right_plate",
    }
    try:
        key = key.lower()
    except AttributeError as err:
        raise TypeError("Invalid key type: {}".format(type(key))) from err
    if key not in valid_keys:
        raise ValueError(
            "Invalid key: {}".format(key)
            + "\nkey must be one of {}".format(valid_keys)
        )

    try:
        extent = extent.lower()
    except AttributeError:
        pass
    if extent == "global":
        extent = (-180.0, 180.0, -90.0, 90.0)
    minx, maxx, miny, maxy = extent

    if shape is not None:
        lons = np.linspace(minx, maxx, shape[1], endpoint=True)
        lats = np.linspace(miny, maxy, shape[0], endpoint=True)
    else:
        lons = np.arange(minx, maxx + resx, resx)
        lats = np.arange(miny, maxy + resy, resy)
    nx = lons.size
    ny = lats.size

    if rotation_model is None:
        if time is not None:
            raise TypeError(
                "Rotation model must be provided if `time` is not `None`"
            )
        rotation_model = pygplates.RotationModel(pygplates.Feature())
        time = 0.0
    features = pygplates.FeaturesFunctionArgument(features).get_features()
    time = float(time)

    reconstructed = []
    pygplates.reconstruct(
        features,
        rotation_model,
        reconstructed,
        time,
    )
    geometries = pygplates_to_shapely(reconstructed)

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

    out = _rasterize(
        shapes=zip(geometries, values),
        out_shape=(ny, nx),
        fill=fill_value,
        dtype=dtype,
        merge_alg=MergeAlg.replace,
        transform=_from_bounds(minx, miny, maxx, maxy, nx, ny),
    )
    if origin == "lower":
        out = np.flipud(out)
    return out


rasterize = rasterise


def _lat_lon_to_vector(lat, lon, degrees=False, threads=1):
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

    if threads == 1:
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)
    else:
        n = lat.size
        step = np.ceil(n / threads).astype(np.int_)
        executor = concurrent.futures.ThreadPoolExecutor(threads)

        def _fill(out_x, out_y, out_z, first, last):
            np.multiply(
                np.cos(lat[first:last]),
                np.cos(lon[first:last]),
                out=out_x[first:last],
            )
            np.multiply(
                np.cos(lat[first:last]),
                np.sin(lon[first:last]),
                out=out_y[first:last],
            )
            np.sin(lat[first:last], out=out_z[first:last])

        futures = {}
        x = np.zeros_like(lat)
        y = np.zeros_like(x)
        z = np.zeros_like(x)
        for i in range(threads):
            args = (
                _fill,
                x,
                y,
                z,
                i * step,
                (i + 1) * step,
            )
            futures[executor.submit(*args)] = i
        concurrent.futures.wait(futures)
        executor.shutdown(False)

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
    threads=1,
):
    """Convert one or more (x, y, z) vectors (on the unit sphere) to
    (lat, lon) coordinate pairs, in degrees or radians. Optionally, use
    more than one thread.
    """
    x = np.atleast_1d(x).flatten()
    y = np.atleast_1d(y).flatten()
    z = np.atleast_1d(z).flatten()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if threads == 1:
            lat = np.arcsin(z)
            lon = np.arctan2(y, x)
            if degrees:
                lat = np.rad2deg(lat)
                lon = np.rad2deg(lon)
        else:
            n = x.size
            step = np.ceil(n / threads).astype(np.int_)
            executor = concurrent.futures.ThreadPoolExecutor(threads)

            def _fill(out_lat, out_lon, first, last, degrees=False):
                if degrees:
                    np.rad2deg(
                        np.arcsin(z[first:last]),
                        out=out_lat[first:last],
                    )
                    np.rad2deg(
                        np.arctan2(
                            y[first:last],
                            x[first:last],
                        ),
                        out=out_lon[first:last],
                    )
                else:
                    np.arcsin(z[first:last], out=out_lat[first:last])
                    np.arctan2(
                        y[first:last],
                        x[first:last],
                        out=out_lon[first:last],
                    )

            futures = {}
            lat = np.zeros_like(x)
            lon = np.zeros_like(lat)
            for i in range(threads):
                args = (
                    _fill,
                    lat,
                    lon,
                    i * step,
                    (i + 1) * step,
                    degrees,
                )
                futures[executor.submit(*args)] = i
            concurrent.futures.wait(futures)
            executor.shutdown(False)

    if lat.size == 1 and not return_array:
        lat = np.atleast_1d(np.squeeze(lat))[0]
        lon = np.atleast_1d(np.squeeze(lon))[0]
        return (lat, lon)

    lat = lat.reshape((-1, 1))
    lon = lon.reshape((-1, 1))
    return lat, lon


def _rotate(vectors, rotation, angle):
    cross = _cross_products
    dot = np.dot

    invalid_dims_err = ValueError(
        "Invalid shapes: {}, {}".format(vectors.shape, rotation.shape)
    )
    vectors = np.atleast_2d(vectors)
    rotation = np.squeeze(rotation)
    if vectors.shape[1] != 3:
        vectors = vectors.T
    if vectors.shape[1] != 3 or rotation.shape != (3,):
        raise invalid_dims_err

    angle = float(angle)

    t1 = np.cos(angle) * vectors
    t2 = np.sin(angle) * cross(rotation, vectors)
    t3 = (
        (1.0 - np.cos(angle))
        * dot(vectors, rotation.reshape((-1, 1))).reshape((-1, 1))
        * vectors
    )
    return t1 + t2 + t3


def _cross_products(a, b):
    """Cross products of a vector and a list of vectors."""
    if a.ndim == 2 and b.ndim == 1:
        return -1.0 * _cross_products(b, a)
    vec = a
    arr = b
    invalid_dims_err = ValueError(
        "Invalid dimensions: {}, {}".format(vec.ndim, arr.ndim)
    )
    if vec.ndim != 1 or arr.ndim != 2:
        raise invalid_dims_err

    if arr.shape[1] != 3:
        arr = arr.T
    if arr.shape[1] != 3:
        raise invalid_dims_err

    out = np.zeros_like(arr)
    out[:, 0] = vec[1] * arr[:, 2] - vec[2] * arr[:, 1]
    out[:, 1] = vec[2] * arr[:, 0] - vec[0] * arr[:, 2]
    out[:, 2] = vec[0] * arr[:, 1] - vec[1] * arr[:, 0]
    return out


class Raster(object):
    """A class for working with raster data. 

    `Raster`'s functionalities inclue interpolating point data on rasters using Scipy’s 
    RegularGridInterpolator, resampling rasters with new X and Y-direction spacings and 
    resizing rasters using new X and Y grid pixel resolutions. NaN-type data in rasters
    can be replaced with the values of their nearest valid neighbours. 

    Attributes
    ----------
    PlateReconstruction_object : object pointer
<<<<<<< HEAD
        A pointer to GPlately's `PlateReconstruction` object and its attributes, like the 
        `rotation_model`, a set of reconstructable `topology_featues` and `static_polygons`
        that belong to a particular plate model. These attributes can be used in the `Points` 
        object if called using “self.PlateReconstruction_object.X”, where X is the attribute.

    filename : str, default=None
        Full string path to netCDF raster file
    OR
    array : ndarray, default=None
        An array with elements that define a grid. The number of rows corresponds to the number of 
        latitudinal points, while the number of columns corresponds to the number of longitudinal 
        points.

    extent : 1D numpy array, default=None
        Four-element array to specify [min lon, max lon, min lat, max lat] extents of any sampling 
        points. If no extents are supplied, full global extent [-180,180,-90,90] is assumed. 
=======
    filename
    array
    extent
    Resample
    data 
    lons
    lats
    method

    Methods
    -------
    __init__(self, filename=None, array=None, extent=None, resample=None)
        Constructs all necessary attributes for the Raster object.
        
    _update(self)
        Allows RegularGridInterpolator attributes ((self.lats, self.lons), self.data, method='linear') and methods 
        (__call__(), or RegularGridInterpolator) to be accessible from the Raster object.
        
    interpolate(self, lons, lats, method='linear', return_indices=False, return_distances=False)
        Sample gridded data on a set of points using interpolation from RegularGridInterpolator.
        
    resample(self, spacingX, spacingY, overwrite=False)
        Resamples the grid using X & Y-spaced lat-lon arrays, meshed with linear interpolation.
        
    resize(self, resX, resY, overwrite=False)
        Resizes the grid with a specific resolution and samples points using linear interpolation.
        
    fill_NaNs(self, overwrite=False)
        Searches for invalid ‘data’ cells containing NaN-type entries and replaces NaNs with the value of the nearest
        valid data cell.
    """
    def __init__(self, PlateReconstruction_object=None, filename=None, array=None, extent=None, resample=None, time=0):
        """Constructs all necessary attributes for the raster object.

        Note: either a str path to a netCDF file OR an ndarray representing a grid must be specified. 

        Parameters
        ----------
        PlateReconstruction_object : object pointer
            Allows for the accessibility of PlateReconstruction object attributes. Namely, PlateReconstruction object 
            attributes rotation_model, topology_featues and static_polygons can be used in the points object if called using
            “self.PlateReconstruction_object.X”, where X is the attribute.

        filename : str, default=None
            Path to netCDF file
        OR
        array : ndarray, default=None
            An array with elements that define a grid. The number of rows corresponds to the number of latitudinal points, while
            the number of columns corresponds to the number of longitudinal points.

        extent : 1D numpy array, default=None
            Four-element array to specify [min lon, max lon, min lat, max lat] extents of any sampling points. If no extents are 
            supplied, full global extent [-180,180,-90,90] is assumed. 

        resample : tuple, default=None
            Optionally resample grid, pass spacing in X and Y direction as a tuple
            e.g. resample=(spacingX, spacingY)

        Returns
        -------
        __init__ generates the following attributes for the raster object:
        data : ndarray
            The grid - either a read netCDF4 file, or the ndarray supplied to __init__.

        extent : 1d array
            The [min lon, max lon, min lat, max lat] extents supplied to __init__. If not supplied, it is taken to be
            [-180,180,-90,90].
>>>>>>> raster-reconstruction

    resample : tuple, default=None
        Optionally resample grid, pass spacing in X and Y direction as a tuple
        e.g. resample=(spacingX, spacingY)

    """
    def __init__(self, PlateReconstruction_object=None, filename=None, array=None, extent=None, resample=None, time=0):
        self.PlateReconstruction_object = PlateReconstruction_object

        # we initialise an empty points object as we do not want to build this before any resampling takes place.
        self.time = float(time)

        if filename is None and array is None:
            raise ValueError("Supply either a filename or numpy array")

        elif filename and array:
            raise ValueError("Supply either a filename or numpy array")

        elif filename is not None:
            self.data, lons, lats = read_netcdf_grid(filename, return_grids=True, resample=resample)
            self.extent = [lons.min(), lons.max(), lats.min(), lats.max()]
            self.lons = lons
            self.lats = lats

        elif array is not None:
            if extent is None:
                extent = [-180,180,-90,90]
            self.data = array
            self.extent = extent
            self.lons = np.linspace(extent[0], extent[1], self.data.shape[1])
            self.lats = np.linspace(extent[2], extent[3], self.data.shape[0])

        self._update()

        if array is not None and resample is not None:
            self.resample(*resample, override=True)


    def _update(self):
        """Stores the RegularGridInterpolator object’s method for sampling gridded data at a set of 
        point coordinates. 

        Allows methods of the Raster object to access grid sampling functionalities. The gridded data 
        used is the “data” attribute - either read from a netCDF4 file, or supplied as an ndarray. 
        Points to sample are either variables of the netCDF4 file, or are generated from the “extent” 
        attribute and scaled to fit the grid “data”.
        """
        # store interpolation object
        interpolator = RegularGridInterpolator((self.lats, self.lons), self.data, method='linear')
        self._interpolator = interpolator


    def interpolate(self, lons, lats, method='linear', return_indices=False, return_distances=False):
        """Interpolate a set of point data (either linearly or through nearest-neighbour methods) 
        onto the gridded data provided to the `Raster` object. 

        Notes
        -----
        If `return_indices` is set to `True`, the indices of raster data used for interpolating the
        points are returned as an array containing two arrays:

        * array [0] is for the raster row coordinate (lat), and 
        * array [1] is for the raster column (lon) coordinate.

        An example output:

            # The first array holds the rows of the raster where point data spatially falls near.
            # The second array holds the columns of the raster where point data spatially falls near.
            sampled_indices = [array([1019, 1019, 1019, ..., 1086, 1086, 1087]), array([2237, 2237, 2237, ...,  983,  983,  983])]


        If `return_distances` is set to `True`, the distances between the raster data used for 
        interpolating the points are returned as an array containing two arrays:

        * array [0] is for the latitudinal component of distance between the raster sampling 
        point and the interpolated point.
        * array [1] is for the longitudinal component of distance between the raster sampling 
        point and the interpolated point.

        An example output:

            # The first array holds the lat-component of the normal dist, while the second array holds the lon-component.
            sampled_dist = [array([5.30689060e-05, 3.47557804e-02, 1.03967049e-01, ..., 3.46526690e-02, 5.77772021e-01, 1.20890767e-01]), 
            array([4.41756600e-04, 2.89440621e-01, 8.66576791e-01, ..., 4.08341107e-01, 3.74526858e-01, 3.40690957e-01])]
    
        Parameters
        ----------
        lons, lats : array
            1d arrays containing the longitudes and latitudes of the points to interpolate onto the
            gridded data. 

        method : str, default=’linear’
            The method of interpolation to perform. Supported are `linear` and `Nearest`. Assumes 
            `linear` interpolation if `None` provided.  

        return_indices : bool, default=False
            Choose whether to return the row and column indices of points on the `grid` used to 
            interpolate the point data. 

        return_distances : bool, default=False
            Choose whether to return the row and column normal distances between interpolated 
            points and neighbouring sampling points.

        Returns
        -------
        data_interp : tuple of ndarrays
            By default, `data_interp` has one ndarray - this holds the values of the grid data 
            where interpolated points lie. If sample point indices and/or distances have been 
            requested (by setting `return_indices` and/or `return_distances`
            to `True`), these are returned as subsequent tuple elements. 

        Raises
        ------
        ValueError
            * Raised if the string method supplied is not `linear` or `nearest`.
            * Raised if the provided lat, lon arrays generate sample points that do not have the 
            same dimensions as the 
            supplied grid. 
            * Raised if the provided lat, lon arrays generate sample points that include any point 
            out of grid bounds. 
            Alerts user which dimension (index) the point is located. 
        """
        interp = self._interpolator
        interp.values = self.data
        data_interp = interp((lats,lons), method=method, return_indices=return_indices, return_distances=return_distances)
        return data_interp


    def resample(self, spacingX, spacingY, overwrite=False):
        """Resample the `grid` passed to the `Raster` object with a new `spacingX` and 
        `spacingY` using linear interpolation.

        Notes
        -----
        Ultimately, `resample` changes the lat-lon resolution of the gridded data. The
        larger the x and y spacings given are, the larger the pixellation of raster data. 

        `resample` creates new latitude and longitude arrays with specified spacings in the
        X and Y directions (`spacingX` and `spacingY`). These arrays are linearly interpolated 
        into a new raster. If `overwrite` is set to `True`, the respaced latitude array, longitude 
        array and raster will overwrite the ones currently attributed to the `Raster` object.

        Parameters
        ----------
        spacingX, spacingY : ndarray
            Specify the spacing in the X and Y directions with which to resample. The larger 
            `spacingX` and `spacingY` are, the larger the raster pixels become (less resolved).
            Note: to keep the size of the raster consistent, set `spacingX = spacingY`; 
            otherwise, if for example `spacingX > spacingY`, the raster will appear stretched 
            longitudinally. 

        overwrite : bool, default=False
            Choose to overwrite the raster (the `self.data` attribute), latitude array 
            (`self.lats`) and longitude array (`self.lons`) currently attributed to the 
            `Raster` object. 

        Returns
        -------
        data : ndarray grid
            A new version of the raster data attributed to the `Raster` object resampled to 
            the given `spacingX` and `spacingY` spacings.
        """
        lons = np.arange(self.extent[0], self.extent[1]+spacingX, spacingX)
        lats = np.arange(self.extent[2], self.extent[3]+spacingY, spacingY)
        lonq, latq = np.meshgrid(lons, lats)

        data = self.interpolate(lonq, latq)
        if overwrite:
            self.data = data
            self.lons = lons
            self.lats = lats
            self._update()

        return data


    def resize(self, resX, resY, overwrite=False):
        """Resize the grid passed to the `Raster` object with a new x and y resolution 
        (`resX` and `resY`) using linear interpolation. 

        Notes
        -----
        Ultimately, `resize` "stretches" a raster in the x and y directions. The larger
        the resolutions in x and y, the more stretched the raster appears in x and y.

        It creates new latitude and longitude arrays with specific resolutions in 
        the X and Y directions (`resX` and `resY`). These arrays are linearly interpolated
        into a new raster. If `overwrite` is set to `True`, the resized latitude, longitude 
        arrays and raster will overwrite the ones currently attributed to the `Raster` object.

        Parameters
        ----------
        resX, resY : ndarray
            Specify the resolutions with which to resize the raster. The larger `resX` is,
            the more longitudinally-stretched the raster becomes. The larger `resY` is, the
            more latitudinally-stretched the raster becomes.

        overwrite : bool, default=False
            Choose to overwrite the raster (the `self.data` attribute), latitude array 
            (`self.lats`) and longitude array (`self.lons`) currently attributed to the 
            `Raster` object. 

        Returns
        -------
        data : meshed ndarray grid
            A new resized raster. If `overwrite` is set to `True`, this raster overwrites the
            one attributed to `data`.
        """
        # construct grid
        lons = np.linspace(self.extent[0], self.extent[1], resX)
        lats = np.linspace(self.extent[2], self.extent[3], resY)
        lonq, latq = np.meshgrid(lons, lats)

        data = self.interpolate(lonq, latq)
        if overwrite:
            self.data = data
            self.lons = lons
            self.lats = lats
            self._update()

        return data


    def fill_NaNs(self, overwrite=False):
        """Search raster for invalid ‘data’ cells containing NaN-type entries replaces them 
        with the value of their nearest valid data cells.

        Parameters
        ---------
        overwrite : bool, default=False
            Choose whether to overwrite the grid currently held in the `data` attribute with
            the filled grid.

        Returns
        --------
        data : ndarray
            An updated grid of data where each invalid cell has been replaced with the v
            alue of its nearest valid neighbour. If `overwrite` is set to `True`, this raster 
            overwrites the one attributed to `data`.
        """
        data = fill_raster(self.data)
        if overwrite:
            self.data = data

        return data


    def save_to_NetCDF4(self, filename):
        """ Saves the grid attributed to the `Raster` object to the given `filename` (including
        the ".nc" extension) in netCDF4 format."""
        write_netcdf_grid(str(filename), self.data, self.extent)


    def reconstruct(self, time):
        rotation_model = self.PlateReconstruction_object.rotation_model
        static_polygons = self.PlateReconstruction_object.static_polygons
        return reconstruct_grid(self.data, static_polygons, rotation_model, from_time=self.time, to_time=float(time), extent=self.extent)


class TimeRaster(Raster):
    """ A class for the temporal manipulation of raster data. To be added soon!
    """
    def __init__(self, PlateReconstruction_object, filename=None, array=None, extent=None, resample=None):

        super(TimeRaster, self).__init__(filename, array, extent, resample)

        self.PlateReconstruction_object = PlateReconstruction_object
