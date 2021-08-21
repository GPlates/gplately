"""A module that uses Scipy interpolation tools for working with raster/grid data.

Gridded data can be sampled at a set of point coordinates using either linear or nearest-neighbour interpolation. 
These grids can also be resampled using X and Y-direction spacing, and can be resized using given X and Y resolutions.
Grids can be searched for invalid, NaN-type data cells. These can be replaced with the values of their nearest valid neighbours. 

Classes
-------
RegularGridInterpolator
Raster
TimeRaster
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator as _RGI
from scipy.ndimage import distance_transform_edt
from .reconstruction import Points as _Points
from .tools import EARTH_RADIUS as _EARTH_RADIUS

def fill_raster(data,invalid=None):
    """Searches grid for invalid ‘data’ cells containing NaN-type entries (as indicated by ‘invalid’), locates the index 
    of the nearest valid data cell and replaces NaNs with the value of the nearest valid data cell.

    Searches a supplied data frame “data” for invalid data cells containing NaN-type entries. If these invalidities have been  
    replaced with a “fill_value” attribute before being passed into fill_raster, there may no longer be invalid cells in the 
    data. If so, the data’s NaN entries are recovered. Where there is an invalid cell in “data”, the index of the closest 
    valid data cell entry is located. This is superimposed onto “data”’s invalid cells and replaces NaNs with the nearest
    valid value.

    Parameters
    ----------
    data : ndarray
        A numpy array enclosing grid data that may have invalid cells (entries of type NaN), or formerly-invalid cells masked 
        with a fill value. 

    invalid : ndarray, optional, default=None
        A boolean-binary array with the same shape as “data”. Entries should be 1 if its corresponding entry in “data” is of 
        type NaN, and 0 if its corresponding entry in “data” is valid. Used to locate the indices of the nearest valid data 
        cells. An optional parameter; by default, this method assumes that “invalid” isn’t provided and will create it if 
        not supplied.

    Returns
    -------
    data : ndarray
        An updated grid of data where each invalid cell has been replaced with the value of its nearest valid neighbour. 
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
    """Reads in a netCDF file and re-aligns its grid, lat and lon variables from -180 to 180 degrees.

    Can optionally resample grid if given required spacing in X and Y direction. Depending on user preference, it can 
    return the grid read from the file, or the grid along with its associated lat, lon arrays.
    
    Parameters
    ----------
    filename : str
        Path to netCDF file
        
    return_grids : bool, default=False
        If set to True, optionally returns lon, lat arrays associated with grid.
        
    resample : tuple, default=None
        Optionally resample grid, pass spacing in X and Y direction as a tuple
        e.g. resample=(spacingX, spacingY)

    Returns
    -------
    cdf_grid_z : array-like
        A numpy array of the grid defined by the supplied netCDF4 file. Can be resampled if given a specific spacing in 
        the X and Y directions. Entries are rescaled using longitudes between -180 and 180 degrees.

    cdf_lon, cdf_lat : array-like
        Numpy arrays encasing longitude and latitude variables belonging to the supplied netCDF4 file. Longitudes are 
        rescaled between -180 and 180 degrees.  
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
    """ Writes grid, latitude and longitude variables to a given netCDF4 file using specified longitudinal and latitudinal
    extents. 

    Latitude and longitude arrays correspond to the size (num of rows and columns respectively) of the given numpy grid and 
    set between specified latitudinal and longitudinal angular extents. The given grid and generated lat,lon arrays are 
    ascribed to a given netCDF4 filename and written as additional variables of the file. 

    Parameters
    ----------
    filename : str
        Path to the netCDF file

    grid : array-like
        An array with elements that define a grid. The number of rows corresponds to the number of latitudinal points, 
        while the number of columns corresponds to the number of longitudinal points.

    extent : 1D numpy array, default=[-180,180,-90,90]
        Four elements must specify the [min lon, max lon, min lat, max lat] with which to constrain the lat and lon 
        variables to write to the netCDF file. If no extents are supplied, full global extent is assumed. 

    Returns
    -------
    cdf_lon, cdf_lat : 1D numpy arrays
        Longitude and latitude variables that have been written to the supplied netCDF4 file. Lengths of these arrays 
        equal the number of cols and rows respectively of the supplied “grid”. Defined between angular extents specified 
        in “extent”.

    cdf_data : array-like
        The supplied grid is ascribed to the given netCDF4 file.
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
    interpolation methods.

    It is a child class of the scipy.interpolate module’s RegularGridInterpolator class. 

    Attributes
    ----------
    points
    values
    method
    bounds_error
    Fill_value

    Methods
    -------
    __init__(self, points, values, method="linear", bounds_error=False, fill_value=np.nan)
        Constructs all necessary attributes for the RegularGridInterpolator object.
    __call__(self, xi, method=None, return_indices=False, return_distances=False)
        Allows the RegularGridInterpolator object to be called as a method.
    """
    def __init__(self, points, values, method="linear", bounds_error=False, fill_value=np.nan):
        """Constructs all necessary attributes for the RegularGridInterpolator object.

        Parameters
        ----------
        points : tuple of 1d arrays
            Each array contains point coordinates (e.g. 2 arrays; 1 for each point’s lat, 1 for each point’s lon).
            Defines the points to sample data with. 
        values : ndarray
            Defines a grid. The number of rows corresponds to the number of latitudinal points, while the number
            of columns corresponds to the number of longitudinal points.
        method : str, default=’linear’
            The method of interpolation to perform. Supported are "linear" and "nearest". Assumes “linear” by default.
        bounds_error : bool, default=false
            Choose whether to return a ValueError and terminate the interpolation if any provided sample points are out
            of grid bounds. By default, it is set to false. In this case, all out-of-bound point values are replaced 
            with the fill_value (defined below) if supplied.
        fill_value : float, default=np.nan
            Used to replace point values that are out of grid bounds, provided that ‘bounds_error’ is false.
        """ 
        super(RegularGridInterpolator, self).__init__(points, values, method, bounds_error, fill_value)

    def __call__(self, xi, method=None, return_indices=False, return_distances=False):
        """Samples gridded data at a set of point coordinates. Uses either linear or nearest-neighbour interpolation methods.

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
    """Samples gridded data at a set of point coordinates. Uses either linear or nearest-neighbour interpolation methods.
    
    Note: if any provided sample points are out of grid bounds and a corresponding error message was suppressed (by 
    specifying bounds_error=False), all out-of-bound point values are replaced with the RegularGridInterpolator object’s
    self.fill_value attribute (if it exists). Terminates otherwise. 

    Parameters
    ----------
    lon, lat : 1d arrays
        Two arrays each specifying the longitude and latitude of sampling points for interpolation.

    grid : ndarray
        An array with elements that define a grid. The number of rows corresponds to the number of latitudinal points, while
        the number of columns corresponds to the number of longitudinal points.

    extent : 1D numpy array, default=[-180,180,-90,90]
        Four-element array to specify the [min lon, max lon, min lat, max lat] with which to constrain lat and lon sampling
        points with respect to the given grid. If no extents are supplied, full global extent is assumed. 

    return_indices : bool, default=False
        Choose whether to return indices of neighbouring sampling points. 

    return_distances : bool, default=False
        Choose whether to return normal distances between interpolated points and neighbouring sampling points.

    method : str, default=’linear’
        The method of interpolation to perform. Supported are "linear" and "nearest". Assumes “linear” by default.

    Returns
    -----
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
    interpolator = RegularGridInterpolator((np.linspace(extent[2], extent[3], grid.shape[0]),
                                            np.linspace(extent[0], extent[1], grid.shape[1])),
                                            grid, method=method)

    return interpolator(np.c_[lat, lon], return_indices=return_indices, return_distances=return_distances)




class Raster(object):
    """A class providing Scipy’s RegularGridInterpolator functionalities for interpolation. 

    Gridded data are sampled at a set of point coordinates using either linear or nearest-neighbour interpolation. 
    These grids can also be resampled using X and Y-direction spacing, and can be resized using X and Y resolutions.
    Grids can be searched for invalid, NaN-type data cells. These can be replaced with the values of their nearest
    valid neighbours. 

    Attributes
    ----------
    PlateReconstruction_object : object pointer
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
    def __init__(self, PlateReconstruction_object=None, filename=None, array=None, extent=None, resample=None):
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

        lons, lats : 1d arrays
            Either the longitude and latitude variables belonging to the netCDF4 file provided, or arrays linearly spaced 
            between the given lon & lat extents to match the dimensions of the grid array.


        The following objects + methods can be accessed in the raster object:
        _update() : method of RegularGridInterpolator
            Stored as _interpolator, this samples the “data” attribute at a set of point coordinates (generated from the 
            attributes “lats” & “lons”). Uses linear interpolation.
        """

        self.PlateReconstruction_object = PlateReconstruction_object

        # we initialise an empty points object as we do not want to build this before any resampling takes place.
        self.points = None

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
            self.data = np.array(array)
            self.extent = extent
            self.lons = np.linspace(extent[0], extent[1], self.data.shape[1])
            self.lats = np.linspace(extent[2], extent[3], self.data.shape[0])

        self._update()

        if array is not None and resample is not None:
            self.resample(*resample, override=True)


    def _update(self):
        """Stores the RegularGridInterpolator object’s method for sampling gridded data at a set of point coordinates. 

        Allows methods of the Raster object to access grid sampling functionalities. The gridded data used is the “data” 
        attribute - either read from a netCDF4 file, or supplied as an ndarray. Points to sample are either variables of the 
        netCDF4 file, or are generated from the “extent” attribute and scaled to fit the grid “data”.
        """
        # store interpolation object
        interpolator = RegularGridInterpolator((self.lats, self.lons), self.data, method='linear')
        self._interpolator = interpolator


    def interpolate(self, lons, lats, method='linear', return_indices=False, return_distances=False):
        """Samples gridded data at a set of point coordinates and interpolates points on grid. Uses either linear or 
        nearest-neighbour interpolation methods.

        Uses the grid stored in the raster object “data” attribute, and samples a series of points generated with the 
        “lons” and “lats” function parameters.
    
        Parameters
        ----------
        lons, lats : ndarray
            Longitudes and latitudes of points to sample the gridded data with. Used to generate the points ndarrays of 
            shape (..., ndim). Should have the same dimension as the grid “data” attribute.

        method : str, default=’linear’
            The method of interpolation to perform. Supported are "linear" and "Nearest". Assumes “linear” interpolation
            if None provided.  

        return_indices : bool, default=False
            Choose whether to return indices of neighbouring sampling points. 

        return_distances : bool, default=False
            Choose whether to return normal distances between interpolated points and neighbouring sampling points.

        Returns
        -------
        data_interp : tuple of ndarrays
            The first ndarray in the output tuple holds the interpolated grid data. If sample point distances and indices are
            required, these are returned as subsequent tuple elements. 

        Raises
        ------
        ValueError
            * Raised if the string method supplied is not “linear” or “nearest”.
            * Raised if the provided lat, lon arrays generate sample points that do not have the same dimensions as the 
            supplied grid. 
            * Raised if the provided lat, lon arrays generate sample points that include any point out of grid bounds. 
            Alerts user which dimension (index) the point is located. 
        """
        interp = self._interpolator
        data_interp = interp((lats,lons), method=method, return_indices=return_indices, return_distances=return_distances)
        return data_interp


    def resample(self, spacingX, spacingY, overwrite=False):
        """Resamples the grid using linear interpolation. New grid overwrites the current grid stored in the “data” attribute.
        Optional: can also resample and overwrite the arrays in the lats and lons attributes and overwrite the interpolation
        object “_update()”.

        Generates latitude and longitude arrays based on a specific spacing in X and Y directions, and the latitude and 
        longitude extents held in the “extent” raster object attribute. These lat-lon arrays are meshed into a set of sample
        points that are linearly interpolated onto the grid currently held in the “data” attribute. This final grid overwrites
        the current “data” grid. If specified by the user, the generated lat-lon arrays can also overwrite the arrays in the
        “lats” and “lons” raster object attributes.

        Parameters
        ----------
        spacingX, spacingY : ndarray
            Specify the spacing in the X and Y directions with which to resample.

        overwrite : bool, default=False
            Choose to also overwrite lons and lats currently stored in the self.lons andself.lats attributes. Doing so will 
            also overwrite the interpolation object. By default, it is false, so only the “data” grid is overwritten in that 
            case. 

        Returns
        -------
        data : meshed ndarray grid
            A new resampled and linearly-interpolated grid stored to the “data” attribute. Overwrites the current grid held
            in “data”.  
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
        """Resizes the grid with a specific resolution and samples points using linear interpolation. New grid overwrites
        the current grid stored in the “data” attribute. Optional: can also resample and overwrite the arrays in the lats
        and lons attributes and overwrite the interpolation object “_update()”.

        Generates latitude and longitude arrays based on a specific resolution in X and Y directions, and the latitude and
        longitude extents held in the “extent” raster object attribute. These lat-lon arrays are meshed into a set of sample
        points that are linearly interpolated onto the grid currently held in the “data” attribute. This final grid 
        overwrites the current “data” grid. If specified by the user, the generated lat-lon arrays can also overwrite the
        arrays in the “lats” and “lons” raster object attributes.

        Parameters
        ----------
        resX, resY : ndarray
            Specify the resolution (the larger, the finer the grid and lat-lon arrays) with which to resize.

        overwrite : bool, default=False
            Choose to also overwrite lons and lats currently stored in the self.lons andself.lats attributes. Doing so will
            also overwrite the interpolation object. By default, it is false, so only the “data” grid is overwritten in 
            that case. 

        Returns
        -------
        data : meshed ndarray grid
            A new resized and linearly-interpolated grid stored to the “data” attribute. Overwrites the current grid held
            in “data”.
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
        """Searches for invalid ‘data’ cells containing NaN-type entries (as indicated by ‘invalid’), locates the index
        of the nearest valid data cell and replaces NaNs with the value of the nearest valid data cell.

        Searches a supplied data frame “data” for invalid data cells containing NaN-type entries. If these invalidities 
        have been replaced with a “fill_value” attribute before being passed into fill_raster, there may no longer be 
        invalid cells in the data. If so, the data’s NaN entries are recovered. Where there is an invalid cell in “data”,
        the index of the closest valid data cell entry is located. This is superimposed onto “data”’s invalid cells and
        replaces NaNs with the nearest valid value.

        Parameters
        ---------
        overwrite : bool, default=False
            Choose whether to overwrite the grid currently held in the “data” raster object attribute with the new grid
            (which will have any NaNs filled).

        Returns
        --------
        data : ndarray
            An updated grid of data where each invalid cell has been replaced with the value of its nearest valid neighbour. 
        """
        data = fill_raster(self.data)
        if overwrite:
            self.data = data

        return data


    def save_to_NetCDF4(self, filename):
        """ Saves file to netCDF4 format"""
        write_netcdf_grid(str(filename), self.data, self.extent)


    def reconstruct(self, to_time, from_time=0, anchor_plate_id=0, **kwargs):

        import stripy

        lonq, latq = np.meshgrid(self.lons, self.lats)
        lonq_ = lonq.ravel()
        latq_ = latq.ravel()

        if self.points is None:
            self.points = _Points(self.PlateReconstruction_object, lonq_, latq_, from_time, anchor_plate_id)

        lons, lats = self.points.reconstruct(to_time, anchor_plate_id=anchor_plate_id, **kwargs)

        # also remove duplicate entries - # BUT this sorts the indices!!
        _, uindex = np.unique(np.c_[lons,lats], return_index=True, axis=0)
        uindex_sorted = sorted(uindex)
        ilons = lons[uindex_sorted]
        ilats = lats[uindex_sorted]
        idata = self.data.flat[uindex_sorted]


        # interpolate onto sphere
        # this is not very elegant - need to work out why stripy is struggling here.
        for i in range(10):
            try:
                mesh = stripy.sTriangulation(np.radians(ilons), np.radians(ilats), tree=True, permute=True)
            except ValueError:
                pass
            else:
                break
        zi, ierr = mesh.interpolate(np.radians(lonq_), np.radians(latq_), idata, order=1)

        # get cell spacing
        dx = np.diff(self.lons).mean()
        dy = np.diff(self.lats).mean()
        dxy = np.hypot(dx,dy)
        rxy = np.radians(dxy)

        # find angular separation / great circle distance between mesh and reconstructed points
        angles, idx = mesh.nearest_vertices(np.radians(lonq_), np.radians(latq_))

        zi[angles.ravel() > rxy] = np.nan
        zi = zi.reshape(lonq.shape)
        return zi





class TimeRaster(Raster):

    def __init__(self, PlateReconstruction_object, filename=None, array=None, extent=None, resample=None):

        super(TimeRaster, self).__init__(filename, array, extent, resample)

        self.PlateReconstruction_object = PlateReconstruction_object