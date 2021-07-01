import numpy as np
from scipy.interpolate import RegularGridInterpolator as _RGI
from scipy.ndimage import distance_transform_edt

def fill_raster(data,invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid')
    by the value of the nearest valid data cell
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
    """
    Read in a netCDF file and re-align from -180 to 180 degrees
    
    Parameters
    ----------
    filename : str
        path to netCDF file
    return_grids : bool
        optionally return lon, lat arrays associated with grid
    resample : tuple
        optionally resample grid, pass spacing in X and Y direction as a tuple
        e.g. resample=(spacingX, spacingY)
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
        
        if cdf_lon_mask.any():
            cdf_grid_z = np.hstack([cdf_grid[:,cdf_lon_mask], cdf_grid[:,~cdf_lon_mask]])
            cdf_lon = np.hstack([cdf_lon[cdf_lon_mask], cdf_lon[~cdf_lon_mask]])
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

    def __init__(self, points, values, method="linear", bounds_error=False, fill_value=np.nan):

        super(RegularGridInterpolator, self).__init__(points, values, method, bounds_error, fill_value)

    def __call__(self, xi, method=None, return_indices=False, return_distances=False):
        """
        Interpolation at coordinates
        Parameters
        ----------
        xi : ndarray of shape (..., ndim)
            The coordinates to sample the gridded data at
        method : str
            The method of interpolation to perform. Supported are "linear" and
            "nearest".
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
    interpolator = RegularGridInterpolator((np.linspace(extent[2], extent[3], grid.shape[0]),
                                            np.linspace(extent[0], extent[1], grid.shape[1])),
                                            grid, method=method)

    return interpolator(np.c_[lat, lon], return_indices=return_indices, return_distances=return_distances)




class Raster(object):

    def __init__(self, filename=None, array=None, extent=None, resample=None):

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


    def _update(self):
        # store interpolation object
        interpolator = RegularGridInterpolator((self.lats, self.lons), self.data, method='linear')
        self._interpolator = interpolator


    def interpolate(self, lons, lats, method='linear', return_indices=False, return_distances=False):
        interp = self._interpolator
        data_interp = interp((lats,lons), method=method, return_indices=return_indices, return_distances=return_distances)
        return data_interp


    def resample(self, spacingX, spacingY, overwrite=False):

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

        data = fill_raster(self.data)
        if overwrite:
            self.data = data

        return data


    def save_to_NetCDF4(filename):
        pass


class TimeRaster(Raster):

    def __init__(self, PlateReconstruction_object, filename=None, array=None, extent=None, resample=None):

        super(TimeRaster, self).__init__(filename, array, extent, resample)

        self.PlateReconstruction_object = PlateReconstruction_object