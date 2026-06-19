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
Tools to read geometry data from input files and output them as `Shapely`
geometries. These geometries can be plotted directly with GPlately's
`PlotTopologies` object.

By default, input files are read with `GeoPandas` and output as a
`geopandas.GeoSeries` object that contains `Shapely` geometries.
If `GeoPandas` is not found on the system, input files are read with
`Shapely` instead and are still returned as `Shapely` geometries.

"""

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false


from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry
import xarray as xr

gpd = None
shpreader = None
try:
    import geopandas as gpd

    USE_GEOPANDAS = True
except ImportError:
    import shapefile as shpreader

    USE_GEOPANDAS = False

__all__ = [
    "get_geometries",
    "get_valid_geometries",
]


def get_geometries(filename, buffer=None):
    """Read a file and return feature geometries.

    If `geopandas` is available, it will be used to read the file,
    returning a `geopandas.GeoSeries`. If `geopandas` is not found,
    only shapefiles can be read, and a list of `shapely` geometries
    will be returned instead of a `geopandas.GeoSeries`.

    Parameters
    ----------
    filename : str
        Path to the file to be read.

    Returns
    -------
    geometries : list or geopandas.GeoSeries
        `shapely` geometries that define the feature geometry held in the
        shapefile.
    """
    if USE_GEOPANDAS:
        return _get_geometries_geopandas(filename, buffer=buffer)
    return _get_geometries_cartopy(filename, buffer=buffer)


def get_valid_geometries(filename):
    """Read a file and return valid feature geometries.

    If `geopandas` is available, it will be used to read the file,
    returning a `geopandas.GeoSeries`. If `geopandas` is not found,
    only shapefiles can be read, and a list of `shapely` geometries
    will be returned instead of a `geopandas.GeoSeries`.

    Parameters
    ----------
    filename : str
        Path to the file to be read.

    Returns
    -------
    geometries : list or geopandas.GeoSeries
        Valid `shapely` geometries that define the feature geometry held in the
        shapefile.
    """
    return get_geometries(filename, buffer=0.0)


def _get_geometries_geopandas(filename, buffer=None):
    def buffer_func(geoms, buffer=None):
        if buffer is not None:
            geoms = geoms.buffer(buffer)
        return geoms

    if gpd is None:
        raise ImportError(
            "Geopandas is not available. Please install geopandas to read files with geopandas."
        )

    if isinstance(filename, gpd.GeoDataFrame):
        return buffer_func(filename.geometry, buffer)
    if isinstance(filename, gpd.GeoSeries):
        return buffer_func(filename, buffer)
    if isinstance(filename, BaseGeometry):
        return buffer_func(gpd.GeoSeries([filename]), buffer)
    try:
        for i in filename:
            if isinstance(i, BaseGeometry):
                # Iterable of geometries
                return buffer_func(gpd.GeoSeries(filename), buffer)
            break
    except TypeError:
        # Not an iterable
        # Since strings are iterable, anything that's not an iterable
        # will probably fail at the next step anyway
        pass
    # If it gets to this line, `filename` should actually be a filename
    gdf = gpd.read_file(filename)
    return buffer_func(gdf.geometry, buffer)


def _get_geometries_cartopy(filename, buffer=None):
    def buffer_func(geoms, buffer=None):
        if buffer is None:
            return list(geoms)
        out = []
        for i in geoms:
            out.append(i.buffer(buffer))
        return out

    if isinstance(filename, BaseGeometry):
        return buffer_func([filename], buffer)
    try:
        for i in filename:
            if isinstance(i, BaseGeometry):
                return buffer_func(filename, buffer)
            break
    except TypeError:
        pass

    if shpreader is None:
        raise ImportError(
            "Shapefile reader is not available. Please install pyshp to read shapefiles without geopandas."
        )

    with shpreader.Reader(filename) as reader:
        shape_records = reader.shapeRecords()
        shapes = [i.shape for i in shape_records]
        geoms = [shape(i.__geo_interface__) for i in shapes]
    return buffer_func(geoms, buffer)


def load_data_array_from_netcdf(filename, var_name=None):
    """Load a data array from a netCDF file.

    Parameters
    ----------
    filename : str
        Path to the netCDF file to be read.
    var_name : str, optional
        The variable name of the raster data in the netCDF file. If not provided, the first variable in the netCDF file will be used.

    Returns
    -------
    data_array : xarray.DataArray
        The data array loaded from the netCDF file.
    """

    try:
        raster_xr = xr.open_dataarray(filename)
    except ValueError:
        # Fallback for NetCDF files with multiple variables.
        dataset = xr.open_dataset(filename, decode_times=False)
        if var_name and var_name in dataset.data_vars:
            raster_xr = dataset[var_name]
        else:
            first_var = next(iter(dataset.data_vars))
            raster_xr = dataset[first_var]
    return raster_xr


def to_geographic_data_array(data_array):
    """Convert a DataArray to a geographic grid format that PyGMT can understand.

    This function ensures that the DataArray has the correct coordinate names and attributes for PyGMT
    to treat it as a geographic grid. Specifically, it sets the 'gmt.gtype' attribute to 1 (indicating a geographic grid)
    and ensures that the coordinate names are 'lat' and 'lon'.

    Parameters
    ----------
    data_array : xarray.DataArray
        The input DataArray with lat/lon coordinates.

    Returns
    -------
    xarray.DataArray
        A DataArray formatted for use with PyGMT, with appropriate attributes set.
    """
    if data_array.ndim != 2:
        raise ValueError("Input data array must be 2-dimensional.")
    ret_da = data_array.copy()

    if (
        "x" in ret_da.dims
        and "y" in ret_da.dims
        or "x" in ret_da.coords
        and "y" in ret_da.coords
    ):
        ret_da = ret_da.rename({"x": "lon", "y": "lat"})
    elif (
        "longitude" in ret_da.dims
        and "latitude" in ret_da.dims
        or "longitude" in ret_da.coords
        and "latitude" in ret_da.coords
    ):
        ret_da = ret_da.rename({"longitude": "lon", "latitude": "lat"})
    else:
        pass  # assume the dims are already 'lon' and 'lat'

    if "lon" not in ret_da.coords or "lat" not in ret_da.coords:
        raise ValueError("Input data array must contain 'lon' and 'lat' coordinates.")

    # Normalize longitudes to [-180, 180) and sort so GMT receives a valid region.
    wrapped_lon = ((ret_da.coords["lon"] + 180.0) % 360.0) - 180.0
    if not wrapped_lon.equals(ret_da.coords["lon"]):
        ret_da = ret_da.assign_coords(lon=wrapped_lon).sortby("lon")

    ret_da.coords["lon"].attrs.update(
        {
            "standard_name": "longitude",
            "long_name": "longitude",
            "units": "degrees_east",
        }
    )

    ret_da.coords["lat"].attrs.update(
        {
            "standard_name": "latitude",
            "long_name": "latitude",
            "units": "degrees_north",
        }
    )

    ret_da.gmt.gtype = 1  # 1 = geographic
    # ret_da.gmt.registration = 0  # 0 = gridline node
    return ret_da
