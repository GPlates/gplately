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

from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

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

    with shpreader.Reader(filename) as reader:
        shape_records = reader.shapeRecords()
        shapes = [i.shape for i in shape_records]
        geoms = [shape(i.__geo_interface__) for i in shapes]
    return buffer_func(geoms, buffer)
