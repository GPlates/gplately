from shapely.geometry.base import BaseGeometry
try:
    import geopandas as gpd

    GEOPANDAS_AVAILABLE = True
except ImportError:
    import cartopy.io.shapereader as shpreader

    GEOPANDAS_AVAILABLE = False

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
        shapefile. Can be plotted directly using
        `gplately.plot.add_geometries`.
    """
    if GEOPANDAS_AVAILABLE:
        return _get_geometries_geopandas(filename, buffer=buffer)
    else:
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
        shapefile. Can be plotted directly using
        `gplately.plot.add_geometries`.
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

    reader = shpreader.Reader(filename)
    return buffer_func(reader.geometries(), buffer)
