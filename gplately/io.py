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
        Shapefile filename.

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
    """Read a shapefile and return valid feature geometries.

    If `geopandas` is available, it will be used to read the file,
    returning a `geopandas.GeoSeries`. If `geopandas` is not found,
    only shapefiles can be read, and a list of `shapely` geometries
    will be returned instead of a `geopandas.GeoSeries`.

    Parameters
    ----------
    filename : str
        Shapefile filename.

    Returns
    -------
    geometries : list or geopandas.GeoSeries
        Valid `shapely` geometries that define the feature geometry held in the
        shapefile. Can be plotted directly using
        `gplately.plot.add_geometries`.
    """
    return get_geometries(filename, buffer=0.0)


def _get_geometries_geopandas(filename, buffer=None):
    gdf = gpd.read_file(filename)
    geoms = gdf.geometry
    if buffer is not None:
        geoms = geoms.buffer(buffer)
    return geoms


def _get_geometries_cartopy(filename, buffer=None):
    reader = shpreader.Reader(filename)
    shp_geom = reader.geometries()
    geometries = []
    for record in shp_geom:
        if buffer is not None:
            geometries.append(record.buffer(0.0))
        else:
            geometries.append(record)
    return geometries
