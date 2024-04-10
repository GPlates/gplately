import logging
import warnings

import geopandas as gpd

from .._utils.feature_utils import shapelify_features
from .._utils.plot_utils import _clean_polygons, _meridian_from_ax

logger = logging.getLogger("gplately")


def get_ridges(self, central_meridian=0.0, tessellate_degrees=1):
    """Create a geopandas.GeoDataFrame object containing geometries of reconstructed ridge lines.

    Notes
    -----
    The `ridges` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
    parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
    either when `PlotTopologies` is first called...

        gplot = gplately.PlotTopologies(..., time=100,...)

    or anytime afterwards, by setting:

        time = 100 #Ma
        gplot.time = time

    ...after which this function can be re-run. Once the `ridges` are reconstructed, they are
    converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

    Returns
    -------
    gdf : instance of <geopandas.GeoDataFrame>
        A pandas.DataFrame that has a column with `ridges` geometry.
    central_meridian : float
        Central meridian around which to perform wrapping; default: 0.0.
    tessellate_degrees : float or None
        If provided, geometries will be tessellated to this resolution prior
        to wrapping.

    Raises
    ------
    ValueError
        If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
        `ridges` to the requested `time` and thus populate the GeoDataFrame.

    """
    if self._time is None:
        raise ValueError(
            "No ridges have been resolved. Set `PlotTopologies.time` to construct ridges."
        )

    if self.ridges is None:
        raise ValueError("No ridge topologies passed to PlotTopologies.")

    ridge_lines = shapelify_features(
        self.ridges,
        central_meridian=central_meridian,
        tessellate_degrees=tessellate_degrees,
    )
    gdf = gpd.GeoDataFrame({"geometry": ridge_lines}, geometry="geometry")
    return gdf


def plot_ridges(self, ax, color="black", **kwargs):
    """Plot reconstructed ridge polylines onto a standard map Projection.

    Notes
    -----
    The `ridges` for plotting are accessed from the `PlotTopologies` object's
    `ridges` attribute. These `ridges` are reconstructed to the `time`
    passed to the `PlotTopologies` object and converted into Shapely polylines.
    The reconstructed `ridges` are plotted onto the GeoAxes or GeoAxesSubplot map
    `ax` using GeoPandas. Map presentation details (e.g. `facecolor`, `edgecolor`, `alpha`…)
    are permitted as keyword arguments.

    Ridge geometries are wrapped to the dateline using
    pyGPlates' [DateLineWrapper](https://www.gplates.org/docs/pygplates/generated/pygplates.datelinewrapper)
    by splitting a polyline into multiple polylines at the dateline. This is to avoid
    horizontal lines being formed between polylines at longitudes of -180 and 180 degrees.
    Point features near the poles (-89 & 89 degree latitude) are also clipped to ensure
    compatibility with Cartopy.

    Parameters
    ----------
    ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
        A subclass of `matplotlib.axes.Axes` which represents a map Projection.
        The map should be set at a particular Cartopy projection.

    color : str, default=’black’
        The colour of the ridge lines. By default, it is set to black.

    **kwargs :
        Keyword arguments for parameters such as `alpha`, etc. for
        plotting ridge geometries.
        See `Matplotlib` keyword arguments
        [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

    Returns
    -------
    ax : instance of <geopandas.GeoDataFrame.plot>
        A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
        with ridge features plotted onto the chosen map projection.
    """
    if not self.plate_reconstruction.topology_features:
        logger.warn(
            "Plate model does not have topology features. Unable to plot_ridges."
        )
        return

    if "transform" in kwargs.keys():
        warnings.warn(
            "'transform' keyword argument is ignored by PlotTopologies",
            UserWarning,
        )
        kwargs.pop("transform")
    tessellate_degrees = kwargs.pop("tessellate_degrees", 1)
    central_meridian = kwargs.pop("central_meridian", None)
    if central_meridian is None:
        central_meridian = _meridian_from_ax(ax)

    gdf = self.get_ridges(
        central_meridian=central_meridian,
        tessellate_degrees=tessellate_degrees,
    )
    if hasattr(ax, "projection"):
        gdf = _clean_polygons(data=gdf, projection=ax.projection)
    else:
        kwargs["transform"] = self.base_projection
    return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)
