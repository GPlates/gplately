import logging, math
from geopandas.geodataframe import GeoDataFrame
import cartopy.crs as ccrs
from .plot_engine import PlotEngine
from ..utils.plot_utils import _clean_polygons, plot_subduction_teeth
from ..tools import EARTH_RADIUS

logger = logging.getLogger("gplately")

DEFAULT_CARTOPY_PROJECTION = ccrs.PlateCarree()


class CartopyPlotEngine(PlotEngine):
    def __init__(self):
        pass

    def plot_geo_data_frame(self, ax_or_fig, gdf: GeoDataFrame, **kwargs):
        if hasattr(ax_or_fig, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax_or_fig.projection)
        else:
            kwargs["transform"] = DEFAULT_CARTOPY_PROJECTION

        return gdf.plot(ax=ax_or_fig, **kwargs)

    def plot_pygplates_features(self, ax_or_fig, features, **kwargs):
        pass

    def plot_subduction_zones(
        self,
        ax_or_fig,
        gdf_subduction_left: GeoDataFrame,
        gdf_subduction_right: GeoDataFrame,
        color="blue",
        **kwargs,
    ):
        if "transform" in kwargs.keys():
            logger.warning(
                "'transform' keyword argument is ignored by CartopyPlotEngine."
            )
            kwargs.pop("transform")

        spacing = kwargs.pop("spacing")
        size = kwargs.pop("size")
        aspect = kwargs.pop("aspect")

        try:
            projection = ax_or_fig.projection
        except AttributeError:
            logger.warning(
                "The ax.projection does not exist. You must set projection to plot Cartopy maps, such as ax = plt.subplot(211, projection=cartopy.crs.PlateCarree())"
            )
            projection = None

        if isinstance(projection, ccrs.PlateCarree):
            spacing = math.degrees(spacing)
        else:
            spacing = spacing * EARTH_RADIUS * 1e3

        if aspect is None:
            aspect = 2.0 / 3.0
        if size is None:
            size = spacing * 0.5

        height = size * aspect

        plot_subduction_teeth(
            gdf_subduction_left,
            size,
            "l",
            height,
            spacing,
            projection=projection,
            ax=ax_or_fig,
            color=color,
            **kwargs,
        )
        plot_subduction_teeth(
            gdf_subduction_right,
            size,
            "r",
            height,
            spacing,
            projection=projection,
            ax=ax_or_fig,
            color=color,
            **kwargs,
        )
