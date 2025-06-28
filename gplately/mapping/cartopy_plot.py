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
import logging
import math

import cartopy.crs as ccrs
from geopandas.geodataframe import GeoDataFrame

from ..tools import EARTH_RADIUS
from ..utils.plot_utils import plot_subduction_teeth
from .plot_engine import PlotEngine

logger = logging.getLogger("gplately")

DEFAULT_CARTOPY_PROJECTION = ccrs.PlateCarree()


class CartopyPlotEngine(PlotEngine):
    """Use Cartopy for map plotting"""

    def __init__(self):
        pass

    def plot_geo_data_frame(self, ax_or_fig, gdf: GeoDataFrame, **kwargs):
        """Use Cartopy to plot geometries in a GeoDataFrame object onto a map

        Parameters
        ----------
        ax_or_fig : cartopy.mpl.geoaxes.GeoAxes
            Cartopy GeoAxes instance
        gdf : GeoDataFrame
            GeoPandas GeoDataFrame object

        """
        if hasattr(ax_or_fig, "projection"):
            if gdf.crs is None:
                gdf.crs = ccrs.PlateCarree()
            gdf = gdf.to_crs(ax_or_fig.projection)
        else:
            kwargs["transform"] = DEFAULT_CARTOPY_PROJECTION

        return gdf.plot(ax=ax_or_fig, **kwargs)

    def plot_pygplates_features(self, ax_or_fig, features, **kwargs):
        """Not implemented yet"""
        pass

    def plot_subduction_zones(
        self,
        ax_or_fig,
        gdf_subduction_left: GeoDataFrame,
        gdf_subduction_right: GeoDataFrame,
        color="blue",
        **kwargs,
    ):
        """Use Cartopy to plot subduction zones with "teeth" onto a map

        Parameters
        ----------
        ax_or_fig : cartopy.mpl.geoaxes.GeoAxes
            Cartopy GeoAxes instance
        gdf_subduction_left : GeoDataFrame
            subduction zone with "left" polarity
        gdf_subduction_right : GeoDataFrame
            subduction zone with "right" polarity
        color : str
            The colour used to fill the "teeth".

        """
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
