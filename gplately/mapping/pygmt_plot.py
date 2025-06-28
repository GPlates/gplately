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
import pygmt
from geopandas.geodataframe import GeoDataFrame

from .plot_engine import PlotEngine

pygmt.config(
    FONT_ANNOT=8,
    FONT_LABEL=8,
    FONT=8,
    MAP_TICK_PEN="0.75p",
    MAP_FRAME_PEN="0.75p",
    MAP_TICK_LENGTH_PRIMARY="4p",
)

# NW's example is at https://gist.github.com/nickywright/f53018a8eda29223cca6f39ab2cfa25d


class PygmtPlotEngine(PlotEngine):
    """Use pygmt for map plotting"""

    def __init__(self):
        pass

    def plot_geo_data_frame(self, ax_or_fig, gdf: GeoDataFrame, **kwargs):
        """Use pygmt to plot geometries in a GeoDataFrame object onto a map

        Parameters
        ----------
        ax_or_fig : pygmt.Figure()
            pygmt Figure object
        gdf : GeoDataFrame
            GeoPandas GeoDataFrame object
        edgecolor : str
            For polygons, it is the border colour. For polylines, it is the line colour.
            Currently, only colour names are tested and officially supported, for example, "red", "blue", etc.
        facecolor : str
            The colour used to fill the polygon.
        fill : str
            GMT "fill" parameter
        pen : str
            GMT "pen" parameter
        style : str
            GMT "style" parameter
        gmtlabel : str
            GMT "label" parameter

        """
        line_color = kwargs.pop("edgecolor", "blue")
        line_width = f"{kwargs.pop('linewidth',0.1)}p"

        fill = kwargs.pop("facecolor", None)
        if fill and fill.lower() == "none":
            fill = None
        fill = kwargs.pop("fill", fill)  # the "fill" parameter override the "facecolor"

        if line_color.lower() == "none":
            # line_width = "0"
            # line_color = fill
            pen = None
        else:
            pen = kwargs.pop("pen", f"{line_width},{line_color}")
        style = kwargs.pop("style", None)
        label = kwargs.pop("gmtlabel", None)

        ax_or_fig.plot(
            data=gdf.geometry,
            pen=pen,
            fill=fill,
            style=style,
            transparency=0,
            label=label,
        )

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
        """Use pygmt to plot subduction zones with "teeth"

        Parameters
        ----------
        ax_or_fig : pygmt.Figure()
            pygmt Figure object
        gdf_subduction_left : GeoDataFrame
            subduction zone with "left" polarity
        gdf_subduction_right : GeoDataFrame
            subduction zone with "right" polarity
        color : str
            The colour used to fill the "teeth".
        gmtlabel : str
            GMT "label" parameter
        """
        label = kwargs.pop("gmtlabel", None)

        ax_or_fig.plot(
            data=gdf_subduction_left,
            pen=f"0.5p,{color}",
            fill=color,
            style="f0.2/0.08+l+t",
            label=label,
        )
        ax_or_fig.plot(
            data=gdf_subduction_right,
            pen=f"0.5p,{color}",
            fill=color,
            style="f0.2/0.08+r+t",
        )
