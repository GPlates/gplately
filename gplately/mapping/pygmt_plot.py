#
#    Copyright (C) 2024 The University of Sydney, Australia
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
from geopandas.geodataframe import GeoDataFrame
import pygmt
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
    def __init__(self):
        pass

    def plot_geo_data_frame(self, ax_or_fig, gdf: GeoDataFrame, **kwargs):
        plot_geo_data_frame(ax_or_fig, gdf, **kwargs)

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
        plot_subduction_zones(
            ax_or_fig, gdf_subduction_left, gdf_subduction_right, color=color, **kwargs
        )


def get_pygmt_basemap_figure(projection="N180/10c", region="d"):
    fig = pygmt.Figure()
    fig.basemap(region=region, projection=projection, frame="lrtb")
    return fig


def plot_subduction_zones(
    fig: pygmt.Figure,
    gdf_subduction_left: GeoDataFrame,
    gdf_subduction_right: GeoDataFrame,
    color="blue",
    **kwargs,
):
    label = kwargs.pop("gmtlabel", None)

    fig.plot(
        data=gdf_subduction_left,
        pen=f"0.5p,{color}",
        fill=color,
        style="f0.2/0.08+l+t",
        label=label,
    )
    fig.plot(
        data=gdf_subduction_right,
        pen=f"0.5p,{color}",
        fill=color,
        style="f0.2/0.08+r+t",
    )


def plot_geo_data_frame(fig: pygmt.Figure, gdf: GeoDataFrame, **kwargs):

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

    fig.plot(
        data=gdf.geometry, pen=pen, fill=fill, style=style, transparency=0, label=label
    )
