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

# ----- parameters for plot
region = "d"
width = 10
projection = "N180/"
x_offset = width + 2

# plate boundary stuff
plateboundary_width = "0.5p"
age_font = "12p,Helvetica,black"
label_font = "12p,Helvetica,black"
label_offset = "j0/-0.5c"
label_position = "TC"


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
    fig.plot(
        data=gdf_subduction_left, pen=f"0.5p,{color}", fill=color, style="f0.2/0.08+l+t"
    )
    fig.plot(
        data=gdf_subduction_right,
        pen=f"0.5p,{color}",
        fill=color,
        style="f0.2/0.08+r+t",
    )


def plot_geo_data_frame(fig: pygmt.Figure, gdf: GeoDataFrame, **kwargs):
    line_width = "0.1p"
    line_color = "blue"

    if "edgecolor" in kwargs.keys():
        if isinstance(kwargs["edgecolor"], str):
            line_color = kwargs["edgecolor"]
        else:
            raise Exception(
                "The edgecolor parameter is not string. Currently, the pygmt plot engine only supports colour name."
            )

    if "linewidth" in kwargs.keys():
        line_width = f"{kwargs['linewidth']}p"

    fill = None
    if "facecolor" in kwargs.keys() and kwargs["facecolor"].lower() != "none":
        fill = f"{kwargs['facecolor']}"

    if line_color.lower() == "none":
        line_width = "0"
        line_color = fill

    if "fill" in kwargs.keys():
        fill = kwargs["fill"]

    if "pen" in kwargs.keys():
        pen = kwargs["pen"]
    else:
        pen = f"{line_width},{line_color}"

    style = None
    if "style" in kwargs.keys():
        style = kwargs["style"]

    label = None
    if "gmtlabel" in kwargs.keys():
        label = kwargs["gmtlabel"]

    fig.plot(
        data=gdf.geometry, pen=pen, fill=fill, style=style, transparency=0, label=label
    )

    """
    fig.plot(data=gdf_coastlines, fill=coastline_color, frame=["xa0", "ya0"],  transparency=0)

    fig.plot(data=gdf_topo_plates.geometry, pen='%s,%s' % (plateboundary_width, plate_colour), frame="lrtb")
    fig.plot(data=gdf_subduction_left, pen='%s,%s' % (plateboundary_width, subduction_zone_colour), fill=subduction_zone_colour, style='f0.2/0.08+l+t')
    fig.plot(data=gdf_subduction_right, pen='%s,%s' % (plateboundary_width, subduction_zone_colour), fill=subduction_zone_colour, style='f0.2/0.08+r+t')
    fig.plot(data=gdf_ridges_transforms, pen='%s,%s' % (plateboundary_width, ridge_colour))
    fig.plot(data=gplot.get_transforms(), pen='%s,%s' % (plateboundary_width, transform_color))

    fig.text(text='gplot.get_transforms(): %s Ma' % age, position=label_position, no_clip=True, font=label_font, offset=label_offset)

    fig.shift_origin(xshift=x_offset)
    fig.basemap(region=region, projection="%s%sc" % (projection, width), frame="lrtb")
    fig.plot(data=gdf_cobs, fill=COB_color, transparency=0, )
    fig.plot(data=gdf_coastlines, fill=coastline_color, frame=["xa0", "ya0"],  transparency=0)

    fig.plot(data=gdf_topo_plates.geometry, pen='%s,%s' % (plateboundary_width, plate_colour), frame="lrtb", label='other plate boundary types')
    fig.plot(data=gdf_subduction_left, pen='%s,%s' % (plateboundary_width, subduction_zone_colour), fill=subduction_zone_colour, style='f0.2/0.08+l+t', label='subduction zones')
    fig.plot(data=gdf_subduction_right, pen='%s,%s' % (plateboundary_width, subduction_zone_colour), fill=subduction_zone_colour, style='f0.2/0.08+r+t')
    fig.plot(data=gdf_ridges_transforms, pen='%s,%s' % (plateboundary_width, ridge_colour), label='ridges and transforms')

    # from gpml: transforms
    fig.plot(data=gdf_topo_transforms, pen='%s,%s' % (plateboundary_width, transform_color), label = 'transforms')
    fig.text(text='FeatureType.gpml_transform: %s Ma' % age, position=label_position, no_clip=True, font=label_font, offset=label_offset)

    fig.legend(position='jBL+o-2.7/0', box="+gwhite+p0.5p")
    """
