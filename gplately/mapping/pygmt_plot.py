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


def plot_geo_data_frame(fig: pygmt.Figure, gdf: GeoDataFrame, **kwargs):
    fig.plot(data=gdf.geometry, pen="0.5p,blue")
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
