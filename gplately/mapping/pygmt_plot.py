#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
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
from pathlib import Path

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false
import pygplates

logger = logging.getLogger("gplately")
try:
    import pygmt

    pygmt.config(
        FONT_ANNOT=8,
        FONT_LABEL=8,
        FONT=8,
        MAP_TICK_PEN="0.75p",
        MAP_FRAME_PEN="0.75p",
        MAP_TICK_LENGTH_PRIMARY="4p",
    )
except:
    logger.error("Failed to import PyGMT. PyGMT requires Python>=3.11.")
from geopandas.geodataframe import GeoDataFrame

from .plot_engine import PlotEngine

# NW's example is at https://gist.github.com/nickywright/f53018a8eda29223cca6f39ab2cfa25d


class PygmtPlotEngine(PlotEngine):
    """Use PyGMT for map plotting."""

    def __init__(self):
        pass

    def plot_geo_data_frame(self, ax_or_fig, gdf: GeoDataFrame, **kwargs):
        """Use PyGMT to plot geometries in a GeoDataFrame object onto a map.

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
            GMT "label" parameter for lines and polygons
        pointlabel : str
            Optional label applied only to point geometries.

        """
        line_color = kwargs.pop("edgecolor", "blue")
        line_width = f"{kwargs.pop('linewidth',0.5)}p"

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
        kwargs.pop("label", None)
        label = kwargs.pop("gmtlabel", None)
        point_label = kwargs.pop("pointlabel", None)

        point_geometries = []
        line_geometries = []
        polygon_geometries = []
        for geometry in gdf.geometry:
            geom_type = getattr(geometry, "geom_type", "")
            if geom_type in ("Point", "MultiPoint"):
                point_geometries.append(geometry)
            elif geom_type in ("LineString", "MultiLineString"):
                line_geometries.append(geometry)
            else:
                polygon_geometries.append(geometry)

        if point_geometries:
            point_xs = []
            point_ys = []
            for geometry in point_geometries:
                geom_type = getattr(geometry, "geom_type", "")
                if geom_type == "Point":
                    point_xs.append(geometry.x)
                    point_ys.append(geometry.y)
                elif geom_type == "MultiPoint":
                    for point in getattr(geometry, "geoms", ()):
                        x = getattr(point, "x", None)
                        y = getattr(point, "y", None)
                        if x is None or y is None:
                            continue
                        point_xs.append(x)
                        point_ys.append(y)

            if point_xs and point_ys:
                point_style = kwargs.pop("pointstyle", "c0.08c")
                point_fill = kwargs.pop("pointfill", None)
                if point_fill in (None, "none"):
                    point_fill = line_color
                point_pen = kwargs.pop("pointpen", None)

                ax_or_fig.plot(
                    x=point_xs,
                    y=point_ys,
                    pen=point_pen,
                    fill=point_fill,
                    style=point_style,
                    transparency=0,
                    label=point_label,
                    **kwargs,
                )

        if line_geometries:
            line_gdf = GeoDataFrame(
                {"geometry": line_geometries}, geometry="geometry", crs=gdf.crs
            )
            ax_or_fig.plot(
                data=line_gdf.geometry,
                pen=pen,
                transparency=0,
                label=label,
                **kwargs,
            )

        if not polygon_geometries:
            return ax_or_fig

        gdf = GeoDataFrame(
            {"geometry": polygon_geometries}, geometry="geometry", crs=gdf.crs
        )

        ax_or_fig.plot(
            data=gdf.geometry,
            pen=pen,
            fill=fill,
            style=style,
            transparency=0,
            label=label,
        )

    def plot_pygplates_features(self, ax_or_fig, features, **kwargs):
        """Use PyGMT to plot one or more pygplates features onto a map.

        Parameters
        ----------
        ax_or_fig : pygmt.Figure()
            pygmt Figure object
        features : pygplates.Feature or list of pygplates.Feature
            One or more pygplates.Feature objects to be plotted.
        edgecolor : str
            For polygons, it is the border colour. For polylines, it is the line colour.
        facecolor : str
            The colour used to fill the polygon.
        fill : str
            GMT "fill" parameter
        pen : str
            GMT "pen" parameter

        Warnings
        --------
        This method will not check features' valid time. It just simply plots all the geometries in the features.
        You need to filter features by valid time yourself before passing them to this method if you want to plot features at a specific time.

        .. seealso::

            Use the class :class:`ValidTimeFilter` for filtering features by valid time.
        """
        from gplately.geometry import pygplates_to_shapely

        if isinstance(features, pygplates.Feature):
            features = [features]

        shapely_geometries = []
        for feature in features:
            geometries = feature.get_all_geometries()  # type: ignore
            if not geometries:
                continue

            for geometry in geometries:
                shapely_geometry = pygplates_to_shapely(geometry)  # type: ignore
                if shapely_geometry is not None:
                    shapely_geometries.append(shapely_geometry)

        gdf = GeoDataFrame(
            {"geometry": shapely_geometries}, geometry="geometry", crs="EPSG:4326"
        )
        if len(gdf) == 0:
            return ax_or_fig

        return self.plot_geo_data_frame(ax_or_fig, gdf, **kwargs)

    def plot_subduction_zones(
        self,
        ax_or_fig,
        gdf_subduction_left: GeoDataFrame,
        gdf_subduction_right: GeoDataFrame,
        color="blue",
        **kwargs,
    ):
        """Use PyGMT to plot subduction zones with "teeth"

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

    def plot_grid(
        self,
        ax_or_fig,
        grid,
        projection=None,
        extent=(-180, 180, -90, 90),
        cmap="gmt/geo",
        nan_transparent=False,
        shading=None,
        **kwargs,
    ):
        """Use PyGMT to plot a grid onto a map.

        Parameters
        ----------
        ax_or_fig : pygmt.Figure()
            A PyGMT Figure object.
        grid : Raster
            A gplately Raster object or 2D array-like grid data.
        projection : str
            Not used currently.
        extent : str or tuple
            (xmin, xmax, ymin, ymax). See details at
            https://www.pygmt.org/dev/tutorials/basics/regions.html
        cmap : str
            A built-in GMT colormaps name or a CPT file path.
        nan_transparent : bool
            If True, NaN values in the grid will be plotted as transparent.
        shading : bool, str, or grid-like, optional
            Apply illumination/hillshading to the grid image. Accepted values are:

            - ``True``: use default shading parameters (equivalent to GMT ``-I+d``).
            - A string such as ``"+a315+ne0.6"`` to pass directly as the GMT ``-I``
              option (azimuth and intensity specification).
            - An ``xarray.DataArray`` or file path pointing to an illumination grid
              computed externally (e.g. via :func:`pygmt.grdgradient`).
            - ``None`` (default): no shading is applied.
        **kwargs :
            Additional keyword arguments.
        """
        from ..grids import Raster
        import xarray as xr  # pyright: ignore[reportMissingImports]

        # we need to convert the grid data to xarray.DataArray for pygmt.grdimage().
        if isinstance(grid, Raster):
            data = xr.DataArray(
                data=grid.data,
                dims=["lat", "lon"],
                coords=dict(
                    lon=(["lon"], grid.lons),
                    lat=(["lat"], grid.lats),
                ),
            )
        else:
            data = xr.DataArray(grid)

        # check exisence if cmap is a CPT file
        if cmap.endswith(".cpt"):
            if not Path(cmap).exists():
                raise FileNotFoundError(f"The CPT file '{cmap}' does not exist.")

        grdimage_kwargs = dict(
            grid=data,
            cmap=cmap,
            region=extent,
            nan_transparent=nan_transparent,
        )
        if shading is not None:
            grdimage_kwargs["shading"] = shading

        ax_or_fig.grdimage(**grdimage_kwargs)
