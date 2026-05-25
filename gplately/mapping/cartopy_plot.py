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
import math

import pygplates

import cartopy.crs as ccrs
from geopandas.geodataframe import GeoDataFrame

from ..grids import Raster

from ..tools import EARTH_RADIUS
from ..utils.plot_utils import plot_subduction_teeth
from .plot_engine import PlotEngine

logger = logging.getLogger("gplately")

DEFAULT_CARTOPY_PROJECTION = ccrs.PlateCarree()


class CartopyPlotEngine(PlotEngine):
    """Use Cartopy for map plotting."""

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
        """Use Cartopy to plot one or more pygplates features onto a map.

        Point-like geometries are rendered with ``scatter`` so marker styling
        behaves as expected, while line/polygon geometries are rendered with
        ``add_geometries``.

        Parameters
        ----------
        ax_or_fig : cartopy.mpl.geoaxes.GeoAxes
            Cartopy GeoAxes instance
        features : pygplates.Feature or list of pygplates.Feature
            One or more pygplates features to plot
        edgecolor : str
            For polygons, it is the border colour. For polylines, it is the line colour.
        facecolor : str
            The colour used to fill the polygon.
        crs : cartopy.crs.Projection
            The coordinate reference system of the input geometries. Default is PlateCarree (lon/lat).
        **kwargs :
            Keyword arguments for plotting the features.
            see Matplotlib's ``plot()`` and ``scatter()`` keyword arguments for line/polygon and point geometries, respectively.

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

        edgecolor = kwargs.pop("edgecolor", "blue")
        facecolor = kwargs.pop("facecolor", "none")
        crs = kwargs.pop("crs", ccrs.PlateCarree())

        for feature in features:
            geometries = feature.get_all_geometries()  # type: ignore
            if not geometries:
                continue

            for geometry in geometries:
                shapely_geometry = pygplates_to_shapely(geometry)  # type: ignore
                if shapely_geometry is None:
                    continue

                geom_type = getattr(shapely_geometry, "geom_type", "")
                if geom_type in ("Point", "MultiPoint"):
                    if geom_type == "Point":
                        coords = getattr(shapely_geometry, "coords", None)
                        if not coords:
                            continue
                        xs = [coords[0][0]]
                        ys = [coords[0][1]]
                    else:
                        geoms = getattr(shapely_geometry, "geoms", ())
                        xs = [
                            getattr(pt, "x", None)
                            for pt in geoms
                            if getattr(pt, "x", None) is not None
                        ]
                        ys = [
                            getattr(pt, "y", None)
                            for pt in geoms
                            if getattr(pt, "y", None) is not None
                        ]
                        if not xs or not ys:
                            continue

                    point_kwargs = kwargs.copy()
                    point_kwargs.pop("transform", None)
                    point_kwargs.pop("color", None)
                    point_kwargs.pop("c", None)

                    # Use small, filled dots for points unless explicitly overridden.
                    point_size = point_kwargs.pop("s", 12)
                    point_color = point_kwargs.pop("facecolors", None)
                    if point_color in (None, "none"):
                        point_color = edgecolor
                    point_edgecolors = point_kwargs.pop("edgecolors", "none")

                    ax_or_fig.scatter(  # type: ignore
                        xs,
                        ys,
                        transform=crs,
                        s=point_size,
                        edgecolors=point_edgecolors,
                        facecolors=point_color,
                        **point_kwargs,
                    )
                    continue

                ax_or_fig.add_geometries(  # type: ignore
                    [shapely_geometry],
                    crs=crs,
                    edgecolor=edgecolor,
                    facecolor=facecolor,
                    **kwargs,
                )
        return ax_or_fig

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

    def plot_grid(
        self, ax_or_fig, grid, projection=None, extent=(-180, 180, -90, 90), **kwargs
    ):
        """Plot a grid onto a map using Cartopy

        Parameters
        ----------
        ax_or_fig : cartopy.mpl.geoaxes.GeoAxes
            Cartopy GeoAxes instance
        grid : 2D array-like
            The grid data to be plotted
        projection : cartopy.crs.Projection
            The projection to use for the grid
        extent : tuple
            The extent of the grid in the form (min_lon, max_lon, min_lat, max_lat)
        **kwargs :
            Keyword arguments for plotting the grid. See Matplotlib's ``imshow()`` keyword arguments
            `here <https://matplotlib.org/3.5.1/api/_as_gen/matplotlib.axes.Axes.imshow.html>`__.

        """
        # Override matplotlib default origin ('upper')
        origin = kwargs.pop("origin", "lower")

        if isinstance(grid, Raster):
            # extract extent and origin
            extent = grid.extent
            origin = grid.origin
            data = grid.data
        else:
            data = grid

        return ax_or_fig.imshow(
            data,
            extent=extent,
            transform=projection,
            origin=origin,
            **kwargs,
        )


def _plot_feature_collection(
    feature_collection: pygplates.FeatureCollection,  # type: ignore
    title: str = "Untitled Feature Collection",
    ax=None,
    figsize=(8, 4),
    projection=ccrs.Robinson(central_longitude=180),
):
    """Helper function to plot a pygplates FeatureCollection using Cartopy.
    Not part of the public API.
    Mostly this function is for testing and debugging purposes,
    and to provide a simple example of how to plot pygplates features with Cartopy.
    """
    import cartopy.crs as ccrs  # type: ignore
    import matplotlib.pyplot as plt  # type: ignore
    from gplately.geometry import pygplates_to_shapely
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import matplotlib.ticker as mticker

    if ax is None:
        fig = plt.figure(figsize=figsize, dpi=72)
        ax = fig.add_subplot(111, projection=projection)

    ax.set_global()  # type: ignore
    # Add gridlines and lat/lon labels
    gl = ax.gridlines(  # type: ignore
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.8,
        color="gray",
        alpha=0.6,
        linestyle="--",
    )

    # Hide labels on top/right if you want cleaner maps
    # Newer Cartopy
    if hasattr(gl, "top_labels"):
        gl.top_labels = False
    if hasattr(gl, "right_labels"):
        gl.right_labels = False

    # Older Cartopy
    if hasattr(gl, "xlabels_top"):
        gl.xlabels_top = False
    if hasattr(gl, "ylabels_right"):
        gl.ylabels_right = False

    # Control tick locations
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 60))
    gl.ylocator = mticker.FixedLocator(range(-90, 91, 30))

    # Nice lon/lat formatting
    gl.xformatter = LongitudeFormatter(number_format=".0f", degree_symbol="°")
    gl.yformatter = LatitudeFormatter(number_format=".0f", degree_symbol="°")

    # Label style
    gl.xlabel_style = {"size": 10}
    gl.ylabel_style = {"size": 10}

    for feature in feature_collection:
        valid_time = feature.get_valid_time(None)  # type: ignore
        if valid_time is not None:
            if valid_time[1] not in [pygplates.GeoTimeInstant.create_distant_future(), 0]:  # type: ignore
                continue  # skip features that are not valid at 0 Ma
        geometries = feature.get_geometries()  # type: ignore
        if geometries:
            for geometry in geometries:
                shapely_geometry = pygplates_to_shapely(geometry)  # type: ignore
                if shapely_geometry is not None:
                    ax.add_geometries(  # type: ignore
                        [shapely_geometry],
                        crs=ccrs.PlateCarree(),
                        edgecolor="blue",
                        facecolor="none",
                    )
    ax.set_title(title)
    # plt.show()
