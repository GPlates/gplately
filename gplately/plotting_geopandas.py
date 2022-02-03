
import numpy as np
import matplotlib.pyplot as plt
from packaging.version import Version

def _flatten_multi_geoms(geoms, prefix="Multi"):
    """
    Returns Series like geoms and index, except that any Multi geometries
    are split into their components and indices are repeated for all component
    in the same Multi geometry.  Maintains 1:1 matching of geometry to value.
    Prefix specifies type of geometry to be flatten. 'Multi' for MultiPoint and similar,
    "Geom" for GeometryCollection.
    Returns
    -------
    components : list of geometry
    component_index : index array
        indices are repeated for all components in the same Multi geometry
    """
    components, component_index, component_type = [], [], []
    
    for ix, geom in enumerate(geoms):
        geom_type = str(geom.type)
        if geom_type.startswith(prefix) and not geom.is_empty:
            for poly in geom:
                components.append(poly)
                component_index.append(ix)
                component_type.append(str(poly.type))
        else:
            components.append(geom)
            component_index.append(ix)
            component_type.append(geom_type)

    return components, np.array(component_index), np.array(component_type)

def _expand_kwargs(kwargs, multiindex):
    """
    Most arguments to the plot functions must be a (single) value, or a sequence
    of values. This function checks each key-value pair in 'kwargs' and expands
    it (in place) to the correct length/formats with help of 'multiindex', unless
    the value appears to already be a valid (single) value for the key.
    """
    import matplotlib
    from matplotlib.colors import is_color_like
    from typing import Iterable

    mpl = Version(matplotlib.__version__)
    if mpl >= Version("3.4") or mpl > Version("3.3.2"):
        # alpha is supported as array argument with matplotlib 3.4+
        scalar_kwargs = ["marker", "path_effects"]
    else:
        scalar_kwargs = ["marker", "alpha", "path_effects"]

    for att, value in kwargs.items():
        if "color" in att:  # color(s), edgecolor(s), facecolor(s)
            if is_color_like(value):
                continue
        elif "linestyle" in att:  # linestyle(s)
            # A single linestyle can be 2-tuple of a number and an iterable.
            if (
                isinstance(value, tuple)
                and len(value) == 2
                and isinstance(value[1], Iterable)
            ):
                continue
        elif att in scalar_kwargs:
            # For these attributes, only a single value is allowed, so never expand.
            continue

        if np.size(value) > 1:
            kwargs[att] = np.take(value, multiindex, axis=0)


def _PolygonPatch(polygon, **kwargs):
    """Constructs a matplotlib patch from a Polygon geometry
    The `kwargs` are those supported by the matplotlib.patches.PathPatch class
    constructor. Returns an instance of matplotlib.patches.PathPatch.
    Example (using Shapely Point and a matplotlib axes)::
        b = shapely.geometry.Point(0, 0).buffer(1.0)
        patch = _PolygonPatch(b, fc='blue', ec='blue', alpha=0.5)
        ax.add_patch(patch)
    GeoPandas originally relied on the descartes package by Sean Gillies
    (BSD license, https://pypi.org/project/descartes) for PolygonPatch, but
    this dependency was removed in favor of the below matplotlib code.
    """
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path

    path = Path.make_compound_path(
        Path(np.asarray(polygon.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in polygon.interiors],
    )
    return PathPatch(path, **kwargs)


def _plot_polygon_collection(
    ax, geoms, values=None, color=None, cmap=None, vmin=None, vmax=None, **kwargs
):
    """
    Plots a collection of Polygon and MultiPolygon geometries to `ax`
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        where shapes will be plotted
    geoms : a sequence of `N` Polygons and/or MultiPolygons (can be mixed)
    values : a sequence of `N` values, optional
        Values will be mapped to colors using vmin/vmax/cmap. They should
        have 1:1 correspondence with the geometries (not their components).
        Otherwise follows `color` / `facecolor` kwargs.
    edgecolor : single color or sequence of `N` colors
        Color for the edge of the polygons
    facecolor : single color or sequence of `N` colors
        Color to fill the polygons. Cannot be used together with `values`.
    color : single color or sequence of `N` colors
        Sets both `edgecolor` and `facecolor`
    **kwargs
        Additional keyword arguments passed to the collection
    Returns
    -------
    collection : matplotlib.collections.Collection that was plotted
    """
    from matplotlib.collections import PatchCollection

    geoms, multiindex, multitype = _flatten_multi_geoms(geoms)
    if values is not None:
        values = np.take(values, multiindex, axis=0)

    # PatchCollection does not accept some kwargs.
    kwargs = {
        att: value
        for att, value in kwargs.items()
        if att not in ["markersize", "marker"]
    }

    # Add to kwargs for easier checking below.
    if color is not None:
        kwargs["color"] = color

    _expand_kwargs(kwargs, multiindex)

    collection = PatchCollection(
        [_PolygonPatch(poly) for poly in geoms if not poly.is_empty], **kwargs
    )

    if values is not None:
        collection.set_array(np.asarray(values))
        collection.set_cmap(cmap)
        if "norm" not in kwargs:
            collection.set_clim(vmin, vmax)

    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection


def _plot_linestring_collection(
    ax, geoms, values=None, color=None, cmap=None, vmin=None, vmax=None, **kwargs
):
    """
    Plots a collection of LineString and MultiLineString geometries to `ax`
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        where shapes will be plotted
    geoms : a sequence of `N` LineStrings and/or MultiLineStrings (can be
            mixed)
    values : a sequence of `N` values, optional
        Values will be mapped to colors using vmin/vmax/cmap. They should
        have 1:1 correspondence with the geometries (not their components).
    color : single color or sequence of `N` colors
        Cannot be used together with `values`.
    Returns
    -------
    collection : matplotlib.collections.Collection that was plotted
    """
    from matplotlib.collections import LineCollection

    geoms, multiindex, multitype = _flatten_multi_geoms(geoms)
    if values is not None:
        values = np.take(values, multiindex, axis=0)

    # LineCollection does not accept some kwargs.
    kwargs = {
        att: value
        for att, value in kwargs.items()
        if att not in ["markersize", "marker"]
    }

    # Add to kwargs for easier checking below.
    if color is not None:
        kwargs["color"] = color

    _expand_kwargs(kwargs, multiindex)

    segments = [np.array(linestring.coords)[:, :2] for linestring in geoms]
    collection = LineCollection(segments, **kwargs)

    if values is not None:
        collection.set_array(np.asarray(values))
        collection.set_cmap(cmap)
        if "norm" not in kwargs:
            collection.set_clim(vmin, vmax)

    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection


def _plot_point_collection(
    ax,
    geoms,
    values=None,
    color=None,
    cmap=None,
    vmin=None,
    vmax=None,
    marker="o",
    markersize=None,
    **kwargs,
):
    """
    Plots a collection of Point and MultiPoint geometries to `ax`
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        where shapes will be plotted
    geoms : sequence of `N` Points or MultiPoints
    values : a sequence of `N` values, optional
        Values mapped to colors using vmin, vmax, and cmap.
        Cannot be specified together with `color`.
    markersize : scalar or array-like, optional
        Size of the markers. Note that under the hood ``scatter`` is
        used, so the specified value will be proportional to the
        area of the marker (size in points^2).
    Returns
    -------
    collection : matplotlib.collections.Collection that was plotted
    """
    if values is not None and color is not None:
        raise ValueError("Can only specify one of 'values' and 'color' kwargs")

    geoms, multiindex, multitype = _flatten_multi_geoms(geoms)
    # values are expanded below as kwargs["c"]

    x = [p.x if not p.is_empty else None for p in geoms]
    y = [p.y if not p.is_empty else None for p in geoms]

    # matplotlib 1.4 does not support c=None, and < 2.0 does not support s=None
    if values is not None:
        kwargs["c"] = values
    if markersize is not None:
        kwargs["s"] = markersize

    # Add to kwargs for easier checking below.
    if color is not None:
        kwargs["color"] = color
    if marker is not None:
        kwargs["marker"] = marker
    _expand_kwargs(kwargs, multiindex)

    if "norm" not in kwargs:
        collection = ax.scatter(x, y, vmin=vmin, vmax=vmax, cmap=cmap, **kwargs)
    else:
        collection = ax.scatter(x, y, cmap=cmap, **kwargs)

    return collection

def plot_series(
    s, cmap=None, color=None, ax=None, figsize=None, aspect="auto", **style_kwds
):
    """
    Plot a GeoSeries.
    Generate a plot of a GeoSeries geometry with matplotlib.
    Parameters
    ----------
    s : Series
        The GeoSeries to be plotted. Currently Polygon,
        MultiPolygon, LineString, MultiLineString and Point
        geometries can be plotted.
    cmap : str (default None)
        The name of a colormap recognized by matplotlib. Any
        colormap will work, but categorical colormaps are
        generally recommended. Examples of useful discrete
        colormaps include:
            tab10, tab20, Accent, Dark2, Paired, Pastel1, Set1, Set2
    color : str (default None)
        If specified, all objects will be colored uniformly.
    ax : matplotlib.pyplot.Artist (default None)
        axes on which to draw the plot
    figsize : pair of floats (default None)
        Size of the resulting matplotlib.figure.Figure. If the argument
        ax is given explicitly, figsize is ignored.
    aspect : 'auto', 'equal', None or float (default 'auto')
        Set aspect of axis. If 'auto', the default aspect for map plots is 'equal'; if
        however data are not projected (coordinates are long/lat), the aspect is by
        default set to 1/cos(s_y * pi/180) with s_y the y coordinate of the middle of
        the GeoSeries (the mean of the y range of bounding box) so that a long/lat
        square appears square in the middle of the plot. This implies an
        Equirectangular projection. If None, the aspect of `ax` won't be changed. It can
        also be set manually (float) as the ratio of y-unit to x-unit.
    **style_kwds : dict
        Color options to be passed on to the actual plot function, such
        as ``edgecolor``, ``facecolor``, ``linewidth``, ``markersize``,
        ``alpha``.
    Returns
    -------
    ax : matplotlib axes instance
    """
    # if cmap is specified, create range of colors based on cmap
    values = None
    if cmap is not None:
        values = np.arange(len(s))
        if hasattr(cmap, "N"):
            values = values % cmap.N
        style_kwds["vmin"] = style_kwds.get("vmin", values.min())
        style_kwds["vmax"] = style_kwds.get("vmax", values.max())

    # decompose GeometryCollections
    geoms, multiindex, multitype = _flatten_multi_geoms(s, prefix="Multi")
    values = np.take(values, multiindex, axis=0) if cmap else None
    # expl_series = geopandas.GeoSeries(geoms)
    expl_series = geoms

    polys = []
    lines = []
    points = []
    
    geom_types = multitype
    poly_idx  = np.asarray((geom_types == "Polygon") | (geom_types == "MultiPolygon"))
    line_idx  = np.asarray((geom_types == "LineString") | (geom_types == "MultiLineString") | (geom_types == "LinearRing"))
    point_idx = np.asarray((geom_types == "Point") | (geom_types == "MultiPoint"))

    if poly_idx.all():
        polys = geoms

    elif line_idx.all():
        lines = geoms

    elif point_idx.all():
        points = geoms

    else:
        # categorise geometries into separate lists
        for ix in np.nonzero(poly_idx):
            polys.append(geoms[ix])

        for ix in np.nonzero(line_idx):
            lines.append(geoms[ix])

        for ix in np.nonzero(point_idx):
            points.append(geoms[ix])

    # plot all Polygons and all MultiPolygon components in the same collection
    if polys:
        # color overrides both face and edgecolor. As we want people to be
        # able to use edgecolor as well, pass color to facecolor
        facecolor = style_kwds.pop("facecolor", None)
        if color is not None:
            facecolor = color

        values_ = values[poly_idx] if cmap else None
        _plot_polygon_collection(
            ax, polys, values_, facecolor=facecolor, cmap=cmap, **style_kwds
        )

    # plot all LineStrings and MultiLineString components in same collection
    if lines:
        values_ = values[line_idx] if cmap else None
        _plot_linestring_collection(
            ax, lines, values_, color=color, cmap=cmap, **style_kwds
        )

    # plot all Points in the same collection
    if points:
        values_ = values[point_idx] if cmap else None
        _plot_point_collection(
            ax, points, values_, color=color, cmap=cmap, **style_kwds
        )

    plt.draw()
    return ax
