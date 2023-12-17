from typing import (
    Any,
    Dict,
    List,
    Sequence,
    Union,
    Optional,
)

import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.axes import Axes
from matplotlib.cbook import CallbackRegistry
from matplotlib.transforms import Affine2D
from shapely.geometry import (
    LineString,
    MultiLineString,
    Point,
    box,
)
from shapely.geometry.base import (
    BaseGeometry,
    BaseMultipartGeometry,
)
from shapely.ops import linemerge

from ._axes_tools import (
    transform_distance_axes as _transform_distance_axes,
)
from ..geometry import (
    explode_geometries as _explode_geometries,
)


class SubductionTeeth:
    """Add subduction zone teeth to a map."""
    def __init__(
        self,
        left: Sequence[Union[LineString, MultiLineString]],
        right: Sequence[Union[LineString, MultiLineString]],
        ax: Optional[Union[Axes, GeoAxes]] = None,
        size: float = 6.0,
        aspect: float = 1.0,
        spacing: Optional[float] = None,
        color="black",
        **kwargs
    ):
        """Add subduction zone teeth to a map.

        Parameters
        ----------
        left, right : sequence of LineString or MultiLineString
            Shapely geometries representing the left- and right-polarity
            subduction zones.

        ax : matplotlib Axes or cartopy GeoAxes, optional
            The axes on which to plot the subduction zone teeth. If not specified,
            will use the current axes.

        size : float, default: 6.0
            Teeth size in points (alias: `markersize`).

        aspect : float, default: 1.0
            Aspect ratio of teeth triangles (height / width).

        spacing : float, optional
            Teeth spacing, in display units (usually inches). The default
            of `None` will choose a value based on the teeth size.

        color : str, default='black'
            The colour of the teeth (`markerfacecolor` and `markeredgecolor`).

        **kwargs :
            Further keyword arguments are passed to `matplotlib.pyplot.plot`.
            See `matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).
        """
        if ax is None:
            ax = plt.gca()
        self._ax = ax
        self._left = _explode_geometries(left)
        self._left_projected = None
        self._right = _explode_geometries(right)
        self._right_projected = None
        self._aspect = float(aspect)
        self._teeth = []

        if spacing is not None:
            spacing = float(spacing)
        self._spacing = spacing

        self._plot_kw = dict(kwargs)
        if "markersize" not in self._plot_kw.keys():
            self._plot_kw["markersize"] = size
        if "facecolor" in self._plot_kw.keys():
            self._plot_kw["markerfacecolor"] = self._plot_kw.pop("facecolor")
        if "edgecolor" in self._plot_kw.keys():
            self._plot_kw["markeredgecolor"] = self._plot_kw.pop("edgecolor")
        if "markerfacecolor" not in self._plot_kw.keys():
            self._plot_kw["markerfacecolor"] = color
        if "markeredgecolor" not in self._plot_kw.keys():
            self._plot_kw["markeredgecolor"] = color

        self._triangle = mpath.Path(
            vertices=[
                (-0.5, 0),
                (0.5, 0),
                (0, self.aspect),
                (-0.5, 0),
            ]
        )

        self.ax.set_xlim(*self.ax.get_xlim())
        self.ax.set_ylim(*self.ax.get_ylim())
        self._draw_teeth()

        self._callbacks = CallbackRegistry()
        self._callback_ids = set()

        def callback_func(ax):
            return self._draw_teeth(ax)

        for event in ("xlim_changed", "ylim_changed"):
            self._callback_ids.add(self._callbacks.connect(event, callback_func))
        self._ax.callbacks = self._callbacks


    def __del__(self):
        for callback_id in self._callback_ids:
            self._callbacks.disconnect(callback_id)
        del self._ax.callbacks


    def _draw_teeth(self, ax=None):
        if ax is None:
            ax = self.ax

        if self._teeth is not None:
            for i in self._teeth:
                i.remove()
        self._teeth = []

        spacing = _transform_distance_axes(self.spacing, self.ax)
        left = self.left_projected
        right = self.right_projected

        domain = (ax.transData.inverted().transform_bbox(ax.bbox))
        domain = box(domain.x0, domain.y0, domain.x1, domain.y1)
        left = domain.intersection(left)
        right = domain.intersection(right)

        for polarity, geometries in zip(
            ("left", "right"),
            (left, right),
        ):
            geometries = _explode_geometries(geometries)
            for geometry in geometries:
                if not isinstance(geometry, BaseGeometry):
                    continue
                if geometry.is_empty:
                    continue

                length = geometry.length
                geom_points = [Point(i) for i in geometry.coords]
                cumlen = np.concatenate(
                    (
                        [0.0],
                        np.cumsum(
                            [
                                geom_points[i + 1].distance(geom_points[i])
                                for i in range(len(geom_points) - 1)
                            ]
                        )
                    )
                )
                for distance in np.arange(spacing, length, spacing):
                    point = Point(geometry.interpolate(distance))
                    after = int(np.where(cumlen >= distance)[0][0])
                    before = after - 1
                    p1 = geom_points[before]
                    p2 = geom_points[after]
                    normal = np.array([p1.y - p2.y, p2.x - p1.x])
                    if polarity == "right":
                        normal *= -1.0
                    angle = np.arctan2(*(normal[::-1]))
                    marker = self._triangle.transformed(
                        Affine2D().rotate_deg(-90).rotate(angle)
                    )
                    p = ax.plot(point.x, point.y, marker=marker, **self.plot_kw)
                    self._teeth.extend(p)


    @staticmethod
    def _get_default_spacing(markersize: float, aspect: float) -> float:
        """Default spacing is approximately 2 times triangle width."""
        width = np.sqrt(2 * markersize / aspect)  # approximately
        width /= 72  # convert from points to inches
        spacing = 2 * width
        return spacing


    @property
    def ax(self):
        return self._ax


    @property
    def figure(self):
        return self.ax.figure


    @property
    def projection(self) -> Optional[ccrs.Projection]:
        if hasattr(self.ax, "projection"):
            return self.ax.projection
        return None


    @property
    def left(self) -> List[LineString]:
        return self._left


    @property
    def left_projected(self) -> List[LineString]:
        if self.projection is None:
            return self.left
        if self._left_projected is None:
            projected = [self.projection.project_geometry(i) for i in self.left]
            projected = [
                i for i in projected
                if isinstance(i, BaseGeometry) and not i.is_empty
            ]
            projected = _explode_geometries(projected)
            if len(projected) == 0:
                self._left_projected = []
                return self._left_projected
            merged = linemerge(projected)
            if hasattr(merged, "geometries"):
                self._left_projected = list(merged.geometries)
            else:
                self._left_projected = [merged]
        return self._left_projected


    @property
    def right(self) -> List[LineString]:
        return self._right


    @property
    def right_projected(self) -> List[LineString]:
        if self.projection is None:
            return self.right
        if self._right_projected is None:
            projected = [self.projection.project_geometry(i) for i in self.right]
            projected = [
                i for i in projected
                if isinstance(i, BaseGeometry) and not i.is_empty
            ]
            projected = _explode_geometries(projected)
            if len(projected) == 0:
                self._right_projected = []
                return self._right_projected
            merged = linemerge(projected)
            if hasattr(merged, "geometries"):
                self._right_projected = list(merged.geometries)
            else:
                self._right_projected = [merged]
        return self._right_projected


    @property
    def spacing(self) -> float:
        if self._spacing is None:
            return self._get_default_spacing(
                self.plot_kw["markersize"],
                self.aspect,
            )
        return self._spacing

    @spacing.setter
    def spacing(self, x: Optional[float]):
        if x is not None:
            x = float(x)
        self._spacing = x
        self._draw_teeth()


    @property
    def aspect(self):
        return self._aspect


    @property
    def plot_kw(self) -> Dict[str, Any]:
        return self._plot_kw


def plot_subduction_teeth(
    left: Sequence[Union[LineString, MultiLineString]],
    right: Sequence[Union[LineString, MultiLineString]],
    ax: Optional[Union[Axes, GeoAxes]] = None,
    size: float = 6.0,
    aspect: float = 1.0,
    spacing: Optional[float] = None,
    color="black",
    **kwargs,
):
    """Add subduction zone teeth to a map.

    Parameters
    ----------
    left, right : sequence of LineString or MultiLineString
        Shapely geometries representing the left- and right-polarity
        subduction zones.

    ax : matplotlib Axes or cartopy GeoAxes, optional
        The axes on which to plot the subduction zone teeth. If not specified,
        will use the current axes.

    size : float, default: 6.0
        Teeth size in points (alias: `markersize`).

    aspect : float, default: 1.0
        Aspect ratio of teeth triangles (height / width).

    spacing : float, optional
        Teeth spacing, in display units (usually inches). The default
        of `None` will choose a value based on the teeth size.

    color : str, default='black'
        The colour of the teeth (`markerfacecolor` and `markeredgecolor`).

    **kwargs :
        Further keyword arguments are passed to `matplotlib.pyplot.plot`.
        See `matplotlib` keyword arguments
        [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

    Returns
    -------
    SubductionTeeth

    Notes
    -----
    This function is equivalent to the newer `SubductionTeeth` class, but is
    kept for backwards compatibility.
    """
    return SubductionTeeth(
        ax=ax,
        left=left,
        right=right,
        spacing=spacing,
        aspect=aspect,
        size=size,
        color=color,
        **kwargs
    )
