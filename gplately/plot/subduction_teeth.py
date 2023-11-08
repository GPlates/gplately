from typing import (
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
    def __init__(
        self,
        ax: Union[Axes, GeoAxes],
        left: Sequence[Union[LineString, MultiLineString]],
        right: Sequence[Union[LineString, MultiLineString]],
        spacing: Optional[float] = None,
        aspect: float = 1.0,
        size: float = 6.0,
        color="black",
        **kwargs
    ):
        """Plot subduction teeth onto a standard map Projection.

        Notes
        -----
        Subduction teeth are tessellated from `PlotTopologies` object attributes `trench_left` and
        `trench_right`, and transformed into Shapely polygons for plotting.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        spacing : float, optional
            Teeth spacing, in display units (usually inches). The default
            of `None` will choose a value based on the area of the plot.

        size : float, default: 7.0
            Teeth size (alias: `markersize`).

        aspect : float, default: 1.0
            Aspect ratio of teeth triangles (height / width).

        color : str, default='black'
            The colour of the teeth (`markerfacecolor` and `markeredgecolor`).

        **kwargs :
            Further keyword arguments are passed to `matplotlib.pyplot.plot`.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).
        """
        self._ax = ax
        self._left = list(left)
        self._left_projected = None
        self._right = list(right)
        self._right_projected = None
        self._aspect = float(aspect)

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

        callbacks = CallbackRegistry()

        def callback_func(ax):
            return self._draw_teeth(ax)

        for event in ("xlim_changed", "ylim_changed"):
            callbacks.connect(event, callback_func)
        ax.callbacks = callbacks


    def _draw_teeth(self, ax=None):
        if ax is None:
            ax = self.ax

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
            if isinstance(geometries, BaseMultipartGeometry):
                geometries = list(geometries.geoms)
            elif isinstance(geometries, BaseGeometry):
                geometries = [geometries]
            for geometry in geometries:
                if not isinstance(geometry, BaseGeometry):
                    continue
                if geometry.is_empty:
                    continue

                length = geometry.length
                tessellated_x = []
                tessellated_y = []
                for distance in np.arange(spacing, length, spacing):
                    point = Point(geometry.interpolate(distance))
                    tessellated_x.append(point.x)
                    tessellated_y.append(point.y)
                tessellated_x = np.array(tessellated_x)
                tessellated_y = np.array(tessellated_y)

                for i in range(len(tessellated_x) - 1):
                    normal_x = tessellated_y[i] - tessellated_y[i + 1]
                    normal_y = tessellated_x[i + 1] - tessellated_x[i]
                    midpoint = np.array((tessellated_x[i], tessellated_y[i]))
                    if polarity == "right":
                        normal_x *= -1
                        normal_y *= -1
                    angle = np.arctan2(normal_y, normal_x)
                    marker = self._triangle.transformed(
                        Affine2D().rotate_deg(-90).rotate(angle)
                    )
                    ax.plot(*midpoint, marker=marker, **self.plot_kw)

    def _get_default_spacing(self, ax=None) -> float:
        if ax is None:
            ax = self.ax
        axes_bbox = ax.get_position()
        fig_bbox = ax.figure.bbox_inches
        axes_width = axes_bbox.width * fig_bbox.width
        axes_height = axes_bbox.height * fig_bbox.height
        # axes_area = axes_width * axes_height
        # spacing = (axes_area / 1000.0) ** 0.5
        spacing = axes_height / 30
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
    def left(self):
        return self._left

    @property
    def left_projected(self) -> List[Union[LineString, MultiLineString]]:
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
            self._left_projected = linemerge(projected)
        return self._left_projected

    @property
    def right(self):
        return self._right

    @property
    def right_projected(self) -> List[Union[LineString, MultiLineString]]:
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
            self._right_projected = linemerge(projected)
        return self._right_projected

    @property
    def spacing(self) -> float:
        if self._spacing is None:
            return self._get_default_spacing()
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
    def plot_kw(self):
        return self._plot_kw


def plot_subduction_teeth(
    left,
    right,
    spacing=None,
    ax=None,
    aspect=1.0,
    size=7.0,
    color="black",
    **kwargs,
):
    """Plot subduction teeth onto a standard map Projection.

    Notes
    -----
    Subduction teeth are tessellated from `PlotTopologies` object attributes `trench_left` and
    `trench_right`, and transformed into Shapely polygons for plotting.

    Parameters
    ----------
    ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
        A subclass of `matplotlib.axes.Axes` which represents a map Projection.
        The map should be set at a particular Cartopy projection.

    spacing : float, optional
        Teeth spacing, in display units (usually inches). The default
        of `None` will choose a value based on the area of the plot.

    size : float, default: 7.0
        Teeth size (alias: `markersize`).

    aspect : float, default: 1.0
        Aspect ratio of teeth triangles (height / width).

    color : str, default='black'
        The colour of the teeth (`markerfacecolor` and `markeredgecolor`).

    **kwargs :
        Further keyword arguments are passed to `matplotlib.pyplot.plot`.
        See `Matplotlib` keyword arguments
        [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).
    """
    if ax is None:
        ax = plt.gca()

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
