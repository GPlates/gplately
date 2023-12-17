import cartopy.crs as ccrs
import numpy as np
from matplotlib.axes import Axes


def meridian_from_ax(ax: Axes) -> float:
    if hasattr(ax, "projection") and isinstance(ax.projection, ccrs.Projection):
        proj = ax.projection
        return meridian_from_projection(projection=proj)
    return 0.0


def meridian_from_projection(projection: ccrs.Projection) -> float:
    x = np.mean(projection.x_limits)
    y = np.mean(projection.y_limits)
    return ccrs.PlateCarree().transform_point(x, y, projection)[0]


def transform_distance_axes(d: float, ax: Axes, inverse=False) -> float:
    axes_bbox = ax.get_position()
    fig_bbox = ax.figure.bbox_inches  # display units (inches)

    axes_width = axes_bbox.width * fig_bbox.width
    axes_height = axes_bbox.height * fig_bbox.height

    # Take mean in case they're different somehow
    xlim = ax.get_xlim()
    x_factor = np.abs(xlim[1] - xlim[0]) / axes_width  # map units per display unit
    ylim = ax.get_ylim()
    y_factor = np.abs(ylim[1] - ylim[0]) / axes_height
    factor = 0.5 * (x_factor + y_factor)
    if inverse:
        return d / factor
    return d * factor
