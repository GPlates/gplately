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

"""
This file contains deprecated functions and classes. They will be removed in future versions of gplately.
Keep the deprecated functions and classes here for a while to avoid breaking existing code that depends on them.
please update your code to use the new functions and classes instead.
"""

import logging
import warnings

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false
import numpy as np
import shapely

logger = logging.getLogger("gplately")


# PlotTopologies.misc_transforms()
def _misc_transforms_impl(self):
    """
    Deprecated! DO NOT USE.
    """

    warnings.warn(
        "Deprecated! The 'misc_transforms' property will be removed in the future GPlately releases. "
        "Update your workflow to use the 'transforms' property instead, "
        "otherwise your workflow will not work with the future GPlately releases.",
        DeprecationWarning,
        stacklevel=4,
    )
    return self._transforms


# PlotTopologies.plot_misc_transforms()
def _plot_misc_transforms_impl(self, ax, color="black", **kwargs):
    """
    Deprecated! DO NOT USE.
    """
    warnings.warn(
        "Deprecated! The 'plot_misc_transforms' function will be removed in the future GPlately releases. "
        "Update your workflow to use the 'plot_transforms' function instead, "
        "otherwise your workflow will not work with the future GPlately releases.",
        DeprecationWarning,
        stacklevel=4,
    )
    self.plot_transforms(ax=ax, color=color, **kwargs)


# PlotTopologies.get_misc_transforms()
def _get_misc_transforms_impl(
    self,
    central_meridian=0.0,
    tessellate_degrees=None,
):
    """
    Deprecated! DO NOT USE.
    """

    warnings.warn(
        "Deprecated! The 'get_misc_transforms' function will be removed in the future GPlately releases. "
        "Update your workflow to use the 'get_transforms' function instead, "
        "otherwise your workflow will not work with the future GPlately releases.",
        DeprecationWarning,
        stacklevel=4,
    )
    return self.get_transforms(
        central_meridian=central_meridian, tessellate_degrees=tessellate_degrees
    )


# PlotTopologies.plot_ridges_and_transforms()
def _plot_ridges_and_transforms_impl(self, ax, color="black", **kwargs):
    """
    Deprecated! DO NOT USE!
    """
    warnings.warn(
        "Deprecated! The 'plot_ridges_and_transforms' function will be removed in the future GPlately releases. "
        "Update your workflow to use the 'plot_ridges' and 'plot_transforms' functions instead, "
        "otherwise your workflow will not work with the future GPlately releases.",
        DeprecationWarning,
        stacklevel=4,
    )
    logger.debug(
        "The 'plot_ridges_and_transforms' function has been changed since GPlately 1.3.0. "
        "You need to check your workflow to make sure the new 'plot_ridges_and_transforms' function still suits your purpose. "
        "In earlier releases of GPlately, we used an algorithm to identify the 'ridges' and 'transforms' within the gpml:MidOceanRidge features. "
        "Unfortunately, the algorithm did not work very well. So we have removed the algorithm and now the 'plot_ridges_and_transforms' function plots all the features "
        "which are labelled as gpml:Transform or gpml:MidOceanRidge in the reconstruction model."
    )  # use logger.debug to make the message less aggressive

    self.plot_ridges(ax, color=color, **kwargs)
    self.plot_transforms(ax, color=color, **kwargs)


# PlotTopologies.plot_plate_id()
def _plot_plate_id_impl(self, *args, **kwargs):
    """Deprecated! DO NOT USE!

    The function name plot_plate_id() is bad and should be changed.
    The new name is plot_plate_polygon_by_id().
    For backward compatibility, we allow users to use the old name in their legcy code for now.
    No new code should call this function.
    """
    # TODO: remove this function
    warnings.warn(
        "The class method plot_plate_id is deprecated and will be removed in the future soon. Use plot_plate_polygon_by_id instead.",
        DeprecationWarning,
        stacklevel=4,
    )
    return self.plot_plate_polygon_by_id(*args, **kwargs)


# PlotTopologies.ridge_transforms()
def _ridge_transforms_impl(self):
    """
    Deprecated! DO NOT USE!
    """

    warnings.warn(
        "Deprecated! DO NOT USE!"
        "The 'ridge_transforms' property will be removed in the future GPlately releases. "
        "Update your workflow to use the 'ridges' and 'transforms' properties instead, "
        "otherwise your workflow will not work with the future GPlately releases.",
        DeprecationWarning,
        stacklevel=4,
    )
    logger.debug(
        "The 'ridge_transforms' property has been changed since GPlately 1.3.0. "
        "You need to check your workflow to make sure the new 'ridge_transforms' property still suits your purpose. "
        "In earlier releases of GPlately, the 'ridge_transforms' property contains only the features "
        "which are labelled as gpml:MidOceanRidge in the reconstruction model. "
        "Now, the 'ridge_transforms' property contains both gpml:Transform and gpml:MidOceanRidge features."
    )
    return self._ridges + self._transforms


# PlotTopologies.get_ridges_and_transforms()
def _get_ridges_and_transforms_impl(self, central_meridian=0.0, tessellate_degrees=1):
    """
    Deprecated! DO NOT USE.
    """
    warnings.warn(
        "Deprecated! The 'get_ridges_and_transforms' function will be removed in the future GPlately releases. "
        "Update your workflow to use the 'get_ridges' and 'get_transforms' functions instead, "
        "otherwise your workflow will not work with the future GPlately releases.",
        DeprecationWarning,
        stacklevel=2,
    )
    logger.debug(
        "The 'get_ridges_and_transforms' function has been changed since GPlately 1.3.0. "
        "You need to check your workflow to make sure the new 'get_ridges_and_transforms' function still suits your purpose. "
        "In earlier releases of GPlately, we used an algorithm to identify the 'ridges' and 'transforms' within the gpml:MidOceanRidge features. "
        "Unfortunately, the algorithm did not work very well. So we have removed the algorithm and now the 'get_ridges_and_transforms' function returns all the features "
        "which are labelled as gpml:MidOceanRidge or gpml:Transform in the reconstruction model."
    )  # use logger.debug to make the message less aggressive

    return self.get_feature(
        self._ridges + self._transforms,
        central_meridian=central_meridian,
        tessellate_degrees=tessellate_degrees,
    )


# PlotTopologies._plot_subduction_teeth_deprecated()
def _plot_subduction_teeth_deprecated(
    self, ax, spacing=0.1, size=2.0, aspect=1, color="black", **kwargs
):
    """Plot subduction teeth on a map.

    Notes
    -----
    Subduction teeth are tessellated from :attr:`gplately.PlotTopologies.trench_left` and
    :attr:`gplately.PlotTopologies.trench_right`, and transformed into Shapely polygons for plotting.

    Parameters
    ----------
    ax :
        Cartopy ax or pygmt figure object.

    spacing : float, default=0.1
        The tessellation threshold (in radians). Parametrises subduction tooth density.
        Triangles are generated only along line segments with distances that exceed the given threshold ``spacing``.

    size : float, default=2.0
        Length of teeth triangle base.

    aspect : float, default=1
        Aspect ratio of teeth triangles. Ratio is 1.0 by default.

    color : str, default='black'
        The colour of the teeth. By default, it is set to black.

    **kwargs :
        Keyword arguments for parameters such as 'alpha', etc. for plotting subduction tooth polygons.
        See `Matplotlib` keyword arguments `here <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html>`__.
    """

    # add Subduction Teeth
    subd_xL, subd_yL = self._tessellate_triangles(
        self.trench_left,
        tesselation_radians=spacing,
        triangle_base_length=size,
        triangle_aspect=-aspect,
    )
    subd_xR, subd_yR = self._tessellate_triangles(
        self.trench_right,
        tesselation_radians=spacing,
        triangle_base_length=size,
        triangle_aspect=aspect,
    )

    teeth = []
    for tX, tY in zip(subd_xL, subd_yL):
        triangle_xy_points = np.c_[tX, tY]
        shp = shapely.geometry.Polygon(triangle_xy_points)
        teeth.append(shp)

    for tX, tY in zip(subd_xR, subd_yR):
        triangle_xy_points = np.c_[tX, tY]
        shp = shapely.geometry.Polygon(triangle_xy_points)
        teeth.append(shp)

    return ax.add_geometries(teeth, crs=self.base_projection, color=color, **kwargs)


# subduction teeth
# PlotTopologies._tessellate_triangles()
def _tessellate_triangles_impl(
    self, features, tesselation_radians, triangle_base_length, triangle_aspect=1.0
):
    """Places subduction teeth along subduction boundary line segments within a MultiLineString shapefile.

    Parameters
    ----------
    shapefilename  : str
        Path to shapefile containing the subduction boundary features

    tesselation_radians : float
        Parametrises subduction teeth density. Triangles are generated only along line segments with distances
        that exceed the given threshold tessellation_radians.

    triangle_base_length : float
        Length of teeth triangle base

    triangle_aspect : float, default=1.0
        Aspect ratio of teeth triangles. Ratio is 1.0 by default.

    Returns
    -------
    X_points : (n,3) array
        X points that define the teeth triangles
    Y_points : (n,3) array
        Y points that define the teeth triangles
    """

    tesselation_degrees = np.degrees(tesselation_radians)
    triangle_pointsX = []
    triangle_pointsY = []

    date_line_wrapper = pygplates.DateLineWrapper()  # type: ignore

    for feature in features:
        cum_distance = 0.0

        for geometry in feature.get_geometries():
            wrapped_lines = date_line_wrapper.wrap(geometry)
            for line in wrapped_lines:
                pts = np.array(
                    [(p.get_longitude(), p.get_latitude()) for p in line.get_points()]
                )

                for p in range(0, len(pts) - 1):
                    A = pts[p]
                    B = pts[p + 1]

                    AB_dist = B - A
                    AB_norm = AB_dist / np.hypot(*AB_dist)
                    cum_distance += np.hypot(*AB_dist)

                    # create a new triangle if cumulative distance is exceeded.
                    if cum_distance >= tesselation_degrees:
                        C = A + triangle_base_length * AB_norm

                        # find normal vector
                        AD_dist = np.array([AB_norm[1], -AB_norm[0]])
                        AD_norm = AD_dist / np.linalg.norm(AD_dist)

                        C0 = A + 0.5 * triangle_base_length * AB_norm

                        # project point along normal vector
                        D = C0 + triangle_base_length * triangle_aspect * AD_norm

                        triangle_pointsX.append([A[0], C[0], D[0]])
                        triangle_pointsY.append([A[1], C[1], D[1]])

                        cum_distance = 0.0

    return np.array(triangle_pointsX), np.array(triangle_pointsY)
