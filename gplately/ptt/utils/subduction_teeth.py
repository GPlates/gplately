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
#   You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

import numpy as np


def tesselate_triangles(
    shapefilename, tesselation_radians, triangle_base_length, triangle_aspect=1.0
):
    """
    Place subduction teeth along line segments within a MultiLineString shapefile

    Parameters
    ----------
    shapefilename: str
        path to shapefile

    tesselation_radians: float
        tesselation radians

    triangle_base_length: float
        length of base

    triangle_aspect: float
        aspect ratio, setting triangle_aspect to -1 reverses the tooth direction

    Returns
    -------
    X_points : (n,3)
        array of triangle x points
    Y_points : (n,3)
        array of triangle y points
    """

    import shapefile

    shp = shapefile.Reader(shapefilename)

    tesselation_degrees = np.degrees(tesselation_radians)
    triangle_pointsX = []
    triangle_pointsY = []

    for i in range(len(shp)):
        pts = np.array(shp.shape(i).points)

        cum_distance = 0.0
        for p in range(len(pts) - 1):

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

    shp.close()
    return np.array(triangle_pointsX), np.array(triangle_pointsY)
