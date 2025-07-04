#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
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
#
# https://sinestesia.co/blog/tutorials/python-icospheres/
#

import math

import numpy as np

from .quaternions import lat_lon_to_cart
from .rotation import distance


def _get_vertices_and_faces():
    """Return hardcoded vertices and faces (the corners of three orthogonal golden planes)."""
    r = (1.0 + math.sqrt(5.0)) / 2.0
    vertices = np.array(
        [
            [-1.0, r, 0.0],
            [1.0, r, 0.0],
            [-1.0, -r, 0.0],
            [1.0, -r, 0.0],
            [0.0, -1.0, r],
            [0.0, 1.0, r],
            [0.0, -1.0, -r],
            [0.0, 1.0, -r],
            [r, 0.0, -1.0],
            [r, 0.0, 1.0],
            [-r, 0.0, -1.0],
            [-r, 0.0, 1.0],
        ],
        dtype=float,
    )

    faces = np.array(
        [
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [5, 4, 9],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ]
    )
    return [normalize(v) for v in vertices], faces


def _get_vertices_and_faces_stripy():
    """Return vertices and faces the same as in Stripy (https://github.com/underworldcode/stripy/blob/f87add7efe4db23165d820bfb5394c83596c18da/stripy/spherical_meshes.py#L34-L47)."""
    mid_lat = np.degrees(np.arctan(0.5))  # 26.56505117707799 degrees
    vertices_stripy = np.array(
        [
            [90, 0.0],
            [mid_lat, 0.0],
            [-mid_lat, 36.0],
            [mid_lat, 72.0],
            [-mid_lat, 108.0],
            [mid_lat, 144.0],
            [-mid_lat, 180.0],
            [mid_lat, -72.0],
            [-mid_lat, -36.0],
            [mid_lat, -144.0],
            [-mid_lat, -108.0],
            [-90, 0.0],
        ]
    )

    return np.array(
        [
            lat_lon_to_cart(math.radians(v[0]), math.radians(v[1]))
            for v in vertices_stripy
        ]
    ), np.array(_find_faces(vertices_stripy))


def _find_neighbours(point, all_points):
    """Find the 5 adjacent vertices for a given vertex of an Icosahedron."""
    min_distance = None
    neighbours = []
    for idx, p in enumerate(all_points):
        dist = distance(
            (math.radians(point[0]), math.radians(point[1])),
            (math.radians(p[0]), math.radians(p[1])),
        )
        if dist > 0:
            if min_distance is not None and math.isclose(
                min_distance, dist, abs_tol=1e-04
            ):
                neighbours.append(idx)
            elif min_distance is None or dist < min_distance:
                min_distance = dist
                neighbours = [idx]
            else:
                pass

    # print(neighbours)
    assert len(neighbours) == 5
    return neighbours


def _find_faces(vertices):
    """Return 20 Icosahedron faces for the given Icosahedron vertices([lat, lon] in degrees)."""
    faces = []
    for idx, p in enumerate(vertices):
        neighbours = _find_neighbours(p, vertices)
        for i in range(len(neighbours)):
            dist_i_p = distance(
                (
                    math.radians(vertices[neighbours[i]][0]),
                    math.radians(vertices[neighbours[i]][1]),
                ),
                (math.radians(p[0]), math.radians(p[1])),
            )
            for j in range(len(neighbours)):
                if i >= j:
                    continue
                else:
                    dist_i_j = distance(
                        (
                            math.radians(vertices[neighbours[i]][0]),
                            math.radians(vertices[neighbours[i]][1]),
                        ),
                        (
                            math.radians(vertices[neighbours[j]][0]),
                            math.radians(vertices[neighbours[j]][1]),
                        ),
                    )
                    if math.isclose(dist_i_j, dist_i_p):
                        faces.append(sorted([idx, neighbours[i], neighbours[j]]))

    unique_faces = []
    for f in faces:
        if f not in unique_faces:
            unique_faces.append(f)
    assert len(unique_faces) == 20
    return unique_faces


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def bisect(vertices, faces, level=3):
    if level == 0:
        return vertices, faces

    new_vertices = vertices
    new_faces = None
    for face in faces:
        idx_v0 = face[0]
        idx_v1 = face[1]
        idx_v2 = face[2]
        v0 = vertices[idx_v0]
        v1 = vertices[idx_v1]
        v2 = vertices[idx_v2]
        v3 = normalize(0.5 * (v0 + v1))
        v4 = normalize(0.5 * (v1 + v2))
        v5 = normalize(0.5 * (v2 + v0))
        new_vertices = np.append(new_vertices, [v3, v4, v5], axis=0)
        v_len = new_vertices.shape[0]
        # print(v_len)
        idx_v3 = v_len - 3
        idx_v4 = v_len - 2
        idx_v5 = v_len - 1
        tmp = np.array(
            [
                [idx_v0, idx_v3, idx_v5],
                [idx_v3, idx_v1, idx_v4],
                [idx_v4, idx_v2, idx_v5],
                [idx_v3, idx_v4, idx_v5],
            ]
        )
        if new_faces is None:
            new_faces = tmp
        else:
            new_faces = np.append(new_faces, tmp, axis=0)
        # print(new_faces.shape)
    return bisect(new_vertices, new_faces, level - 1)


def get_mesh(level=5, use_stripy_icosahedron=False):
    """Return the Icospheres mesh."""
    if use_stripy_icosahedron:
        vertices, faces = _get_vertices_and_faces_stripy()
    else:
        vertices, faces = _get_vertices_and_faces()
    return bisect(vertices, faces, level)


def xyz2lonlat(x, y, z):
    """
    x = R * cos(lat) * cos(lon)
    y = R * cos(lat) * sin(lon)
    z = R *sin(lat)
    lat = asin(z / R)
    lon = atan2(y, x)
    """
    lat = np.arcsin(z)
    # print(z,lat, x, y)
    lon = np.arctan2(y, x)
    return np.rad2deg(lon), np.rad2deg(lat)
