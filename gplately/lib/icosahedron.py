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


vertices = [normalize(v) for v in vertices]


def get_mesh(level=5):
    """
    call this function to get the mesh
    """
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


if __name__ == "__main__":
    seen = set()

    vertices_0, faces_0 = get_mesh(6)
    print(vertices_0.shape, faces_0.shape)  # type: ignore

    with open("sphere_mesh.gmt", "w+") as f:
        for v in vertices_0:
            lon, lat = xyz2lonlat(v[0], v[1], v[2])
            line = f"{lon:0.2f} {lat:0.2f}\n"
            if line in seen:
                continue
            f.write(line)
            seen.add(line)
    print("done")
