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

import math


def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = math.sqrt(mag2)
        v = tuple(n / mag for n in v)
    return v


def quat_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z


def quat_conjugate(q):
    q = normalize(q)
    w, x, y, z = q
    return (w, -x, -y, -z)


def quat_vec_mult(q1, v1):
    v1 = normalize(v1)
    q2 = (0.0,) + v1
    return quat_mult(quat_mult(q1, q2), quat_conjugate(q1))[1:]


def axis_angle_to_quat(v, theta):
    v = normalize(v)
    x, y, z = v
    theta /= 2
    w = math.cos(theta)
    x = x * math.sin(theta)
    y = y * math.sin(theta)
    z = z * math.sin(theta)
    return w, x, y, z


def quat_to_axis_angle(quat):
    w, v = quat[0], quat[1:]
    theta = math.acos(w) * 2.0
    return normalize(v), theta


def lat_lon_to_cart(lat, lon):
    x = math.cos(lat) * math.cos(lon)
    y = math.cos(lat) * math.sin(lon)
    z = math.sin(lat)
    return x, y, z


def cart_to_lat_lon(x, y, z):
    lat = math.asin(z)
    lon = math.atan2(y, x)
    return lat, lon


def test_rotate(point, axis, angle):
    """test rotate point

    :param point: point to be rotated (in lat,lon)
    :param axis: rotation axis (in lat,lon)
    :param angle: roation angle (in degrees)

    :returns: new coordinates of the point (in lat,lon)
    """
    v = lat_lon_to_cart(math.radians(point[0]), math.radians(point[1]))
    axis = lat_lon_to_cart(math.radians(axis[0]), math.radians(axis[1]))
    quat = axis_angle_to_quat(axis, math.radians(angle))
    ret = quat_vec_mult(quat, v)
    ret_lat_lon = cart_to_lat_lon(ret[0], ret[1], ret[2])
    return math.degrees(ret_lat_lon[0]), math.degrees(ret_lat_lon[1])


if __name__ == "__main__":
    print(test_rotate([12, 134], [90, 0], 34))
