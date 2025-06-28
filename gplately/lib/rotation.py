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

from . import quaternions

#
# a collection of rotation related functions
#


def cross(a, b):
    """cross product of two vectors defined only in three-dimensional space
    https://www.mathsisfun.com/algebra/vectors-cross-product.html
    https://en.wikipedia.org/wiki/Cross_product#Computing

    :param a: 3D vector, such as [1,2,3]
    :param b: 3D vector, such as [234]

    :returns: cross product 3D vector
    """

    c = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
    return c


def dot(a, b):
    """dot product of two vectors
    https://www.mathsisfun.com/algebra/vectors-dot-product.html

    :param a: a vector, such as [1,2,3,4,5] (must be the same length with b)
    :param b: a vector, such as [6,7,8,9,10] (must be the same length with a)

    :returns: dot product (a number)
    """
    return sum(i * j for i, j in zip(a, b))


def find_axis_and_angle(point_a, point_b):
    """given two points point_a and point_b in (lat, lon) format in radians,
    return axis(lat, lon) and angle in radians

    :param point_a: point a (lat, lon) in radians
    :param point_b: point b (lat, lon) in radians

    :returns: axis(lat, lon) in radians; angle in radians
    """
    # get 3D vectors
    a_v = quaternions.lat_lon_to_cart(point_a[0], point_a[1])
    b_v = quaternions.lat_lon_to_cart(point_b[0], point_b[1])

    ab_cross = quaternions.normalize(cross(a_v, b_v))
    ab_dot = dot(a_v, b_v)

    axis = quaternions.cart_to_lat_lon(*ab_cross)  # in radians
    angle = math.acos(ab_dot)  # in radians

    return axis, angle


def rotate(point, axis, angle):
    """rotate a point by axis and angle

    :param point: (lat, lon) in radians
    :param axis: (lat, lon) in radians
    :param angle: in radians

    :returns: new point(lat, lon) in radians
    """
    v = quaternions.lat_lon_to_cart(point[0], point[1])
    axis = quaternions.lat_lon_to_cart(axis[0], axis[1])
    quat = quaternions.axis_angle_to_quat(axis, angle)
    ret = quaternions.quat_vec_mult(quat, v)
    ret_lat_lon = quaternions.cart_to_lat_lon(ret[0], ret[1], ret[2])
    return (ret_lat_lon[0], ret_lat_lon[1])


def interp_two_points(point_a, point_b, num=10):
    """interpolate between two points(in radians) along the great circle
    return a list of interpolated points (lat, lon) in radians

    :param point_a: point a (lat, lon) in radians
    :param point_b: point b (lat, lon) in radians
    :param num: the number of segments

    :returns: a list of interpolated points (lat, lon) in radians
        the return list has num+1 points(num segments)
    """
    # make the two input points are not the same
    assert not (
        math.isclose(point_a[0], point_b[0]) and math.isclose(point_a[1], point_b[1])
    )
    ret = [point_a]
    axis, angle = find_axis_and_angle(point_a, point_b)
    for angle_ in [(angle / num) * i for i in range(1, num + 1)]:
        lat, lon = rotate(point_a, axis, angle_)
        ret.append((lat, lon))

    return ret


def sample_between_two_points(point_a, point_b, distances_to_point_a):
    """sample points between two points along the great circle according to a list of distances to point a

    :param point_a: point a (lat, lon) in radians
    :param point_b: point b (lat, lon) in radians
    :param distances_to_point_a: a list of distances(in radians) to point a.
        the distance must less than the distance between point a and point b.
        invalid distances will produce Nones in the return list

    :returns: a list of sampled points (lat, lon) in radians.
        The return list has the same length with the list of distances_to_point_a.
        point_a and point_b will not be included.

    """
    # make sure the two input points are not the same
    assert not (
        math.isclose(point_a[0], point_b[0]) and math.isclose(point_a[1], point_b[1])
    )
    ret = []
    axis, angle = find_axis_and_angle(point_a, point_b)
    for distance in distances_to_point_a:
        if distance < 0 or distance > angle:
            ret.append(None)
        else:
            lat, lon = rotate(point_a, axis, distance)
            ret.append((lat, lon))

    return ret


def distance(point_a, point_b, earth_radius=6378.14):
    """calculate the distance between to points (lat, lon) in radians
    return the distance in km

    :param point_a: point a (lat, lon) in radians
    :param point_b: point b (lat, lon) in radians
    :param earth_radius: earth radius

    :returns: distance in km
    """
    if math.isclose(point_a[0], point_b[0]) and math.isclose(point_a[1], point_b[1]):
        return 0
    _, angle = find_axis_and_angle(point_a, point_b)
    return radian_to_km(angle, earth_radius)


def radian_distance(point_a, point_b):
    """calculate the distance between to points (lat, lon) in radians
    return the distance in radians

    :param point_a: point a (lat, lon) in radians
    :param point_b: point b (lat, lon) in radians

    :returns: distance in radians
    """
    if math.isclose(point_a[0], point_b[0]) and math.isclose(point_a[1], point_b[1]):
        return 0
    _, angle = find_axis_and_angle(point_a, point_b)
    return angle


def km_to_radian(km, earth_radius=6378.14):
    return km / earth_radius


def radian_to_km(radian, earth_radius=6378.14):
    return radian * earth_radius


if __name__ == "__main__":
    # test
    # 1. find axis and angle between two points a and b
    # 2. use the axis and angle to rotate point a
    # 3. check if the rotation ends up at point b
    point_a = tuple((math.radians(i) for i in (45, 20)))  # lat, lon
    point_b = tuple((math.radians(i) for i in (-45, 120)))

    axis, angle = find_axis_and_angle(point_a, point_b)

    print(axis, angle)

    point_b_new = rotate(point_a, axis, angle)

    # check if the two points are the same
    print(tuple(round(math.degrees(i), 2) for i in point_b_new))
    print(tuple(round(math.degrees(i), 2) for i in point_b))

    # intepolate between point a and b
    points_in_between = interp_two_points(point_a, point_b)
    with open("interpolation.gmt", "w+") as f:
        f.write(f"{math.degrees(axis[1]):.2f} {math.degrees(axis[0]):.2f}\n")
        f.write(f"{math.degrees(axis[0]):.2f} {math.degrees(axis[1]):.2f}\n")
        for p in points_in_between:
            f.write(f"{math.degrees(p[1]):.2f} {math.degrees(p[0]):.2f}\n")

    print(distance(point_a, point_b))
    print("test finished!")
