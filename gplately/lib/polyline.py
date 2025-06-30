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

from . import rotation


def _find_next_point(line, distance, strict=False, accum_length=0):
    """recursive function to find the next point along the polyline

    :param line: a list of (lat, lon) in radians
    :param distance: the distance to the next point in radians
    :param strict: in strict mode, we guarantee the distance. sacrifice the end point if necessary
    :accum_length: parameter to pass down the accumulated length for this recursive function

    :returns: (next point, new line)
        next point (lat,lon) in radians. None, if no more point
        new line -- the line after taking out this point
    """
    if len(line) < 2:
        return None, []

    start_point = line[0]
    second_point = line[1]

    dd = rotation.radian_distance(start_point, second_point)

    if math.isclose(dd + accum_length, distance):
        return second_point, line[1:]
    elif (dd + accum_length) > distance:
        pp = rotation.sample_between_two_points(
            start_point, second_point, [distance - accum_length]
        )[0]
        assert pp
        return pp, [pp] + line[1:]
    else:
        if len(line[1:]) > 1:
            return _find_next_point(
                line[1:], distance, strict=strict, accum_length=dd + accum_length
            )
        else:
            if strict:
                # in strict mode, we guarantee the distance. It means we will sacrifice the end point if necessary
                return None, []
            else:
                return second_point, line[1:]  # reached the last point


def discretize_polyline(line, distance, strict=False):
    """travel along a polyline; take a sampling point for each `distance`; return a list of the sampling points.

    :param line: a list of (lat, lon) in radians
    :param distance: the distance of the adjacent sampling points in radians
    :param strict: in strict mode, we guarantee the distance, sacrifice the end point if necessary.
        in unstrict mode, we always keep the end point

    :returns: the new line
        a list of (lat, lon) in radians.
        If the input line is shorter than distance,
            in strict mode, return [];
            in unstrict mode, return [start_point, end_point]
    """
    new_line = [line[0]]

    ll = line

    while True:
        pp, ll = _find_next_point(ll, distance, strict=strict)
        if pp:
            new_line.append(pp)
        else:
            break

    if len(new_line) > 1:
        return new_line
    else:
        return []
