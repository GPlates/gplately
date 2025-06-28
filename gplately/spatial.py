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

"""This sub-module contains spatial tools for calculating distances on the Earth."""

import numpy as np
import pygplates
from scipy.spatial import cKDTree as _KDTree

EARTH_RADIUS = pygplates.Earth.mean_radius_in_kms


def lonlat2xyz(lon, lat, degrees=True):
    """Convert lon / lat (radians) for spherical triangulation into Cartesian (x,y,z) coordinates on the unit sphere.

    Parameters
    ----------
    lon, lat : lists
        Longitudes and latitudes of feature points in radians.

    Returns
    -------
    xs, ys, zs : lists
        Cartesian coordinates of each feature point in all 3 dimensions.
    """
    if degrees:
        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)
    cosphi = np.cos(lat)
    xs = cosphi * np.cos(lon)
    ys = cosphi * np.sin(lon)
    zs = np.sin(lat)
    return xs, ys, zs


def xyz2lonlat(x, y, z, validate=False, degrees=True):
    """Converts Cartesian (x,y,z) representation of points (on the unit sphere) for spherical triangulation into
    lon / lat (radians).

    Note: No check is made here that (x,y,z) are unit vectors - it is assumed.

    Parameters
    ----------
    x, y, z : lists
        Cartesian coordinates of each feature point in all 3 dimensions.

    Returns
    -------
    lon, lat : lists
        Longitudes and latitudes of feature points in radians.

    Notes
    -----
    No check is made here that (x,y,z) are unit vectors, unless validate=True is specified
    """
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    if validate:
        mags = np.sqrt(x**2 + y**2 + z**2)
        ones = np.full_like(mags, 1)
        if not np.all(np.equal(mags, ones)):
            raise ValueError("All (x, y, z) must be unit vectors")
    lons = np.arctan2(y, x)
    lats = np.arcsin(z)
    if degrees:
        lons = np.rad2deg(lons)
        lats = np.rad2deg(lats)
    if lons.size == 1:
        lons = np.atleast_1d(np.squeeze(lons))[0]
    if lats.size == 1:
        lats = np.atleast_1d(np.squeeze(lats))[0]
    return lons, lats


def haversine_distance(lon1, lon2, lat1, lat2, degrees=True):
    """Computes the Haversine distance (the shortest distance on the surface of an ideal spherical Earth) between two
    points given their latitudes and longitudes.

    Sources
    -------
    https://en.wikipedia.org/wiki/Haversine_formula
    https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters

    Parameters
    ----------
    lon1, lon2 : float
        Longitudes of both points

    lat1, lat2 : float
        Latitudes of both points

    Returns
    -------
    d : float
        The Haversine distance in metres.

    Notes
    -----
    Default behaviour assumes values in degrees; for radians specify degrees=False
    """
    if degrees:
        dLat = np.deg2rad(lat2) - np.deg2rad(lat1)
        dLon = np.deg2rad(lon2) - np.deg2rad(lon1)
    else:
        dLat = lat2 - lat1
        dLon = lon2 - lon1
    a = (
        np.sin(dLat / 2) ** 2
        + np.cos(lat1 * np.pi / 180)
        * np.cos(lat2 * np.pi / 180)
        * np.sin(dLon / 2) ** 2
    )
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    d = EARTH_RADIUS * c
    return d * 1000


def cartesian_distance(
    lon1, lon2, lat1, lat2, degrees=True, k=1, return_neighbours=False
):

    from .tools import EARTH_RADIUS, lonlat2xyz

    x1, y1, z1 = lonlat2xyz(lon1, lat1, degrees)
    x2, y2, z2 = lonlat2xyz(lon2, lat2, degrees)

    xyz1 = np.c_[x1, y1, z1]
    xyz2 = np.c_[x2, y2, z2]

    tree = _KDTree(xyz1)
    dist, neighbours = tree.query(xyz2, k=k)
    dist *= EARTH_RADIUS
    if return_neighbours:
        return dist, neighbours
    else:
        return dist


def great_circle_distance(
    lon1, lon2, lat1, lat2, degrees=True, k=1, return_neighbours=False
):

    from .tools import EARTH_RADIUS, lonlat2xyz

    x1, y1, z1 = lonlat2xyz(lon1, lat1, degrees)
    x2, y2, z2 = lonlat2xyz(lon2, lat2, degrees)

    xyz1 = np.c_[x1, y1, z1]
    xyz2 = np.c_[x2, y2, z2]

    tree = _KDTree(xyz1)
    dist, neighbours = tree.query(xyz2, k=k)

    if k == 1:
        neighbours_ = neighbours.reshape(-1, 1)
    else:
        neighbours_ = neighbours

    ## Now find the angular separation / great circle distance: dlatlon
    xyz1_neighbours = xyz1[neighbours_].transpose(0, 2, 1)
    xyz2_ext = np.repeat(xyz2, k, axis=1).reshape(xyz1_neighbours.shape)

    angles = np.arccos((xyz2_ext * xyz1_neighbours).sum(axis=1))
    angles *= EARTH_RADIUS

    if return_neighbours:
        return angles, neighbours
    else:
        return angles


def geocentric_radius(lat, degrees=True):
    """Calculates the latitude-dependent radius of an ellipsoid Earth.

    Parameters
    ----------
    lat : float
        The geodetic latitude at which to calculate the Earth's radius
    degrees : bool, default=True
        Specify whether the given latitude is in degrees.

    Returns
    -------
    earth_radius : float
        The Earth's geocentric radius (in metres) at the given geodetic latitude.
    """
    if degrees:
        rlat = np.radians(lat)
    else:
        rlat = lat

    coslat = np.cos(rlat)
    sinlat = np.sin(rlat)
    r1 = 6384.4e3
    r2 = 6352.8e3
    num = (r1**2 * coslat) ** 2 + (r2**2 * sinlat) ** 2
    den = (r1 * coslat) ** 2 + (r2 * sinlat) ** 2
    earth_radius = np.sqrt(num / den)
    return earth_radius
