from itertools import product

import numpy as np
import pytest
from gplately import EARTH_RADIUS
from gplately.geometry import *
from shapely.geometry import (
    Point,
    LineString,
    Polygon,
)

from conftest import (
    test_geometry_n_points as N_POINTS,
    test_geometry_origins as origins,
    test_geometry_radii as RADII,
)


@pytest.mark.parametrize("origin,radius", product(origins, RADII))
def test_polygon_conversion(origin, radius, n_points=N_POINTS):
    lons, lats = _get_circle(origin, radius, n=n_points)
    circle = Polygon(zip(lons, lats))
    converted = PolygonOnSphere.from_shapely(circle)
    reconverted = converted.to_shapely(
        central_meridian=converted.get_centroid().to_lat_lon()[1],
    )
    assert np.allclose(
        circle.exterior.coords,
        reconverted.exterior.coords,
    ), "Polygon coordinates not equal"


@pytest.mark.parametrize("origin", origins)
def test_polygon_areas(origin, radii=RADII, n_points=N_POINTS):
    areas = [_get_circle_areas(r, origin=origin) for r in radii]
    geometric_areas, pygplates_areas = zip(*areas)
    assert np.allclose(
        geometric_areas,
        pygplates_areas,
    ), "Polygon areas not equal"


@pytest.mark.parametrize("origin", origins)
def test_polygon_perimeters(origin, radii=RADII, n_points=N_POINTS):
    perimeters = [_get_circle_perimeters(r, origin=origin) for r in radii]
    geometric_perimeters, pygplates_perimeters = zip(*perimeters)
    assert np.allclose(
        geometric_perimeters,
        pygplates_perimeters,
    ), "Polygon perimeters not equal"


def _get_circle_perimeters(radius, origin=(0, 0), n=N_POINTS):
    lons, lats = _get_circle(origin, radius, n=n)
    circle = Polygon(zip(lons, lats))
    converted = PolygonOnSphere.from_shapely(circle)

    pygplates_perimeter = converted.get_arc_length() * EARTH_RADIUS
    geometric_perimeter = (
        2 * np.pi
        * EARTH_RADIUS * np.sin(radius / EARTH_RADIUS)
    )
    return geometric_perimeter, pygplates_perimeter


def _get_circle_areas(radius, origin=(0, 0), n=N_POINTS):
    lons, lats = _get_circle(origin, radius, n=n)
    circle = Polygon(zip(lons, lats))
    converted = PolygonOnSphere.from_shapely(circle)

    pygplates_area = converted.get_area() * (EARTH_RADIUS ** 2)
    geometric_area = (
        2 * np.pi *
        (EARTH_RADIUS ** 2)
        * (1 - np.cos(radius / EARTH_RADIUS))
    )
    return geometric_area, pygplates_area


def _point_from_origin_distance_azimuth(
    origin,
    distance,
    azimuth,
    radius=EARTH_RADIUS,
):
    lon1, lat1 = origin
    lon1 = np.deg2rad(lon1)
    lat1 = np.deg2rad(lat1)
    azimuth = np.deg2rad(azimuth)

    b = distance / radius
    a = np.arccos(
        np.cos(b) * np.cos(0.5 * np.pi - lat1)
        + np.sin(0.5 * np.pi - lat1) * np.sin(b) * np.cos(azimuth)
    )
    B = np.arcsin(
        np.sin(b) * np.sin(azimuth)
        / np.sin(a)
    )
    lat2 = np.rad2deg(0.5 * np.pi - a)
    lon2 = np.rad2deg(B + lon1)

    return lon2, lat2


def _get_circle(origin, radius, earth_radius=EARTH_RADIUS, n=N_POINTS):
    azimuth = np.linspace(0, 360.0, n, endpoint=False)
    lons, lats = _point_from_origin_distance_azimuth(
        origin,
        radius,
        azimuth,
        earth_radius,
    )
    return lons, lats
