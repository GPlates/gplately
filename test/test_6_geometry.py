from itertools import product

import numpy as np
import pygplates
import pytest
from gplately import EARTH_RADIUS
from gplately.geometry import *
from shapely.geometry import (
    LineString,
    MultiPoint,
    Point,
    Polygon,
)
from shapely.geometry.base import BaseMultipartGeometry

from conftest import (
    test_geometry_n_points as N_POINTS,
    test_geometry_origins as origins,
    test_geometry_radii as RADII,
    test_geometry_azimuths as azimuths,
)


@pytest.mark.parametrize("lat", (-100, 95))
def test_invalid_point(lat):
    p = Point(0, lat)
    with pytest.raises(pygplates.InvalidLatLonError):
        PointOnSphere.from_shapely(p)


def test_multipoint_conversion(n_points=N_POINTS):
    n_lats = int(np.sqrt(n_points / 2))
    point_lons = np.linspace(-179, 179, n_lats * 2)
    point_lats = np.linspace(-89, 89, n_lats)
    grid_lons, grid_lats = np.meshgrid(point_lons, point_lats)

    mp_shapely = MultiPoint(
        np.column_stack(
            [np.ravel(i) for i in (grid_lons, grid_lats)]
        )
    )
    mp_pgp = MultiPointOnSphere.from_shapely(mp_shapely)

    a1 = np.row_stack([i.coords for i in mp_shapely.geoms])
    a2 = np.fliplr(mp_pgp.to_lat_lon_array())
    assert np.allclose(
        a1,
        a2,
    ), "Multi-point coordinates not equal"


@pytest.mark.parametrize(
    "origin,azimuth,distance",
    product(origins, azimuths, RADII),
)
def test_polyline_conversion(origin, azimuth, distance, n_points=N_POINTS):
    lons, lats = _point_from_origin_distance_azimuth(
        origin,
        np.linspace(0, distance, n_points),
        azimuth,
    )
    line = LineString(zip(lons, lats))
    converted = PolylineOnSphere.from_shapely(line)
    reconverted = converted.to_shapely(
        central_meridian=converted.get_centroid().to_lat_lon()[1],
    )

    # Ensure consistency
    a1 = np.array(line.coords)
    inds1 = np.where(a1[:, 0] < -180)
    a1[inds1, 0] += 360
    inds2 = np.where(a1[:, 0] > 180)
    a1[inds2, 0] -= 360

    a2 = np.array(reconverted.coords)
    inds3 = np.where(a2[:, 0] < -180)
    a2[inds3, 0] += 360
    inds4 = np.where(a2[:, 0] > 180)
    a2[inds4, 0] -= 360

    assert np.allclose(a1, a2), "Polyline coordinates not equal"


@pytest.mark.parametrize(
    "origin,azimuth,distance",
    product(origins, azimuths, RADII),
)
def test_polyline_length(origin, azimuth, distance, n_points=N_POINTS):
    lons, lats = _point_from_origin_distance_azimuth(
        origin,
        np.linspace(0, distance, n_points),
        azimuth,
    )
    line = LineString(zip(lons, lats))
    converted = PolylineOnSphere.from_shapely(line)
    assert np.allclose(
        distance,
        converted.get_arc_length() * EARTH_RADIUS,
    ), "Polyline length not consistent"


@pytest.mark.parametrize("origin,radius", product(origins, RADII))
def test_polygon_conversion(origin, radius, n_points=N_POINTS):
    lons, lats = _get_circle(origin, radius, n=n_points)
    circle = Polygon(zip(lons, lats))
    converted = PolygonOnSphere.from_shapely(circle)
    reconverted = converted.to_shapely(
        central_meridian=converted.get_interior_centroid().to_lat_lon()[1],
    )
    assert np.allclose(
        circle.exterior.coords,
        reconverted.exterior.coords,
    ), "Polygon coordinates not equal"


def test_polygon_splitting(n_points=N_POINTS):
    origin = (-179, 0)
    radius = 1000
    lons, lats = _get_circle(origin, radius, n=n_points)
    circle = Polygon(zip(lons, lats))
    converted = PolygonOnSphere.from_shapely(circle)
    split = converted.to_shapely(central_meridian=0)
    unsplit = converted.to_shapely(
        central_meridian=converted.get_interior_centroid().to_lat_lon()[1]
    )

    assert isinstance(
        split,
        BaseMultipartGeometry,
    ), "Polygon splitting failed (incorrect output type: {})".format(type(split))
    assert isinstance(
        unsplit,
        Polygon,
    ), "Polygon conversion failed (incorrect output type: {})".format(
        type(unsplit)
    )

    assert (
        len(split.geoms) == 2
    ), "Polygon splitting failed (incorrect number of outputs: {})".format(
        len(split.geoms)
    )

    assert (
        np.allclose(split.area, unsplit.area)
    ), "Polygon splitting area mismatch"


@pytest.mark.parametrize("origin", origins)
def test_polygon_areas(origin, radii=RADII, n_points=N_POINTS):
    areas = [_get_circle_areas(r, origin=origin, n=n_points) for r in radii]
    geometric_areas, pygplates_areas = zip(*areas)
    assert np.allclose(
        geometric_areas,
        pygplates_areas,
    ), "Polygon areas not equal"


@pytest.mark.parametrize("origin", origins)
def test_polygon_perimeters(origin, radii=RADII, n_points=N_POINTS):
    perimeters = [_get_circle_perimeters(r, origin=origin, n=n_points) for r in radii]
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
