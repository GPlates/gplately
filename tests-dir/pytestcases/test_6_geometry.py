from itertools import product
from typing import cast

import numpy as np
import pygplates
import pytest
from conftest import logger
from conftest import test_geometry_azimuths as azimuths
from conftest import test_geometry_n_points as N_POINTS
from conftest import test_geometry_origins as origins
from conftest import test_geometry_radii as RADII
from shapely.geometry import LineString, MultiPoint, Point, Polygon
from shapely.geometry.base import BaseMultipartGeometry

from gplately import EARTH_RADIUS
from gplately.geometry import *

logger.info(__name__)


@pytest.mark.parametrize("lat", (-100, 95))
def test_invalid_point(lat):
    p = Point(0, lat)
    with pytest.raises(pygplates.InvalidLatLonError):
        shapely_to_pygplates(p)


def test_multipoint_conversion(n_points=N_POINTS):
    n_lats = int(np.sqrt(n_points / 2))
    point_lons = np.linspace(-179, 179, n_lats * 2)
    point_lats = np.linspace(-89, 89, n_lats)
    grid_lons, grid_lats = np.meshgrid(point_lons, point_lats)

    mp_shapely = MultiPoint(
        np.column_stack([np.ravel(i) for i in (grid_lons, grid_lats)])
    )
    mp_pgp = cast(pygplates.MultiPointOnSphere, shapely_to_pygplates(mp_shapely))

    a1 = np.vstack([i.coords for i in mp_shapely.geoms])
    a2 = np.fliplr(mp_pgp.to_lat_lon_array())
    assert np.allclose(
        a1,
        a2,
    ), "Multi-point coordinates not equal"


@pytest.mark.parametrize(
    "origin,azimuth,distance",
    list(product(origins, azimuths, RADII)),
)
def test_polyline_conversion(origin, azimuth, distance, n_points=N_POINTS):
    lons, lats = _point_from_origin_distance_azimuth(
        origin,
        np.linspace(0, distance, n_points),
        azimuth,
    )
    line = LineString(zip(lons, lats))
    converted = shapely_to_pygplates(line)
    assert isinstance(
        converted, pygplates.PolylineOnSphere
    ), f"Expected a PolylineOnSphere, got {type(converted)}"
    converted = cast(pygplates.PolylineOnSphere, converted)

    reconverted = pygplates_to_shapely(
        converted,
        central_meridian=converted.get_centroid().to_lat_lon()[1],
    )
    assert isinstance(
        reconverted, LineString
    ), f"Expected a LineString, got {type(reconverted)}"

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
    list(product(origins, azimuths, RADII)),
)
def test_polyline_length(origin, azimuth, distance, n_points=N_POINTS):
    lons, lats = _point_from_origin_distance_azimuth(
        origin,
        np.linspace(0, distance, n_points),
        azimuth,
    )
    line = LineString(zip(lons, lats))
    converted = cast(pygplates.PolylineOnSphere, shapely_to_pygplates(line))
    assert np.allclose(
        distance,
        converted.get_arc_length() * EARTH_RADIUS,
    ), "Polyline length not consistent"


@pytest.mark.parametrize("origin,radius", list(product(origins, RADII)))
def test_polygon_conversion(origin, radius, n_points=N_POINTS):
    lons, lats = _get_circle(origin, radius, n=n_points)
    circle = Polygon(zip(lons, lats))
    converted = cast(pygplates.PolygonOnSphere, shapely_to_pygplates(circle))
    reconverted = pygplates_to_shapely(
        converted,
        central_meridian=converted.get_interior_centroid().to_lat_lon()[1],
    )
    assert isinstance(
        reconverted, Polygon
    ), f"Expected a Polygon, got {type(reconverted)}"
    assert np.allclose(
        circle.exterior.coords,
        reconverted.exterior.coords,
    ), "Polygon coordinates not equal"


def test_polygon_splitting(n_points=N_POINTS):
    origin = (-179, 0)
    radius = 1000
    lons, lats = _get_circle(origin, radius, n=n_points)
    circle = Polygon(zip(lons, lats))
    converted = cast(pygplates.PolygonOnSphere, shapely_to_pygplates(circle))
    split = pygplates_to_shapely(converted, central_meridian=0)
    unsplit = pygplates_to_shapely(
        converted,
        central_meridian=converted.get_interior_centroid().to_lat_lon()[1],
    )

    assert isinstance(
        split,
        BaseMultipartGeometry,
    ), "Polygon splitting failed (incorrect output type: {})".format(type(split))
    assert isinstance(
        unsplit,
        Polygon,
    ), "Polygon conversion failed (incorrect output type: {})".format(type(unsplit))

    assert (
        len(split.geoms) == 2
    ), "Polygon splitting failed (incorrect number of outputs: {})".format(
        len(split.geoms)
    )

    assert np.allclose(split.area, unsplit.area), "Polygon splitting area mismatch"


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


# ---------------------------------------------------------------------------
# `to_shapely` / `from_shapely`: direct, unwrapped conversions (no
# antimeridian splitting, tessellation, or validation), unlike
# `pygplates_to_shapely` / `shapely_to_pygplates` above.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("lon,lat", [(0, 0), (179.5, -89.9), (-155.4696, 19.8202)])
def test_to_shapely_point(lon, lat):
    result = to_shapely(pygplates.PointOnSphere(lat, lon))
    assert isinstance(result, Point)
    assert np.allclose(result.coords[0], (lon, lat))


@pytest.mark.parametrize("lon,lat", [(0, 0), (179.5, -89.9), (-155.4696, 19.8202)])
def test_from_shapely_point(lon, lat):
    result = from_shapely(Point(lon, lat))
    assert isinstance(result, pygplates.PointOnSphere)
    assert np.allclose(result.to_lat_lon(), (lat, lon))


def test_to_shapely_multipoint(n_points=N_POINTS):
    lons = np.linspace(-170, 170, n_points)
    lats = np.linspace(-80, 80, n_points)
    result = to_shapely(pygplates.MultiPointOnSphere(list(zip(lats, lons))))
    assert isinstance(result, MultiPoint)
    assert np.allclose(
        np.array([point.coords[0] for point in result.geoms]),
        np.column_stack((lons, lats)),
    )


def test_from_shapely_multipoint(n_points=N_POINTS):
    lons = np.linspace(-170, 170, n_points)
    lats = np.linspace(-80, 80, n_points)
    result = from_shapely(MultiPoint(np.column_stack((lons, lats))))
    assert isinstance(result, pygplates.MultiPointOnSphere)
    assert np.allclose(
        np.fliplr(result.to_lat_lon_array()),
        np.column_stack((lons, lats)),
    )


@pytest.mark.parametrize(
    "origin,azimuth,distance",
    list(product(origins, azimuths, RADII)),
)
def test_to_shapely_polyline(origin, azimuth, distance, n_points=N_POINTS):
    lons, lats = _point_from_origin_distance_azimuth(
        origin,
        np.linspace(0, distance, n_points),
        azimuth,
    )
    result = to_shapely(pygplates.PolylineOnSphere(list(zip(lats, lons))))
    assert isinstance(result, LineString)
    assert np.allclose(
        np.array(result.coords), np.column_stack((_normalize_lon(lons), lats))
    )


@pytest.mark.parametrize(
    "origin,azimuth,distance",
    list(product(origins, azimuths, RADII)),
)
def test_from_shapely_polyline(origin, azimuth, distance, n_points=N_POINTS):
    lons, lats = _point_from_origin_distance_azimuth(
        origin,
        np.linspace(0, distance, n_points),
        azimuth,
    )
    line = LineString(np.column_stack((lons, lats)))
    result = from_shapely(line)
    assert isinstance(result, pygplates.PolylineOnSphere)
    assert np.allclose(
        np.fliplr(result.to_lat_lon_array()),
        np.column_stack((_normalize_lon(lons), lats)),
    )


@pytest.mark.parametrize("origin,radius", list(product(origins, RADII)))
def test_to_shapely_polygon(origin, radius, n_points=N_POINTS):
    lons, lats = _get_circle(origin, radius, n=n_points)
    result = to_shapely(pygplates.PolygonOnSphere(list(zip(lats, lons))))
    assert isinstance(result, Polygon)
    assert np.allclose(
        np.array(result.exterior.coords)[:-1],
        np.column_stack((_normalize_lon(lons), lats)),
    )


@pytest.mark.parametrize("origin,radius", list(product(origins, RADII)))
def test_from_shapely_polygon(origin, radius, n_points=N_POINTS):
    lons, lats = _get_circle(origin, radius, n=n_points)
    polygon = Polygon(np.column_stack((lons, lats)))
    result = from_shapely(polygon)
    assert isinstance(result, pygplates.PolygonOnSphere)
    assert np.allclose(
        np.fliplr(
            np.array(
                [point.to_lat_lon() for point in result.get_exterior_ring_points()]
            )
        ),
        np.column_stack((_normalize_lon(lons), lats)),
    )


def test_polygon_with_holes_round_trip():
    exterior = list(zip(*_get_circle((0, 0), 2000)))
    hole = list(zip(*_get_circle((0, 0), 500)))
    polygon = Polygon(exterior, holes=[hole])

    converted = from_shapely(polygon)
    assert converted.get_number_of_interior_rings() == 1

    round_tripped = to_shapely(converted)
    assert isinstance(round_tripped, Polygon)
    assert len(round_tripped.interiors) == 1
    assert np.allclose(polygon.exterior.coords, round_tripped.exterior.coords)
    assert np.allclose(polygon.interiors[0].coords, round_tripped.interiors[0].coords)


@pytest.mark.parametrize("invalid", ["not a geometry", 42, None])
def test_to_shapely_invalid_type(invalid):
    with pytest.raises(TypeError):
        to_shapely(invalid)


@pytest.mark.parametrize("invalid", ["not a geometry", 42, None])
def test_from_shapely_invalid_type(invalid):
    with pytest.raises(TypeError):
        from_shapely(invalid)


def _get_circle_perimeters(radius, origin=(0, 0), n=N_POINTS):
    lons, lats = _get_circle(origin, radius, n=n)
    circle = Polygon(zip(lons, lats))
    converted = cast(pygplates.PolygonOnSphere, shapely_to_pygplates(circle))

    pygplates_perimeter = converted.get_arc_length() * EARTH_RADIUS
    geometric_perimeter = 2 * np.pi * EARTH_RADIUS * np.sin(radius / EARTH_RADIUS)
    return geometric_perimeter, pygplates_perimeter


def _get_circle_areas(radius, origin=(0, 0), n=N_POINTS):
    lons, lats = _get_circle(origin, radius, n=n)
    circle = Polygon(zip(lons, lats))
    converted = cast(pygplates.PolygonOnSphere, shapely_to_pygplates(circle))

    pygplates_area = converted.get_area() * (EARTH_RADIUS**2)
    geometric_area = 2 * np.pi * (EARTH_RADIUS**2) * (1 - np.cos(radius / EARTH_RADIUS))
    return geometric_area, pygplates_area


def _normalize_lon(lon):
    """Wrap longitude values into the canonical (-180, 180] range used
    internally by pygplates, so raw expected values (which may fall outside
    that range, e.g. from `_point_from_origin_distance_azimuth`) can be
    compared against pygplates output."""
    lon = np.asarray(lon, dtype=float)
    return ((lon + 180.0) % 360.0) - 180.0


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
    B = np.arcsin(np.sin(b) * np.sin(azimuth) / np.sin(a))
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
