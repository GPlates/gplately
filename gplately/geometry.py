#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
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

"""This module contains tools for handling geometries,
such as converting PyGPlates geometries to Shapely geometries (and vice versa)."""

import numpy as np
import pygplates
from shapely.geometry import LinearRing as _LinearRing
from shapely.geometry import LineString as _LineString
from shapely.geometry import MultiLineString as _MultiLineString
from shapely.geometry import MultiPoint as _MultiPoint
from shapely.geometry import MultiPolygon as _MultiPolygon
from shapely.geometry import Point as _Point
from shapely.geometry import Polygon as _Polygon
from shapely.geometry.base import BaseGeometry as _BaseGeometry
from shapely.geometry.base import BaseMultipartGeometry as _BaseMultipartGeometry

from pygplates import (
    GeometryOnSphere as _GeometryOnSphere,
    PointOnSphere as _PointOnSphere,
    MultiPointOnSphere as _MultiPointOnSphere,
    PolylineOnSphere as _PolylineOnSphere,
    PolygonOnSphere as _PolygonOnSphere,
)

__all__ = [
    "pygplates_to_shapely",
    "shapely_to_pygplates",
    "wrap_geometries",
    "to_shapely",
    "from_shapely",
]


def from_shapely(geometry: _BaseGeometry) -> _GeometryOnSphere:
    """Convert one Shapely geometry to a PyGPlates geometry.

    Parameters
    ----------
    geometry : shapely.geometry.base.BaseGeometry
        A `Point`, `MultiPoint`, `LineString`, or `Polygon`.

    Returns
    -------
    pygplates.GeometryOnSphere
        A `PointOnSphere`, `MultiPointOnSphere`, `PolylineOnSphere`, or `PolygonOnSphere`, respectively.
    """
    if isinstance(geometry, _Point):
        return _PointOnSphere(geometry.y, geometry.x)

    if isinstance(geometry, _MultiPoint):
        coords = np.array([point.coords[0] for point in geometry.geoms])
        return _MultiPointOnSphere(np.fliplr(coords))

    if isinstance(geometry, _LineString):
        return _PolylineOnSphere(np.fliplr(np.array(geometry.coords)))

    if isinstance(geometry, _Polygon):
        if geometry.exterior is None:
            raise AttributeError("Polygon geometry has no exterior")
        exterior = np.fliplr(np.array(geometry.exterior.coords)[:-1])
        interiors = [
            np.fliplr(np.array(interior.coords)[:-1]) for interior in geometry.interiors
        ]
        return _PolygonOnSphere(exterior, interiors)

    raise TypeError("Invalid geometry type: " + str(type(geometry)))


def to_shapely(geometry: _GeometryOnSphere) -> _BaseGeometry:
    """Convert one PyGPlates geometry to a Shapely geometry.

    A direct, unwrapped conversion: coordinates are used as-is, with no
    antimeridian wrapping, tessellation, or validation. For those, use
    `pygplates_to_shapely` instead.

    Parameters
    ----------
    geometry : pygplates.GeometryOnSphere
        A `PointOnSphere`, `MultiPointOnSphere`, `PolylineOnSphere`, or `PolygonOnSphere`.

    Returns
    -------
    shapely.geometry.base.BaseGeometry
        A `Point`, `MultiPoint`, `LineString`, or `Polygon`, respectively.
    """
    if isinstance(geometry, _PointOnSphere):
        lat, lon = geometry.to_lat_lon()
        return _Point(lon, lat)

    if isinstance(geometry, _MultiPointOnSphere):
        return _MultiPoint(np.fliplr(geometry.to_lat_lon_array()))

    if isinstance(geometry, _PolylineOnSphere):
        return _LineString(np.fliplr(geometry.to_lat_lon_array()))

    if isinstance(geometry, _PolygonOnSphere):
        exterior = np.fliplr(
            np.array(
                [point.to_lat_lon() for point in geometry.get_exterior_ring_points()]
            )
        )
        interiors = [
            np.fliplr(
                np.array(
                    [
                        point.to_lat_lon()
                        for point in geometry.get_interior_ring_points(i)
                    ]
                )
            )
            for i in range(geometry.get_number_of_interior_rings())
        ]
        return _Polygon(exterior, interiors)

    raise TypeError("Invalid geometry type: " + str(type(geometry)))


def pygplates_to_shapely(
    geometry,
    central_meridian=0.0,
    tessellate_degrees=None,
    validate=False,
    force_ccw=False,
    explode=False,
):
    """Convert one or more PyGPlates or GPlately geometries to Shapely format.

    Parameters
    ----------
    geometry : pygplates.GeometryOnSphere or pygplates.LatLonPoint or list
        The geometry or geometries to convert.
    central_meridian : float, default: 0.0
        The central meridian around which to wrap geometries;
        geometries will be split at the antimeridian.
    tessellate_degrees : float, optional
        If provided, the geometry will be tessellated to this resolution prior to conversion.
    validate : bool, default: False
        Attempt to ensure output geometry is valid by applying a buffer of 0.
    force_ccw : bool, default: False
        Ensure the coordinates of the output geometry are counter-clockwise(only applies to polygons).
    explode : bool, default: False
        Convert multi-part output geometries to multiple single-part geometries.

    Returns
    -------
    output_geometry : shapely.geometry.base.BaseGeometry or list
        Converted Shapely geometry or geometries.

    Notes
    -----
    If a single input geometry was passed, `output_geometry` will be a
    subclass of `shapely.geometry.base.BaseGeometry`. Otherwise,
    `output_geometry` will be a list of the same length as the input.

    Input geometries that were split while wrapping around
    `central_meridian` will produce multi-part output geometries, unless `explode=True` is specified.

    Input geometry types are converted as follows:
        - `PointOnSphere` or `LatLonPoint`:
            `Point`
        - `MultiPointOnSphere`:
            `MultiPoint`
        - `PolylineOnSphere`:
            `LineString` or
            `MultiLineString`
        - `PolygonOnSphere`:
            `Polygon` or
            `MultiPolygon`
    """
    if _contains_pygplates_geometries(geometry):
        return [
            pygplates_to_shapely(
                i,
                central_meridian=central_meridian,
                tessellate_degrees=tessellate_degrees,
                validate=validate,
                force_ccw=force_ccw,
            )
            for i in geometry
        ]

    if isinstance(geometry, pygplates.LatLonPoint):
        geometry = geometry.to_point_on_sphere()
    if isinstance(geometry, pygplates.ReconstructedFeatureGeometry):
        geometry = geometry.get_reconstructed_geometry()
    if isinstance(
        geometry,
        (
            pygplates.ResolvedTopologicalLine,
            pygplates.ResolvedTopologicalBoundary,
            pygplates.ResolvedTopologicalNetwork,
        ),
    ):
        geometry = geometry.get_resolved_geometry()
    if not isinstance(geometry, pygplates.GeometryOnSphere):
        raise TypeError("Invalid geometry type: " + str(type(geometry)))

    wrapper = pygplates.DateLineWrapper(central_meridian=central_meridian)
    wrapped = wrapper.wrap(geometry, tessellate_degrees=tessellate_degrees)

    if isinstance(wrapped, pygplates.LatLonPoint):
        return _Point(wrapped.to_lat_lon()[::-1])
    if isinstance(wrapped, pygplates.DateLineWrapper.LatLonMultiPoint):
        points = wrapped.get_points()
        return _MultiPoint([i.to_lat_lon()[::-1] for i in points])

    # pygplates.DateLineWrapper wraps correctly to the dateline but Cartopy can sometimes
    # move a point at central_meridian+180 to central_meridian-180 (or vice versa).
    # To correct this we move points on the dateline slightly inwards (inside the map projection).
    #
    # Empirically this can go down to about 1e-11 (fails if 1e-12).
    # But we leave enough of headroom for errors to accumulate.
    dateline_clip_threshold = 1e-8

    output_geoms = []
    output_type = None
    for i in wrapped:
        if isinstance(i, pygplates.DateLineWrapper.LatLonPolyline):
            tmp = np.array([j.to_lat_lon()[::-1] for j in i.get_points()])
            # Clip near the dateline.
            tmp[:, 0] = np.clip(
                tmp[:, 0],
                central_meridian - (180 - dateline_clip_threshold),
                central_meridian + (180 - dateline_clip_threshold),
            )
            tmp = _LineString(tmp)
            output_geoms.append(tmp)
            output_type = _MultiLineString
        elif isinstance(i, pygplates.DateLineWrapper.LatLonPolygon):
            tmp = np.array([j.to_lat_lon()[::-1] for j in i.get_exterior_points()])
            # Clip near the dateline.
            tmp[:, 0] = np.clip(
                tmp[:, 0],
                central_meridian - (180 - dateline_clip_threshold),
                central_meridian + (180 - dateline_clip_threshold),
            )
            # tmp[:,1] = np.clip(tmp[:,1], -89, 89) # clip polygons near poles
            tmp = _Polygon(tmp)
            if force_ccw and tmp.exterior is not None and not tmp.exterior.is_ccw:
                tmp = _Polygon(list(tmp.exterior.coords)[::-1])
                # tmp.exterior.coords = list(tmp.exterior.coords)[::-1]
            if validate:
                tmp = tmp.buffer(0.0)
            # this is for pole-clipped polygons turned into MultiPolygons
            if isinstance(tmp, _MultiPolygon):
                # for geom in list(tmp):
                for geom in tmp.geoms:
                    output_geoms.append(geom)
            else:
                output_geoms.append(tmp)
            output_type = _MultiPolygon
        else:
            raise TypeError(
                "Unrecognised output from `pygplates.DateLineWrapper.wrap`: "
                + str(type(i))
            )
    if output_type is None:
        raise TypeError(
            "Unrecognised output from `pygplates.DateLineWrapper.wrap`: "
            + str(type(wrapped[0]))
        )
    # Empty geometries can sometimes occur by this point, causing nearly all
    # subsequent geometric operations to fail
    output_geoms = [i for i in output_geoms if not i.is_empty]
    if force_ccw:
        output_geoms = [_ensure_ccw(i) for i in output_geoms]
    if len(output_geoms) == 1:
        return output_geoms[0]
    if explode:
        return output_geoms
    return output_type(output_geoms)  # pyright: ignore[reportArgumentType]


def shapely_to_pygplates(geometry):
    """Convert one or more Shapely geometries to gplately format.

    Parameters
    ----------
    geometry : shapely.geometry.base.BaseGeometry or list
        The geometry or geometries to convert.

    Returns
    -------
    output_geometry : GeometryOnSphere or list
        Converted gplately geometry or geometries.

    Notes
    -----
    If a single input geometry was passed, `output_geometry` will be a subclass of `GeometryOnSphere`.
    Otherwise, `output_geometry` will be a list of `GeometryOnSphere`, of the same length as the input.

    Input geometry types are converted as follows:
        - `Point`: `PointOnSphere`
        - `MultiPoint`: `MultiPointOnSphere`
        - `LineString`: `PolylineOnSphere`
        - `LinearRing` or `Polygon`: `PolygonOnSphere`

    Multi-part input geometry types other than `MultiPoint` will be treated
    as an iterable of their component single-part geometries.
    """
    pygplates_conversion = {
        _Point: _PointOnSphere,
        _MultiPoint: _MultiPointOnSphere,
        _LineString: _PolylineOnSphere,
        _LinearRing: _PolygonOnSphere,
        _Polygon: _PolygonOnSphere,
    }

    if isinstance(geometry, _BaseMultipartGeometry) and not isinstance(
        geometry, _MultiPoint
    ):
        return [shapely_to_pygplates(i) for i in geometry.geoms]
    if _contains_shapely_geometries(geometry):
        # Recursively convert all elements in iterable of geometries
        out = []
        for i in geometry:  # pyright: ignore[reportGeneralTypeIssues]
            tmp = shapely_to_pygplates(i)
            if isinstance(tmp, pygplates.GeometryOnSphere):
                # Output is a single geometry
                out.append(tmp)
            else:
                # Output should be a list of geometries
                out.extend(tmp)
        return out

    for input_type in pygplates_conversion:
        if isinstance(geometry, input_type):
            output_type = pygplates_conversion[input_type]
            break
    else:
        raise TypeError("Invalid geometry type: " + str(type(geometry)))
    if isinstance(geometry, _MultiPoint):
        coords = np.array([i.coords for i in geometry.geoms]).squeeze()
    elif isinstance(geometry, _Polygon):
        if geometry.exterior is None:
            raise AttributeError("Polygon geometry has no exterior")
        coords = np.array(geometry.exterior.coords).squeeze()[:-1, ...]
    elif hasattr(geometry, "coords"):
        coords = np.array(geometry.coords).squeeze()
    else:
        raise TypeError("Invalid geometry type: " + str(type(geometry)))
    if coords.ndim > 1:
        coords = np.fliplr(coords)
    else:
        coords = np.flip(coords)
    return output_type(coords)


def wrap_geometries(
    geometries,
    central_meridian=0.0,
    tessellate_degrees=None,
    validate=False,
    force_ccw=False,
    explode=False,
):
    """Wrap one or more Shapely geometries around a central meridian.

    Wrapped geometries will be split at the antimeridian.

    Parameters
    ----------
    geometry : shapely.geometry.base.BaseGeometry or list
        The geometry or geometries to wrap.
    central_meridian : float, default: 0.0
        The central meridian around which to wrap geometries;
        geometries will be split at the antimeridian.
    tessellate_degrees : float, optional
        If provided, the geometry will be tessellated to this resolution prior to wrapping.
    validate : bool, default: False
        Attempt to ensure output geometry is valid by applying a buffer of 0.
    force_ccw : bool, default: False
        Ensure the coordinates of the output geometry are counter-clockwise(only applies to polygons).
    explode : bool, default: False
        Convert multi-part output geometries to multiple single-part geometries.

    Returns
    -------
    output_geometries : shapely.geometry.base.BaseGeometry or list
        Wrapped Shapely geometry or geometries.

    Notes
    -----
    If a single input geometry was passed, `output_geometry` will be a
    subclass of `shapely.geometry.base.BaseGeometry`.
    Otherwise, `output_geometry` will be a list of the same length as the input, unless `explode=True` is specified.

    Input geometries that were split while wrapping around
    `central_meridian` will produce multi-part output geometries, unless `explode=True` is specified.
    """
    if isinstance(geometries, _BaseGeometry):
        return _wrap_geometry(
            geometry=geometries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
            validate=validate,
            force_ccw=force_ccw,
            explode=explode,
        )
    else:
        out = []
        for i in geometries:
            tmp = _wrap_geometry(
                geometry=i,
                central_meridian=central_meridian,
                tessellate_degrees=tessellate_degrees,
                validate=validate,
                force_ccw=force_ccw,
                explode=explode,
            )
            if isinstance(tmp, _BaseGeometry):
                out.append(tmp)
            else:
                out.extend(tmp)  # pyright: ignore[reportArgumentType]
        return out


def geographic_circle(lon, lat, radius_deg, n_points=360):
    """
    Return (lons, lats) tracing a great-circle buffer.

    A great-circle buffer around a point is the set of all points that are exactly
    the same great-circle distance away from it — which traces out a circle on the sphere's surface.

    Parameters
    ----------
    lon, lat    : centre coordinates in degrees
    radius_deg  : radius in great-circle degrees (1° ≈ 111 km)
    n_points    : number of vertices

    Returns
    -------
    lons, lats : arrays of longitudes and latitudes of the circle vertices, in degrees
    """
    azimuths = np.linspace(0, 360, n_points, endpoint=False)
    lons, lats = [], []
    lat_r = np.radians(lat)
    d_r = np.radians(radius_deg)

    for az in azimuths:
        az_r = np.radians(az)
        out_lat = np.degrees(
            np.arcsin(
                np.sin(lat_r) * np.cos(d_r) + np.cos(lat_r) * np.sin(d_r) * np.cos(az_r)
            )
        )
        out_lon = lon + np.degrees(
            np.arctan2(
                np.sin(az_r) * np.sin(d_r) * np.cos(lat_r),
                np.cos(d_r) - np.sin(lat_r) * np.sin(np.radians(out_lat)),
            )
        )
        lons.append(out_lon)
        lats.append(out_lat)

    lons.append(lons[0])
    lats.append(lats[0])
    return np.array(lons), np.array(lats)


def _wrap_geometry(
    geometry,
    central_meridian=0.0,
    tessellate_degrees=None,
    validate=False,
    force_ccw=False,
    explode=False,
):
    if not isinstance(geometry, _BaseGeometry):
        raise TypeError("Invalid geometry type: " + str(type(geometry)))

    converted = shapely_to_pygplates(geometry)
    if isinstance(converted, pygplates.GeometryOnSphere):
        out = [converted]
    else:
        out = []
        for i in converted:
            if isinstance(i, pygplates.GeometryOnSphere):
                out.append(i)
            else:
                out.extend(i)
    if explode:
        out_tmp = []
        for i in out:
            tmp = pygplates_to_shapely(
                geometry=i,
                central_meridian=central_meridian,
                tessellate_degrees=tessellate_degrees,
                validate=validate,
                force_ccw=force_ccw,
                explode=explode,
            )
            if isinstance(tmp, _BaseGeometry):
                out_tmp.append(tmp)
            else:
                out_tmp.extend(tmp)
        out = out_tmp
    else:
        out = [
            pygplates_to_shapely(
                geometry=i,
                central_meridian=central_meridian,
                tessellate_degrees=tessellate_degrees,
                validate=validate,
                force_ccw=force_ccw,
                explode=explode,
            )
            for i in out
        ]
    if len(out) == 1:
        return out[0]
    if explode:
        return out
    if isinstance(geometry, (_Point, _MultiPoint)):
        return _MultiPoint(out)  # pyright: ignore[reportArgumentType]
    if isinstance(geometry, (_LineString, _MultiLineString)):
        return _MultiLineString(out)  # pyright: ignore[reportArgumentType]
    if isinstance(geometry, (_LinearRing, _Polygon, _MultiPolygon)):
        return _MultiPolygon(out)  # pyright: ignore[reportArgumentType]


def _contains_shapely_geometries(i):
    """Check if input is an iterable containing only Shapely geometries."""
    if isinstance(i, _BaseGeometry):
        return False
    try:
        # Check all elements in i are Shapely geometries
        for j in i:
            if not isinstance(j, _BaseGeometry):
                break
        else:
            return True
    except TypeError:  # i is not iterable
        pass
    return False


def _is_pygplates_geometry(geom):
    return isinstance(
        geom,
        (
            pygplates.GeometryOnSphere,
            pygplates.LatLonPoint,
            pygplates.ReconstructedFeatureGeometry,
            pygplates.ResolvedTopologicalLine,
            pygplates.ResolvedTopologicalBoundary,
            pygplates.ResolvedTopologicalNetwork,
        ),
    )


def _contains_pygplates_geometries(i):
    """Check if input is an iterable containing only PyGPlates geometries."""
    if _is_pygplates_geometry(i):
        return False
    try:
        # Check all elements in i are PyGPlates geometries
        for j in i:
            if not _is_pygplates_geometry(j):
                break
        else:
            return True
    except TypeError:  # i is not iterable
        pass
    return False


def _ensure_ccw(geometry):
    if (
        isinstance(geometry, _Polygon)
        and geometry.exterior is not None
        and not geometry.exterior.is_ccw
    ):
        return _Polygon(list(geometry.exterior.coords)[::-1])
    return geometry
