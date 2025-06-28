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

"""This sub-module contains tools for converting PyGPlates or GPlately geometries to Shapely geometries for mapping (and vice versa).

Supported PyGPlates geometries inherit from the following classes:

* [pygplates.GeometryOnSphere](https://www.gplates.org/docs/pygplates/generated/pygplates.geometryonsphere): This
class has the following derived GeometryOnSphere classes:
    * [pygplates.PointOnSphere](https://www.gplates.org/docs/pygplates/generated/pygplates.pointonsphere#pygplates.PointOnSphere)
    * [pygplates.MultiPointOnSphere](https://www.gplates.org/docs/pygplates/generated/pygplates.multipointonsphere#pygplates.MultiPointOnSphere)
    * [pygplates.PolylineOnSphere](https://www.gplates.org/docs/pygplates/generated/pygplates.polylineonsphere#pygplates.PolylineOnSphere)
    * [pygplates.PolygonOnSphere](https://www.gplates.org/docs/pygplates/generated/pygplates.polygononsphere#pygplates.PolygonOnSphere)


* [pygplates.LatLonPoint](https://www.gplates.org/docs/pygplates/generated/pygplates.latlonpoint)
* [pygplates.ReconstructedFeatureGeometry](https://www.gplates.org/docs/pygplates/generated/pygplates.reconstructedfeaturegeometry)
* [pygplates.ResolvedTopologicalLine](https://www.gplates.org/docs/pygplates/generated/pygplates.resolvedtopologicalline)
* [pygplates.ResolvedTopologicalBoundary](https://www.gplates.org/docs/pygplates/generated/pygplates.resolvedtopologicalboundary)
* [pygplates.ResolvedTopologicalNetwork](https://www.gplates.org/docs/pygplates/generated/pygplates.resolvedtopologicalnetwork)

Note: GPlately geometries derive from the `GeometryOnSphere` and `pygplates.GeometryOnSphere` base classes.

Supported Shapely geometric objects include:

* __Point__: a single point in 2D space with coordinate tuple (x,y) or 3D space with coordinate tuple (x,y,z).
* __LineString__: a sequence of points joined together to form a line (a list of point coordinate tuples).
* __Polygon__: a sequence of points joined together to form the outer ring of a filled area, or a hole (a list of at least
three point coordinate tuples).

Also supported are collections of geometric objects, such as:

* __MultiPoint__: a list of __Point__ objects (a list of point coordinate tuples).
* __MultiLineString__: a list of __LineString__ objects (a list containing lists of point coordinate tuples).
* __MultiPolygon__:  a list of __Polygon__ objects (a list containing lists of point coordinate tuples that define exterior rings and/or holes).


__Converting PyGPlates geometries into Shapely geometries__ involves:

* __wrapping geometries at the dateline__: this involves splitting a polygon, MultiPolygon, line segment or MultiLine segment between
connecting points at the dateline. This is to ensure the geometry's points are joined along the short path rather than
the long path horizontally across the 2D map projection display.
* __ordering geometries counter-clockwise__

Input PyGPlates geometries are converted to the following Shapely geometries:

- `PointOnSphere` or `LatLonPoint`: `Point`
- `MultiPointOnSphere`: `MultiPoint`
- `PolylineOnSphere`: `LineString` or `MultiLineString`
- `PolygonOnSphere`: `Polygon` or `MultiPolygon`

__Converting Shapely geometries into PyGPlates geometries__:
Input Shapely geometries are converted to the following PyGPlates geometries:

- `Point`: `PointOnSphere`
- `MultiPoint`: `MultiPointOnSphere`
- `LineString`: `PolylineOnSphere`
- `LinearRing` or `Polygon`: `PolygonOnSphere`

"""

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

__all__ = [
    "GeometryOnSphere",
    "LatLonPoint",
    "MultiPointOnSphere",
    "PointOnSphere",
    "PolygonOnSphere",
    "PolylineOnSphere",
    "pygplates_to_shapely",
    "shapely_to_pygplates",
    "wrap_geometries",
]


class GeometryOnSphere(pygplates.GeometryOnSphere):
    """Class to mix in `to_shapely` method to all GPlately geometry classes.

    All GPlately geometry classes inherit from this class, in addition
    to their PyGPlates base class.
    """

    def to_shapely(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
        validate=False,
        force_ccw=False,
        explode=False,
    ):
        """Convert to Shapely geometry.

        See Also
        --------
        pygplates_to_shapely : Equivalent function.
        """
        return pygplates_to_shapely(
            self,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
            validate=validate,
            force_ccw=force_ccw,
            explode=explode,
        )

    @classmethod
    def from_shapely(cls, geom):
        converted = shapely_to_pygplates(geom)
        return cls(converted)


class PointOnSphere(pygplates.PointOnSphere, GeometryOnSphere):
    """GPlately equivalent of `pygplates.PointOnSphere`, incorporating
    `to_shapely` method
    """

    pass


class MultiPointOnSphere(pygplates.MultiPointOnSphere, GeometryOnSphere):
    """GPlately equivalent of `pygplates.MultiPointOnSphere`, incorporating
    `to_shapely` method
    """

    pass


class PolylineOnSphere(pygplates.PolylineOnSphere, GeometryOnSphere):
    """GPlately equivalent of `pygplates.PolylineOnSphere`, incorporating
    `to_shapely` method
    """

    pass


class PolygonOnSphere(pygplates.PolygonOnSphere, GeometryOnSphere):
    """GPlately equivalent of `pygplates.PolygonOnSphere`, incorporating
    `to_shapely` method
    """

    pass


class LatLonPoint(pygplates.LatLonPoint):
    """GPlately equivalent of `pygplates.LatLonPoint`, incorporating
    `to_shapely` method
    """

    def to_shapely(self, central_meridian=0.0, tessellate_degrees=None):
        return pygplates_to_shapely(
            self,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )


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
    return output_type(output_geoms)


def _ensure_ccw(geometry):
    if (
        isinstance(geometry, _Polygon)
        and geometry.exterior is not None
        and not geometry.exterior.is_ccw
    ):
        return _Polygon(list(geometry.exterior.coords)[::-1])
    return geometry


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
        _Point: PointOnSphere,
        _MultiPoint: MultiPointOnSphere,
        _LineString: PolylineOnSphere,
        _LinearRing: PolygonOnSphere,
        _Polygon: PolygonOnSphere,
    }

    if isinstance(geometry, _BaseMultipartGeometry) and not isinstance(
        geometry, _MultiPoint
    ):
        return [shapely_to_pygplates(i) for i in geometry.geoms]
    if _contains_shapely_geometries(geometry):
        # Recursively convert all elements in iterable of geometries
        out = []
        for i in geometry:
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
                out.extend(tmp)
        return out


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
        return _MultiPoint(out)
    if isinstance(geometry, (_LineString, _MultiLineString)):
        return _MultiLineString(out)
    if isinstance(geometry, (_LinearRing, _Polygon, _MultiPolygon)):
        return _MultiPolygon(out)


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


__pdoc__ = {
    "PointOnSphere": PointOnSphere.__doc__,
    "PolygonOnSphere": PolygonOnSphere.__doc__,
    "LatLonPoint": LatLonPoint.__doc__,
    "MultiPointOnSphere": MultiPointOnSphere.__doc__,
    "PolylineOnSphere": PolylineOnSphere.__doc__,
}
