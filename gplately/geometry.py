import numpy as np
import pygplates
from shapely.geometry import (
    LinearRing as _LinearRing,
    LineString as _LineString,
    MultiLineString as _MultiLineString,
    MultiPoint as _MultiPoint,
    MultiPolygon as _MultiPolygon,
    Point as _Point,
    Polygon as _Polygon,
)
from shapely.geometry.base import (
    BaseGeometry as _BaseGeometry,
    BaseMultipartGeometry as _BaseMultipartGeometry,
)

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
    area_threshold=1e-6
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
        If provided, the geometry will be tessellated to this
        resolution prior to conversion.
    validate : bool, default: False
        Attempt to ensure output geometry is valid by applying a buffer of 0.
    force_ccw : bool, default: False
        Ensure the coordinates of the output geometry are counter-clockwise
        (only applies to polygons).
    explode : bool, default: False
        Convert multi-part output geometries to multiple single-part
        geometries.

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
    `central_meridian` will produce multi-part output geometries, unless
    `explode=True` is specified.

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
    if not isinstance(geometry, pygplates.GeometryOnSphere):
        raise TypeError("Invalid geometry type: " + str(type(geometry)))

    wrapper = pygplates.DateLineWrapper(central_meridian=central_meridian)
    wrapped = wrapper.wrap(geometry, tessellate_degrees=tessellate_degrees)

    if isinstance(wrapped, pygplates.LatLonPoint):
        return _Point(wrapped.to_lat_lon()[::-1])
    if isinstance(wrapped, pygplates.DateLineWrapper.LatLonMultiPoint):
        points = wrapped.get_points()
        return _MultiPoint([i.to_lat_lon()[::-1] for i in points])

    output_geoms = []
    output_type = None
    for i in wrapped:
        if isinstance(i, pygplates.DateLineWrapper.LatLonPolyline):
            tmp = _LineString([j.to_lat_lon()[::-1] for j in i.get_points()])
            output_geoms.append(tmp)
            output_type = _MultiLineString
        elif isinstance(i, pygplates.DateLineWrapper.LatLonPolygon):
            tmp = np.array([j.to_lat_lon()[::-1] for j in i.get_exterior_points()])
            # tmp[:,1] = np.clip(tmp[:,1], -89, 89) # clip polygons near poles
            tmp = _Polygon(tmp)
            if (
                force_ccw
                and tmp.exterior is not None
                and not tmp.exterior.is_ccw
            ):
                tmp = _Polygon(list(tmp.exterior.coords)[::-1])
                # tmp.exterior.coords = list(tmp.exterior.coords)[::-1]
            if validate:
                tmp = tmp.buffer(0.0)
            # this is for pole-clipped polygons turned into MultiPolygons
            if isinstance(tmp, _MultiPolygon):
                for geom in list(tmp):
                    if geom.area > area_threshold:
                        output_geoms.append(geom)
            else:
                if tmp.area > area_threshold:
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
    """Convert one or more Shapely geometries to PyGPlates format.

    Parameters
    ----------
    geometry : shapely.geometry.base.BaseGeometry or list
        The geometry or geometries to convert.

    Returns
    -------
    output_geometry : pygplates.GeometryOnSphere or list
        Converted PyGPlates geometry or geometries.

    Notes
    -----
    If a single input geometry was passed, `output_geometry` will be a
    subclass of `pygplates.GeometryOnSphere`. Otherwise, `output_geometry`
    will be a list of `pygplates.GeometryOnSphere`, of the same length as
    the input.

    Input geometry types are converted as follows:
        - `Point`: `PointOnSphere`
        - `MultiPoint`: `MultiPointOnSphere`
        - `LineString`: `PolylineOnSphere`
        - `LinearRing` or `Polygon`:
            `PolygonOnSphere`

    Multi-part input geometry types other than `MultiPoint` will be treated
    as an iterable of their component single-part geometries.
    """
    pygplates_conversion = {
        _Point: pygplates.PointOnSphere,
        _MultiPoint: pygplates.MultiPointOnSphere,
        _LineString: pygplates.PolylineOnSphere,
        _LinearRing: pygplates.PolygonOnSphere,
        _Polygon: pygplates.PolygonOnSphere,
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
        If provided, the geometry will be tessellated to this
        resolution prior to wrapping.
    validate : bool, default: False
        Attempt to ensure output geometry is valid by applying a buffer of 0.
    force_ccw : bool, default: False
        Ensure the coordinates of the output geometry are counter-clockwise
        (only applies to polygons).
    explode : bool, default: False
        Convert multi-part output geometries to multiple single-part
        geometries.

    Returns
    -------
    output_geometries : shapely.geometry.base.BaseGeometry or list
        Wrapped Shapely geometry or geometries.

    Notes
    -----
    If a single input geometry was passed, `output_geometry` will be a
    subclass of `shapely.geometry.base.BaseGeometry`. Otherwise,
    `output_geometry` will be a list of the same length as the input,
    unless `explode=True` is specified.

    Input geometries that were split while wrapping around
    `central_meridian` will produce multi-part output geometries, unless
    `explode=True` is specified.
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


def _contains_pygplates_geometries(i):
    """Check if input is an iterable containing only PyGPlates geometries."""
    is_pygplates_geometry = lambda x: isinstance(
        x, (pygplates.GeometryOnSphere, pygplates.LatLonPoint)
    )
    if is_pygplates_geometry(i):
        return False
    try:
        # Check all elements in i are PyGPlates geometries
        for j in i:
            if not is_pygplates_geometry(j):
                break
        else:
            return True
    except TypeError:  # i is not iterable
        pass
    return False
