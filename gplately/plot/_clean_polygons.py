import cartopy.crs as ccrs
import geopandas as gpd
import numpy as np
from shapely.geometry import (
    LineString,
    MultiPolygon,
    Point,
    Polygon,
    box,
)
from shapely.geometry.base import BaseMultipartGeometry
from shapely.ops import linemerge, substring

from ._axes_tools import meridian_from_projection


def clean_polygons(data, projection):
    data = gpd.GeoDataFrame(data)
    data = data.explode(ignore_index=True)

    if data.crs is None:
        data.crs = ccrs.PlateCarree()

    if isinstance(
        projection,
        (
            ccrs._RectangularProjection,
            ccrs._WarpedRectangularProjection,
        ),
    ):
        central_longitude = meridian_from_projection(projection)
        dx = 1.0e-3
        dy = 5.0e-2
        rects = (
            box(
                central_longitude - 180,
                -90,
                central_longitude - 180 + dx,
                90,
            ),
            box(
                central_longitude + 180 - dx,
                -90,
                central_longitude + 180,
                90,
            ),
            box(
                central_longitude - 180,
                -90 - dy * 0.5,
                central_longitude + 180,
                -90 + dy * 0.5,
            ),
            box(
                central_longitude - 180,
                90 - dy * 0.5,
                central_longitude + 180,
                90 + dy * 0.5,
            ),
        )
        rects = gpd.GeoDataFrame(
            {"geometry": rects},
            geometry="geometry",
            crs=ccrs.PlateCarree(),
        )
        data = data.overlay(rects, how="difference")

    projected = data.to_crs(projection)

    # If no [Multi]Polygons, return projected data
    for geom in projected.geometry:
        if isinstance(geom, (Polygon, MultiPolygon)):
            break
    else:
        return projected

    proj_width = np.abs(projection.x_limits[1] - projection.x_limits[0])
    proj_height = np.abs(projection.y_limits[1] - projection.y_limits[0])
    min_distance = np.mean((proj_width, proj_height)) * 1.0e-4

    boundary = projection.boundary
    if np.array(boundary.coords).shape[1] == 3:
        boundary = type(boundary)(np.array(boundary.coords)[:, :2])
    return fill_all_edges(projected, boundary, min_distance=min_distance)


def fill_all_edges(data, boundary, min_distance=None):
    data = gpd.GeoDataFrame(data).explode(ignore_index=True)

    def drop_func(geom):
        if hasattr(geom, "exterior"):
            geom = geom.exterior
        coords = np.array(geom.coords)
        return np.all(np.abs(coords) == np.inf)

    to_drop = data.geometry.apply(drop_func)
    data = (data[~to_drop]).copy()

    def filt_func(geom):
        if hasattr(geom, "exterior"):
            geom = geom.exterior
        coords = np.array(geom.coords)
        return np.any(np.abs(coords) == np.inf) or (
            min_distance is not None and geom.distance(boundary) < min_distance
        )

    to_fix = data.index[data.geometry.apply(filt_func)]
    for index in to_fix:
        fixed = fill_edge_polygon(
            data.geometry.at[index],
            boundary,
            min_distance=min_distance,
        )
        data.geometry.at[index] = fixed
    return data


def fill_edge_polygon(geometry, boundary, min_distance=None):
    if isinstance(geometry, BaseMultipartGeometry):
        return type(geometry)(
            [fill_edge_polygon(i, boundary, min_distance) for i in geometry.geoms]
        )
    if not isinstance(geometry, Polygon):
        geometry = Polygon(geometry)
    coords = np.array(geometry.exterior.coords)

    segments_list = []
    segment = []
    for x, y in coords:
        if (np.abs(x) == np.inf or np.abs(y) == np.inf) or (
            min_distance is not None and boundary.distance(Point(x, y)) <= min_distance
        ):
            if len(segment) > 1:
                segments_list.append(segment)
                segment = []
            continue
        segment.append((x, y))
    if len(segments_list) == 0:
        return geometry
    segments_list = [LineString(i) for i in segments_list]

    out = []
    for i in range(-1, len(segments_list) - 1):
        segment_before = segments_list[i]
        point_before = Point(segment_before.coords[-1])

        segment_after = segments_list[i + 1]
        point_after = Point(segment_after.coords[0])

        d0 = boundary.project(point_before, normalized=True)
        d1 = boundary.project(point_after, normalized=True)
        boundary_segment = substring(boundary, d0, d1, normalized=True)

        if boundary_segment.length > 0.5 * boundary.length:
            if d1 > d0:
                seg0 = substring(boundary, d0, 0, normalized=True)
                seg1 = substring(boundary, 1, d1, normalized=True)
            else:
                seg0 = substring(boundary, d0, 1, normalized=True)
                seg1 = substring(boundary, 0, d1, normalized=True)
            boundary_segment = linemerge([seg0, seg1])

        if i == -1:
            out.append(segment_before)
        out.append(boundary_segment)
        if i != len(segments_list) - 2:
            out.append(segment_after)

    return Polygon(np.vstack([i.coords for i in out]))
