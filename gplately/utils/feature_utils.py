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

import pygplates

from ..geometry import pygplates_to_shapely


def shapelify_features(features, central_meridian=0.0, tessellate_degrees=None):
    """Generate Shapely `MultiPolygon` or `MultiLineString` geometries from reconstructed feature polygons.

    Notes
    -----
    Some Shapely polygons generated by `shapelify_features` cut longitudes of 180
    or -180 degrees. These features may appear unclosed at the dateline, so Shapely
    "closes" these polygons by connecting any of their open ends with lines. These
    lines may manifest on GeoAxes plots as horizontal lines that span the entire
    global extent. To prevent this, `shapelify_features` uses pyGPlates'
    [DateLineWrapper](https://www.gplates.org/docs/pygplates/generated/pygplates.datelinewrapper)
    to split a feature polygon into multiple closed polygons if it happens to cut the antimeridian.
    Another measure taken to ensure features are valid is to order exterior coordinates
    of Shapely polygons anti-clockwise.

    Parameters
    ----------
    features : iterable of <pygplates.Feature>, <ReconstructedFeatureGeometry> or <GeometryOnSphere>
        Iterable containing reconstructed polygon features.
    central_meridian : float
        Central meridian around which to perform wrapping; default: 0.0.
    tessellate_degrees : float or None
        If provided, geometries will be tessellated to this resolution prior to wrapping.

    Returns
    -------
    all_geometries : list of `shapely.geometry.BaseGeometry`
        Shapely geometries converted from the given reconstructed features. Any
        geometries at the dateline are split.

    See Also
    --------
    geometry.pygplates_to_shapely : convert PyGPlates geometry objects to
    Shapely geometries.
    """
    if isinstance(
        features,
        (
            pygplates.Feature,
            pygplates.ReconstructedFeatureGeometry,
            pygplates.GeometryOnSphere,
        ),
    ):
        features = [features]

    geometries = []
    for feature in features:
        if isinstance(feature, pygplates.Feature):
            geometries.extend(feature.get_all_geometries())
        elif isinstance(feature, pygplates.ReconstructedFeatureGeometry):
            geometries.append(feature.get_reconstructed_geometry())
        elif isinstance(feature, (pygplates.GeometryOnSphere, pygplates.LatLonPoint)):
            geometries.append(feature)
        elif isinstance(feature, pygplates.DateLineWrapper.LatLonMultiPoint):
            geometries.append(
                pygplates.MultiPointOnSphere(
                    [i.to_lat_lon() for i in feature.get_points()]
                )
            )
        elif isinstance(feature, pygplates.DateLineWrapper.LatLonPolyline):
            geometries.append(pygplates.PolylineOnSphere(feature.get_points()))
        elif isinstance(feature, pygplates.DateLineWrapper.LatLonPolygon):
            geometries.append(
                pygplates.PolygonOnSphere(
                    [i.to_lat_lon() for i in feature.get_exterior_points()]
                )
            )

    return [
        pygplates_to_shapely(
            i,
            force_ccw=True,
            validate=True,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
            explode=False,
        )
        for i in geometries
    ]
