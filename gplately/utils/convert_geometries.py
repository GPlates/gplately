#
#    Copyright (C) 2024 The University of Sydney, Australia
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

# work in progress!!


def convert_polylines_to_polygons(feature_collection: pygplates.FeatureCollection):
    """convert all polylines  in a given feature collection to polygons

    Parameters
    ----------
    feature_collection: pygplates.FeatureCollection
        all the polylines in this feature collection will be replaced with polygons
    """
    for feature in feature_collection:
        geometry = feature.get_geometry()
        polygon = pygplates.PolygonOnSphere(geometry)
        feature.set_geometry(polygon)


def convert_polygons_to_polylines(feature_collection: pygplates.FeatureCollection):
    """similar to the function above"""
    # TODO:
    # get the exterior ring of the polygon and convert it to polyline?
    pass
