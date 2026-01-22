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


def convert_polylines_to_polygons_within_feature(
    feature: pygplates.Feature,
    only_convert_closed_polylines=True,
    verify_information_model=pygplates.VerifyInformationModel.yes,
):
    """TODO: describe what this funciton do and explain the parameters"""
    # First step: gather all geometry property names (ignore other properties).
    geometry_property_names = []
    for property in feature:
        # If it's a geometry property then add the property name to the list.
        if property.get_value().get_geometry():
            geometry_property_names.append(property.get_name())

    # Second step: convert polylines to polygons (don't convert other geometry types).
    for geometry_property_name in geometry_property_names:
        # Get all geometries with the current property name.
        #
        # NOTE: I'm specifying the property name rather than relying on the default.
        geometries_with_property_name = feature.get_geometries(geometry_property_name)

        # Among the current geometries, convert any polylines to polygons (leave the other types alone).
        for geometry_index, geometry in enumerate(geometries_with_property_name):
            if isinstance(geometry, pygplates.PolylineOnSphere):
                # As a bonus, only convert polyline if it is closed (ie, its first and last vertices are equal).
                # Otherwise converting to a polygon will be dodgy (as Nicky pointed out).
                # if only_convert_closed_polylines has been set to False by users, we assume the users know what they are doing and respect the users
                if not only_convert_closed_polylines or geometry[0] == geometry[-1]:
                    polygon = pygplates.PolygonOnSphere(geometry)
                    geometries_with_property_name[geometry_index] = polygon

        # Set all geometries with the current property name.
        #
        # NOTE: I'm specifying the property name rather than relying on the default.
        feature.set_geometry(
            geometries_with_property_name,
            geometry_property_name,
            verify_information_model=verify_information_model,
        )


def convert_polylines_to_polygons(
    feature_collection: pygplates.FeatureCollection,
    only_convert_closed_polylines=True,
    verify_information_model=pygplates.VerifyInformationModel.yes,
):
    """convert all polylines in a given feature collection to polygons

    Parameters
    ----------
    feature_collection: pygplates.FeatureCollection
        all the polylines in this feature collection will be replaced with polygons
    """
    for feature in feature_collection:
        convert_polylines_to_polygons_within_feature(
            feature,
            only_convert_closed_polylines=only_convert_closed_polylines,
            verify_information_model=verify_information_model,
        )


def convert_polygons_to_polylines(feature_collection: pygplates.FeatureCollection):
    """similar to the function above"""
    # TODO:
    # get the exterior ring of the polygon and convert it to polyline?
    # Bianca, let's try to do this function
    # step 1: create a new function "convert_polygons_to_polylines_within_feature"
    #         the new function should be very similar to "convert_polylines_to_polygons_within_feature"
    #         You may copy the code from "convert_polylines_to_polygons_within_feature" to the new function and make changes to the copy
    # step 2: get the exterior ring and interior rings from the polygon and convert them to PolylineOnSphere
    #         John may have better ideas on this?
    #         see the links below for relavant pygplates doc
    # https://www.gplates.org/docs/pygplates/generated/pygplates.polygononsphere#pygplates.PolygonOnSphere.get_exterior_ring_points
    # https://www.gplates.org/docs/pygplates/generated/pygplates.polygononsphere#pygplates.PolygonOnSphere.get_interior_ring_points
    # https://www.gplates.org/docs/pygplates/generated/pygplates.polylineonsphere#pygplates.PolylineOnSphere.__init__
    # step 3: try to call this function in test_convert_geometries.py

    #
    # Please feel free to ask me or John if you got stuck at some point
    # mostly, this is a training task. John and I are obliged to provide support.
    #
    pass
