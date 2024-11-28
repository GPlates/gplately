#!/usr/bin/env python3

import pygplates

import gplately

print(gplately.__file__)

from gplately.utils import convert_geometries

# Bianca, finish this file to test the new functions

if __name__ == "__main__":

    filename = "Global_230-0Ma_GK07_AREPS_COB_Terranes.gpml"
    feature_collection = pygplates.FeatureCollection(filename)
    geometry_types = []
    for feature in feature_collection:
        geoms = feature.get_all_geometries()
        for geometry in geoms:
            type_name = type(geometry).__name__
            if type_name not in geometry_types:
                geometry_types.append(type_name)

    print("feature types before conversion: ", geometry_types)

    convert_geometries.convert_polylines_to_polygons(feature_collection)

    # now let's see if the geometries have been converted
    geometry_types = []
    for feature in feature_collection:
        geoms = feature.get_all_geometries()
        for geometry in geoms:
            type_name = type(geometry).__name__
            # print(type_name)
            if type_name == "PolylineOnSphere":
                print(geoms)
            if type_name not in geometry_types:
                geometry_types.append(type_name)

    print("feature types after conversion: ", geometry_types)

    # Write the modified feature collection to output file.
    feature_collection.write("converted-cob-muller2016.gpmlz")
