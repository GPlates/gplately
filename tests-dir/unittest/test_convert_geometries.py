#!/usr/bin/env python3

import pygplates

import gplately

print(gplately.__file__)

from gplately.utils import convert_geometries

# Bianca, finish this file to test the new functions

if __name__ == "__main__":

    filename = "replace this with the COB file name in muller2016"
    feature_collection = pygplates.FeatureCollection(filename)
    convert_geometries.convert_polylines_to_polygons(feature_collection)

    # now let's see if the geometries have been converted
    for feature in feature_collection:
        geometry = feature.get_geometry()
        print(geometry)

    # Write the modified feature collection to output file.
    pygplates.FeatureCollection(feature_collection).write(
        "converted-cob-muller2016.gpmlz"
    )
