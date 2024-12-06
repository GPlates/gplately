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

    # Bianca, let's create a .gpml file and put all PolylineOnSphere geometries into the .gpml file
    # we would like to see what these PolylineOnSphere geometries look like
    # step 1: create a new FeatureCollection, something like "new_feature_collection = pygplates.FeatureCollection()"
    # step 2: loop through the original feature collection, something like "for feature in feature_collection:"
    # step 3: for each feature, get all geometries, something like "geoms = feature.get_all_geometries()"
    # step 4: loop through all geometries with an inner loop, something like "for geometry in geoms:"
    # step 5: check if the geometry is a PolylineOnSphere, something like "if isinstance(geometry, pygplates.PolylineOnSphere):"
    # step 6: for each PolylineOnSphere geometries, create a new feature, something like "new_feature  = pygplates.Feature()"
    # step 7: put the PolylineOnSphere geometry into the new feature, something like "new_feature.set_geometry(geometry)"
    # step 8: add the new feature into the new feature collection, something like "new_feature_collection.add(new_feature)"
    # step 9: save the new feature collection to a .gpml file
    # step 10: open the .gpml file in GPlates desktop software. https://www.earthbyte.org/download-gplates-2-5/

    #############
    new_feature_collection = pygplates.FeatureCollection()

    for feature in feature_collection:
        geoms = feature.get_all_geometries()
        for geometry in geoms:
            if isinstance(geometry, pygplates.PolylineOnSphere):
                new_feature = pygplates.Feature()
                new_feature.set_geometry(geometry)
                new_feature_collection.add(new_feature)

    new_feature_collection.write("polylines-in-muller2016.gpml")
    #############

    convert_geometries.convert_polylines_to_polygons(
        feature_collection,
        only_convert_closed_polylines=False,
        verify_information_model=pygplates.VerifyInformationModel.no,
    )

    # now let's see if the geometries have been converted
    geometry_types = []
    for feature in feature_collection:
        geoms = feature.get_all_geometries()
        for geometry in geoms:
            type_name = type(geometry).__name__
            # print(type_name)
            # if type_name == "PolylineOnSphere":
            #    print(geoms)
            if type_name not in geometry_types:
                geometry_types.append(type_name)

    print("feature types after conversion: ", geometry_types)

    # Write the modified feature collection to output file.
    feature_collection.write("converted-cob-muller2016.gpmlz")
