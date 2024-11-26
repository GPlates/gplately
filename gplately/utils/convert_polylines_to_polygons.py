import os.path
import pygplates


input_filenames = [
    'Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons_2019_v1.shp'
]

output_filename_append = '_polygon'

for input_filename in input_filenames:
    # Convert polylines to polygons in the features.
    feature_collection = pygplates.FeatureCollection(input_filename)
    for feature in feature_collection:
        # polygons = [pygplates.PolygonOnSphere(geometry) for geometry in feature.get_geometries()]
        # changed feature get statement to 'all_geometries' - since feature has issues getting different
        # types of geometries
        for geometry in feature.get_all_geometries():
            polygons = pygplates.PolygonOnSphere(geometry)
            feature.set_geometry(polygons)
    
    # Get the output (polygon) filename from the input (polyline) filename.
    input_filename_root, input_filename_ext = os.path.splitext(input_filename)
    output_filename = ''.join((input_filename_root, output_filename_append, input_filename_ext))
    
    # Write the modified feature collection to output file.
    pygplates.FeatureCollection(feature_collection).write(output_filename)