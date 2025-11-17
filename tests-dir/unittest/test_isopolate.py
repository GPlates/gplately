from __future__ import print_function
import pygplates as pgp
import isopolate
import sys

# Read files and reconstruction time from bash script
ReconTime = int(sys.argv[1])
RotFile = str(sys.argv[2])
RidgeFile = str(sys.argv[3])
IsochronFile = str(sys.argv[4])
IsoCOBFile = str(sys.argv[5])
# Date = str(sys.argv[6])

print("     Current reconstruction time is", ReconTime, "Ma")

# This section gets the features that have not been subducted at the desired time

# -- For the isochron file
featureCollection_iso = pgp.FeatureCollectionFileFormatRegistry()
features_iso = featureCollection_iso.read(IsochronFile)
IsochronFile_subset = pgp.FeatureCollection()

for feature in features_iso:
    if feature.get_valid_time()[1] <= ReconTime:
        IsochronFile_subset.add(feature)

# -- For the ridge file
featureCollection_ri = pgp.FeatureCollectionFileFormatRegistry()
features_ri = featureCollection_ri.read(RidgeFile)
RidgeFile_subset = pgp.FeatureCollection()

for feature in features_ri:
    if feature.get_valid_time()[1] <= ReconTime:
        RidgeFile_subset.add(feature)

# -- For the isocob file
featureCollection_isocob = pgp.FeatureCollectionFileFormatRegistry()
features_isocob = featureCollection_isocob.read(IsoCOBFile)
IsoCOBFile_subset = pgp.FeatureCollection()

for feature in features_isocob:
    if feature.get_valid_time()[1] <= ReconTime:
        IsoCOBFile_subset.add(feature)

rotation_model = pgp.RotationModel(RotFile)

# This is where isopolate is actually called - note that these one line interpolates the isochrons,
# but does not reconstruct them

# The spacing between lines (in radians) at which to generate interpolated isochrons.
# Note that 'tessellate_threshold_radians' is the spacing ALONG lines.
#interval_spacing_radians = isopolate.DEFAULT_INTERPOLATE_RESOLUTION_RADIANS  # Default spacing (0.1 degrees).
interval_spacing_radians = 0.002

# Keywords arguments for the interpolate_isochrons() function.
interpolate_isochrons_kwargs = { }
interpolate_isochrons_kwargs['print_debug_output'] = 1
interpolate_isochrons_kwargs['tessellate_threshold_radians'] = interval_spacing_radians  # Make spacing the same along *and* between lines.
# This is where we enable/disable outputing of various scalar types.
# This determines what gets written to the '.xy' file.
# If they're all False then only the geometry (lat, lon) is written (to '.xy' file).
interpolate_isochrons_kwargs['output_scalar_age'] = True
interpolate_isochrons_kwargs['output_scalar_spreading_rate'] = True
interpolate_isochrons_kwargs['output_scalar_spreading_asymmetry'] = True
interpolate_isochrons_kwargs['output_scalar_full_spreading_rate'] = True
interpolate_isochrons_kwargs['output_scalar_spreading_direction'] = True
interpolate_isochrons_kwargs['output_scalar_spreading_obliquity'] = True


output_features = isopolate.interpolate_isochrons(
        rotation_model,
        [RidgeFile_subset, IsochronFile_subset, IsoCOBFile_subset],
        interval_spacing_radians,
        # Unpack the keyword arguments dict into keyword arguments...
        **interpolate_isochrons_kwargs)

# The scalar types exported to the '.xy' file - depends on 'interpolate_isochrons_kwargs' settings above.
output_scalar_types = []
if interpolate_isochrons_kwargs['output_scalar_age']:
    output_scalar_types.append(pgp.ScalarType.create_gpml('Age'))
if interpolate_isochrons_kwargs['output_scalar_spreading_rate']:
    output_scalar_types.append(pgp.ScalarType.create_gpml('SpreadingRate'))
if interpolate_isochrons_kwargs['output_scalar_spreading_asymmetry']:
    output_scalar_types.append(pgp.ScalarType.create_gpml('SpreadingAsymmetry'))
if interpolate_isochrons_kwargs['output_scalar_full_spreading_rate']:
    output_scalar_types.append(pgp.ScalarType.create_gpml('FullSpreadingRate'))
if interpolate_isochrons_kwargs['output_scalar_spreading_direction']:
    output_scalar_types.append(pgp.ScalarType.create_gpml('SpreadingDirection'))
if interpolate_isochrons_kwargs['output_scalar_spreading_obliquity']:
    output_scalar_types.append(pgp.ScalarType.create_gpml('SpreadingObliquity'))

# We are outputting scalar coverages to '.xy' format.
#
# We handle this as a special case so we can write out the scalar values
# after the xy (lat/lon) values.
#
# We write the scalar types in the same order as they appear in 'output_scalar_types'.
#
isopolate.write_coverage_features_to_xy_file(
        'InterpolatedIsochrons_' + str(ReconTime) + 'Ma.xy',
        output_features,
        output_scalar_types,
        rotation_model,
        reconstruction_time=ReconTime,
        print_debug_output=1)

# the output features only exist in memory, these lines save to a file
OutFile = './InterpolatedIsochrons.gpml'
file_registry = pgp.FeatureCollectionFileFormatRegistry()
file_registry.write(pgp.FeatureCollection(output_features), OutFile)

# Here we do the last step that the isopolate command-line script does, reconstructing the interpolated
# isochrons to the selected time. The output file format is picked up from the file name
# extension (eg .gmt, .shp). The interpolated isochrons are created in a 'tmp' file
# print("     Creating interpolated isochrons at", ReconTime, "Ma")
#FNAME = './tmp/InterpolatedIsochrons_'+str(ReconTime)+'Ma.shp'
# FNAME = 'InterpolatedIsochrons_'+str(ReconTime)+'Ma.shp'
# pgp.reconstruct(OutFile, RotFile, FNAME, ReconTime, 0)
# print("     Creating interpolated isochrons at", ReconTime, "Ma  --- done!")
