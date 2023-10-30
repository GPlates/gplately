
"""
    Copyright (C) 2015 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

from __future__ import print_function
import argparse
import codecs
import sys
import os.path
import re

import pygplates

DEFAULT_OUTPUT_FILENAME_EXTENSION = 'gpml'


# Class to hold a feature's type, list of property name/value pairs and list of geometry points.
class _FeatureData(object):
    def __init__(self, line_number):
        self.line_number = line_number
        self.type = None # Default to 'unclassified' feature.
        self.properties = []
        self.points_data = []


def _read_feature_metadata(feature_data, line):
    
    # See if line looks like "name=value".
    line_data = line.split('=')
    if len(line_data) != 2:
        # It's not a metadata line (does not have a single '=' char).
        # So just ignore it.
        return
    
    name = line_data[0].strip().lower() # Convert to lower case.
    value = line_data[1].strip()
    
    # Convert 'unicode' back to UTF8-encoded 'str' since
    # pyGPlates revision 12 (and below) cannot accept 'unicode'.
    value = value.encode('utf-8')
    
    # See if specifying feature type or a single feature property.
    if name == 'featuretype':
        feature_data.type = pygplates.FeatureType.create_gpml(value)
    elif name == 'name':
        property_name = pygplates.PropertyName.create_gml('name')
        property_value = pygplates.XsString(value)
        feature_data.properties.append((property_name, property_value))
    elif name == 'description':
        property_name = pygplates.PropertyName.create_gml('description')
        property_value = pygplates.XsString(value)
        feature_data.properties.append((property_name, property_value))
    elif name == 'reconstructionplateid':
        property_name = pygplates.PropertyName.create_gpml('reconstructionPlateId')
        try:
            property_value = pygplates.GpmlPlateId(int(value))
            feature_data.properties.append((property_name, property_value))
        except ValueError:
            print(u'Line {0}: Ignoring feature property - property name "{1}" expects an integer value.'
                    .format(feature_data.line_number, property_name.to_qualified_string()), file=sys.stderr)
    elif name == 'validtime':
        property_name = pygplates.PropertyName.create_gml('validTime')
        try:
            valid_time_match = re.match(r'\s*\(\s*(\S*)\s*,\s*(\S*)\s*\)\s*', value)
            if not valid_time_match:
                raise ValueError('Not a 2-tuple')
            begin_time = float(valid_time_match.group(1))
            end_time = float(valid_time_match.group(2))
            property_value = pygplates.GmlTimePeriod(begin_time, end_time)
            feature_data.properties.append((property_name, property_value))
        except ValueError:
            print(u'Line {0}: Ignoring feature property - property name "{1}" expects "(begin_time, end_time)".'
                    .format(feature_data.line_number, property_name.to_qualified_string()), file=sys.stderr)
        except pygplates.GmlTimePeriodBeginTimeLaterThanEndTimeError:
            print(u'Line {0}: Ignoring feature property - property name "{1}" requires begin time to be older than end time".'
                    .format(feature_data.line_number, property_name.to_qualified_string()), file=sys.stderr)
    else:
        print(u'Line {0}: Ignoring feature property - "{1}" is not a recognised feature property name.'
                .format(feature_data.line_number, name), file=sys.stderr)


def _create_feature(feature_data, multi_geometry_type, lon_lat_order, scalar_types):
    
    # Create a new feature.
    try:
        feature = pygplates.Feature(feature_data.type)
    except pygplates.InformationModelError:
        print(u'Line {0}: Ignoring feature - "{1}" is not a recognised feature type.'
                .format(feature_data.line_number, feature_data.type.to_qualified_string()), file=sys.stderr)
        return
    
    # Add any feature properties.
    for property_name, property_value in feature_data.properties:
        try:
            feature.add(property_name, property_value)
        except pygplates.InformationModelError:
            print(u'Line {0}: Ignoring feature property - "{1}" is not a feature property name supported by feature type "{2}".'
                    .format(feature_data.line_number, property_name.to_qualified_string(), feature_data.type.to_qualified_string()), file=sys.stderr)
    
    points = []
    
    # Store points as (lat, lon) tuples (ie, in lat/lon order) since pygplates uses that order.
    if lon_lat_order:
        for feature_point_data in feature_data.points_data:
            point = (feature_point_data[1], feature_point_data[0])
            points.append(point)
    else:
        for feature_point_data in feature_data.points_data:
            point = feature_point_data[0], feature_point_data[1]
            points.append(point)
    
    # If only one point then store as a point, otherwise the requested multi-geometry type.
    if len(points) == 1:
        geometry = pygplates.PointOnSphere(points[0])
    else:
        geometry = multi_geometry_type(points)
    
    if scalar_types:
        # Store the scalar values in a dict using scalar type as key.
        scalar_coverages = {}
        for scalar_type_index, scalar_type in enumerate(scalar_types):
            scalars = [feature_point_data[scalar_type_index + 2]
                    for feature_point_data in feature_data.points_data]
            scalar_coverages[scalar_type] = scalars
        
        # Save geometry as a coverage.
        coverage = (geometry, scalar_coverages)
        feature.set_geometry(coverage)
    else:
        feature.set_geometry(geometry)
    
    return feature


def import_geometry_from_xy_file(input_geometry_filename, multi_geometry_type=pygplates.PolylineOnSphere, lon_lat_order=True, scalar_types=None):
    """Parse ascii file containing latitude/longitude geometries.
    
    Returns a list of imported features, or None on parse error.
    
    If an input ascii file contains '>' separators at the beginning of lines then these indicate a new geometry
    (there can be multiple consecutive '>' lines with arbitrary text, or feature metadata, to delineate two geometries).
    If there is only one 'lon lat' line between '>' markers then a point geometry is created.
    Multiple consecutive 'lon lat' lines are combined into a single polyline (or multipoint/polygon depending on multi-geometry type specified).
    For example, the following is one point geometry and one line geometry:
    
    >
    > This is a point geometry...
    0 0
    > This is a polyline (or multipoint/polygon) geometry...
    0 0
    0 10
    
    Furthermore, each line starting with '>' can optionally contain arbitrary text, or feature metadata.
    Feature metadata can be the type of feature or a feature property.
    The following example shows all currently supported metadata...
    
    >
    > An arbitrary text line does not contain the '=' character.
    > Whereas each metadata text line looks like 'name=value'...
    >
    > FeatureType = Coastline
    > Name = Australia
    > Description = Australian coastline
    > ReconstructionPlateId = 801
    > ValidTime = (4540, -inf)
    141.6863 -12.5592
    ...
    
    ...where all metadata is optional.
    If 'FeatureType' is not specified then it defaults to an 'unclassified' feature.
    Note that 'ValidTime' supports 'inf' and '-inf' which represent distant past and future respectively.
    
    If there are no '>' markers (ie, not actually a GMT '.xy' file) then each point will be a separate geometry
    (instead of one large polyline geometry), except if the multipoint type is 'pygplates.MultiPointOnSphere'
    (which results in one large multipoint).
    This treats the ascii file simply as a list of points.
    For example, the following is three point geometries:
    
    0 0
    0 10
    0 20
    
    By default a GMT '.xy' file has lon/lat order.
    Specifying 'lon_lat_order' as False changes this to lat/lon order.
    
    Specify 'scalar_types' to convert extra scalar values to coverage data.
    If 'scalar_types' is specified it should be a sequence of 'pygplates.ScalarType' instances.
    Each line in the file containing a point has a longitude and a latitude value as a minimum.
    Extra values are scalar values associated with each point.
    This option specifies the scalar types of each extra field.
    For example:
    
    import_geometry_from_xy_file(
        'geometries.xy',
        scalar_types = [
            pygplates.ScalarType.create_gpml('SpreadingRate'), 
            pygplates.ScalarType.create_gpml('SpreadingDirection'),
            pygplates.ScalarType.create_gpml('SpreadingAsymmetry')])
    
    ...could represent three extra scalar values for each point (in that order).
    So a line containing:
    
    10 20 15 80 48
    
    ...would have longitude=10, latitude=20, SpreadingRate=15, SpreadingDirection=80 and SpreadingAsymmetry=48.
    """
    
    # At minimum each point has a longitude and latitude.
    num_scalars = 2
    if scalar_types:
        num_scalars += len(scalar_types)
    
    num_lines_with_more_scalars_than_expected = 0
    
    features = []
    
    feature_data = _FeatureData(1) # First feature on line 1
    encountered_marker = False # Encountered a line starting with '>'.
    
    # Assume file is encoded as UTF8 (which includes basic 7-bit ascii).
    with codecs.open(input_geometry_filename, 'r', 'utf-8') as geometry_file:
        for line_number, line in enumerate(geometry_file):
            # Make line number 1-based instead of 0-based.
            line_number = line_number + 1
            
            line = line.strip()
            
            # See if line begins with '>'.
            if line.startswith('>'):
                # If have points then generate the previous feature.
                if feature_data.points_data:
                    feature = _create_feature(feature_data, multi_geometry_type, lon_lat_order, scalar_types)
                    if feature:
                        features.append(feature)
                    
                    # Reset for next feature.
                    feature_data = _FeatureData(line_number)

                # Attempt to read the feature type or a single feature property (one property per line).
                _read_feature_metadata(feature_data, line.lstrip('>'))
                
                encountered_marker = True
                
                # Skip to next line.
                continue
            
            point_data = line.split()
            
            num_scalars_in_point = len(point_data)
            if num_scalars_in_point < num_scalars:
                # If just a line containing white-space then skip to next line.
                if num_scalars_in_point == 0:
                    continue
                
                print(u'Line {0}: Ignoring file "{1}" - point has fewer than {2} white-space separated numbers.'
                        .format(line_number, input_geometry_filename, num_scalars), file=sys.stderr)
                return
            
            # If we have more scalars than expected then we'll print a warning message before returning.
            if num_scalars_in_point > num_scalars:
                num_lines_with_more_scalars_than_expected += 1
            
            try:
                # Convert strings to numbers.
                feature_point_data = [float(scalar_string) for scalar_string in point_data]
                feature_data.points_data.append(feature_point_data)
            except ValueError:
                print(u'Line {0}: Ignoring file "{1}" - cannot read point longitude and latitude values '
                    '(and optional scalar values).'
                        .format(line_number, input_geometry_filename), file=sys.stderr)
                return
    
    # If have points then either create the last feature (if encountered '>' markers) or
    # create a feature for each point in the input file.
    if feature_data.points_data:
        if encountered_marker: # last feature
            feature = _create_feature(feature_data, multi_geometry_type, lon_lat_order, scalar_types)
            if feature:
                features.append(feature)
        else:
            # There were no '>' markers.
            # Treat it like a non-GMT file that is just a list of points.
            # If the multi-geometry type is a multipoint then output a single multipoint,
            # otherwise output multiple points.
            if multi_geometry_type == pygplates.MultiPointOnSphere:
                # Create a single multipoint feature.
                feature = _create_feature(feature_data, multi_geometry_type, lon_lat_order, scalar_types)
                if feature:
                    features.append(feature)
            else:
                # Create a feature for each point.
                line_number = 1 # Make line number 1-based instead of 0-based.
                for feature_point_data in feature_data.points_data:
                    point_feature_data = _FeatureData(line_number)
                    point_feature_data.points_data = [feature_point_data]
                    
                    feature = _create_feature(point_feature_data, multi_geometry_type, lon_lat_order, scalar_types)
                    if feature:
                        features.append(feature)
                    
                    line_number = line_number + 1
    
    # Print a warning if any lines contained points with more numbers than expected.
    if num_lines_with_more_scalars_than_expected:
        print(u'Warning: There were {0} points in file {1} containing greater than {2} white-space '
            'separated numbers - extra numbers were ignored.'
                .format(num_lines_with_more_scalars_than_expected, input_geometry_filename, num_scalars),
            file=sys.stderr)
    
    return features


if __name__ == "__main__":

    # Check the imported pygplates version.
    required_version = pygplates.Version(7)
    if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < required_version:
        print(u'{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
                os.path.basename(__file__), pygplates.Version.get_imported_version(), required_version),
            file=sys.stderr)
        sys.exit(1)
    
    lat_lon_order_option = ('-l', '--lat_lon_order')
    multipoint_option = ('-m', '--multipoint')
    polygon_option = ('-p', '--polygon')
    scalar_coverages_option = ('-s', '--scalar_coverages')
    
    __description__ = \
    """Converts geometry in one or more input ascii files (such as '.xy' files) to output files suitable for loading into GPlates.
    
    If an input ascii file contains '>' separators at the beginning of lines then these indicate a new geometry
    (there can be multiple consecutive '>' lines with arbitrary text, or feature metadata, to delineate two geometries).
    If there is only one 'lon lat' line between '>' markers then a point geometry is created.
    Multiple consecutive 'lon lat' lines are combined into a single polyline or a single multipoint (if the {1} option is specified)
    or a single polygon (if the {2} option is specified).
    For example, the following is one point geometry and one line geometry:
    
    >
    > This is a point geometry...
    0 0
    > This is a polyline (or multipoint/polygon) geometry...
    0 0
    0 10
    
    Furthermore, each line starting with '>' can optionally contain arbitrary text, or feature metadata.
    Feature metadata can be the type of feature or a feature property.
    The following example shows all currently supported metadata...
    
    >
    > An arbitrary text line does not contain the '=' character.
    > Whereas each metadata text line looks like 'name=value'...
    >
    > FeatureType = Coastline
    > Name = Australia
    > Description = Australian coastline
    > ReconstructionPlateId = 801
    > ValidTime = (4540, -inf)
    141.6863 -12.5592
    ...
    
    ...where all metadata is optional.
    If 'FeatureType' is not specified then it defaults to an 'unclassified' feature.
    Note that 'ValidTime' supports 'inf' and '-inf' which represent distant past and future respectively.
    
    If there are no '>' markers (ie, not actually a GMT '.xy' file) then each point will be a separate geometry
    (instead of one large polyline geometry), except if the {1} option is specified (which results in one large multipoint).
    This treats the ascii file simply as a list of points.
    For example, the following is three point geometries:
    
    0 0
    0 10
    0 20
    
    The name of each output file is the associated input filename with a different filename extension.
    For example:
       'data/test.xy' -> 'data/test.{0}'
    
    By default a GMT '.xy' file has lon/lat order.
    This can be changed to lat/lon order by specifying the {3} option.
    
    Specify the {4} option to convert extra scalar values to coverage data.
    If specified it should be a sequence of 'pygplates.ScalarType' instances.
    Each line in the file containing a point has a longitude and a latitude value as a minimum.
    Extra values are scalar values associated with each point.
    This option specifies the scalar types of each extra field.
    For example:
    
    import_geometry_from_xy_file(
        'geometries.xy',
        scalar_types = [
            pygplates.ScalarType.create_gpml('SpreadingRate'), 
            pygplates.ScalarType.create_gpml('SpreadingDirection'),
            pygplates.ScalarType.create_gpml('SpreadingAsymmetry')])
    
    ...could represent three extra scalar values for each point (in that order).
    So a line containing:
    
    10 20 15 80 48
    
    ...would have longitude=10, latitude=20, SpreadingRate=15, SpreadingDirection=80 and SpreadingAsymmetry=48.


    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -e shp -- input1.xy input2.xy
     """.format(DEFAULT_OUTPUT_FILENAME_EXTENSION, multipoint_option, polygon_option, lat_lon_order_option, scalar_coverages_option)

    # The command-line parser.
    parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-e', '--output_filename_extension', type=str,
            default='{0}'.format(DEFAULT_OUTPUT_FILENAME_EXTENSION),
            help="The filename extension of each input filename is changed to get each output filename "
                "- the default extension is '{0}' - supported extensions include 'gpml', 'shp', 'gmt', 'dat'."
                .format(DEFAULT_OUTPUT_FILENAME_EXTENSION))
    
    parser.add_argument(lat_lon_order_option[0], lat_lon_order_option[1],
            action='store_true',
            help="By default a GMT '.xy' file stores each point in lon/lat order - specifying this option "
                "interprets each point in lat/lon order instead - default order is lon/lat.")

    geometry_option_group = parser.add_mutually_exclusive_group()
    geometry_option_group.add_argument(multipoint_option[0], multipoint_option[1],
            action='store_true',
            help="Create multipoint, instead of polyline, geometries (when using '>' markers). "
                "Also, if there are no '>' markers, then create one multipoint geometry instead multiple point geometries.")
    geometry_option_group.add_argument(polygon_option[0], polygon_option[1],
            action='store_true',
            help="Create polygon, instead of polyline, geometries (when using '>' markers).")
    
    parser.add_argument(scalar_coverages_option[0], scalar_coverages_option[1], type=str, nargs='+',
            metavar='scalar_type',
            help='Convert extra scalar values to coverage data. '
                'NOTE: Scalar coverage data is only saved to ".gpml" format. '
                'Each line in the file containing a point has a longitude and a latitude value as a minimum. '
                'Extra values are scalar values associated with each point. '
                'This option specifies the scalar types of each extra field. '
                'For example SpreadingRate SpreadingDirection SpreadingAsymmetry could represent three '
                'extra scalar values for each point (in that order) - so a line containing '
                '"10 20 15 80 48" would have longitude=10, latitude=20, SpreadingRate=15, '
                'SpreadingDirection=80 and SpreadingAsymmetry=48.')
    
    def unicode_filename(value_string):
        try:
            # Filename uses the system encoding - decode from 'str' to 'unicode'.
            filename = value_string.decode(sys.getfilesystemencoding())
        except UnicodeDecodeError:
            raise argparse.ArgumentTypeError("Unable to convert filename %s to unicode" % value_string)
        
        return filename
    
    parser.add_argument('input_filenames', type=unicode_filename, nargs='+',
            metavar='input_filename',
            help='The ascii input files containing the geometry in latitude/longitude coordinates.')
    
    # Parse command-line options.
    args = parser.parse_args()
    
    # Determine whether to use multipoints or polylines or polygons.
    if args.multipoint:
        multi_geometry_type = pygplates.MultiPointOnSphere
    elif args.polygon:
        multi_geometry_type = pygplates.PolygonOnSphere
    else:
        multi_geometry_type = pygplates.PolylineOnSphere
    
    # Import the input geometry files.
    lon_lat_order = not args.lat_lon_order
    if args.scalar_coverages:
        # Convert each string to a pygplates.ScalarType.
        scalar_types = [pygplates.ScalarType.create_gpml(scalar_coverage)
                for scalar_coverage in args.scalar_coverages]
    else:
        scalar_types = None
    imported_geometry_feature_collections = [
            import_geometry_from_xy_file(input_filename, multi_geometry_type, lon_lat_order, scalar_types)
            for input_filename in args.input_filenames]
    
    # Write the imported geometry feature collections to disk.
    for feature_collection_index in range(len(imported_geometry_feature_collections)):
        imported_geometry_feature_collection = imported_geometry_feature_collections[feature_collection_index]
        
        # Skip files that failed to import.
        if imported_geometry_feature_collection is None:
            continue

        # Each output filename is the input filename with a different extension.
        input_filename = args.input_filenames[feature_collection_index]
        filename_root, filename_ext = os.path.splitext(input_filename)
        output_filename = ''.join((filename_root, '.', args.output_filename_extension))
        
        # Convert output filename from 'unicode' to UTF8-encoded 'str' since pyGPlates versions
        # prior to revision 13 cannot accept 'unicode' strings.
        output_filename = output_filename.encode('utf-8')
        
        # Write the output file.
        pygplates.FeatureCollection(imported_geometry_feature_collection).write(output_filename)
