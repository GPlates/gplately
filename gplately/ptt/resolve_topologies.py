
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
import math
from . import separate_ridge_transform_segments
import sys
import os.path
import pygplates

# PyGPlates version 26 supports default anchor plate ID in pygplates.RotationModel.__init__().
USING_PYGPLATES_VERSION_GREATER_EQUAL_26 = (hasattr(pygplates, 'Version') and pygplates.Version.get_imported_version() >= pygplates.Version(26))
# PyGPlates version 29 supports passing a RotationModel to RotationModel.__init__
USING_PYGPLATES_VERSION_GREATER_EQUAL_29 = (hasattr(pygplates, 'Version') and pygplates.Version.get_imported_version() >= pygplates.Version(29))


DEFAULT_OUTPUT_FILENAME_PREFIX = 'topology_'
DEFAULT_OUTPUT_FILENAME_EXTENSION = 'shp'


def resolve_topologies(
        rotation_features_or_model,
        topology_features,
        time,
        output_filename_prefix,
        output_filename_extension,
        transform_segment_deviation_in_radians = separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS,
        anchor_plate_id=None):
    """
    Resolves topologies at specified time and saves (to separate files) the resolved topologies, and their boundary sections as subduction zones,
    mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges).
    
    rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).
    
    topology_features: Topology feature collection(s), or list of features, or filename(s) or any combination of those.
    
    time: Reconstruction time to resolved topologies.
    
    transform_segment_deviation_in_radians: How much a mid-ocean ridge segment can deviate from the stage pole before
                                            it's considered a transform segment (in radians).

    anchor_plate_id: Anchor plate ID (optional, int).
    
    Writes output files containing the following features...
            - resolved topology features (topological plates and networks)
            - ridge and transform boundary sections (resolved features)
            - ridge boundary sections (resolved features)
            - transform boundary sections (resolved features)
            - subduction boundary sections (resolved features)
            - left subduction boundary sections (resolved features)
            - right subduction boundary sections (resolved features)
            - other boundary sections (resolved features) that are not subduction zones or mid-ocean ridges (ridge/transform)
    """
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = _parse_rotation_args(
        rotation_features_or_model,
        anchor_plate_id=anchor_plate_id,
    )

    # Get topology features (could include filenames).    
    topology_features = pygplates.FeaturesFunctionArgument(topology_features).get_features()
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    (resolved_topology_features,
     ridge_transform_boundary_section_features,
     ridge_boundary_section_features,
     transform_boundary_section_features,
     subduction_boundary_section_features,
     left_subduction_boundary_section_features,
     right_subduction_boundary_section_features,
     other_boundary_section_features) = resolve_topologies_into_features(
            rotation_model,
            topology_features,
            time,
            transform_segment_deviation_in_radians,
            anchor_plate_id=anchor_plate_id)
    
    # Write each list of features to a separate file.
    _write_resolved_topologies(
        time,
        output_filename_prefix,
        output_filename_extension,
        resolved_topology_features,
        ridge_transform_boundary_section_features,
        ridge_boundary_section_features,
        transform_boundary_section_features,
        subduction_boundary_section_features,
        left_subduction_boundary_section_features,
        right_subduction_boundary_section_features,
        other_boundary_section_features)


def resolve_topologies_into_features(
        rotation_features_or_model,
        topology_features,
        time,
        transform_segment_deviation_in_radians = separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS,
        anchor_plate_id=None):
    """
    Resolves topologies at specified time and returns resolved topologies and their boundary sections as subduction zones,
    mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges).
    
    rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).
    
    topology_features: Topology feature collection(s), or list of features, or filename(s) or any combination of those.
    
    time: Reconstruction time to resolved topologies.
    
    transform_segment_deviation_in_radians: How much a mid-ocean ridge segment can deviate from the stage pole before
                                            it's considered a transform segment (in radians).

    anchor_plate_id: Anchor plate ID (optional, int).
    
    Returns: A tuple containing the following lists...
            - resolved topology features (topological plates and networks)
            - ridge and transform boundary sections (resolved features)
            - ridge boundary sections (resolved features)
            - transform boundary sections (resolved features)
            - subduction boundary sections (resolved features)
            - left subduction boundary sections (resolved features)
            - right subduction boundary sections (resolved features)
            - other boundary sections (resolved features) that are not subduction zones or mid-ocean ridges (ridge/transform)
    """
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = _parse_rotation_args(
        rotation_features_or_model,
        anchor_plate_id=anchor_plate_id,
    )
    time = float(time)

    # Turn topology data into a list of features (if not already).
    topology_features = pygplates.FeaturesFunctionArgument(topology_features)
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(
        topology_features.get_features(),
        rotation_model,
        resolved_topologies,
        time,
        shared_boundary_sections,
        anchor_plate_id=anchor_plate_id,
    )

    # We'll create a feature for each boundary polygon feature and each type of
    # resolved topological section feature we find.
    resolved_topology_features = []
    ridge_transform_boundary_section_features = []
    ridge_boundary_section_features = []
    transform_boundary_section_features = []
    subduction_boundary_section_features = []
    left_subduction_boundary_section_features = []
    right_subduction_boundary_section_features = []
    other_boundary_section_features = []

    # Iterate over the resolved topologies.
    for resolved_topology in resolved_topologies:
        resolved_topology_features.append(resolved_topology.get_resolved_feature())

    # Iterate over the shared boundary sections.
    for shared_boundary_section in shared_boundary_sections:
        
        # Get all the geometries of the current boundary section.
        boundary_section_features = [shared_sub_segment.get_resolved_feature()
                for shared_sub_segment in shared_boundary_section.get_shared_sub_segments()]
        
        # Add the feature to the correct list depending on feature type, etc.
        if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.gpml_subduction_zone:
            
            # Put all subduction zones in one collection/file.
            subduction_boundary_section_features.extend(boundary_section_features)
            
            # Also put subduction zones in left/right collection/file.
            polarity_property = shared_boundary_section.get_feature().get(
                    pygplates.PropertyName.create_gpml('subductionPolarity'))
            if polarity_property:
                polarity = polarity_property.get_value().get_content()
                if polarity == 'Left':
                    left_subduction_boundary_section_features.extend(boundary_section_features)
                elif polarity == 'Right':
                    right_subduction_boundary_section_features.extend(boundary_section_features)

        elif shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.gpml_mid_ocean_ridge:
            ridge_transform_boundary_section_features.extend(boundary_section_features)
            
            # Find the stage rotation of the MOR in the frame of reference of its reconstructed geometry at the current 'time'.
            # The stage pole can then be directly geometrically compared to the reconstructed spreading geometry.
            spreading_stage_rotation = separate_ridge_transform_segments.get_stage_rotation_for_reconstructed_geometry(
                    shared_boundary_section.get_feature(), rotation_model, time)
            if spreading_stage_rotation:
                boundary_section_geometries = [shared_sub_segment.get_resolved_geometry()
                        for shared_sub_segment in shared_boundary_section.get_shared_sub_segments()]
                for boundary_section_index, boundary_section_feature in enumerate(boundary_section_features):
                    # Split into ridge and transform segments.
                    ridge_and_transform_segment_geometries = separate_ridge_transform_segments.separate_geometry_into_ridges_and_transforms(
                            spreading_stage_rotation,
                            boundary_section_geometries[boundary_section_index],
                            transform_segment_deviation_in_radians)
                    if ridge_and_transform_segment_geometries:
                        ridge_sub_segment_geometries, transform_sub_segment_geometries = ridge_and_transform_segment_geometries

                        if ridge_sub_segment_geometries:
                            ridge_boundary_section_feature = boundary_section_feature.clone()
                            ridge_boundary_section_feature.set_geometry(ridge_sub_segment_geometries)
                            ridge_boundary_section_features.append(ridge_boundary_section_feature)

                        if transform_sub_segment_geometries:
                            transform_boundary_section_feature = boundary_section_feature.clone()
                            transform_boundary_section_feature.set_geometry(transform_sub_segment_geometries)
                            transform_boundary_section_features.append(transform_boundary_section_feature)

                    # else skip shared sub segment - it's not a polyline (or polygon).
                
            # else most likely MOR doesn't have left/right plate IDs or reconstruction/conjugate plate IDs,
            # so don't split into ridge and transform segments.

        else:
            other_boundary_section_features.extend(boundary_section_features)

    return (resolved_topology_features,
            ridge_transform_boundary_section_features,
            ridge_boundary_section_features,
            transform_boundary_section_features,
            subduction_boundary_section_features,
            left_subduction_boundary_section_features,
            right_subduction_boundary_section_features,
            other_boundary_section_features)


def find_total_boundary_length_in_kms(
        ridge_transform_boundary_section_features,
        ridge_boundary_section_features,
        transform_boundary_section_features,
        subduction_boundary_section_features,
        left_subduction_boundary_section_features,
        right_subduction_boundary_section_features,
        other_boundary_section_features):
    """
    Find the total length (in kms) of resolved topology boundary sections (in the various categories).
    """

    def find_total_length(features):
        total_length_features = 0.0
        for feature in features:
            for geometry in feature.get_geometries():
                total_length_features += geometry.get_arc_length()
        # Return total lengths in kms (converted from radians).
        return  total_length_features * pygplates.Earth.mean_radius_in_kms

    # Return total lengths in kms (converted from radians).
    return (find_total_length(ridge_transform_boundary_section_features),
            find_total_length(ridge_boundary_section_features),
            find_total_length(transform_boundary_section_features),
            find_total_length(subduction_boundary_section_features),
            find_total_length(left_subduction_boundary_section_features),
            find_total_length(right_subduction_boundary_section_features),
            find_total_length(other_boundary_section_features))


def _write_resolved_topologies(
        time,
        output_filename_prefix,
        output_filename_extension,
        resolved_topology_features,
        ridge_transform_boundary_section_features,
        ridge_boundary_section_features,
        transform_boundary_section_features,
        subduction_boundary_section_features,
        left_subduction_boundary_section_features,
        right_subduction_boundary_section_features,
        other_boundary_section_features):

    if resolved_topology_features:
        # Put the features in a feature collection so we can write them to a file.
        resolved_topology_feature_collection = pygplates.FeatureCollection(resolved_topology_features)
        resolved_topology_features_filename = '{0}boundary_polygons_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, time, output_filename_extension)
        resolved_topology_feature_collection.write(resolved_topology_features_filename)
        
    if ridge_transform_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        ridge_transform_boundary_section_feature_collection = pygplates.FeatureCollection(ridge_transform_boundary_section_features)
        ridge_transform_boundary_section_features_filename = '{0}ridge_transform_boundaries_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, time, output_filename_extension)
        ridge_transform_boundary_section_feature_collection.write(ridge_transform_boundary_section_features_filename)
        
    if ridge_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        ridge_boundary_section_feature_collection = pygplates.FeatureCollection(ridge_boundary_section_features)
        ridge_boundary_section_features_filename = '{0}ridge_boundaries_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, time, output_filename_extension)
        ridge_boundary_section_feature_collection.write(ridge_boundary_section_features_filename)
        
    if transform_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        transform_boundary_section_feature_collection = pygplates.FeatureCollection(transform_boundary_section_features)
        transform_boundary_section_features_filename = '{0}transform_boundaries_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, time, output_filename_extension)
        transform_boundary_section_feature_collection.write(transform_boundary_section_features_filename)
        
    if subduction_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        subduction_boundary_section_feature_collection = pygplates.FeatureCollection(subduction_boundary_section_features)
        subduction_boundary_section_features_filename = '{0}subduction_boundaries_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, time, output_filename_extension)
        subduction_boundary_section_feature_collection.write(subduction_boundary_section_features_filename)
        
    if left_subduction_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        left_subduction_boundary_section_feature_collection = pygplates.FeatureCollection(left_subduction_boundary_section_features)
        left_subduction_boundary_section_features_filename = '{0}subduction_boundaries_sL_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, time, output_filename_extension)
        left_subduction_boundary_section_feature_collection.write(left_subduction_boundary_section_features_filename)
        
    if right_subduction_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        right_subduction_boundary_section_feature_collection = pygplates.FeatureCollection(right_subduction_boundary_section_features)
        right_subduction_boundary_section_features_filename = '{0}subduction_boundaries_sR_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, time, output_filename_extension)
        right_subduction_boundary_section_feature_collection.write(right_subduction_boundary_section_features_filename)

    if other_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        other_boundary_section_feature_collection = pygplates.FeatureCollection(other_boundary_section_features)
        other_boundary_section_features_filename = '{0}other_boundaries_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, time, output_filename_extension)
        other_boundary_section_feature_collection.write(other_boundary_section_features_filename)


def _parse_rotation_args(rotation_features_or_model, anchor_plate_id=None):
    if isinstance(rotation_features_or_model, pygplates.RotationModel):
        if anchor_plate_id is None:
            # Keep the original anchor plate
            return rotation_features_or_model
        if not USING_PYGPLATES_VERSION_GREATER_EQUAL_29:
            # Cannot change anchor plate
            raise RuntimeError(
                "Using pygplates version {} ".format(
                    pygplates.Version.get_imported_version()
                )
                + "but version 0.29 or greater is required "
                + "to specify a new anchor plate ID"
            )
        # Return RotationModel with new anchor plate
        return pygplates.RotationModel(
            rotation_features_or_model,
            default_anchor_plate_id=anchor_plate_id,
        )

    # Input is not a RotationModel
    if anchor_plate_id is not None:
        if anchor_plate_id == 0:
            # Zero is the default for RotationModel, which means it works
            # even in versions older than 0.26
            return pygplates.RotationModel(rotation_features_or_model)
        if not USING_PYGPLATES_VERSION_GREATER_EQUAL_26:
            raise RuntimeError(
                "Using pygplates version {} ".format(
                    pygplates.Version.get_imported_version()
                )
                + "but version 0.26 or greater is required "
                + "to specify a default anchor plate ID"
            )
        return pygplates.RotationModel(
            rotation_features_or_model,
            default_anchor_plate_id=anchor_plate_id,
        )
    # Use default anchor plate ID
    return pygplates.RotationModel(rotation_features_or_model)


if __name__ == "__main__":

    # Check the imported pygplates version.
    required_version = pygplates.Version(9)
    if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < required_version:
        print('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
                os.path.basename(__file__), pygplates.Version.get_imported_version(), required_version),
            file=sys.stderr)
        sys.exit(1)


    __description__ = \
    """Resolve topological plate polygons (and deforming networks) and saves (to separate files) the resolved topologies, and their
    boundary sections as subduction zones, mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges).

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -r rotations1.rot rotations2.rot -m topologies1.gpml topologies2.gpml -t 10 -- topology_"""

    # The command-line parser.
    parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-r', '--rotation_filenames', type=str, nargs='+', required=True,
            metavar='rotation_filename', help='One or more rotation files.')
    parser.add_argument('-m', '--topology_filenames', type=str, nargs='+', required=True,
            metavar='topology_filename', help='One or more topology files.')
    parser.add_argument('-a', '--anchor', type=int, default=0,
            dest='anchor_plate_id',
            help='Anchor plate id used for reconstructing. Defaults to zero.')
    
    # User must specify a sequence of reconstruction times or a time range (but not both).
    reconstruction_times_argument_group = parser.add_mutually_exclusive_group(required=True)
    reconstruction_times_argument_group.add_argument('-t', '--reconstruction_times', type=int, nargs='+',
            metavar='reconstruction_time',
            help='One or more times at which to reconstruct/resolve topologies.')
    reconstruction_times_argument_group.add_argument('-i', '--reconstruction_time_range', type=int, nargs=2,
            metavar=('young_time', 'old_time'),
            help='The time range (in Ma) from young time to old time over which to reconstruct/resolve topologies (in increments of 1Myr).')
    
    parser.add_argument('-d', '--transform_segment_deviation_degrees', type=float,
            default='{0}'.format(separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_DEGREES),
            help="How many degrees a mid-ocean ridge spreading segment can deviate from its stage pole before it's considered a transform segment - "
                "default is '{0}'".format(separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_DEGREES))
    
    parser.add_argument('-e', '--output_filename_extension', type=str,
            default='{0}'.format(DEFAULT_OUTPUT_FILENAME_EXTENSION),
            help="The filename extension of the output files containing the resolved topological boundaries and sections "
                "- the default extension is '{0}' - supported extensions include 'shp', 'gmt' and 'xy'."
                .format(DEFAULT_OUTPUT_FILENAME_EXTENSION))
    
    parser.add_argument(
        '-nb', '--no_output_boundaries', action='store_true',
        help='Do not write geometries of resolved topologies and their boundaries to files. '
             'This is most useful when used with the "-l" option to generate only the text file containing boundary lengths. '
             'By default geometry files are written.')
    
    parser.add_argument(
        '-l', '--output_boundary_lengths', action='store_true',
        help='Also generate a text file containing the total boundary lengths (in kms) of the resolved topology boundary sections '
             'for the requested reconstruction times. By default no text file is generated.')
    
    parser.add_argument('output_filename_prefix', type=str, nargs='?',
            default='{0}'.format(DEFAULT_OUTPUT_FILENAME_PREFIX),
            help="The prefix of the output files containing the resolved topological boundaries and sections "
                "- the default prefix is '{0}'".format(DEFAULT_OUTPUT_FILENAME_PREFIX))
    
    
    # Parse command-line options.
    args = parser.parse_args()

    if args.reconstruction_times:
        reconstruction_times = args.reconstruction_times
    else:
        if args.reconstruction_time_range[0] > args.reconstruction_time_range[1]:
            raise argparse.ArgumentTypeError("First (young) value in reconstruction time range is greater than second (old) value")
        reconstruction_times = list(range(args.reconstruction_time_range[0], args.reconstruction_time_range[1] + 1))

    # Create a rotation model from rotation files.
    #    
    # Note: We only special case 'anchor_plate_id != 0' since that is the uncommon and the 'default_anchor_plate_id'
    #       argument of 'pygplates.RotationModel.__init__()' was only added in pyGPlates revision 26.
    #       If 'anchor_plate_id == 0' then any pyGPlates revision will work.
    if args.anchor_plate_id != 0:
        if not USING_PYGPLATES_VERSION_GREATER_EQUAL_26:
            raise RuntimeError('Using pygplates version {0} but version 0.26 or greater is required for non-zero anchor plate IDs'.format(
                    pygplates.Version.get_imported_version()))
        rotation_model = pygplates.RotationModel(args.rotation_filenames, default_anchor_plate_id=args.anchor_plate_id)
    else:
        rotation_model = pygplates.RotationModel(args.rotation_filenames)
    
    # Read topology features from topology files.
    topology_features = []
    for topology_filename in args.topology_filenames:
        for topology_feature in pygplates.FeatureCollection(topology_filename):
            # FIXME: Temporary fix to avoid getting OGR GMT/Shapefile error "Mismatch in field names..." and
            # missing geometries when saving resolved topologies/sections to GMT/Shapefile.
            # It's caused by the OGR writer inside pyglates trying to write out features with different
            # shapefiles attribute field (key) names to the same file. We get around this by removing
            # all shapefile attributes.
            topology_feature.remove(pygplates.PropertyName.gpml_shapefile_attributes)
            topology_features.append(topology_feature)
    
    if args.output_boundary_lengths:
        # Create the text file for writing.
        boundary_lengths_filename = '{0}boundary_lengths.txt'.format(args.output_filename_prefix)
        boundary_lengths_file = open(boundary_lengths_filename, mode='w')
        # Write the header line (showing what each column represents).
        boundary_lengths_file.write(
            '# time(Ma)'
            ' ridge-transform-boundary-sections(km)'
            ' ridge-boundary-sections(km)'
            ' transform-boundary-sections(km)'
            ' subduction-boundary-sections(km)'
            ' left-subduction-boundary-sections(km)'
            ' right-subduction-boundary-sections(km)'
            ' other-boundary-sections(km)\n')

    
    for reconstruction_time in reconstruction_times:
   
        # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
        # We generate both the resolved topology boundaries and the boundary sections between them.
        (resolved_topology_features,
         ridge_transform_boundary_section_features,
         ridge_boundary_section_features,
         transform_boundary_section_features,
         subduction_boundary_section_features,
         left_subduction_boundary_section_features,
         right_subduction_boundary_section_features,
         other_boundary_section_features) = resolve_topologies_into_features(
                rotation_model,
                topology_features,
                reconstruction_time,
                math.radians(args.transform_segment_deviation_degrees))
        
        # Write each list of features to a separate file.
        if not args.no_output_boundaries:
            _write_resolved_topologies(
                reconstruction_time,
                args.output_filename_prefix,
                args.output_filename_extension,
                resolved_topology_features,
                ridge_transform_boundary_section_features,
                ridge_boundary_section_features,
                transform_boundary_section_features,
                subduction_boundary_section_features,
                left_subduction_boundary_section_features,
                right_subduction_boundary_section_features,
                other_boundary_section_features)
        
        # Write the total boundary lengths to a text file (if requested).
        if args.output_boundary_lengths:
            (total_length_kms_ridge_transform_boundary_section_features,
             total_length_kms_ridge_boundary_section_features,
             total_length_kms_transform_boundary_section_features,
             total_length_kms_subduction_boundary_section_features,
             total_length_kms_left_subduction_boundary_section_features,
             total_length_kms_right_subduction_boundary_section_features,
             total_length_kms_other_boundary_section_features) = find_total_boundary_length_in_kms(
                    ridge_transform_boundary_section_features,
                    ridge_boundary_section_features,
                    transform_boundary_section_features,
                    subduction_boundary_section_features,
                    left_subduction_boundary_section_features,
                    right_subduction_boundary_section_features,
                    other_boundary_section_features)

            boundary_lengths_file.write(
                    '{0} {1} {2} {3} {4} {5} {6} {7}\n'.format(
                        reconstruction_time,
                        total_length_kms_ridge_transform_boundary_section_features,
                        total_length_kms_ridge_boundary_section_features,
                        total_length_kms_transform_boundary_section_features,
                        total_length_kms_subduction_boundary_section_features,
                        total_length_kms_left_subduction_boundary_section_features,
                        total_length_kms_right_subduction_boundary_section_features,
                        total_length_kms_other_boundary_section_features))
            
