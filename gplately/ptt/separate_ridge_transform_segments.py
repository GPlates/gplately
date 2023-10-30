#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Copyright (C) 2017 The University of Sydney, Australia
    
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


#################################################################################################################
# Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments
# based on each segment’s alignment with the geometry’s stage pole at its time of appearance.
#
# Source code is based on:
#   http://www.gplates.org/docs/pygplates/sample-code/pygplates_split_isochron_into_ridges_and_transforms.html
#
#################################################################################################################


from __future__ import print_function
import math
import os.path
import pygplates
import sys


# How much a segment can deviate from the stage pole before it's considered a transform segment.
DEFAULT_TRANSFORM_SEGMENT_DEVIATION_DEGREES = 45   # An even 45 degrees split
DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS = math.radians(DEFAULT_TRANSFORM_SEGMENT_DEVIATION_DEGREES)


def separate_features_into_ridges_and_transforms(
        rotation_features_or_model,
        spreading_features,
        spreading_feature_types = None,
        transform_segment_deviation_in_radians = DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS):
    """
    Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments based
    on each segment’s alignment with the geometry’s stage pole at its time of appearance.
    
    rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).
    
    spreading_features: Spreading feature collection(s), or list of features, or filename(s) or any combination of those.
    
    spreading_feature_types: Only spreading features with a feature type contained in this list are considered.
                             If None then all spreading features are considered.
    
    transform_segment_deviation_in_radians: How much a segment can deviate from the stage pole before
                                            it's considered a transform segment (in radians).
    
    Returns: The separated ridge and transform features respectively, of type
             2-tuple (list of pygplates.Feature, list of pygplates.Feature).
    """
    
    # Turn rotation data into a RotationModel (if not already).
    #
    # OPTIMISATION:
    # We will be reconstructing (and reverse reconstructing) mid-ocean ridges in groups, where ridges in each
    # group have the same time-of-appearance. They need to have the same time of appearance because, in pyGPlates,
    # version 3 half-stage rotations start spreading from the time-of-appearance in 10 My intervals
    # (ie, the 10My intervals are 'begin_time', 'begin_time-10', begin_time-20', ..., reconstruction_time).
    # And because the mid-ocean ridges with the same time-of-appearance also have the same time intervals
    # they'll reuse the cached reconstruction trees (one cached tree per time interval).
    # This will avoid a lot of wasted time recreating these trees if the cache is continually flushed
    # (eg, by mixing mid-ocean ridges with different appearance times).
    #
    # A cache size of 100 is enough to go back to 1,000Ma (100 entries * 10My interval).
    rotation_model = pygplates.RotationModel(rotation_features_or_model, 100)
    
    # Turn spreading feature data into a list of features (if not already).
    spreading_features = pygplates.FeaturesFunctionArgument(spreading_features)
    
    # Gather all spreading features with the same begin time (time-of-appearance) into groups.
    #
    # This is an optimisation that enables reconstructing multiple mid-ocean ridges with the same begin time together.
    # This can make a *big* difference to the running time (see note above regarding rotation model cache size).
    spreading_features_grouped_by_begin_time = {}
    
    # Iterate over all geometries in spreading features.
    for spreading_feature in spreading_features.get_features():
        
        # Filter spreading feature types if requested.
        if (spreading_feature_types and
            spreading_feature.get_feature_type() not in spreading_feature_types):
            continue
        
        begin_time, _ = spreading_feature.get_valid_time()
        
        # Add to list of spreading features with same begin time.
        if begin_time not in spreading_features_grouped_by_begin_time:
            spreading_features_grouped_by_begin_time[begin_time] = []
        spreading_features_grouped_by_begin_time[begin_time].append(spreading_feature)
    
    # The separated ridge/transform segment features.
    # Both types of segment feature will have the same feature type as the feature they are extracted from.
    ridge_segment_features = []
    transform_segment_features = []
    
    # Iterate over groups of spreading features with the same begin time (time-of-appearance).
    for begin_time, spreading_features_with_begin_time in spreading_features_grouped_by_begin_time.items():
        
        # Reconstruct the spreading features to their common birth time.
        reconstructed_spreading_features = []
        pygplates.reconstruct(
                spreading_features_with_begin_time,
                rotation_model,
                reconstructed_spreading_features,
                begin_time,
                group_with_feature=True)
        
        ridge_segment_features_with_begin_time = []
        transform_segment_features_with_begin_time = []
        
        # Iterate over reconstructed spreading features.
        for spreading_feature, reconstructed_spreading_geometries in reconstructed_spreading_features:
            
            # Find the stage rotation of the spreading feature in the frame of reference of its
            # geometry at its birth time. The stage pole can then be directly geometrically compared
            # to the reconstructed spreading geometry.
            stage_rotation = get_stage_rotation_for_reconstructed_geometry(spreading_feature, rotation_model, begin_time)
            if not stage_rotation:
                # Skip current feature - it's not a spreading feature.
                continue
            
            ridge_segment_geometries = []
            transform_segment_geometries = []
            
            # A feature usually has a single geometry but it could have more - iterate over them all.
            for reconstructed_spreading_geometry in reconstructed_spreading_geometries:
                ridge_and_transform_segment_geometries = separate_geometry_into_ridges_and_transforms(
                        stage_rotation,
                        reconstructed_spreading_geometry.get_reconstructed_geometry(),
                        transform_segment_deviation_in_radians)
                if ridge_and_transform_segment_geometries:
                    ridge_segment_geometries.extend(ridge_and_transform_segment_geometries[0])
                    transform_segment_geometries.extend(ridge_and_transform_segment_geometries[1])
            
            # Put all ridge segment geometries into one feature and transform segment geometries into another.
            if ridge_segment_geometries:
                ridge_segment_feature = spreading_feature.clone()
                ridge_segment_feature.set_geometry(ridge_segment_geometries)
                ridge_segment_features_with_begin_time.append(ridge_segment_feature)
            if transform_segment_geometries:
                transform_segment_feature = spreading_feature.clone()
                transform_segment_feature.set_geometry(transform_segment_geometries)
                transform_segment_features_with_begin_time.append(transform_segment_feature)
        
        # Reverse reconstruct the segmented spreading features from their common birth time.
        #
        # Each new feature needs to be reverse reconstructed from birth time to present day because the
        # geometries are reconstructed (but need to be stored in present day positions within the feature).
        pygplates.reverse_reconstruct(
                (ridge_segment_features_with_begin_time, transform_segment_features_with_begin_time),
                rotation_model,
                begin_time)
        
        ridge_segment_features.extend(ridge_segment_features_with_begin_time)
        transform_segment_features.extend(transform_segment_features_with_begin_time)
    
    return ridge_segment_features, transform_segment_features


def separate_geometry_into_ridges_and_transforms(
        stage_rotation,
        geometry_at_spreading_time,
        transform_segment_deviation_in_radians = DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS):
    """
    Split the geometry of an isochron or mid-ocean ridge (at a time when there is spreading) into ridge and
    transform segments based on each segment’s alignment with the geometry’s stage pole at its time of appearance.
    
    For isochrons the geometry should be at its time of appearance (ie, when formed at mid-ocean ridge).
    For mid-ocean ridges the geometry can be any time when the ridge is actively spreading.
    
    stage_rotation: The stage rotation that can be applied to the geometry at the spreading time.
                    NOTE: It must have already had transforms to and from the stage pole reference frame applied.
                    In other words, if you get the stage pole from it, using 'get_euler_pole_and_angle()', then it
                    should be the stage pole in the frame of reference of the geometry at the spreading time.
    
    geometry_at_spreading_time: The polyline (or polygon) at the spreading time.
    
    transform_segment_deviation_in_radians: How much a segment can deviate from the stage pole before
                                            it's considered a transform segment (in radians).
    
    Returns: The separated ridge and transform geometries respectively, of type
             2-tuple (list of pygplates.Polyline, list of pygplates.Polyline).
             Returns None if 'geometry_at_spreading_time' is not a polyline (or polygon).
    """
    
    # Iterate over the segments of the geometry.
    # Note that we're assuming the geometry is a polyline (or polygon) - otherwise we ignore the geometry.
    try:
        segments = geometry_at_spreading_time.get_segments()
    except AttributeError:
        return
    
    # Get the stage pole of the stage rotation.
    # Note that the stage rotation is already in frame of reference of the geometry at the spreading time.
    stage_pole, _ = stage_rotation.get_euler_pole_and_angle()
    
    ridge_segment_geometries = []
    transform_segment_geometries = []
    
    # Points of contiguous segments belonging either to a ridge or transform.
    contiguous_segment_points = []
    is_transform = False
    
    for segment in segments:
        # Ignore zero length segments - they don't have a direction.
        if segment.is_zero_length():
            continue
        
        # Get the point in the middle of the segment and its tangential direction.
        segment_midpoint = segment.get_arc_point(0.5)
        segment_direction_at_midpoint = segment.get_arc_direction(0.5)
        
        # Get the direction from the segment midpoint to the stage pole.
        # This is the tangential direction at the start of an arc from the segment
        # midpoint to the stage pole (the zero parameter indicates the arc start point
        # which is the segment midpoint).
        segment_to_stage_pole_direction = pygplates.GreatCircleArc(
                segment_midpoint, stage_pole).get_arc_direction(0)
        
        # The angle that the segment deviates from the stage pole direction.
        deviation_of_segment_direction_from_stage_pole = pygplates.Vector3D.angle_between(
                segment_direction_at_midpoint, segment_to_stage_pole_direction)
        
        # When comparing the deviation angle we need to consider the case where the two
        # direction vectors are aligned but pointing in opposite directions.
        if (deviation_of_segment_direction_from_stage_pole < transform_segment_deviation_in_radians or
            deviation_of_segment_direction_from_stage_pole > math.pi - transform_segment_deviation_in_radians):
            
            # If switching from transform to ridge.
            if is_transform:
                # Emit transform polyline (if any) and restart contiguous segments.
                if contiguous_segment_points:
                    transform_segment_geometries.append(pygplates.PolylineOnSphere(contiguous_segment_points))
                    contiguous_segment_points = []
            
            is_transform = False
            
        else: # transform
            
            # If switching from ridge to transform.
            if not is_transform:
                # Emit ridge polyline (if any) and restart contiguous segments.
                if contiguous_segment_points:
                    ridge_segment_geometries.append(pygplates.PolylineOnSphere(contiguous_segment_points))
                    contiguous_segment_points = []
            
            is_transform = True
        
        # Add segment start point if first segment.
        if not contiguous_segment_points:
            contiguous_segment_points.append(segment.get_start_point())
        
        # Add segment end point.
        contiguous_segment_points.append(segment.get_end_point())
    
    # Emit last ridge or transform polyline (if any).
    if contiguous_segment_points:
        if is_transform:
            transform_segment_geometries.append(pygplates.PolylineOnSphere(contiguous_segment_points))
        else:
            ridge_segment_geometries.append(pygplates.PolylineOnSphere(contiguous_segment_points))
    
    return ridge_segment_geometries, transform_segment_geometries


def get_stage_rotation_for_reconstructed_geometry(
        spreading_feature,
        rotation_model,
        spreading_time = None):
    """
    Find the stage rotation of the spreading feature in the frame of reference of its geometry at the spreading time.
    The stage pole can then be directly geometrically compared to the reconstructed spreading geometry.
    
    spreading_feature: Can be a feature with half-stage rotation reconstruction (using left/right plate IDs)
                       or a regular feature with a conjugate plate ID.
                       An example of the former is a mid-ocean ridge, and of the latter an isochron.
    
    rotation_model: Rotation model of type pygplates.RotationModel.
    
    spreading_time: A time at which spreading is happening.
                    For isochrons this should be its time of appearance (ie, when formed at mid-ocean ridge).
                    For mid-ocean ridges this can be any time when the ridge is actively spreading.
                    Defaults to the time of appearance of 'spreading_feature'.
    
    Returns: The stage rotation that can be applied to the geometry at the spreading time.
             NOTE: It has already had transforms to and from the stage pole reference frame applied.
             So if you get the stage pole from it, using 'get_euler_pole_and_angle()', then it
             will be the stage pole in the frame of reference of the geometry at the spreading time.
             
             Returns None if 'spreading_feature' does not satisfy requirements of a spreading feature.
             (ie, have left/right plate IDs or reconstruction/conjugate plate IDs, and
             have spreading time not in distant past or future, and
             have non-zero stage rotation from 'spreading_time + 1' to 'spreading_time').
    """
    
    # If the spreading time is not specified then default to the feature's time of appearance.
    if spreading_time is None:
        spreading_time, _ = spreading_feature.get_valid_time()
    
    # Spreading time must not be distant past or future.
    if (pygplates.GeoTimeInstant(spreading_time).is_distant_past() or
        pygplates.GeoTimeInstant(spreading_time).is_distant_future()):
        return
    
    # Reconstructing either by plate ID or by half stage rotation.
    if spreading_feature.get_reconstruction_method() == 'ByPlateId':
        
        # See if spreading feature has reconstruction and conjugate plate ids.
        reconstruction_and_conjugate_plate_ids = _get_reconstruction_and_conjugate_plate_ids(spreading_feature)
        if not reconstruction_and_conjugate_plate_ids:
            # Spreading feature has no reconstruction/conjugate plate pair.
            return
        
        reconstruction_plate_id, conjugate_plate_id = reconstruction_and_conjugate_plate_ids
        
        #
        # In order to compare spreading geometries with the pole of the stage rotation at spreading time
        # we need to transform either (1) present day spreading geometries, or (2) geometries reconstructed
        # to spreading time, into the reference frame of the stage rotation (so can compare to stage pole).
        #
        # To help us decide this we start by writing the equation for a regular feature (with a conjugate plate)...
        #
        #   geometry_reconstructed = R(0->t, A->Recon) * geometry_present_day
        #                          = R(0->t, A->Conj) * R(0->t, Conj->Recon) * geometry_present_day
        #                          = R(0->t, A->Conj) * R(t+1->t, Conj->Recon) * R(0->t+1, Conj->Recon) * geometry_present_day
        #
        # ...where 'Recon' is reconstruction plate ID and 'Conj' is conjugate plate ID.
        #
        # We want to transform the spreading geometry into the stage pole reference frame at time 't=spreading_time'.
        # The easiest way to do this is to transform 'geometry_reconstructed' instead of 'geometry_present_day' since
        # it's easier to get into the reference frame of the 'R(t+1->t, Conj->Recon)' rotation
        # which is the stage rotation we're interested in when 't=spreading_time'.
        # Rearranging the above equation we get...
        #
        #   geometry_present_day = inverse[R(0->spreading_time+1, Conj->Recon)]
        #                          * inverse[R(spreading_time+1->spreading_time, Conj->Recon)]
        #                          * inverse[R(0->spreading_time, A->Conj)]
        #                          * geometry_reconstructed
        #                        = inverse[R(0->spreading_time+1, Conj->Recon)]
        #                          * inverse[R(spreading_time+1->spreading_time, Conj->Recon)]
        #                          * geometry_in_stage_pole_reference_frame
        #
        #   geometry_in_stage_pole_reference_frame = inverse[R(0->spreading_time, A->Conj)] * geometry_reconstructed
        #
        # ...where 'geometry_in_stage_pole_reference_frame' is in the stage pole reference frame because it gets rotated by
        # the stage pole rotation 'inverse[R(spreading_time+1->spreading_time, Conj->Recon)]' which differs from
        # 'R(spreading_time+1->spreading_time, Conj->Recon)' only in angle (has negated angle but pole remains the same).
        #
        # So to get reconstructed spreading geometry in the stage pole reference frame we reverse rotate 'geometry_reconstructed'
        # by 'inverse[R(0->spreading_time, A->Conj)]'.
        #
        # So to apply the stage rotation to the reconstructed spreading geometry we rotate it in stage pole reference frame,
        # then apply stage rotation and then rotate back from the stage pole reference frame...
        #
        #   stage_rotate_geometry_reconstructed = R(0->spreading_time, A->Conj)
        #                                         * R(spreading_time+1->spreading_time, Conj->Recon)
        #                                         * inverse[R(0->spreading_time, A->Conj)]
        #                                         * geometry_reconstructed
        #
        # For more detail see:
        #   http://www.gplates.org/docs/pygplates/sample-code/pygplates_split_isochron_into_ridges_and_transforms.html
        #
        stage_rotation = rotation_model.get_rotation(spreading_time, reconstruction_plate_id, spreading_time + 1, conjugate_plate_id)
        if stage_rotation.represents_identity_rotation():
            return
        from_stage_pole_reference_frame = rotation_model.get_rotation(spreading_time, conjugate_plate_id)
        to_stage_pole_reference_frame = from_stage_pole_reference_frame.get_inverse()
        stage_rotation = from_stage_pole_reference_frame * stage_rotation * to_stage_pole_reference_frame
    
    else: # Reconstruction is by half stage rotation...
        
        # See if spreading feature has left and right plate ids (it should).
        left_and_right_plate_ids = _get_left_and_right_plate_ids(spreading_feature)
        if not left_and_right_plate_ids:
            # Spreading feature has no left/right plate pair.
            return
        
        left_plate_id, right_plate_id = left_and_right_plate_ids
        
        #
        # In order to compare spreading geometries with the pole of the stage rotation at birth time
        # we need to transform either (1) present day spreading geometries, or (2) geometries reconstructed
        # to birth time, into the reference frame of the stage rotation (so can compare to stage pole).
        #
        # To help us decide this we start by writing the equation for a mid-ocean ridge (MOR)...
        #
        #   geometry_reconstructed = R(0->t, A->MOR) * geometry_present_day
        #                          = R(0->t, A->Left) * R(0->t, Left->MOR) * geometry_present_day
        #                          = R(0->t, A->Left) * spread(ts->t, Left->Right) * geometry_present_day
        #
        # ...where 'MOR' is not a plate ID, which is why we do half-spreading (or asymmetric spreading) of
        # right plate relative to left plate. The function 'spread()' usually splits the time interval from
        # spreading start time 'ts' to time 't' into N stages and accumulates spreading over those N stages...
        #
        #   geometry_reconstructed = R(0->t, A->Left) * spread(ts->t, Left->Right) * geometry_present_day
        #                          = R(0->t, A->Left)
        #                            * spread(t[N-1]->t, Left->Right) * spread(t[N-2]->t[N-1], Left->Right) * ... * spread(t1->t2, Left->Right) * spread(ts->t1, Left->Right)
        #                            * geometry_present_day
        #
        # ...in GPlates the "gpml:ReconstructionMethodEnumeration" property currently supports 'HalfStageRotation' versions 1, 2 and 3.
        # They only differ in the spreading start time 'ts' and the number of stages N.
        # Version 1 has 'ts=0' and 'N=1'.
        # Version 2 has 'ts=0' and 'N>1'.
        # Version 3 has 'ts=spreading_time' and 'N>1'.
        #
        # We want to transform the spreading geometry into the stage pole reference frame at time 't=spreading_time'.
        # The easiest way to do this is to transform 'geometry_reconstructed' instead of 'geometry_present_day' since
        # it's easier to get into the reference frame of the 'spread(t[N-1]->t, Left->Right)' rotation
        # which is the stage rotation we're interested in when 't=spreading_time'.
        # Rearranging the above equation we get...
        #
        #   geometry_present_day = inverse[spread(ts->t1, Left->Right)] * ... * inverse[spread(t[N-1]->spreading_time, Left->Right)]
        #                          * inverse[R(0->t, A->Left)]
        #                          * geometry_reconstructed
        #                        = inverse[spread(ts->t1, Left->Right)] * ... * inverse[spread(t[N-1]->spreading_time, Left->Right)]
        #                          * geometry_in_stage_pole_reference_frame
        #
        #   geometry_in_stage_pole_reference_frame = inverse[R(0->spreading_time, A->Left)] * geometry_reconstructed
        #
        # ...where 'geometry_in_stage_pole_reference_frame' is in the stage pole reference frame because it gets rotated by
        # the stage pole rotation 'inverse[spread(t[N-1]->spreading_time, Left->Right)]' which differs from
        # 'spread(t[N-1]->spreading_time, Left->Right)' only in angle (has negated angle but pole remains the same).
        #
        # So to get reconstructed spreading geometry in the stage pole reference frame we reverse rotate 'geometry_reconstructed' by 'inverse[R(0->spreading_time, A->Left)]'.
        #
        # So to apply the stage rotation to the reconstructed spreading geometry we rotate it in stage pole reference frame,
        # then apply stage rotation and then rotate back from the stage pole reference frame...
        #
        #   stage_rotate_geometry_reconstructed = R(0->spreading_time, A->Left)
        #                                         * R(spreading_time+1->spreading_time, Left->Right)
        #                                         * inverse[R(0->spreading_time, A->Left)]
        #                                         * geometry_reconstructed
        #
        stage_rotation = rotation_model.get_rotation(spreading_time, right_plate_id, spreading_time + 1, left_plate_id)
        if stage_rotation.represents_identity_rotation():
            return
        from_stage_pole_reference_frame = rotation_model.get_rotation(spreading_time, left_plate_id)
        to_stage_pole_reference_frame = from_stage_pole_reference_frame.get_inverse()
        stage_rotation = from_stage_pole_reference_frame * stage_rotation * to_stage_pole_reference_frame
    
    return stage_rotation

#
# Private function.
#
# Returns a tuple of left and right plate ids from a 'feature', or None if not found.
#
def _get_left_and_right_plate_ids(feature):
    # Get left and right plate ids (if any).
    left_plate_id = feature.get_left_plate(None)
    right_plate_id = feature.get_right_plate(None)
    if left_plate_id is not None and right_plate_id is not None:
        return left_plate_id, right_plate_id


#
# Private function.
#
# Returns a tuple of reconstruction and conjugate plate ids from a 'feature',
# otherwise looks for plate/conjugate ids if feature came from a PLATES data file,
# otherwise returns None.
#
def _get_reconstruction_and_conjugate_plate_ids(feature):
    reconstruction_plate_id = feature.get_reconstruction_plate_id(None)
    conjugate_plate_id = feature.get_conjugate_plate_id(None)
    # If missing either then attempt to get reconstruction/conjugate from the 'gpml:OldPlatesHeader' property.
    if reconstruction_plate_id is not None and conjugate_plate_id is not None:
        return reconstruction_plate_id, conjugate_plate_id
    
    gpml_old_plates_header = feature.get_value(pygplates.PropertyName.create_gpml('oldPlatesHeader'))
    if gpml_old_plates_header:
        try:
            reconstruction_plate_id = gpml_old_plates_header.get_plate_id_number()
            conjugate_plate_id = gpml_old_plates_header.get_conjugate_plate_id_number()
            return reconstruction_plate_id, conjugate_plate_id
        except AttributeError:
            # The property value type did not match the property name.
            # This indicates the data does not conform to the GPlates Geological Information Model (GPGIM).
            pass




if __name__ == '__main__':
    
    import argparse
    
    
    DEFAULT_OUTPUT_RIDGES_FILENAME_SUFFIX = '_ridges'
    DEFAULT_OUTPUT_TRANSFORMS_FILENAME_SUFFIX = '_transforms'

    __description__ = \
    """Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments.
    
    The splitting is based on each segment's alignment with the geometry's stage pole at its time of appearance.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -r rotations.rot -d 45 -s _ridges -t _transforms -- spreading_features.gpml
    """
    
    # The command-line parser.
    parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-r', '--rotation_filenames', type=str, nargs='+', required=True,
            metavar='rotation_filename', help='One or more rotation files.')
    
    parser.add_argument('-s', '--output_ridges_filename_suffix', type=str,
            default='{0}'.format(DEFAULT_OUTPUT_RIDGES_FILENAME_SUFFIX),
            help="The suffix to append to each input filename to get each output ridges filename - "
                "the default suffix is '{0}'".format(DEFAULT_OUTPUT_RIDGES_FILENAME_SUFFIX))
    parser.add_argument('-t', '--output_transforms_filename_suffix', type=str,
            default='{0}'.format(DEFAULT_OUTPUT_TRANSFORMS_FILENAME_SUFFIX),
            help="The suffix to append to each input filename to get each output transforms filename - "
                "the default suffix is '{0}'".format(DEFAULT_OUTPUT_TRANSFORMS_FILENAME_SUFFIX))
    
    parser.add_argument('-d', '--transform_segment_deviation_degrees', type=float,
            default='{0}'.format(DEFAULT_TRANSFORM_SEGMENT_DEVIATION_DEGREES),
            help="How many degrees a spreading segment can deviate from the stage pole before it's considered a transform segment - "
                "default is '{0}'".format(DEFAULT_TRANSFORM_SEGMENT_DEVIATION_DEGREES))
    
    parser.add_argument('-f', '--spreading_feature_types', type=str, nargs='+',
            metavar='spreading_feature_type',
            help='The feature type(s) to split into ridge/transform segments. '
                'All other feature types will be ignored (and not end up in separated ridge/transform output files). '
                'The format should match the format of '
                'http://www.gplates.org/docs/pygplates/generated/pygplates.FeatureType.html#pygplates.FeatureType.get_name . '
                'For example, mid-ocean ridges are specified as MidOceanRidge (without the gpml: prefix). '
                'Defaults to splitting all features (although features that are not spreading are ignored).')
    
    parser.add_argument('input_filenames', type=str, nargs='+',
            metavar='input_filename',
            help='One or more input filenames (original files).')
    
    # Parse command-line options.
    args = parser.parse_args()
    
    # Convert strings into feature types.
    # For example, 'MidOceanRidge' into pygplates.FeatureType.create_gpml('MidOceanRidge')
    if args.spreading_feature_types:
        args.spreading_feature_types = [pygplates.FeatureType.create_gpml(feature_type)
                for feature_type in args.spreading_feature_types]
    
    for input_filename in args.input_filenames:
        ridge_features, transform_features = separate_features_into_ridges_and_transforms(
                args.rotation_filenames,
                args.input_filenames,
                args.spreading_feature_types,
                math.radians(args.transform_segment_deviation_degrees))
        
        # Each output filename is the input filename with a suffix appended.
        filename_root, filename_ext = os.path.splitext(input_filename)
        output_ridges_filename = ''.join((filename_root, args.output_ridges_filename_suffix, filename_ext))
        pygplates.FeatureCollection(ridge_features).write(output_ridges_filename)
        
        # Each output filename is the input filename with a suffix appended.
        filename_root, filename_ext = os.path.splitext(input_filename)
        output_transforms_filename = ''.join((filename_root, args.output_transforms_filename_suffix, filename_ext))
        pygplates.FeatureCollection(transform_features).write(output_transforms_filename)
