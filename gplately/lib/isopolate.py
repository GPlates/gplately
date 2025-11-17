
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


##############################################################################
# You can either:
#  (1) Import this module ('import isopolate') into your script and call the
#      'interpolate_isochrons()' function in your script, or
#  (2) Run this script with command-line options:
#        python isopolate.py ...
#
# For (1) see the docstring of interpolate_isochrons() for more details.
# For (2) run this script as 'python isopolate.py -h' to see more options.
##############################################################################


from __future__ import print_function
import argparse
import math
import sys
import os.path
import pygplates


# Required pygplates version.
PYGPLATES_VERSION_REQUIRED = pygplates.Version(8)


# The default interval spacing between interpolated isochrons.
DEFAULT_INTERPOLATE_RESOLUTION_DEGREES = 0.1
DEFAULT_INTERPOLATE_RESOLUTION_RADIANS = math.radians(DEFAULT_INTERPOLATE_RESOLUTION_DEGREES)

# The default value for the minimum latitude overlap when detecting if
# we can interpolate between two polylines.
DEFAULT_MINIMUM_LATITUDE_OVERLAP_DEGREES = 1
DEFAULT_MINIMUM_LATITUDE_OVERLAP_RADIANS = math.radians(DEFAULT_MINIMUM_LATITUDE_OVERLAP_DEGREES)

# Default values for the extra range of non-overlapping latitudes at the North and South
# (in stage rotation pole reference frame).
# For isochrons...
DEFAULT_MAXIMUM_ISOCHRON_LATITUDE_NON_OVERLAP_DEGREES = 5
DEFAULT_MAXIMUM_ISOCHRON_LATITUDE_NON_OVERLAP_RADIANS = math.radians(DEFAULT_MAXIMUM_ISOCHRON_LATITUDE_NON_OVERLAP_DEGREES)
# For COBs...
DEFAULT_MAXIMUM_COB_LATITUDE_NON_OVERLAP_DEGREES = 1
DEFAULT_MAXIMUM_COB_LATITUDE_NON_OVERLAP_RADIANS = math.radians(DEFAULT_MAXIMUM_COB_LATITUDE_NON_OVERLAP_DEGREES)

# The maximum distance between any corresponding pair of points (same latitude), in degrees,
# above which two polylines will not be interpolated.
DEFAULT_MAXIMUM_DISTANCE_THRESHOLD_DEGREES = 30
DEFAULT_MAXIMUM_DISTANCE_THRESHOLD_RADIANS = math.radians(DEFAULT_MAXIMUM_DISTANCE_THRESHOLD_DEGREES)

# The default value for the maximum difference in age between two features to be interpolated.
# If exceeded then no new features are interpolated between them.
DEFAULT_MAXIMUM_AGE_DIFFERENCE = 50

# The default value to join polylines if their end points are within a threshold distance of each other.
DEFAULT_JOIN_POLYLINES_THRESHOLD_DEGREES = 3.0
DEFAULT_JOIN_POLYLINES_THRESHOLD_RADIANS = math.radians(DEFAULT_JOIN_POLYLINES_THRESHOLD_DEGREES)


# Scalar coverages.
SCALAR_COVERAGE_SPREADING_ASYMMETRY = 'SpreadingAsymmetry'
SCALAR_COVERAGE_SPREADING_RATE = 'SpreadingRate'
SCALAR_COVERAGE_FULL_SPREADING_RATE = 'FullSpreadingRate'
SCALAR_COVERAGE_SPREADING_DIRECTION = 'SpreadingDirection'
SCALAR_COVERAGE_SPREADING_OBLIQUITY = 'SpreadingObliquity'
SCALAR_COVERAGE_AGE = 'Age'

# Scalar types.
SCALAR_TYPE_SPREADING_ASYMMETRY = pygplates.ScalarType.create_gpml(SCALAR_COVERAGE_SPREADING_ASYMMETRY)
SCALAR_TYPE_SPREADING_RATE = pygplates.ScalarType.create_gpml(SCALAR_COVERAGE_SPREADING_RATE)
SCALAR_TYPE_FULL_SPREADING_RATE = pygplates.ScalarType.create_gpml(SCALAR_COVERAGE_FULL_SPREADING_RATE)
SCALAR_TYPE_SPREADING_DIRECTION = pygplates.ScalarType.create_gpml(SCALAR_COVERAGE_SPREADING_DIRECTION)
SCALAR_TYPE_SPREADING_OBLIQUITY = pygplates.ScalarType.create_gpml(SCALAR_COVERAGE_SPREADING_OBLIQUITY)
SCALAR_TYPE_AGE = pygplates.ScalarType.create_gpml(SCALAR_COVERAGE_AGE)

# Known scalar coverages.
KNOWN_SCALAR_COVERAGES = set([
        SCALAR_COVERAGE_SPREADING_ASYMMETRY,
        SCALAR_COVERAGE_SPREADING_RATE,
        SCALAR_COVERAGE_FULL_SPREADING_RATE,
        SCALAR_COVERAGE_SPREADING_DIRECTION,
        SCALAR_COVERAGE_SPREADING_OBLIQUITY,
        SCALAR_COVERAGE_AGE])


COB_FEATURE_TYPE = pygplates.FeatureType.create_gpml('PassiveContinentalBoundary')
ISOCHRON_FEATURE_TYPE = pygplates.FeatureType.create_gpml('Isochron')
MID_OCEAN_RIDGE_FEATURE_TYPE = pygplates.FeatureType.create_gpml('MidOceanRidge')


#
# Private function.
#
# Returns a tuple of left and right plate ids from a 'feature', or None if not found.
# If feature does not have left/right plate ids then looks for reconstruction/conjugate plate ids,
# otherwise looks for plate/conjugate ids if feature came from a PLATES data file, otherwise returns None.
#
def _get_left_right_plate_ids(feature):
    # Get left and right plate ids (if any).
    left_plate_id = feature.get_left_plate(None)
    right_plate_id = feature.get_right_plate(None)
    # If missing either then attempt to get reconstruction/conjugate plate id.
    if left_plate_id is not None and right_plate_id is not None:
        return left_plate_id, right_plate_id
    
    left_plate_id = feature.get_reconstruction_plate_id(None)
    right_plate_id = feature.get_conjugate_plate_id(None)
    # If missing either then attempt to get reconstruction/conjugate from the 'gpml:OldPlatesHeader' property.
    if left_plate_id is not None and right_plate_id is not None:
        return left_plate_id, right_plate_id
    
    gpml_old_plates_header = feature.get_value(pygplates.PropertyName.create_gpml('oldPlatesHeader'))
    
    if gpml_old_plates_header:
        try:
            left_plate_id = gpml_old_plates_header.get_plate_id_number()
            right_plate_id = gpml_old_plates_header.get_conjugate_plate_id_number()
            return left_plate_id, right_plate_id
        except AttributeError:
            # The property value type did not match the property name.
            # This indicates the data does not conform to the GPlates Geological Information Model (GPGIM).
            pass


#
# Private function.
#
# Remove from 'features' those features that do not have conjugate plate ids in the 'plate_pairs' list.
# 
def _filter_plate_pairs(features, plate_pairs):
    feature_index = 0
    while feature_index < len(features):
        feature = features[feature_index]
        
        # Get left and right plate ids (in one form or another), otherwise remove the current feature.
        left_right_plate_ids = _get_left_right_plate_ids(feature)
        if left_right_plate_ids:
            left_plate_id, right_plate_id = left_right_plate_ids
            
            if (left_plate_id, right_plate_id) in plate_pairs:
                # Keep current feature and move to next.
                feature_index += 1
                continue
            
            if (right_plate_id, left_plate_id) in plate_pairs:
                # Keep current feature and move to next.
                feature_index += 1
                continue
        
        # Remove current feature.
        del features[feature_index]


#
# Private function.
#
# Removes duplicate geometries from the list 'features'.
# This can reduce the number of features in the list.
# 
def _remove_features_with_duplicate_geometries(features):
    if len(features) == 1:
        return
    
    #
    # Do N^2 search over pairs of features to test for duplicates.
    #
    
    # Iterate over all features except last feature.
    # Using len() since some features are removed during iteration.
    feature1_index = 0
    while feature1_index < len(features) - 1:
        feature1 = features[feature1_index]

        # Get all geometries in 'feature1'.
        feature1_geoms = feature1.get_all_geometries()
        
        # Iterate over the remaining features (after 'feature1').
        # Using len() since some features are removed during iteration.
        feature2_index = feature1_index + 1
        while feature2_index < len(features):
            feature2 = features[feature2_index]
            
            # Get all geometries in 'feature2'.
            feature2_geoms = feature2.get_all_geometries()
            
            # Compare the geometries of feature1 and feature2.
            removed = False
            for feature1_geom_index, feature1_geom in enumerate(feature1_geoms):
                feature2_geom_index = 0
                # Using len() since some geometries are removed during iteration.
                while feature2_geom_index < len(feature2_geoms):
                    feature2_geom = feature2_geoms[feature2_geom_index]
                    
                    # Test for duplicate geometries.
                    if feature1_geom == feature2_geom:
                        # Remove the feature 2 geometry (from the list).
                        del feature2_geoms[feature2_geom_index]
                        feature2_geom_index -= 1
                        removed = True
                    
                    feature2_geom_index += 1
            
            if removed:
                # Replace 'feature2's geometries with the reduced set of geometries, or
                # remove 'feature2' from the list if all its geometries have been removed.
                if feature2_geoms:
                    feature2.set_geometry(feature2_geoms)
                else:
                    del features[feature2_index]
                    feature2_index -= 1
            
            feature2_index += 1
        
        feature1_index += 1


#
# Private function.
#
# Joins adjacent geometries in the 'features' list and returns a new joined list.
# 
def _join_adjacent_features(features, distance_threshold_radians, print_debug_output):
    geometries = []
    begin_time = None
    end_time = float('inf')
    
    # Extract all the geometries from the features and expand time range to include all features.
    for feature in features:
        feature_begin_time, feature_end_time = feature.get_valid_time()
        
        # Some features erroneously have end time earlier than begin time.
        # Ignore these features.
        if feature_end_time > feature_begin_time:
            if print_debug_output >= 1:
                print(" ...skipping '{0}' "
                        "feature (id: {1}): begin time {2} later than end time {3}."
                                .format(feature.get_feature_type(), feature.get_feature_id(),
                                        feature_begin_time, feature_end_time))
            continue
        
        # Expand time range to include the current feature.
        # Note: All features should have the same begin time but can have differing end times.
        if begin_time is None:
            begin_time = feature_begin_time
        if feature_end_time < end_time:
            end_time = feature_end_time
        
        # Get all geometries in 'feature'.
        geometries.extend(feature.get_all_geometries())
                    
    # Join polylines that have end points closer than a distance threshold.
    # Print debug message if not all geometries are polylines.
    try:
        joined_geometries = pygplates.PolylineOnSphere.join(
                geometries, distance_threshold_radians, pygplates.PolylineConversion.raise_if_non_polyline)
    except pygplates.GeometryTypeError:
        if print_debug_output >= 1:
            print(' ...skipping non-polyline geometries at time {0}'.format(begin_time))
        # Try again but this time just ignore any geometry that's not a polyline.
        joined_geometries = pygplates.PolylineOnSphere.join(
                geometries, distance_threshold_radians, pygplates.PolylineConversion.ignore_non_polyline)
    
    #print('joined begin time {0}, plate ids {1}'.format(begin_time, _get_left_right_plate_ids(feature)))
    
    # Potentially reduced set of joined features.
    joined_features = []
    
    for joined_geometry in joined_geometries:
        # Clone the properties in the first feature.
        joined_feature = features[0].clone()
        # Modify the time range and geometry.
        joined_feature.set_valid_time(begin_time, end_time)
        joined_feature.set_geometry(joined_geometry)
        
        joined_features.append(joined_feature)
    
    return joined_features


#
# Private function.
#
# Calculate a spreading asymmetry for each point in the young/old polyline.
#
# Both polylines should have the same number of points and they are expected to be
# from the output of 'pygplates.PolylineOnSphere.rotation_interpolate()' which
# keeps matching points latitude aligned (in coordinate system of stage rotation),
# except for the non-overlapping latitude regions at the ends.
#
# The stage rotation pole should be in the present day reference frame (since young/old polylines are present day).
#
# Also 'stage_pole_angle_radians' must not be zero (identity rotation).
# 
def _calc_spreading_asymmetries(
        young_polyline,
        old_polyline,
        stage_rotation_pole_for_present_day_geometry,
        stage_pole_angle_radians):
    
    asymmetries = []
    
    # Store pole as a pygplates.Vector3D instead of a pygplates.PointOnSphere.
    stage_rotation_vector_for_present_day_geometry = pygplates.Vector3D(stage_rotation_pole_for_present_day_geometry.to_xyz())
    abs_stage_pole_angle_radians = abs(stage_pole_angle_radians)
    
    for point_index in range(len(young_polyline)):
        
        # The young point and the stage pole lie in a plane (also through origin of globe).
        # We find the vector perpendicular to the young point and stage pole (this will help us below).
        # 
        # The same applies to the old point (it has a plane and perpendicular vector).
        #
        # The angle between the young and old planes can then be calculated from the great circle arc
        # angle between the young and the old perpendicular vectors. That is because the rotation
        # about the stage pole vector at these perpendicular positions is a great circle (not a small circle).
        #
        # That angle represents the rotation angle about the stage pole that moves the young point
        # onto the old points plane. In most cases the young and old points are latitude aligned
        # and hence the young point will rotate on top of the old point but this won't be the
        # case for the non-overlapping latitude regions at the ends of the polylines due to
        # pygplates.PolylineOnSphere.rotation_interpolate().
        
        young_vector_perp_pole = pygplates.Vector3D.cross(
                pygplates.Vector3D.cross(stage_rotation_vector_for_present_day_geometry, young_polyline[point_index].to_xyz()),
                stage_rotation_vector_for_present_day_geometry)
        
        old_vector_perp_pole = pygplates.Vector3D.cross(
                pygplates.Vector3D.cross(stage_rotation_vector_for_present_day_geometry, old_polyline[point_index].to_xyz()),
                stage_rotation_vector_for_present_day_geometry)
        
        try:
            angle_between_young_and_old_points = pygplates.Vector3D.angle_between(
                    young_vector_perp_pole,
                    old_vector_perp_pole)
        except pygplates.UnableToNormaliseZeroVectorError:
            # Young or old point is coincident with the stage rotation pole so just set asymmetry to zero.
            asymmetries.append(0)
            continue
        
        # Asymmetry is in the range [-1,1].
        #
        # The "min(..., 1)" clamps to the range [0,1].
        # The "2 * () - 1" maps the range [0,1] to [-1,1].
        asymmetry = 2 * min(angle_between_young_and_old_points / abs_stage_pole_angle_radians, 1) - 1
        asymmetries.append(asymmetry)
    
    return asymmetries


#
# Private function.
#
# Adds requested scalar coverages (if any) to an interpolated polyline and tessellates it (if requested).
#
# Any scalar types in 'output_scalar_types' are calculated at each point in the interpolated polyline
# and the returned geometry is then a coverage (ie, a geometry and scalar values).
# The length of 'spreading_asymmetries' must match the number of points in 'interpolated_polyline'.
# If 'output_scalar_types' is empty then 'spreading_asymmetries' is ignored.
#
# If 'tessellate_threshold_radians' is None then no tessellation is performed.
#
# The stage rotation should be in the present day reference frame (since interpolated polyline is present day).
# 
def _add_tessellated_scalar_coverages(
        interpolated_polyline,
        tessellate_threshold_radians,
        output_scalar_types,
        spreading_asymmetries,
        stage_rotation_for_present_day_geometry,
        inverse_present_day_rotation,
        interpolated_time,
        stage_rotation_time_interval,
        rotation_model,
        interpolated_isochron_plate_id):
    
    # If not outputting scalar coverages then just tessellate (if requested) and return.
    if not output_scalar_types:
    
        # Tessellate if requested.
        if tessellate_threshold_radians:
            return interpolated_polyline.to_tessellated(tessellate_threshold_radians)
        
        return interpolated_polyline
    
    
    scalar_coverages = {}
                        
    # Tessellate the interpolated polyline if requested.
    if tessellate_threshold_radians:
        tessellated_points = []
        # We also need to tessellate the spreading asymmetries since its length
        # must match the number of tessellated points.
        tessellated_spreading_asymmetries = []
        
        # Iterate over all polyline segments.
        segments = interpolated_polyline.get_segments()
        for segment_index in range(len(segments)):
            segment = segments[segment_index]
            
            tessellated_segment_points = segment.to_tessellated(tessellate_threshold_radians)
            # Copy all points except the last one (it's a duplicate of the start point of next segment).
            tessellated_points.extend(tessellated_segment_points[:-1])
            
            # Add the spreading asymmetry associated with the segment's start point.
            segment_start_asymmetry = spreading_asymmetries[segment_index]
            tessellated_spreading_asymmetries.append(segment_start_asymmetry)
            # Tessellate the spreading asymmetries except the last one
            # (it's a duplicate of the asymmetry of the start point of next segment).
            num_segment_sub_intervals = len(tessellated_segment_points) - 1
            # Do we even need to tessellate ?
            if num_segment_sub_intervals > 1:
                inv_num_segment_sub_intervals = 1.0 / num_segment_sub_intervals
                segment_end_asymmetry = spreading_asymmetries[segment_index + 1]
                # Generate the asymmetries at the segment internal points (excludes segment end points).
                for interval_index in range(1, num_segment_sub_intervals):
                    interpolate_ratio = interval_index * inv_num_segment_sub_intervals
                    tessellated_spreading_asymmetries.append(
                            (1 - interpolate_ratio) * segment_start_asymmetry +
                            interpolate_ratio * segment_end_asymmetry)
        
        # Add the last point/asymmetry.
        tessellated_points.append(interpolated_polyline[-1])
        tessellated_spreading_asymmetries.append(spreading_asymmetries[-1])
        
        # Replace the interpolated polyline with a tessellated version.
        interpolated_polyline = pygplates.PolylineOnSphere(tessellated_points)
        
        # Replace the spreading asymmetries with the tessellated versions.
        spreading_asymmetries = tessellated_spreading_asymmetries
    
    # Output spreading asymmetry scalar coverage (if requested) - per-point spreading asymmetry values.
    if SCALAR_TYPE_SPREADING_ASYMMETRY in output_scalar_types:
        scalar_coverages[SCALAR_TYPE_SPREADING_ASYMMETRY] = spreading_asymmetries
    
    # Output spreading rate/direction/obliquity scalar coverage (if requested) - per-point spreading rate/direction values.
    if (SCALAR_TYPE_SPREADING_RATE in output_scalar_types or
        SCALAR_TYPE_FULL_SPREADING_RATE in output_scalar_types or
        SCALAR_TYPE_SPREADING_DIRECTION in output_scalar_types or
        SCALAR_TYPE_SPREADING_OBLIQUITY in output_scalar_types):
        
        interpolated_points = interpolated_polyline.get_points()
        
        # Calculate full spreading velocities at the points of the interpolated/tessellated polyline.
        #
        # First we calculate velocities at present day using the present day polyline point locations and the
        # stage rotation in the present day reference frame.
        full_spreading_velocity_vectors_at_present_day = pygplates.calculate_velocities(
                interpolated_points,
                stage_rotation_for_present_day_geometry,
                stage_rotation_time_interval)
        #
        # Then reconstruct the points and velocities from present day to their reconstructed positions/directions at the
        # interpolated isochron's birth time (time of crustal accretion) since we're recording a snapshot of crustal accretion.
        #
        # A small issue is there might be a non-zero finite rotation at present day, so we cannot assume the un-reconstructed
        # isochron is the same as the isochron reconstructed to present day. We're working with isochron geometry at present day
        # so we need to reverse reconstruct it to its un-reconstructed position 'inverse[R(0, A->Plate)] * geometry_present_day'
        # before we can reconstruct it to its begin time (time of crustal accretion).
        #
        #   geometry_moving_plate = R(0->begin_time, A->Plate) * geometry_present_day
        #                         = R(begin_time, A->Plate) * inverse[R(0, A->Plate)] * geometry_present_day
        #
        # ...where '0->t' means reconstructing from its "present day" (not un-reconstructed) position to time 't'.
        #
        # Note: We can't just calculate 'R(0->begin_time, A->Plate)' in one call to 'rotation_model.get_rotation(interpolated_time, interpolated_isochron_plate_id, 0)'
        # because explicitly setting the 'from_time' argument to '0' results in pyGPlates (versions < 0.27) assuming that the present day rotation is zero
        # (and so it essentially just calculates 'R(begin_time, A->Plate)' instead of 'R(begin_time, A->Plate) * inverse[R(0, A->Plate)]').
        # So to ensure the same results for pyGPlates versions before and after 0.27, we'll exclude the 'from_time' argument for 'R(begin_time, A->Plate)'
        # so versions >= 0.27 don't calculate 'R(0->begin_time, A->Plate)' and explicitly include 'inverse[R(0, A->Plate)]' so that versions < 0.27 are accounted for.
        # 
        #
        interpolated_isochron_reconstruction = rotation_model.get_rotation(interpolated_time, interpolated_isochron_plate_id) * inverse_present_day_rotation
        full_spreading_velocity_vectors = [interpolated_isochron_reconstruction * velocity
            for velocity in full_spreading_velocity_vectors_at_present_day]
        reconstructed_interpolated_points = [interpolated_isochron_reconstruction * point
            for point in interpolated_points]

        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
        full_spreading_velocities = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                reconstructed_interpolated_points,
                full_spreading_velocity_vectors)
        
        # Extract the 'magnitude' element of each velocity tuple (if requested) and adjust for asymmetry.
        if SCALAR_TYPE_SPREADING_RATE in output_scalar_types:
            # Convert asymmetry from range [-1,1] to [0,1] and multiply by the full spreading velocity.
            spreading_rates = [0.5 * (spreading_asymmetries[point_index] + 1) * full_spreading_velocities[point_index][0]
                    for point_index in range(len(full_spreading_velocities))]
            scalar_coverages[SCALAR_TYPE_SPREADING_RATE] = spreading_rates
        
        # Extract the 'magnitude' element of each velocity tuple (if requested).
        if SCALAR_TYPE_FULL_SPREADING_RATE in output_scalar_types:
            full_spreading_rates = [velocity[0] for velocity in full_spreading_velocities]
            scalar_coverages[SCALAR_TYPE_FULL_SPREADING_RATE] = full_spreading_rates
        
        # Extract the 'azimuth' element of each velocity tuple (if requested).
        if SCALAR_TYPE_SPREADING_DIRECTION in output_scalar_types:
            spreading_directions = []
            for velocity in full_spreading_velocities:
                spreading_direction = math.degrees(velocity[1])
                # Since spreading at crustal accretion is in two opposite directions,
                # we pick the direction with the smaller azimuth (an azimuth is in the range [0, 360]).
                # For example, for an azimuth of 225 and its opposite (225 - 180 = 45) we choose 45.
                if spreading_direction > 180:
                    spreading_direction -= 180
                spreading_directions.append(spreading_direction)
            scalar_coverages[SCALAR_TYPE_SPREADING_DIRECTION] = spreading_directions
        
        # Calculate deviation of spreading direction from isochron normal.
        if SCALAR_TYPE_SPREADING_OBLIQUITY in output_scalar_types:
            spreading_obliquities = []

            # First calculate a normal for each arc segment (between two adjacent points)
            # in the interpolated polyline.
            segment_normals = []
            for segment in interpolated_polyline.get_segments():
                if segment.is_zero_length():
                    segment_normal = None
                else:
                    segment_normal = segment.get_great_circle_normal()
                
                segment_normals.append(segment_normal)
            
            # For each isochron point calculate a spreading obliquity for its two adjoining arc segment normals and then
            # average the two obliquities (if both are ridges, ie, less than 45 degrees), or take the one obliquity that
            # is a ridge (if other is transform), or both are transforms so hard-wire obliquity to 90 degrees.
            num_points = len(interpolated_points)
            for point_index in range(num_points):
                # Normal of segment before current point.
                if point_index > 0:
                    prev_normal = segment_normals[point_index - 1]
                else:
                    prev_normal = None
                
                # Normal of segment after current point.
                if point_index < num_points - 1:
                    next_normal = segment_normals[point_index]
                else:
                    next_normal = None
                
                #
                # The angle between the velocity and normal vectors is in the range [0, 180].
                #
                # Obliquity is how much the small circle of spreading rotation pole (-/+ velocity_vector)
                # deviates from the perpendicular line at the isochron point (-/+ normal).
                # This is the minimum deviation of 'velocity_vector' and '-velocity_vector' from 'normal' and '-normal'.
                # Obliquity angle is in range [0, 90].
                #

                velocity_vector = full_spreading_velocity_vectors_at_present_day[point_index]

                # Calculate spreading obliquity for previous normal. If there is no previous normal or
                # obliquity is larger than 45 degrees then set to None.
                if prev_normal:
                    # Range [0,180].
                    spreading_obliquity_for_prev_normal = math.degrees(pygplates.Vector3D.angle_between(
                        velocity_vector,
                        prev_normal))
                    # Range [0,90].
                    if spreading_obliquity_for_prev_normal > 90:
                        spreading_obliquity_for_prev_normal = 180 - spreading_obliquity_for_prev_normal
                    # Exclude transform sections (obliquity larger than 45 degrees).
                    if spreading_obliquity_for_prev_normal > 45:
                        spreading_obliquity_for_prev_normal = None
                else:
                    spreading_obliquity_for_prev_normal = None
                
                # Calculate spreading obliquity for next normal. If there is no next normal or
                # obliquity is larger than 45 degrees then set to None.
                if next_normal:
                    # Range [0,180].
                    spreading_obliquity_for_next_normal = math.degrees(pygplates.Vector3D.angle_between(
                        velocity_vector,
                        next_normal))
                    # Range [0,90].
                    if spreading_obliquity_for_next_normal > 90:
                        spreading_obliquity_for_next_normal = 180 - spreading_obliquity_for_next_normal
                    # Exclude transform sections (obliquity larger than 45 degrees).
                    if spreading_obliquity_for_next_normal > 45:
                        spreading_obliquity_for_next_normal = None
                else:
                    spreading_obliquity_for_next_normal = None
                
                if (spreading_obliquity_for_prev_normal is not None and
                    spreading_obliquity_for_next_normal is not None):
                    # Both previous and next normals have ridge-like spreading obliquities.
                    spreading_obliquity = 0.5 * (spreading_obliquity_for_prev_normal + spreading_obliquity_for_next_normal)
                elif spreading_obliquity_for_prev_normal is not None:
                    # Only previous normal has ridge-like spreading obliquity.
                    spreading_obliquity = spreading_obliquity_for_prev_normal
                elif spreading_obliquity_for_next_normal is not None:
                    # Only next normal has ridge-like spreading obliquity.
                    spreading_obliquity = spreading_obliquity_for_next_normal
                else:
                    # Is a transform, so hard-wire to 90 degrees.
                    spreading_obliquity = 90

                spreading_obliquities.append(spreading_obliquity)
            scalar_coverages[SCALAR_TYPE_SPREADING_OBLIQUITY] = spreading_obliquities
    
    # Output age scalar coverage (if requested) - per-point age values.
    if SCALAR_TYPE_AGE in output_scalar_types:
        # Each point in the polyline has the same age.
        scalar_coverages[SCALAR_TYPE_AGE] = [interpolated_time] * len(interpolated_polyline)
    
    # Convert the geometry to a coverage which is a (geometry, scalar-coverages-dict) tuple.
    return (interpolated_polyline, scalar_coverages)


#
# Private function.
#
# Write a GMT ".xy" feature containing an interpolated isochron polyline and associated
# per-point scalar values (one point and scalar values per line in the file).
# 
def _write_coverage_to_xy_file(xy_file, coverage, output_scalar_types):
    
    coverage_geometry, coverage_scalars = coverage
    coverage_points = coverage_geometry.get_points()
    
    # Get the points as a list of latitudes and a list of longitudes.
    lat_lon_points = [point.to_lat_lon() for point in coverage_points]
    latitudes, longitudes = zip(*lat_lon_points)

    # A list of lists of scalar values.
    # The first two lists are the longitudes and latitudes.
    scalars_list = [longitudes, latitudes]
    for output_scalar_type in output_scalar_types:
        scalars = coverage_scalars[output_scalar_type]
        # If asymmetries then convert from range [-1,1] to [0,100] when writing to '.xy'.
        if output_scalar_type == SCALAR_TYPE_SPREADING_ASYMMETRY:
            scalars = [50 * (scalar + 1) for scalar in scalars]
        
        scalars_list.append(scalars)
    
    # Determine the output line format string (which is '{0} {1} {2} ...').
    line_format_string = ' '.join(['{{{0}}}'.format(i) for i in range(len(scalars_list))])
    line_format_string += '\n'
    
    # Write the start of the feature.
    xy_file.write('>\n')
    
    # Iterate over the points/scalars and write each point to a separate line.
    for point_index in range(len(coverage_points)):
        # The x, y, z, ... scalar values to write to the current line of the file.
        line_scalars = [scalars[point_index] for scalars in scalars_list]
        
        xy_file.write(line_format_string.format(*line_scalars))


def write_coverage_features_to_xy_file(
        xy_filename,
        coverage_features,
        output_scalar_types,
        rotation_model,
        anchor_plate_id=0,
        reconstruction_time=None,
        print_debug_output=0):
    """write_coverage_features_to_xy_file(xy_filename, coverage_features, output_scalar_types, rotation_model, anchor_plate_id, reconstruction_time, print_debug_output)
    Write a GMT ".xy" ascii file containing the interpolated isochron polylines and associated per-point scalar values (one point and scalar values per line in the file).
    
    :param xy_filename: filename of '.xy' file
    :type xy_filename: str
    :param coverage_features: the features
    :type coverage_features: a sequence of pygplates.Feature
    :param output_scalar_types: list of scalar types (must match those of coverages in 'coverage_features')
    :type output_scalar_types: list of pygplates.ScalarType
    :param rotation_model: the rotation model
    :type rotation_model: pygplates.RotationModel
    :param anchor_plate_id: the closeness threshold, in radians, used to determine if two geometries are adjacent
    :type anchor_plate_id: int - defaults to zero
    :param reconstruction_time: the reconstruction time (or None if not reconstructing)
    :type reconstruction_time: float or None - defaults to None
    :param print_debug_output: the level at which to print debug output - zero means no debug output, \
            one or more means debug output
    :type print_debug_output: int - defaults to zero
    """
    
    if reconstruction_time is None:
        if print_debug_output >= 1:
            print('Write features...')
        
        with open(xy_filename, 'w') as xy_file:
            for feature in coverage_features:
                
                coverage = feature.get_geometry(coverage_return=pygplates.CoverageReturn.geometry_and_scalars)
                if not coverage:
                    # Shouldn't be able to get here since all feature should have coverages as geometries.
                    continue
                
                _write_coverage_to_xy_file(xy_file, coverage, output_scalar_types)
    
    else:
        if print_debug_output >= 1:
            print('Exporting reconstructed features...')
        
        # Reconstruct the features.
        reconstructed_feature_geometries = []
        pygplates.reconstruct(coverage_features, rotation_model, reconstructed_feature_geometries, reconstruction_time, anchor_plate_id)
        
        with open(xy_filename, 'w') as xy_file:
            for reconstructed_feature_geometry in reconstructed_feature_geometries:
                
                # Get the coverage from the feature (contains scalars values and present-day geometry).
                coverage = reconstructed_feature_geometry.get_feature().get_geometry(coverage_return=pygplates.CoverageReturn.geometry_and_scalars)
                if not coverage:
                    # Shouldn't be able to get here since all feature should have coverages as geometries.
                    continue
                
                # Replace the present day geometry with the reconstructed geometry.
                # The coverage is a tuple of (geometry, scalar-values-dict).
                coverage = (reconstructed_feature_geometry.get_reconstructed_geometry(), coverage[1])
                
                _write_coverage_to_xy_file(xy_file, coverage, output_scalar_types)


def join_adjacent_features_with_same_type_and_begin_time_and_plate_ids(features, distance_threshold_radians, print_debug_output):
    """join_adjacent_features_with_same_type_and_begin_time_and_plate_ids(features, distance_threshold_radians)
    Joins features that have adjacent geometries and have the same feature type, begin time and plate ids.
    
    :param features: the features
    :type features: a sequence of pygplates.Feature
    :param distance_threshold_radians: the closeness threshold, in radians, used to determine if two geometries are adjacent
    :type distance_threshold_radians: float
    :param print_debug_output: the level at which to print debug output (0, 1, 2 or 3) - zero means no debug output
    :type print_debug_output: int - defaults to zero
    """
    
    if print_debug_output >= 1:
        print('Join features...')
    
    # Group features by (begin time, reconstruction plate id, conjugate plate id, feature type).
    join_feature_groups = {}
    for feature in features:
        begin_time, end_time = feature.get_valid_time()
        
        # Get left and right plate ids (in one form or another), otherwise skip the current feature.
        left_right_plate_ids = _get_left_right_plate_ids(feature)
        if not left_right_plate_ids:
            if print_debug_output >= 1:
                print(" ...skipping '{0}' "
                        "feature (id: {1}): No left/right or reconstruction/conjugate plate IDs and no "
                        "'gpml:OldPlatesHeader' property."
                                .format(feature.get_feature_type(), feature.get_feature_id()))
            continue
        
        left_plate_id, right_plate_id = left_right_plate_ids
        
        join_feature_group = join_feature_groups.setdefault(
                (begin_time, left_plate_id, right_plate_id, feature.get_feature_type()),
                [])
        join_feature_group.append(feature)
    
    joined_features = []
    
    # Iterate over the groups (lists) of features to join and attempt to join.
    for join_feature_group in join_feature_groups.values():
        # Join any adjacent features in the group - this potentially modifies 'join_feature_group'.
        _remove_features_with_duplicate_geometries(join_feature_group)
        joined_features.extend(
                _join_adjacent_features(join_feature_group, distance_threshold_radians, print_debug_output))
    
    return joined_features


def interpolate_isochrons(
        rotation_model,
        isochron_cob_ridge_features,
        # Note: can be a distance spacing in radians, or a list of interpolation times (see docstring below)...
        interpolate=DEFAULT_INTERPOLATE_RESOLUTION_RADIANS,
        minimum_latitude_overlap_radians=DEFAULT_MINIMUM_LATITUDE_OVERLAP_RADIANS,
        maximum_isochron_latitude_non_overlap_radians=DEFAULT_MAXIMUM_ISOCHRON_LATITUDE_NON_OVERLAP_RADIANS,
        maximum_cob_latitude_non_overlap_radians=DEFAULT_MAXIMUM_COB_LATITUDE_NON_OVERLAP_RADIANS,
        maximum_distance_threshold_radians=DEFAULT_MAXIMUM_DISTANCE_THRESHOLD_RADIANS,
        maximum_age_difference=DEFAULT_MAXIMUM_AGE_DIFFERENCE,
        join_polylines_threshold_radians=DEFAULT_JOIN_POLYLINES_THRESHOLD_RADIANS,
        **kwargs):
    """interpolate_isochrons(rotation_model, isochron_cob_ridge_features, \
            [interpolate=DEFAULT_INTERPOLATE_RESOLUTION_RADIANS], \
            [minimum_latitude_overlap_radians=DEFAULT_MINIMUM_LATITUDE_OVERLAP_RADIANS], \
            [maximum_isochron_latitude_non_overlap_radians=DEFAULT_MAXIMUM_ISOCHRON_LATITUDE_NON_OVERLAP_RADIANS], \
            [maximum_cob_latitude_non_overlap_radians=DEFAULT_MAXIMUM_COB_LATITUDE_NON_OVERLAP_RADIANS], \
            [maximum_distance_threshold_radians=DEFAULT_MAXIMUM_DISTANCE_THRESHOLD_RADIANS], \
            [maximum_age_difference=DEFAULT_MAXIMUM_AGE_DIFFERENCE], \
            [join_polylines_threshold_radians=DEFAULT_JOIN_POLYLINES_THRESHOLD_RADIANS], \
            [\\*\\*kwargs]) \
    
    Generate interpolated isochrons using isochrons, continental-oceanic boundaries and present-day ridges.
    
    
    :param rotation_model: the rotation model
    :type rotation_model: pygplates.RotationModel
    
    :param isochron_cob_ridge_features: input features containing isochrons, continental-oceanic \
           boundaries and ridge boundaries
    :type isochron_cob_ridge_features: pygplates.FeatureCollection, or string (filename), or pygplates.Feature, \
          or sequence of pygplates.Feature, or sequence of any combination of those four types
    
    :param interpolate: if a single number then *interpolate* is the interval spacing, in radians, \
           between input features at which to generate interpolated isochrons - otherwise if a sequence of \
           numbers (eg, list or tuple) then *interpolate* is the sequence of times at which to generate \
           interpolated isochrons - by default it is the single number DEFAULT_INTERPOLATE_RESOLUTION_RADIANS
    :type interpolate: float or list of float
    
    :param minimum_latitude_overlap_radians: required amount of latitude overlap of two input feature \
           polyline geometries in order for them to be interpolated
    :type minimum_latitude_overlap_radians: float
    
    :param maximum_isochron_latitude_non_overlap_radians: allowed non-overlapping latitude region \
           when the older of two input feature geometries being interpolated belongs to an isochron feature
    :type maximum_isochron_latitude_non_overlap_radians: float
    
    :param maximum_cob_latitude_non_overlap_radians: allowed non-overlapping latitude region \
           when the older of two input feature geometries being interpolated belongs to a COB feature
    :type maximum_cob_latitude_non_overlap_radians: float
    
    :param maximum_distance_threshold_radians: maximum distance (in radians) between two input \
           feature geometries above which they will not be interpolated
    :type maximum_distance_threshold_radians: float
    
    :param maximum_age_difference: maximum difference in age, in My, between two input feature \
           geometries above which they will not be interpolated
    :type maximum_age_difference: float
    
    :param join_polylines_threshold_radians: the closeness threshold, in radians, used to determine if \
           two polyline input feature geometries are adjacent and hence can be joined into a single geometry
    :type join_polylines_threshold_radians: float
    
    :returns: the interpolated features
    :rtype: list of pygplates.Feature
    :raises: OpenFileForReadingError if any file is not readable (when filenames specified)
    :raises: FileFormatNotSupportedError if any file format (identified by the filename \
             extensions) does not support reading (when filenames specified)
    
    
    The following optional keyword arguments are supported by *kwargs*:

    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | Name                                | Type  | Default | Description                                                                 |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | tessellate_threshold_radians        | float | None    | The maximum tessellation angle, in radians, between adjacent points in      |
    |                                     |       |         | each interpolated polyline.                                                 |
    |                                     |       |         | Note that this is *along* each polyline (not between polylines).            |
    |                                     |       |         | If not specified (None) then there is no tessellation.                      |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | plate_pairs                         | list  | None    | Optional list of plate pairs to restrict output to.                         |
    |                                     |       |         | This is a list of (int, int) tuples.                                        |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | print_debug_output                  | int   | 0       | The level at which to print debug output (0, 1, 2 or 3).                    |
    |                                     |       |         | Zero means no debug output.                                                 |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | output_scalar_spreading_asymmetry   | bool  | False   | Store spreading asymmetry at each point in each interpolated isochron.      |
    |                                     |       |         | This is in the range [-1,1] where 0 represents half-stage rotation,         |
    |                                     |       |         | -1 represents spreading only on the conjugate flank and                     |
    |                                     |       |         | +1 represents spreading only on the flank containing the isochron.          |
    |                                     |       |         | This will stored under the 'gpml:SpreadingAsymmetry' scalar type.           |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | output_scalar_spreading_rate        | bool  | False   | Store spreading rate at each point in each interpolated isochron.           |
    |                                     |       |         | This is the spreading rate relative to the mid-ocean ridge.                 |
    |                                     |       |         | This velocity magnitude is in Kms/My.                                       |
    |                                     |       |         | This will stored under the 'gpml:SpreadingRate' scalar type.                |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | output_scalar_full_spreading_rate   | bool  | False   | Store *full* spreading rate at each point in each interpolated isochron.    |
    |                                     |       |         | This is the spreading rate relative to the conjugate plate.                 |
    |                                     |       |         | This velocity magnitude is in Kms/My.                                       |
    |                                     |       |         | This will stored under the 'gpml:SpreadingRate' scalar type.                |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | output_scalar_spreading_direction   | bool  | False   | Store spreading direction at each point in each interpolated isochron.      |
    |                                     |       |         | This is an angle (in degrees) clockwise (East-wise) from North (0 to 180).  |
    |                                     |       |         | The lowest azimuth of the two opposite-pointing directions of spreading     |
    |                                     |       |         | at the time of crustal accretion.                                           |
    |                                     |       |         | This will stored under the 'gpml:SpreadingDirection' scalar type.           |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | output_scalar_spreading_obliquity   | bool  | False   | Store spreading obliquity at each point in each interpolated isochron.      |
    |                                     |       |         | This is an angle from 0 to 90 degrees.                                      |
    |                                     |       |         | The amount that the small circle of spreading rotation pole deviates from   |
    |                                     |       |         | the perpendicular line at the ridge point at the time of crustal accretion. |
    |                                     |       |         | This will stored under the 'gpml:SpreadingObliquity' scalar type.           |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    | output_scalar_age                   | bool  | False   | Store age at each point in each interpolated isochron.                      |
    |                                     |       |         | This is the interpolated isochron's age (in Ma). It is constant for all     |
    |                                     |       |         | points within an isochron and is equal to the isochron's begin time.        |
    |                                     |       |         | This will stored under the 'gpml:Age' scalar type.                          |
    +-------------------------------------+-------+---------+-----------------------------------------------------------------------------+
    
    
    The algorithms used in this function very closely match those used in the C++ program Intertec.
    
    If *interpolate* is a single number then features are interpolated at regularly spaced intervals
    of *interpolate* radians between input features. Most of these interpolated features are isochrons
    with the boundaries being continental-oceanic boundaries (COB).
    
    If *interpolate* is a sequence of numbers then features are interpolated at geological times
    specified by these numbers. All of these interpolated features will be isochrons.
    
    For more details on *minimum_latitude_overlap_radians*, *maximum_isochron_latitude_non_overlap_radians*,
    *maximum_cob_latitude_non_overlap_radians* and *maximum_distance_threshold_radians*
    please see the pygplates API documentation for the function pygplates.PolylineOnSphere.rotation_interpolate.
    
    To interpolate isochrons with a spacing of 2 minutes (with a minimum required latitude
    overlap of 1 degree and with an allowed latitude non-overlap of up to 3 degrees for isochrons
    and up to 1 degree for COBs and with other parameters assuming default values):
    ::
    
        interpolated_isochrons = interpolate_isochrons(
                pygplates.RotationModel('rotations.rot'),
                ('cobs.gpml', 'isochrons.gpml', 'ridges.gpml'),
                math.radians(2 / 60.0),
                math.radians(1),
                math.radians(3),
                math.radians(1),
                print_debug_output = 1)
    
    ...where *print_debug_output* is one of the optional keyword arguments shown in the *kwargs* table above.
    
    To interpolate isochrons at times between 0Ma and 100Ma at 10My intervals (with a minimum
    required latitude overlap of 1 degree and with an allowed latitude non-overlap of up to 3 degrees
    for isochrons and up to 1 degree for COBs and with other parameters assuming default values):
    ::
    
        interpolated_isochrons = interpolate_isochrons(
                pygplates.RotationModel('rotations.rot'),
                ('cobs.gpml', 'isochrons.gpml', 'ridges.gpml'),
                range(0, 101, 10),
                math.radians(1),
                math.radians(3),
                math.radians(1),
                print_debug_output = 1)
    """
    
    # Check the imported pygplates version.
    if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
        raise RuntimeError('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
                os.path.basename(__file__), pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))
    
    # Keyword arguments.
    #
    tessellate_threshold_radians = kwargs.pop('tessellate_threshold_radians', None)
    plate_pairs = kwargs.pop('plate_pairs', None)
    print_debug_output = kwargs.pop('print_debug_output', 0)
    # Scalar coverages to output.
    output_scalar_types = set()
    if kwargs.pop('output_scalar_spreading_asymmetry', False):
        output_scalar_types.add(SCALAR_TYPE_SPREADING_ASYMMETRY)
    if kwargs.pop('output_scalar_spreading_rate', False):
        output_scalar_types.add(SCALAR_TYPE_SPREADING_RATE)
    if kwargs.pop('output_scalar_full_spreading_rate', False):
        output_scalar_types.add(SCALAR_TYPE_FULL_SPREADING_RATE)
    if kwargs.pop('output_scalar_spreading_direction', False):
        output_scalar_types.add(SCALAR_TYPE_SPREADING_DIRECTION)
    if kwargs.pop('output_scalar_spreading_obliquity', False):
        output_scalar_types.add(SCALAR_TYPE_SPREADING_OBLIQUITY)
    if kwargs.pop('output_scalar_age', False):
        output_scalar_types.add(SCALAR_TYPE_AGE)
    # Raise error if any unknown keyword arguments.
    if kwargs:
        raise TypeError('Unknown keyword arguments {0}'.format(kwargs.keys()))
    
    # 'interpolate' is either a sequence of times or an interval spacing (in radians).
    if hasattr(interpolate, '__iter__'):
        # Sort sequence of interpolation times to make it easier to find input features
        # whose times overlap them.
        interpolation_times = sorted(float(time) for time in interpolate)
        using_interpolation_times = True
    else:
        interpolate_interval_radians = float(interpolate)
        using_interpolation_times = False
    
    # Turn function argument into something more convenient for extracting features.
    isochron_cob_ridge_features = pygplates.FeaturesFunctionArgument(isochron_cob_ridge_features)
    
    if print_debug_output >= 1:
        print('Read features...')
    isochron_cob_ridge_features = isochron_cob_ridge_features.get_features()
    
    # Filter the features if requested.
    if plate_pairs is not None:
        _filter_plate_pairs(isochron_cob_ridge_features, plate_pairs)
    
    joined_isochron_cob_ridge_features = join_adjacent_features_with_same_type_and_begin_time_and_plate_ids(
            isochron_cob_ridge_features,
            join_polylines_threshold_radians,
            print_debug_output)
    
    if print_debug_output >= 1:
        if using_interpolation_times:
            print('Interpolate isochrons at geological times {0}...'.format(interpolation_times))
        else:
            print('Interpolate isochrons at {0} minute distance intervals...'.format(
                    60 * math.degrees(interpolate_interval_radians)))
    
    # Group features by conjugate plate ID pairs (plate pair flank).
    feature_groups = {}
    for feature in joined_isochron_cob_ridge_features:
        # Get left and right plate ids (in one form or another), otherwise skip the current feature.
        left_right_plate_ids = _get_left_right_plate_ids(feature)
        if not left_right_plate_ids:
            # We shouldn't get here since join_adjacent_features_with_same_type_and_begin_time_and_plate_ids()
            # will have removed these features for us.
            if print_debug_output >= 1:
                print(" ...skipping '{0}' feature (id: {1}): No left/right or "
                        "reconstruction/conjugate plate IDs and no 'gpml:OldPlatesHeader' property."
                                .format(feature.get_feature_type(), feature.get_feature_id()))
            continue
        
        left_plate_id, right_plate_id = left_right_plate_ids

        # Add feature to (plate pair) group.
        feature_groups.setdefault((left_plate_id, right_plate_id), []).append(feature)
        
        # Add mid-ocean ridges to both flanks (plate pairs) since they are not unique to a
        # particular flank (like isochrons and COBs are).
        if feature.get_feature_type() == MID_OCEAN_RIDGE_FEATURE_TYPE:
            feature_groups.setdefault((right_plate_id, left_plate_id), []).append(feature)
    
    # Sort the features in each (plate pair) group by begin time.
    # If begin times are same then sort such that COB features go last since they
    # terminate the sequence of isochrons on the current flank.
    for feature_group in feature_groups.values():
        feature_group.sort(key = lambda feature: (feature.get_valid_time()[0], feature.get_feature_type()==COB_FEATURE_TYPE))
    
    interpolated_features = []
    
    # Iterate over the feature groups (plate pair flanks).
    for (left_plate_id, right_plate_id), feature_group in feature_groups.items():
        if print_debug_output >= 1:
            print(' Interpolating left plate {0}, right plate {1}'.format(left_plate_id, right_plate_id))
            if print_debug_output >= 2:
                if len(feature_group) == 1:
                    print('  ...skipping flank - only one feature')
        
        # Some features have non-zero finite rotations at present day (0Ma).
        # To account for this we reconstruct features to present day before interpolating and
        # reverse reconstruct when storing interpolated lines in new features.
        #
        # Interpolated isochrons have a reconstruction plate ID equal to the left plate ID.
        # We need to use the plate ID assigned to the interpolated isochrons since they'll use that in turn to reconstruct.
        inverse_present_day_rotation_for_interpolated_isochrons = rotation_model.get_rotation(0, left_plate_id).get_inverse()
        
        # Iterate over the flank's features progressing from oldest time to youngest time.
        # Note that the features have already been sorted by begin time.
        for old_feature_index in range(len(feature_group) - 1, 0, -1): # Note: excludes index 0
            old_feature = feature_group[old_feature_index]
            old_feature_begin_time, old_feature_end_time = old_feature.get_valid_time()
            old_time = old_feature_begin_time
            
            # If we're interpolating at specific geological times
            # rather than at regular distance intervals...
            if using_interpolation_times:
                # If the youngest interpolation time is larger than 'old_time' then
                # we've finished with the current flank.
                if (not interpolation_times or
                    interpolation_times[0] > old_time):
                    break
            
            if print_debug_output >= 2:
                print('  Old time: {0}'.format(old_time))
            
            # Get all geometries in 'old_feature'.
            old_feature_geoms = old_feature.get_all_geometries()
            #
            # NOTE: Some features have non-zero finite rotations at present day (0Ma).
            # To account for this we reconstruct the feature to present day.
            # This will not change the present day geometry position for features with zero present day rotations.
            #
            # We can't just use the left plate ID because each mid-ocean ridge is added to its two adjacent flanks
            # with swapped left/right plates - so we should use its reconstruction plate ID instead.
            present_day_rotation_for_old_feature = rotation_model.get_rotation(0, old_feature.get_reconstruction_plate_id())
            # Note that we don't use 'pygplates.reconstruct' since that excludes features that don't exist at present day.
            old_feature_geoms = [present_day_rotation_for_old_feature * geom for geom in old_feature_geoms]
            
            # Find the next feature that has a smaller time value and is not a COB.
            young_feature_index = old_feature_index - 1
            while young_feature_index >= 0:
                young_feature = feature_group[young_feature_index]
                if young_feature.get_feature_type() != COB_FEATURE_TYPE:
                    young_feature_begin_time, young_feature_end_time = young_feature.get_valid_time()
                    if young_feature_begin_time < old_feature_begin_time:
                        break
                young_feature_index -= 1
            
            if young_feature_index < 0:
                if print_debug_output >= 2:
                    print('   ...skipping old feature - unable to find a younger, non-COB feature')
                continue
            
            old_feature_and_young_feature_overlapped = False
            young_time = None
            
            # Iterate over all remaining young_feature's until we find an overlap.
            # All young_feature's with the same begin time and that overlap old_feature are interpolated
            # with old_feature - this is in case those young_feature's didn't get joined together.
            for young_feature_index in range(young_feature_index, -1, -1): # Note: includes index 0
                young_feature = feature_group[young_feature_index]
                young_feature_begin_time, young_feature_end_time = young_feature.get_valid_time()
                
                # If the next young_feature has a different begin time (than the previous young_feature) then:
                #   1) If old_feature and young_feature have overlapped then we're finished, otherwise
                #   2) Try a new subset of young_feature's that have the next smaller begin time.
                if young_feature_begin_time != young_time:
                    if old_feature_and_young_feature_overlapped:
                        break
                    
                    # Note: we'll always get here on the first iteration on the loop.
                    
                    young_time = young_feature_begin_time
                    
                    # The end time for all interpolated isochrons - set to longer lasting of young and old.
                    interpolation_end_time = min(young_feature_end_time, old_feature_end_time)
                    
                    if print_debug_output >= 2:
                        print('   Young time: {0}'.format(young_time))
                    
                    # Skip to the next old feature if the age difference is too large.
                    if old_time - young_time > maximum_age_difference:
                        if print_debug_output >= 2:
                            print('    ...skipping old/young feature pair - age difference {0} exceeded maximum {1}'.format(
                                    old_time - young_time, maximum_age_difference))
                        break
                    
                    # Calculate the stage rotation from time 'old_time' to time 'young_time' and
                    # of the left plate (moving plate) relative to the right plate (fixed plate).
                    #
                    # We are interpolating polylines belonging to the left plate so the moving plate
                    # is the left plate (and the fixed plate is the right plate).
                    #
                    # NOTE: We also set the anchor plate to the fixed plate in case there is no
                    # plate circuit path from the default anchor plate (zero) to either the moving
                    # or fixed plate (but where a path still exists between moving and fixed plates).
                    stage_rotation = rotation_model.get_rotation(
                            young_time, left_plate_id, old_time, anchor_plate_id=right_plate_id)
                    # If rotation not found (or is identity rotation - has no axis) then skip interpolation altogether.
                    if pygplates.FiniteRotation.represents_identity_rotation(stage_rotation):
                        if print_debug_output >= 2:
                            print('    ...skipping old/young feature pair - identity stage rotation')
                        break
                    stage_rotation_pole, stage_pole_angle_radians = stage_rotation.get_euler_pole_and_angle()
                    
                    # Present day geometries need to be rotated, relative to the right plate, to the 'from' time
                    # of the above stage rotation (which is 'old_time') so that they can then be rotated by
                    # the stage rotation. To avoid having to rotate the present day geometries into this stage pole
                    # reference frame we can instead apply the inverse rotation to the stage pole itself.
                    #
                    #   geometry_moving_plate = R(0->young, A->Left) * geometry_present_day
                    #                         = R(0->young, A->Right) * R(0->young, Right->Left) * geometry_present_day
                    #                         = R(0->young, A->Right) * R(old->young, Right->Left) * R(0->old, Right->Left) * geometry_present_day
                    #
                    # ...so the present day geometry needs to be rotated by 'R(0->old, Right->Left)' before the stage rotation
                    # 'R(old->young, Right->Left)' can be applied to it. Alternatively we can reverse rotate the stage pole
                    # (of the stage rotation) by 'inverse[R(0->old, Right->Left)]' so that it can be applied to the present day geometry.
                    #
                    # For more detail see:
                    #   http://www.gplates.org/docs/pygplates/sample-code/pygplates_split_isochron_into_ridges_and_transforms.html
                    #
                    # However, there's a slight complication because some plates have non-zero finite rotations at present day (0Ma).
                    # This means the geometry stored in the feature and the geometry reconstructed to present day are in different.
                    # Normally the geometry in the feature should be in its present day position, but we'll be lenient since that's not always the case.
                    #
                    #   geometry_present_day = R(0, A->Left) * geometry_in_feature
                    #    geometry_in_feature = inverse[R(0, A->Left)] * geometry_present_day
                    #
                    # ...where 'R(0, A->Left)' means reconstructing from its "un-reconstructed" position (ie, not present day) to present day.
                    #
                    # So now the above equation looks like:
                    #
                    #   geometry_moving_plate = R(young, A->Left) * geometry_in_feature
                    #                         = R(young, A->Right) * R(young, Right->Left) * geometry_in_feature
                    #                         = R(young, A->Right) * R(old->young, Right->Left) * R(old, Right->Left) * geometry_in_feature
                    #                         = R(young, A->Right) * R(old->young, Right->Left) * R(old, Right->Left) * inverse[R(0, A->Left)] * geometry_present_day
                    #
                    # ...so the present day geometry needs to be rotated by 'R(old, Right->Left) * inverse[R(0, A->Left)]' before the stage rotation
                    # 'R(old->young, Right->Left)' can be applied to it. Alternatively we can reverse rotate the stage pole (of the stage rotation) by
                    # the inverse of that so it can be applied to the present day geometry.
                    #
                    to_stage_pole_reference_frame_for_present_day_geometry = rotation_model.get_rotation(
                            old_time, left_plate_id, anchor_plate_id=right_plate_id) * inverse_present_day_rotation_for_interpolated_isochrons
                    from_stage_pole_reference_frame_for_present_day_geometry = to_stage_pole_reference_frame_for_present_day_geometry.get_inverse()
                    stage_rotation_pole_for_present_day_geometry = from_stage_pole_reference_frame_for_present_day_geometry * stage_rotation_pole
                    if print_debug_output >= 3:
                        print('    Stage pole (lat,lon) for present day geometry: {0}'.format(stage_rotation_pole_for_present_day_geometry.to_lat_lon()))

                    #
                    # Also the stage rotation gets applied to present day isochron geometries when calculating velocities (well, before being reconstructed
                    # to the crustal accretion time) so it needs to work with present day geometries.
                    # Applied to an interpolated isochron with time of appearance 'begin_time' (in range [old, young]):
                    #
                    #   geometry_moving_plate = R(begin_time, A->Left) * geometry_in_feature
                    #                         = R(begin_time, A->Right) * R(begin_time, Right->Left) * geometry_in_feature
                    #                         = R(begin_time, A->Right) * R(old->begin_time, Right->Left) * R(old, Right->Left) * geometry_in_feature
                    #                         = R(begin_time, A->Right) * R(young->begin_time, Right->Left) * R(old->young, Right->Left) * R(old, Right->Left) * geometry_in_feature
                    #                         = R(begin_time, A->Right) * R(young->begin_time, Right->Left) * R(old->young, Right->Left) * R(old, Right->Left) * inverse[R(0, A->Left)] * geometry_present_day
                    #
                    # ...so we need to rotate present day geometry by 'R(old, Right->Left) * inverse[R(0, A->Left)]' before the stage rotation 'R(old->young, Right->Left)'
                    # can be applied to it. Then the stage rotation is applied. And finally we reverse-rotate back again.
                    # These three operations can be combined into a single stage rotation:
                    #
                    #   stage_rotation_for_present_day = inverse[R(old, Right->Left) * inverse[R(0, A->Left)]] * stage_rotation * R(old, Right->Left) * inverse[R(0, A->Left)]
                    #
                    stage_rotation_for_present_day_geometry = from_stage_pole_reference_frame_for_present_day_geometry * stage_rotation * to_stage_pole_reference_frame_for_present_day_geometry
                     
                    if using_interpolation_times:
                        # Create a list of interpolation ratios where:
                        #      0 <= (interpolation_time - young_time) / (old_time - young_time) < 1
                        #
                        # ...corresponding to interpolation time between 'old_time' and 'young_time'...
                        #      young_time <= interpolation_time < old_time
                        #
                        # Note: the interpolation times are sorted from youngest to oldest.
                        # Note: the list can be empty in which case pygplates just checks that
                        # interpolation is possible (ie, latitude overlap and distance threshold).
                        interpolation_ratios = []
                        # Search forward for youngest interpolation time satisfying the condition.
                        for young_time_index in range(len(interpolation_times)):
                            if interpolation_times[young_time_index] >= young_time:
                                # Search backward for oldest interpolation time satisfying the condition.
                                for old_time_index in range(len(interpolation_times),young_time_index,-1):
                                    if interpolation_times[old_time_index-1] < old_time:
                                        # Map range [young_time, old_time] to range [0, 1].
                                        # Note that we won't get a divide-by-zero error because we've guaranteed
                                        # 'young_time < old_time' above (just before entering the young feature loop).
                                        interp_denom = 1.0 / (old_time - young_time)
                                        interpolation_ratios.extend(
                                                # Use min/max in case numerical tolerance puts outside [0,1] range...
                                                min(max((time - young_time) * interp_denom, 0), 1)
                                                for time in interpolation_times[young_time_index:old_time_index])
                                        break
                                break
                    
                if young_feature.get_feature_type() != COB_FEATURE_TYPE:
                    # Get all geometries in 'young_feature'.
                    young_feature_geoms = young_feature.get_all_geometries()
                    #
                    # NOTE: Some features have non-zero finite rotations at present day (0Ma).
                    # To account for this we reconstruct the feature to present day.
                    # This will not change the present day geometry position for features with zero present day rotations.
                    #
                    # We can't just use the left plate ID because each mid-ocean ridge is added to its two adjacent flanks
                    # with swapped left/right plates - so we should use its reconstruction plate ID instead.
                    present_day_rotation_for_young_feature = rotation_model.get_rotation(0, young_feature.get_reconstruction_plate_id())
                    # Note that we don't use 'pygplates.reconstruct' since that excludes features that don't exist at present day.
                    young_feature_geoms = [present_day_rotation_for_young_feature * geom for geom in young_feature_geoms]
                    
                    # Iterate over all old_feature/young_feature geometry pairs.
                    for old_feature_geom in old_feature_geoms:
                        for young_feature_geom in young_feature_geoms:
                            # Test whether we're interpolating:
                            #  1) at specific geological times, or
                            #  2) at regular distance intervals.
                            if using_interpolation_times:
                                # To pass to 'pygplates.PolylineOnSphere.rotation_interpolate'.
                                #
                                # Note: this can be an empty list in which case pygplates just checks
                                # that interpolation is possible (ie, latitude overlap and distance threshold).
                                
                                # If we're going to calculate spreading asymmetries (ie, if we're going to
                                # output scalar coverages) then we also need interpolated polylines at the
                                # young and old times (versions of the original young and old polylines but
                                # with the same number of latitude-aligned points as the interpolated polylines).
                                if output_scalar_types:
                                    # Copy list and add 0 to the front and 1 to the back.
                                    pygplates_interpolate_parameter = interpolation_ratios[:]
                                    pygplates_interpolate_parameter.insert(0, 0)
                                    pygplates_interpolate_parameter.append(1)
                                else:
                                    pygplates_interpolate_parameter = interpolation_ratios
                                
                            else:
                                # To pass to 'pygplates.PolylineOnSphere.rotation_interpolate'.
                                #
                                # Note that we don't need to add for two extra polylines when calculating
                                # spreading asymmetries because we will always get interpolated polylines
                                # at the young and old times anyway.
                                pygplates_interpolate_parameter = interpolate_interval_radians
                            
                            # The amount of latitude non-overlapping to allow is based on the feature type.
                            if old_feature.get_feature_type() == COB_FEATURE_TYPE:
                                maximum_latitude_non_overlap_radians = maximum_cob_latitude_non_overlap_radians
                            else:
                                maximum_latitude_non_overlap_radians = maximum_isochron_latitude_non_overlap_radians
                            
                            # Interpolate isochrons between old_feature and young_feature if they overlap.
                            #
                            # Note: Some differences in overlap compared to the original Intertec program:
                            # Original Intertec tested for latitude overlap of 3 degrees if there
                            # was more than one geometry with the same begin time (otherwise it did no
                            # overlap testing). But we test for overlap in all cases (also our default overlap
                            # is less than 3 degrees). Also we filter out geometries that are smaller
                            # than the minimum overlap whereas original Intertec did not.
                            try:
                                interpolated_polylines = pygplates.PolylineOnSphere.rotation_interpolate(
                                        young_feature_geom,
                                        old_feature_geom,
                                        stage_rotation_pole_for_present_day_geometry,
                                        pygplates_interpolate_parameter,
                                        minimum_latitude_overlap_radians,
                                        maximum_latitude_non_overlap_radians,
                                        maximum_distance_threshold_radians,
                                        # Flatten longitude overlaps using the young feature's geometry...
                                        pygplates.FlattenLongitudeOverlaps.use_from,
                                        pygplates.PolylineConversion.raise_if_non_polyline)
                            except pygplates.GeometryTypeError:
                                if print_debug_output >= 2:
                                    print('    ...skipping old/young geometry pair - geometry(s) not polyline(s)')
                                # If either geometry is not a polyline then skip the interpolation.
                                continue
                            
                            # If the polylines do not overlap in latitude, or are too far away, then skip.
                            #
                            # Note: We explicitly test for None since an empty list, which also evaluates
                            # to false, means (for 'using_interpolation_times') that the polylines passed
                            # but the list is empty because no interpolation times where between
                            # 'old_time' and 'young_time' (but we still needed to go through the same
                            # code path as 'not using_interpolation_times' since this affects
                            # 'old_feature_and_young_feature_overlapped' which determines which
                            # old and young features to subsequently interpolate).
                            if interpolated_polylines is None:
                                if print_debug_output >= 2:
                                    print('    ...skipping old/young geometry pair - geometries do not '
                                            'overlap sufficiently or are too far apart')
                                continue
                            
                            old_feature_and_young_feature_overlapped = True
                                
                            # Calculate spreading asymmetries if we're going to output scalar coverages
                            # because most scalar types rely on the asymmetry.
                            if output_scalar_types:
                                # There's always at least two polylines.
                                # These are versions of the original young and old polylines but with
                                # the same number of latitude-aligned points as the interpolated polylines.
                                # Hence the number of asymmetries will match the number of interpolated points.
                                spreading_asymmetries = _calc_spreading_asymmetries(
                                        interpolated_polylines[0],
                                        interpolated_polylines[-1],
                                        stage_rotation_pole_for_present_day_geometry,
                                        stage_pole_angle_radians)
                                
                                # If we added two extra interpolated polylines (young and old) to
                                # calculate spreading asymmetries then we can remove them now.
                                # Note that we only needed to add two extra polylines when
                                # using interpolation times because we already had them when
                                # doing regular interval interpolation.
                                if using_interpolation_times:
                                    # Remove first and last interpolated polylines.
                                    del interpolated_polylines[0]
                                    del interpolated_polylines[-1]
                                
                            else:
                                spreading_asymmetries = None
                            
                            # Test whether we're interpolating:
                            #  1) at specific geological times, or
                            #  2) at regular distance intervals.
                            if using_interpolation_times:
                                
                                # Lengths of 'interpolation_ratios' and
                                # 'interpolated_polylines' should be equal.
                                for interpolation_index, interpolation_ratio in enumerate(interpolation_ratios):
                                    # The interpolated polyline associated with the interpolation time...
                                    interpolated_polyline = interpolated_polylines[interpolation_index]
                                    interpolation_time = young_time + interpolation_ratio * (old_time - young_time)
                                    
                                    # Tessellate and add scalar coverages if requested.
                                    interpolated_polyline = _add_tessellated_scalar_coverages(
                                            interpolated_polyline,
                                            tessellate_threshold_radians,
                                            output_scalar_types,
                                            spreading_asymmetries,
                                            stage_rotation_for_present_day_geometry,
                                            inverse_present_day_rotation_for_interpolated_isochrons,
                                            interpolation_time,
                                            old_time - young_time,
                                            rotation_model,
                                            left_plate_id)
                                    
                                    # Some features had non-zero finite rotations at present day (0Ma).
                                    # To account for this we reconstructed the feature to present day.
                                    # However we now need to reverse these non-zero finite rotations at present day since
                                    # the geometries are stored in features - this is so when they get reconstructed to 0Ma
                                    # they will have the correct position again. We do this after adding scalar coverages since
                                    # they should use the actual (reconstructed) 0Ma positions.
                                    if output_scalar_types:
                                        coverage_polyline, coverage_scalars = interpolated_polyline
                                        present_day_feature_geometry = (inverse_present_day_rotation_for_interpolated_isochrons * coverage_polyline, coverage_scalars)
                                    else:
                                        present_day_feature_geometry = inverse_present_day_rotation_for_interpolated_isochrons * interpolated_polyline
    
                                    # Finally we can create an interpolated isochron.
                                    interpolated_feature = pygplates.Feature.create_reconstructable_feature(
                                            ISOCHRON_FEATURE_TYPE,
                                            present_day_feature_geometry,
                                            valid_time=(interpolation_time, interpolation_end_time),
                                            reconstruction_plate_id=left_plate_id,
                                            conjugate_plate_id=right_plate_id)
                                        
                                    interpolated_features.append(interpolated_feature)
                                
                                if print_debug_output >= 2:
                                    # Only print debug message if there's at least one interpolation.
                                    # Zero interpolations just means we passed the overlap test and
                                    # thus changed our path through this algorithm.
                                    if interpolated_polylines:
                                        print('    ...successful interpolation - number of interpolations {0}'.format(
                                                len(interpolated_polylines)))
                                
                            else: # ...we're interpolating at regular distance intervals...
                                
                                # Number of interpolation intervals is at least one (since always have at least two polylines).
                                num_interpolation_intervals = len(interpolated_polylines) - 1
                                interpolation_time_interval = (old_time - young_time) / num_interpolation_intervals
                                
                                # If the old feature is a COB then create an interpolated feature for it,
                                # otherwise don't (by removing the oldest interpolated polyline).
                                # The COB is typically the oldest feature and hence terminates the sequence
                                # of interpolated isochrons on the current flank.
                                if old_feature.get_feature_type() != COB_FEATURE_TYPE:
                                    del interpolated_polylines[-1]
                                
                                # Create an isochron feature for each interpolated polyline.
                                for interpolate_index, interpolated_polyline in enumerate(interpolated_polylines):
                                    # If the oldest feature is being included then it's always a COB.
                                    if interpolate_index == num_interpolation_intervals:
                                        interpolated_feature_type = COB_FEATURE_TYPE
                                    else:
                                        interpolated_feature_type = ISOCHRON_FEATURE_TYPE
                                    interpolation_time = young_time + interpolate_index * interpolation_time_interval
                                    
                                    # Tessellate and add scalar coverages if requested.
                                    interpolated_polyline = _add_tessellated_scalar_coverages(
                                            interpolated_polyline,
                                            tessellate_threshold_radians,
                                            output_scalar_types,
                                            spreading_asymmetries,
                                            stage_rotation_for_present_day_geometry,
                                            inverse_present_day_rotation_for_interpolated_isochrons,
                                            interpolation_time,
                                            old_time - young_time,
                                            rotation_model,
                                            left_plate_id)
                                    
                                    # Some features had non-zero finite rotations at present day (0Ma).
                                    # To account for this we reconstructed the feature to present day.
                                    # However we now need to reverse these non-zero finite rotations at present day since
                                    # the geometries are stored in features - this is so when they get reconstructed to 0Ma
                                    # they will have the correct position again. We do this after adding scalar coverages since
                                    # they should use the actual (reconstructed) 0Ma positions.
                                    if output_scalar_types:
                                        coverage_polyline, coverage_scalars = interpolated_polyline
                                        present_day_feature_geometry = (inverse_present_day_rotation_for_interpolated_isochrons * coverage_polyline, coverage_scalars)
                                    else:
                                        present_day_feature_geometry = inverse_present_day_rotation_for_interpolated_isochrons * interpolated_polyline
                                    
                                    # Finally we can create an interpolated isochron.
                                    interpolated_feature = pygplates.Feature.create_reconstructable_feature(
                                            interpolated_feature_type,
                                            present_day_feature_geometry,
                                            valid_time=(interpolation_time, interpolation_end_time),
                                            reconstruction_plate_id=left_plate_id,
                                            conjugate_plate_id=right_plate_id)
                                    
                                    interpolated_features.append(interpolated_feature)
                                
                                if print_debug_output >= 2:
                                    print('    ...successful interpolation - number of interpolations {0}'.format(
                                            len(interpolated_polylines)))
    
    if print_debug_output >= 1:
        print('Generate isochron names...')
    
    # Sort the interpolated features by left plate id, right plate id and begin time (in that order).
    interpolated_features.sort(
            key = lambda feature: (_get_left_right_plate_ids(feature), feature.get_valid_time()[0]))
    
    # Produce interpolated feature names based on the sort order.
    pygplates_imported_version_string = str(pygplates.Version.get_imported_version())
    for interpolated_feature_index, interpolated_feature in enumerate(interpolated_features):
        interpolated_feature.set_name(
                '{0} Produced by pygplates revision {1}'.format(
                        interpolated_feature_index + 1,
                        pygplates_imported_version_string))
    
    return interpolated_features


if __name__ == "__main__":

    # Check the imported pygplates version.
    if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
        print('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
                os.path.basename(__file__), pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED),
            file=sys.stderr)
        sys.exit(1)


    __description__ = \
    """Generates interpolated isochrons using continental-oceanic boundaries, isochrons and present-day ridges.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -r rotations.rot -f cobs.gpml isochrons.gpml ridges.gpml -- output.gpml"""

    # The command-line parser.
    parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-r', '--rotation_filenames', type=str, nargs='+', required=True,
            metavar='rotation_filename', help='One or more rotation files.')
    parser.add_argument('-f', '--isochron_cob_ridge_filenames', type=str, nargs='+', required=True,
            metavar='isochron_cob_ridge_filename',
            help='One or more files containing continental-oceanic boundaries, isochrons and ridge boundaries.')
    parser.add_argument('-e', '--export', type=float,
            dest='export_time',
            metavar='reconstruction_time',
            help='Reconstruct the interpolated isochrons to the specified time and export to the '
                'output file (export formats ".shp", ".gmt" or ".xy"). If option not specified then '
                'the unreconstructed features are saved to the output file (feature formats ".gpml", ".dat", etc).')
    parser.add_argument('-s', '--scalar_coverages', type=str, nargs='+',
            metavar='scalar_type',
            help='Calculate scalar values at each point in each interpolated isochron in the output data. '
                'NOTE: Only applies when saving to ".gpml" format or saving/exporting to ".xy" format. '
                'Supported types of scalar data are "{0}", "{1}", "{2}", "{3}", "{4}" and "{5}". '
                'When saving/exporting to ".xy" format, the order of these types determines the order in which the '
                'associated scalar values are written to each line of ".xy" file (eg, specifing "{5} {1}" generates '
                'a line such as '
                '"10 20 40 5" where longitude=10, latitude=20, age=40Ma and spreading_rate=5Kms/My). '
                'When saving to ".gpml" format, the order is not important and the scalar values can be '
                'visualised/coloured in the latest GPlates. '
                '"{0}" is in the range [-1,1] (when saving to ".gpml") or [0,100] (when saving/exporting to ".xy") where '
                '-1 represents spreading only on the conjugate flank and +1 represents spreading only on the flank containing the isochron. '
                '"{1}" is the spreading rate relative to the mid-ocean ridge (taking into account asymmetry) in Kms/My. '
                '"{2}" is the spreading rate relative to the conjugate plate in Kms/My. '
                '"{3}" is an angle clockwise (East-wise) from North (0 to 180 degrees) of the lowest azimuth of the '
                'two opposite-pointing directions of spreading at the time of crustal accretion. '
                '"{4}" is an angle from 0 to 90 degrees that the small circle of spreading rotation pole deviates from '
                'the perpendicular line at the ridge point at the time of crustal accretion. And '
                '"{5}" is the interpolated isochron age in Ma.'.format(
                        SCALAR_COVERAGE_SPREADING_ASYMMETRY,
                        SCALAR_COVERAGE_SPREADING_RATE,
                        SCALAR_COVERAGE_FULL_SPREADING_RATE,
                        SCALAR_COVERAGE_SPREADING_DIRECTION,
                        SCALAR_COVERAGE_SPREADING_OBLIQUITY,
                        SCALAR_COVERAGE_AGE))
    parser.add_argument('-a', '--anchor', type=int, default=0,
            dest='anchor_plate_id',
            help='Anchor plate id used for reconstructing during export (if using "-e" option). Defaults to zero.')
    
    # Action to parse a list of plate pairs.
    class PlatePairsAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            # Should be an even number of numbers.
            if len(values) % 2 != 0:
                parser.error('the plate pairs must be supplied as *pairs* of numbers')
            
            try:
                # Convert strings to integers.
                integer_values = list(map(int, values))
                list_of_plate_pair_tuples = zip(integer_values[::2], integer_values[1::2])
            except ValueError:
                raise argparse.ArgumentTypeError("encountered a plate id that is not an integer")

            setattr(namespace, self.dest, list_of_plate_pair_tuples)
    
    parser.add_argument('-p', '--plate_pairs', nargs='+', action=PlatePairsAction,
            metavar='plate_id conjugate_plate_id',
            help='Zero or more conjugate plate pairs. Can be used to limit the output. '
                'If not specified then all plate pairs are processed.')

    parser.add_argument('-d', '--debug_output', type=int, default=0,
            dest='print_debug_output',
            help='Level at which to print debug output (0, 1, 2 or 3) - zero means no debug output - defaults to zero.')
        
    def parse_positive_number(value_string):
        try:
            value = float(value_string)
        except ValueError:
            raise argparse.ArgumentTypeError("%s is not a number" % value_string)
        
        if value < 0:
            raise argparse.ArgumentTypeError("%g is not a positive number" % value)
        
        return value
    
    # Can specify only one of '-i', '-l' or '-t'.
    interpolate_group = parser.add_mutually_exclusive_group()
    # Note: We don't document the '-i' option anymore since the '-l' option replaces it because the
    # former is in minutes and the latter is in degrees (which is what the other options use).
    # The '-i' option could eventually get deprecated.
    interpolate_group.add_argument('-i', '--interpolate_resolution_minutes',
            type=parse_positive_number,
            help=argparse.SUPPRESS)
    interpolate_group.add_argument('-l', '--interpolate_resolution_degrees',
            type=parse_positive_number,
            help='Interval spacing, in degrees, between interpolated isochrons (cannot be specified '
                'with "-t" option) - defaults to {0} if neither "-l" or "-t" options are specified.'.format(
                        DEFAULT_INTERPOLATE_RESOLUTION_DEGREES))
    interpolate_group.add_argument('-t', '--interpolate_times', type=float, nargs='+',
            metavar='INTERPOLATE_TIME',
            help='Interpolation times at which to generate isochrons (cannot be specified with "-l" option).')
    
    parser.add_argument('-m', '--tessellate_threshold_degrees', type=float,
            help='Each interpolated isochron is tessellated such that its adjacent points are separated no '
                'further than this threshold (in degrees) - if not specified then there is no tessellation.')
    
    parser.add_argument('--minimum_latitude_overlap_degrees', type=float,
            default=DEFAULT_MINIMUM_LATITUDE_OVERLAP_DEGREES,
            help='Required amount of latitude (where stage pole is North pole) overlap, in degrees, '
                'of two feature geometries for interpolation to occur - defaults to {0}.'.format(
                        DEFAULT_MINIMUM_LATITUDE_OVERLAP_DEGREES))
    parser.add_argument('--maximum_isochron_latitude_non_overlap_degrees', type=float,
            default=DEFAULT_MAXIMUM_ISOCHRON_LATITUDE_NON_OVERLAP_DEGREES,
            help='Allowed amount of non-overlapping latitude region, in degrees, when the older of two '
                'feature geometries being interpolated belongs to an isochron feature - defaults to {0}.'.format(
                        DEFAULT_MAXIMUM_ISOCHRON_LATITUDE_NON_OVERLAP_DEGREES))
    parser.add_argument('--maximum_cob_latitude_non_overlap_degrees', type=float,
            default=DEFAULT_MAXIMUM_COB_LATITUDE_NON_OVERLAP_DEGREES,
            help='Allowed amount of non-overlapping latitude region, in degrees, when the older of two '
                'feature geometries being interpolated belongs to a COB feature - defaults to {0}.'.format(
                        DEFAULT_MAXIMUM_COB_LATITUDE_NON_OVERLAP_DEGREES))
    parser.add_argument('--maximum_distance_threshold_degrees', type=float,
            default=DEFAULT_MAXIMUM_DISTANCE_THRESHOLD_DEGREES,
            help='Maximum distance, in degrees, between two feature geometries above which they will '
                'not be interpolated - defaults to {0}.'.format(
                        DEFAULT_MAXIMUM_DISTANCE_THRESHOLD_DEGREES))
    parser.add_argument('--maximum_age_difference', type=float,
            default=DEFAULT_MAXIMUM_AGE_DIFFERENCE,
            help='Maximum difference in age, in My, between two feature geometries above which they will '
                'not be interpolated - defaults to {0}.'.format(
                        DEFAULT_MAXIMUM_AGE_DIFFERENCE))
    parser.add_argument('--join_polylines_threshold_degrees', type=float,
            default=DEFAULT_JOIN_POLYLINES_THRESHOLD_DEGREES,
            help='Closeness threshold, in degrees, used to determine if two polyline feature geometries '
                'are adjacent and hence can be joined into a single geometry - defaults to {0}.'.format(
                        DEFAULT_JOIN_POLYLINES_THRESHOLD_DEGREES))
    
    
    parser.add_argument('output_filename', type=str, help='The output file containing the interpolated features.')
    
    
    # Parse command-line options.
    args = parser.parse_args()
    
    # Determine whether to interpolate using distance intervals or a list of interpolation times.
    if args.interpolate_times:
        interpolate = args.interpolate_times
    elif args.interpolate_resolution_minutes:
        # Convert from minutes to radians.
        interpolate = math.radians(args.interpolate_resolution_minutes / 60.0)
    elif args.interpolate_resolution_degrees:
        # Convert from minutes to radians.
        interpolate = math.radians(args.interpolate_resolution_degrees)
    else: # Default interpolation...
        interpolate = DEFAULT_INTERPOLATE_RESOLUTION_RADIANS
    
    # Keywords arguments for the interpolate_isochrons() function.
    interpolate_isochrons_kwargs = { }
    
    # Set some keyword arguments.
    interpolate_isochrons_kwargs['plate_pairs'] = args.plate_pairs
    interpolate_isochrons_kwargs['print_debug_output'] = args.print_debug_output
    interpolate_isochrons_kwargs['tessellate_threshold_radians'] = (
            math.radians(args.tessellate_threshold_degrees) # convert from degrees to radians
                if args.tessellate_threshold_degrees else None)
    
    # Scalar coverage output.
    if args.scalar_coverages:
        # Raise error if any unknown scalar coverage types.
        unknown_scalar_coverages = [scalar_coverage for scalar_coverage in args.scalar_coverages
                if scalar_coverage not in KNOWN_SCALAR_COVERAGES]
        if unknown_scalar_coverages:
            raise TypeError('Unknown scalar coverages {0}'.format(unknown_scalar_coverages))
        
        # Set scalar coverage keyword arguments.
        interpolate_isochrons_kwargs['output_scalar_spreading_asymmetry'] = (SCALAR_COVERAGE_SPREADING_ASYMMETRY in args.scalar_coverages)
        interpolate_isochrons_kwargs['output_scalar_spreading_rate'] = (SCALAR_COVERAGE_SPREADING_RATE in args.scalar_coverages)
        interpolate_isochrons_kwargs['output_scalar_full_spreading_rate'] = (SCALAR_COVERAGE_FULL_SPREADING_RATE in args.scalar_coverages)
        interpolate_isochrons_kwargs['output_scalar_spreading_direction'] = (SCALAR_COVERAGE_SPREADING_DIRECTION in args.scalar_coverages)
        interpolate_isochrons_kwargs['output_scalar_spreading_obliquity'] = (SCALAR_COVERAGE_SPREADING_OBLIQUITY in args.scalar_coverages)
        interpolate_isochrons_kwargs['output_scalar_age'] = (SCALAR_COVERAGE_AGE in args.scalar_coverages)
    
    if args.print_debug_output >= 1:
        print('Read rotations...')
    rotation_model = pygplates.RotationModel(args.rotation_filenames)
    
    # Generate a list of interpolated isochrons.
    output_features = interpolate_isochrons(
            rotation_model,
            args.isochron_cob_ridge_filenames,
            interpolate,
            math.radians(args.minimum_latitude_overlap_degrees), # convert from degrees to radians
            math.radians(args.maximum_isochron_latitude_non_overlap_degrees), # convert from degrees to radians
            math.radians(args.maximum_cob_latitude_non_overlap_degrees), # convert from degrees to radians
            math.radians(args.maximum_distance_threshold_degrees), # convert from degrees to radians
            args.maximum_age_difference,
            math.radians(args.join_polylines_threshold_degrees), # convert from degrees to radians
            # Unpack the keyword arguments dict into keyword arguments...
            **interpolate_isochrons_kwargs)
    
    if (args.scalar_coverages and
        os.path.splitext(args.output_filename)[1] == '.xy'):
        # We are outputting scalar coverages to '.xy' format.
        #
        # We handle this as a special case so we can write out the scalar values
        # after the xy (lat/lon) values.
        
        # We write the scalar types in the same order as they appear on the command-line.
        output_scalar_types = [pygplates.ScalarType.create_gpml(scalar_coverage)
                for scalar_coverage in args.scalar_coverages]
        write_coverage_features_to_xy_file(
                args.output_filename,
                output_features,
                output_scalar_types,
                rotation_model,
                args.anchor_plate_id,
                args.export_time,
                args.print_debug_output)
    
    elif args.export_time is not None:
        if args.print_debug_output >= 1:
            print('Exporting reconstructed features...')
        pygplates.reconstruct(output_features, rotation_model, args.output_filename, args.export_time, args.anchor_plate_id)
        
    else:
        if args.print_debug_output >= 1:
            print('Write features...')
        
        # Write the output features to disk.
        output_feature_collection = pygplates.FeatureCollection(output_features)
        output_feature_collection.write(args.output_filename)
