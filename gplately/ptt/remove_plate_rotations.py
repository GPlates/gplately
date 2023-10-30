
"""
    Copyright (C) 2019 The University of Sydney, Australia
    
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


##################################################################################################
# Remove one or more plate IDs from a rotation model (consisting of one or more rotation files). #
##################################################################################################


from __future__ import print_function
import sys
import math
import pygplates


# Required pygplates version.
# Need pygplates.RotationModel to clone rotation features by default (revision 12).
# Note that since revision 25 cloning is no longer necessary.
PYGPLATES_VERSION_REQUIRED = pygplates.Version(12)


def remove_plates(
        rotation_feature_collections,
        plate_ids,
        accuracy_parameters=None):
    # Docstring in numpydoc format...
    """Remove one or more plate IDs from a rotation model (consisting of one or more rotation feature collections).
    
    Any rotations with a fixed plate referencing one of the removed plates will be adjusted such that
    the rotation model effectively remains unchanged.
    
    Optional accuracy threshold parameters can be specified to ensure the rotation model after removing
    plate rotations is very similar to the rotation model before removal.
    
    The results are returned as a list of pygplates.FeatureCollection (one per input rotation feature collection).
    
    Ensure you specify all input rotation feature collections that contain the plate IDs to be removed (either as a moving or fixed plate ID).
    
    Parameters
    ----------
    rotation_feature_collections : sequence of (str, or sequence of pygplates.Feature, or pygplates.FeatureCollection, or pygplates.Feature)
        A sequence of rotation feature collections.
        Each collection in the sequence can be a rotation filename, or a sequence (eg, list of tuple) or features, or
        a feature collection, or even a single feature.
    plate_ids : sequence of int
        Plate IDs to remove from rotation model.
    accuracy_parameters: tuple of (float, float, bool), optional
        First parameter is the threshold rotation accuracy (in degrees), the second parameter is the threshold time interval and
        the third parameter is True if insert poles should have times that are integer multiples of the threshold time interval.
        The first parameter is used to compare the latitude, longitude and angle of two rotations before and after removing a plate rotation.
        If any of those three parameters differ by more than the rotation accuracy (in degrees) then
        samples at times mid-way between samples are inserted to ensure before/after accuracy of rotations.
        This mid-way adaptive bisection is repeated (when there is inaccuracy) until the interval between samples
        becomes smaller than the second parameter (threshold time interval).
        Rotation threshold is in degrees and threshold time interval is in millions of years.
    
    Returns
    -------
    list of pygplates.FeatureCollection
        The (potentially) modified feature collections.
        Returned list is same length as ``rotation_feature_collections``.
    """
    
    # Convert each feature collection into a list of features so we can more easily remove features
    # and insert features at arbitrary locations within each feature collection (for example removing
    # a plate sequence and replacing it with a sequence with the same plate ID).
    rotation_feature_collections = [list(pygplates.FeatureCollection(rotation_feature_collection))
        for rotation_feature_collection in rotation_feature_collections]
    
    # Iterate over the plates to be removed and remove each one separately.
    for remove_plate_id in plate_ids:
        
        # Rotation model before any modifications to rotation features.
        #
        # Previously we only created one instance of RotationModel in this function (to serve all removed plates).
        # However we need to create a new RotationModel instance for each plate ID being removed, in
        # case any rotation sequence referencing the removed plate (as its fixed plate) has a time range
        # further into the past than the removed plate sequence. This is a subtle issue to do with
        # not having a plate circuit *through the removed plate* for times older than supported by
        # the removed plate sequence - in this case RotationModel would just return an identity rotation.
        # However if we create a new RotationModel *without* the removed plate sequence then we
        # avoid this issue altogether.
        #
        # Note that RotationModel clones the current rotation features (by default), so any subsequent
        # feature modifications (in this loop iteration) should not affect it.
        # However the RotationModel in the next loop iteration will be affected of course.
        # UPDATE: Since pygplates revision 25 cloning is no longer necessary (and has been deprecated).
        rotation_model = pygplates.RotationModel(rotation_feature_collections)
        
        # Rotation sequences with the current remove plate ID as the *moving* plate ID.
        # Each sequence will have a different *fixed* plate ID.
        remove_plate_sequences = []
        for rotation_feature_collection in rotation_feature_collections:
            rotation_feature_index = 0
            while rotation_feature_index < len(rotation_feature_collection):
                rotation_feature = rotation_feature_collection[rotation_feature_index]
                total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
                if total_reconstruction_pole:
                    fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                    if moving_plate_id == remove_plate_id:
                        sample_times = [pygplates.GeoTimeInstant(sample.get_time())
                            for sample in rotation_sequence.get_enabled_time_samples()]
                        if sample_times:
                            remove_plate_sequences.append(
                                (fixed_plate_id, sample_times))
                        # Remove plate sequences whose moving plate is the current remove plate.
                        # Note that this won't affect 'rotation_model' (since it used a cloned version of all features).
                        del rotation_feature_collection[rotation_feature_index]
                        rotation_feature_index -= 1
                
                rotation_feature_index += 1
        
        # Sort the remove plate sequences in time order.
        # This helps out below, to find the max sample time over the removed sequences (all having same *moving* plate ID).
        #
        # Easiest way to do this is to sort based on the first time sample in each sequence
        # (since each sequence should already be sorted internally).
        remove_plate_sequences.sort(key = lambda sequence: sequence[1][0])
        
        # Find those sequences that need adjustment due to the plate removal.
        # These are sequences whose *fixed* plate is the plate currently being removed.
        for rotation_feature_collection in rotation_feature_collections:
            rotation_feature_index = 0
            while rotation_feature_index < len(rotation_feature_collection):
                rotation_feature = rotation_feature_collection[rotation_feature_index]
                total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
                if total_reconstruction_pole:
                    fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                    if fixed_plate_id == remove_plate_id:
                        child_remove_plate_id = moving_plate_id
                        child_remove_plate_rotation_feature = rotation_feature
                        child_remove_plate_samples = rotation_sequence.get_time_samples()
                        
                        child_remove_plate_sample_times = [pygplates.GeoTimeInstant(sample.get_time())
                            for sample in child_remove_plate_samples]
                        child_remove_plate_min_sample_time = child_remove_plate_sample_times[0]
                        child_remove_plate_max_sample_time = child_remove_plate_sample_times[-1]
                        
                        # Iterate over the removed sequences whose moving plate matched the current plate being removed.
                        for remove_plate_sequence_index, (parent_remove_plate_id, remove_plate_sample_times) in enumerate(remove_plate_sequences):
                            remove_plate_min_sample_time = remove_plate_sample_times[0]
                            remove_plate_max_sample_time = remove_plate_sample_times[-1]
                            
                            # Find the time overlap of the removed sequence and the (child) sequence requiring modification.
                            min_sample_time = max(remove_plate_min_sample_time, child_remove_plate_min_sample_time)
                            if remove_plate_sequence_index == len(remove_plate_sequences) - 1:
                                # We want the last remove plate sequence to go back to the child max sample time.
                                # If it doesn't go that far back then we will artificially extend the remove plate sequence that far back.
                                max_sample_time = child_remove_plate_max_sample_time
                            else:
                                max_sample_time = min(remove_plate_max_sample_time, child_remove_plate_max_sample_time)
                            
                            # Note that the remove sequences are ordered by time (ie, first sequence should start at 0Ma, etc).
                            # The two sequences must overlap.
                            # Note that this excludes the case where the min of one sequence equals the max of the other (or max and min).
                            if min_sample_time < max_sample_time:
                                sample_times = []
                                # Find those sample times of the child sequence within the overlap range.
                                for child_remove_plate_sample_time in child_remove_plate_sample_times:
                                    if (child_remove_plate_sample_time >= min_sample_time and
                                        child_remove_plate_sample_time <= max_sample_time):
                                        sample_times.append(child_remove_plate_sample_time)
                                # Find those sample times of the remove sequence within the overlap range.
                                # Also avoiding duplicating sample times (times already in the child sequence).
                                for remove_plate_sample_time in remove_plate_sample_times:
                                    # Only add the sample time if it's not already in the list.
                                    if (remove_plate_sample_time not in child_remove_plate_sample_times and
                                        remove_plate_sample_time >= min_sample_time and
                                        remove_plate_sample_time <= max_sample_time):
                                        sample_times.append(remove_plate_sample_time)
                                # Need to sort the sample times (since they're likely interleaved between remove and child sequences).
                                sample_times.sort()
                                
                                # Gather the rotation samples from the child's moving plate to the removed plate's fixed plate.
                                parent_to_child_rotation_samples = _merge_rotation_samples(
                                    rotation_model,
                                    child_remove_plate_id,
                                    remove_plate_id,
                                    parent_remove_plate_id,
                                    child_remove_plate_samples,
                                    child_remove_plate_sample_times,
                                    sample_times,
                                    remove_plate_max_sample_time)
                                
                                # Insert new samples at times where the difference between original and new rotation models exceeds a threshold.
                                if accuracy_parameters is not None:
                                    threshold_rotation_accuracy_degrees, threshold_time_interval, use_uniform_accuracy_times = accuracy_parameters
                                    _ensure_sequence_accuracy(
                                        rotation_model,
                                        parent_to_child_rotation_samples,
                                        child_remove_plate_id,
                                        remove_plate_id,
                                        parent_remove_plate_id,
                                        remove_plate_max_sample_time,
                                        threshold_rotation_accuracy_degrees,
                                        threshold_time_interval,
                                        use_uniform_accuracy_times)
                                
                                # Create a new rotation sequence.
                                parent_to_child_rotation_feature = pygplates.Feature.create_total_reconstruction_sequence(
                                    parent_remove_plate_id,
                                    child_remove_plate_id,
                                    pygplates.GpmlIrregularSampling(parent_to_child_rotation_samples),
                                    child_remove_plate_rotation_feature.get_name(None),         # Note: specifying None avoids a pygplates crash in revs < 20
                                    child_remove_plate_rotation_feature.get_description(None))  # Note: specifying None avoids a pygplates crash in revs < 20
                                
                                # Insert the new rotation feature to the current location in the feature collection.
                                # This is better than adding to the end of the collection and thus reordering the order of rotation sequences
                                # in the output collection/file (making it harder to visually find it in a text editor).
                                # Also note that this won't affect 'rotation_model' (since it used a cloned version of all features).
                                rotation_feature_collection.insert(rotation_feature_index, parent_to_child_rotation_feature)
                                rotation_feature_index += 1
                        
                        # The original rotation feature will no longer be needed because we remove plate sequences
                        # whose fixed plate is the current remove plate.
                        # We would have added one or more sequences above to replace it though.
                        # Also note that this won't affect 'rotation_model' (since it used a cloned version of all features).
                        del rotation_feature_collection[rotation_feature_index]
                        rotation_feature_index -= 1
                        
                rotation_feature_index += 1
        
        # Note that we don't join rotation sequences having the same moving/fixed plates.
        # However they will show up as a 'duplicate geo-time' warning when loading into GPlates.
        # TODO: Remove duplicate geo-times and join the offending rotation sequences.
        #
        # Details: It's possible that a sequence having a crossover (really two sequences with same moving
        # plate but different fixed plates) can have one of its fixed plates removed and hence replaced
        # by the fixed plate of the removed sequence. In this situation the original crossover sequence
        # (really two sequences) could now have the same fixed plate ID, and since it also has the
        # same moving plate ID it should really be one sequence.
    
    # Return our (potentially) modified feature collections as a list of pygplates.FeatureCollection.
    return [pygplates.FeatureCollection(rotation_feature_collection)
        for rotation_feature_collection in rotation_feature_collections]


def _merge_rotation_samples(
        rotation_model,
        child_remove_plate_id,
        remove_plate_id,
        parent_remove_plate_id,
        child_remove_plate_samples,
        child_remove_plate_sample_times,
        sample_times,
        remove_plate_max_sample_time):
    """Gather the rotation samples from the child's moving plate to the removed plate's fixed plate."""
    
    # Gather the rotation samples from the child's moving plate to the removed plate's fixed plate.
    parent_to_child_rotation_samples = []
    for sample_time in sample_times:
        # Rotation from parent to remove-plate is now replaced by rotation from
        # (parent of remove-plate) to (child of remove-plate).
        #
        #   R(0->t,parent_plate->child_plate) = R(0->t,parent_plate->remove_plate) * R(0->t,remove_plate->child_plate)
        #
        # Also note that below we set fixed plate as the 'anchor_plate_id' argument to pygplates.RotationModel.get_rotation(),
        # not the 'fixed_plate_id' argument, in case there is no plate circuit path to anchor plate 000.
        # This also means the user doesn't have to load all rotations in the model,
        # only those that have the remove plate IDs as a moving or fixed plate.
        
        if sample_time > remove_plate_max_sample_time:
            # The time span of the (oldest) removed plate sequence is too short, so extend its oldest sample further into
            # the past (ie, assume a constant rotation). We do this by calculating R(0->t,parent_plate->remove_plate)
            # at the max sample time of remove plate sequence instead of the current sample time.
            #
            # Note that the remove sequences are ordered by time (ie, first sequence should start at 0Ma, etc).
            # So 'max_sample_time' should be the oldest time of the oldest removed plate sequence.
            
            # R(0->t,parent_plate->remove_plate)
            parent_to_remove_rotation = rotation_model.get_rotation(
                # Note the time is 'remove_plate_max_sample_time' and not 'sample_time'...
                remove_plate_max_sample_time, remove_plate_id, anchor_plate_id=parent_remove_plate_id)
            # R(0->t,remove_plate->child_plate)
            remove_to_child_rotation = rotation_model.get_rotation(
                sample_time, child_remove_plate_id, anchor_plate_id=remove_plate_id)
            
            parent_to_child_rotation = parent_to_remove_rotation * remove_to_child_rotation
        else:
            # Note that here we don't actually need to compose rotations as in the above equation because both rotations
            # are at the same (sample) time so we can just get pygplates.RotationModel.get_rotation() to compose them for us.
            #
            # R(0->t,parent_plate->child_plate)
            parent_to_child_rotation = rotation_model.get_rotation(
                sample_time, child_remove_plate_id, anchor_plate_id=parent_remove_plate_id)
        
        # If the sample time corresponds to an existing sample then use its description,
        # otherwise create a new sample description noting that the new sample is due to
        # the removal of a specific fixed plate.
        if sample_time in child_remove_plate_sample_times:
            child_remove_plate_sample = child_remove_plate_samples[child_remove_plate_sample_times.index(sample_time)]
            sample_description = child_remove_plate_sample.get_description()
        else:
            sample_description = 'Removed fixed plate {0}'.format(remove_plate_id)
        
        parent_to_child_rotation_samples.append(
            pygplates.GpmlTimeSample(
                pygplates.GpmlFiniteRotation(parent_to_child_rotation),
                sample_time,
                sample_description))
    
    return parent_to_child_rotation_samples


def _ensure_sequence_accuracy(
        rotation_model,
        parent_to_child_rotation_samples,
        child_remove_plate_id,
        remove_plate_id,
        parent_remove_plate_id,
        remove_plate_max_sample_time,
        threshold_rotation_accuracy_degrees,
        threshold_time_interval,
        insert_poles_at_integer_multiples_of_time_interval):
    """Insert new samples at times where the difference between original and new rotation models exceeds a threshold."""
    
    num_original_samples = len(parent_to_child_rotation_samples)

    sample_pair_stack = []
    
    # Add the stage rotation intervals to the stack for later processing.
    for sample_index in range(num_original_samples-1):
        sample1, sample2 = parent_to_child_rotation_samples[sample_index], parent_to_child_rotation_samples[sample_index + 1]
        sample_time1, sample_time2 = sample1.get_time(), sample2.get_time()
        if pygplates.GeoTimeInstant(sample_time2 - sample_time1) > threshold_time_interval:
            sample_pair_stack.append((sample1, sample2))
    
    # Process the stage rotation intervals on the stack until it is empty.
    while sample_pair_stack:
        sample1, sample2 = sample_pair_stack.pop()
        sample_time1, sample_time2 = sample1.get_time(), sample2.get_time()
        
        mid_sample_time = 0.5 * (sample_time1 + sample_time2)
        
        if insert_poles_at_integer_multiples_of_time_interval:
            # Round to the nearest uniformly spaced interval.
            interpolated_sample_time = threshold_time_interval * math.floor((mid_sample_time / threshold_time_interval) + 0.5)
            if interpolated_sample_time > mid_sample_time:
                if interpolated_sample_time >= pygplates.GeoTimeInstant(sample_time2):
                    # We rounded up and the time was greater-or-equal to the end sample time so subtract one time interval.
                    # This is guaranteed to remain within the start/end range since that range should exceed the time interval.
                    interpolated_sample_time -= threshold_time_interval
            else:
                if interpolated_sample_time <= pygplates.GeoTimeInstant(sample_time1):
                    # We rounded down and the time was less-or-equal to the start sample time so add one time interval.
                    # This is guaranteed to remain within the start/end range since that range should exceed the time interval.
                    interpolated_sample_time += threshold_time_interval
                
        else:
            # Just use the sample midway between 'sample1' and 'sample2'.
            interpolated_sample_time = mid_sample_time
        
        interpolated_sample = _create_accurate_sample(
            rotation_model,
            interpolated_sample_time,
            sample1,
            sample2,
            child_remove_plate_id,
            remove_plate_id,
            parent_remove_plate_id,
            remove_plate_max_sample_time,
            threshold_rotation_accuracy_degrees)
        
        if interpolated_sample:
            parent_to_child_rotation_samples.append(interpolated_sample)
            
            # Recurse if the time interval between start sample and the interpolated sample exceeds threshold interval.
            if pygplates.GeoTimeInstant(interpolated_sample_time - sample_time1) > threshold_time_interval:
                sample_pair_stack.append((sample1, interpolated_sample))
            
            # Recurse if the time interval between the interpolated sample and end sample exceeds threshold interval.
            if pygplates.GeoTimeInstant(sample_time2 - interpolated_sample_time) > threshold_time_interval:
                sample_pair_stack.append((interpolated_sample, sample2))
    
    # Sort the sample by time (if we added any new samples).
    if len(parent_to_child_rotation_samples) > num_original_samples:
        parent_to_child_rotation_samples.sort(key = lambda sample: sample.get_time())


def _create_accurate_sample(
        rotation_model,
        interpolated_sample_time,
        sample1,
        sample2,
        child_remove_plate_id,
        remove_plate_id,
        parent_remove_plate_id,
        remove_plate_max_sample_time,
        threshold_rotation_accuracy_degrees):
    """Create a new accurate interpolated sample if the difference between original and new rotation models exceeds a threshold, otherwise None."""
    
    # Find the *original* rotation from parent plate to child plate (through removed plate).
    #
    # R(0->t,parent_plate->remove_plate)
    parent_to_remove_rotation = rotation_model.get_rotation(
        # If the time span of the (oldest) removed plate sequence is too short then extend its oldest rotation to the interpolated-sample time...
        min(interpolated_sample_time, remove_plate_max_sample_time), remove_plate_id, anchor_plate_id=parent_remove_plate_id)
    # R(0->t,remove_plate->child_plate)
    remove_to_child_rotation = rotation_model.get_rotation(
        interpolated_sample_time, child_remove_plate_id, anchor_plate_id=remove_plate_id)
    original_parent_to_child_rotation = parent_to_remove_rotation * remove_to_child_rotation
    
    # Find the *new* rotation from parent plate to child plate (through removed plate).
    #
    # This interpolates the newly calculated samples (that go directly from parent to child, ie, not via removed plate).
    new_parent_to_child_rotation = pygplates.FiniteRotation.interpolate(
        sample1.get_value().get_finite_rotation(),
        sample2.get_value().get_finite_rotation(),
        sample1.get_time(),
        sample2.get_time(),
        interpolated_sample_time)
    
    # If original and new parent-to-child rotations differ too much then add a new (accurate) sample at the interpolated-sample time.
    interpolated_sample = None
    if not pygplates.FiniteRotation.are_equal(
        original_parent_to_child_rotation,
        new_parent_to_child_rotation,
        threshold_rotation_accuracy_degrees):
        
        interpolated_sample = pygplates.GpmlTimeSample(
            pygplates.GpmlFiniteRotation(original_parent_to_child_rotation),
            interpolated_sample_time,
            'Inserted pole to improve accuracy after removing fixed plate {0}'.format(remove_plate_id))
    
    return interpolated_sample


if __name__ == '__main__':
    
    import os.path
    
    
    # Check the imported pygplates version.
    if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
        print('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
                os.path.basename(__file__), pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED),
            file=sys.stderr)
        sys.exit(1)
    
    
    import argparse

    # Action to parse a tuple of accuracy parameters.
    class ArgParseAccuracyAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            # Need two numbers (rotation threshold and threshold time interval).
            if len(values) != 2:
                parser.error('accuracy must be specified as two numbers (rotation threshold and threshold time interval)')
            
            try:
                # Convert strings to float.
                threshold_rotation_accuracy_degrees = float(values[0])
                threshold_time_interval = float(values[1])
            except ValueError:
                raise argparse.ArgumentTypeError("encountered a rotation threshold and threshold time interval that is not a number")
            
            if threshold_rotation_accuracy_degrees < 0 or threshold_rotation_accuracy_degrees > 90:
                parser.error('rotation threshold must be in the range [0, 90]')
            if threshold_time_interval <= 0:
                parser.error('threshold time interval must be positive')
            
            setattr(namespace, self.dest, (threshold_rotation_accuracy_degrees, threshold_time_interval))
    
    
    def main():
    
        __description__ = \
    """Remove one or more plate IDs from a rotation model (consisting of one or more rotation files).
    
    Any rotations with a fixed plate referencing one of the removed plates will be adjusted such that
    the rotation model effectively remains unchanged.
    
    Optional accuracy threshold parameters can be specified to ensure the rotation model after removing
    plate rotations is very similar to the rotation model before removal.
    
    Ensure you specify all input rotation files that contain the plate IDs to be removed (either as a moving or fixed plate ID).
    
    The results are written back to the input rotation files unless an output filename prefix is provided.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -p 70 4 3 1 -o removed_ref_frames_ -- rotations.rot
     """

        # The command-line parser.
        parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
        
        parser.add_argument('-p', '--plates', type=int, nargs='+', required=True,
                metavar='remove_plate_ID',
                dest='plate_ids',
                help='Plate IDs of one or more plates to remove.')
        
        parser.add_argument(
            '-a', '--accuracy', nargs=2, action=ArgParseAccuracyAction,
            metavar=('threshold_rotation_accuracy_degrees', 'threshold_time_interval_My'),
            help='Optional accuracy parameters. '
                 'If specified then the first parameter is the threshold rotation accuracy (in degrees) and '
                 'the second parameter is the threshold time interval. '
                 'The first parameter is used to compare the latitude, longitude and angle of two rotations before and '
                 'after removing a plate rotation. If any of those three parameters differ by more than the rotation accuracy (in degrees) then '
                 'samples at times mid-way between samples are inserted to ensure before/after accuracy of rotations. '
                 'This mid-way adaptive bisection is repeated (when there is inaccuracy) until the interval between samples '
                 'becomes smaller than the second parameter (threshold time interval). '
                 'Rotation threshold is in degrees and threshold time interval is in millions of years.')
        
        parser.add_argument('-u', '--use_uniform_accuracy_times', action="store_true",
                help='If specified then rotation poles inserted for accuracy (according to "-a" option) will be restricted to times '
                     'that are integer multiples of the threshold time interval (specified in the "-a" option).')
        
        parser.add_argument('-o', '--output_filename_prefix', type=str,
                metavar='output_filename_prefix',
                help='Optional output filename prefix. If one is provided then an output rotation file '
                    'is created for each input rotation file by prefixing the input filenames. '
                    'If no filename prefix is provided then the input files are overwritten.')
        
        parser.add_argument('input_rotation_filenames', type=str, nargs='+',
                metavar='input_rotation_filename',
                help='One or more rotation files of a rotation model.')
        
        # Parse command-line options.
        args = parser.parse_args()
        
        # Initialise accuracy parameters (if used).
        accuracy_parameters = None
        if args.accuracy:
            accuracy_parameters = args.accuracy[0], args.accuracy[1], args.use_uniform_accuracy_times
        
        # Read the input rotation feature collections.
        input_rotation_feature_collections = [pygplates.FeatureCollection(input_rotation_filename)
                for input_rotation_filename in args.input_rotation_filenames]
        
        # Remove plate rotations.
        output_rotation_feature_collections = remove_plates(
                input_rotation_feature_collections,
                args.plate_ids,
                accuracy_parameters)
        
        # Write the modified rotation feature collections to disk.
        for rotation_feature_collection_index in range(len(output_rotation_feature_collections)):
            output_rotation_feature_collection = output_rotation_feature_collections[rotation_feature_collection_index]

            # Each output filename is the input filename with an optional prefix prepended.
            input_rotation_filename = args.input_rotation_filenames[rotation_feature_collection_index]
            if args.output_filename_prefix:
                dir, file_basename = os.path.split(input_rotation_filename)
                output_rotation_filename = os.path.join(dir, '{0}{1}'.format(args.output_filename_prefix, file_basename))
            else:
                output_rotation_filename = input_rotation_filename
            
            output_rotation_feature_collection.write(output_rotation_filename)
        
        sys.exit(0)
    
    import traceback
    
    try:
        main()
        sys.exit(0)
    except Exception as exc:
        print('ERROR: {0}'.format(exc), file=sys.stderr)
        # Uncomment this to print traceback to location of raised exception.
        # traceback.print_exc()
        sys.exit(1)
