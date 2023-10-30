
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
# Various rotation utilities:                                                                    #
#   * Calculating stage rotations between consecutive finite rotations in a plate pair.          #
##################################################################################################


from __future__ import print_function
import sys
import math
import pygplates

#
# Python 2 and 3 compatibility.
#
# Iterating over a dict.
try:
    dict.iteritems
except AttributeError:
    # Python 3
    def itervalues(d):
        return iter(d.values())
    def iteritems(d):
        return iter(d.items())
    def listvalues(d):
        return list(d.values())
    def listitems(d):
        return list(d.items())
else:
    # Python 2
    def itervalues(d):
        return d.itervalues()
    def iteritems(d):
        return d.iteritems()
    def listvalues(d):
        return d.values()
    def listitems(d):
        return d.items()


# Required pygplates version.
# Let's not go back too far - require at least the 2nd public pygplates release.
PYGPLATES_VERSION_REQUIRED = pygplates.Version(18)


def extract_plate_pair_stage_rotations(
        rotation_feature_collections,
        plate_pair_filter=None):
    # Docstring in numpydoc format...
    """Calculate stage rotations between consecutive finite rotations in each specified plate pair.
    
    Parameters
    ----------
    rotation_feature_collections : sequence of (str, or sequence of pygplates.Feature, or pygplates.FeatureCollection, or pygplates.Feature)
        A sequence of rotation feature collections.
        Each collection in the sequence can be a rotation filename, or a sequence (eg, list of tuple) or features, or
        a feature collection, or even a single feature.
    plate_pair_filter : Filter function accepting accepting 3 arguments (moving_plate_id, fixed_plate_id, rotation_sequence), or sequence of 2-tuple (moving_plate_id, fixed_plate_id), optional
        Optional filtering of plate pairs to apply operation to.
        Filter function (callable) accepting 3 arguments (), or
        a sequence of (moving, fixed) plate pairs to limit operation to.
    
    Returns
    -------
    list of pygplates.FeatureCollection
        The modified feature collections.
        Returned list is same length as ``rotation_feature_collections``.
    
    Notes
    -----
    The results are returned as a list of pygplates.FeatureCollection (one per input rotation feature collection).
    
    Note that only the rotation features satisfying ``plate_pair_filter`` (if specified) are returned.
    So if a returned feature collection is empty then it means all of its features were filtered out.
    """
    
    if plate_pair_filter is None:
        plate_pair_filter = _all_filter
    # If caller specified a sequence of moving/fixed plate pairs then use them, otherwise it's a filter function (callable).
    elif hasattr(plate_pair_filter, '__iter__'):
        plate_pair_filter = _is_in_plate_pair_sequence(plate_pair_filter)
    # else ...'plate_pair_filter' is a callable...
    
    output_rotation_feature_collections = []
    for rotation_feature_collection in rotation_feature_collections:
        # Create an empty output rotation feature collection for each input rotation feature collection.
        # We'll only add output features when an input feature is modified.
        output_rotation_feature_collection = []
        output_rotation_feature_collections.append(output_rotation_feature_collection)
        
        for rotation_feature in rotation_feature_collection:
            # Get the rotation feature information.
            total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
            if not total_reconstruction_pole:
                # Not a rotation feature.
                continue
        
            fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
            # We're only interested in rotation features with matching moving/fixed plate IDs.
            if not plate_pair_filter(fixed_plate_id, moving_plate_id, rotation_sequence):
                continue
            
            # Clone the input feature before we start making modifications.
            # Otherwise we'll be modifying the input rotation feature (which the caller might not expect).
            output_rotation_feature = rotation_feature.clone()
            
            _, _, output_rotation_sequence = output_rotation_feature.get_total_reconstruction_pole()
            
            # Get the enabled rotation samples - ignore the disabled samples.
            output_enabled_rotation_samples = output_rotation_sequence.get_enabled_time_samples()
            if not output_enabled_rotation_samples:
                # No time samples are enabled.
                continue
            
            prev_finite_rotation = output_enabled_rotation_samples[0].get_value().get_finite_rotation()
            # Replace first finite rotation with identity rotation.
            # This is the stage rotation of first finite rotation (which has no previous finite rotation).
            output_enabled_rotation_samples[0].get_value().set_finite_rotation(pygplates.FiniteRotation())
            
            for rotation_sample_index in range(1, len(output_enabled_rotation_samples)):
                time_sample = output_enabled_rotation_samples[rotation_sample_index]
                finite_rotation = time_sample.get_value().get_finite_rotation()
                
                # The finite rotation at current time tc is composed of finite rotation at
                # previous (younger) time tp and stage rotation between them:
                #   R(0->tc) = R(tp->tc) * R(0->tp)
                # The stage rotation from tp->tc:
                #   R(tp->tc) = R(0->tc) * inverse[R(0->tp)]
                stage_rotation = finite_rotation * prev_finite_rotation.get_inverse()
                
                # Replace current finite rotation with stage rotation.
                time_sample.get_value().set_finite_rotation(stage_rotation)
                
                prev_finite_rotation = finite_rotation
            
            output_rotation_feature_collection.append(output_rotation_feature)
    
    # Return our output feature collections as a list of pygplates.FeatureCollection.
    return [pygplates.FeatureCollection(rotation_feature_collection)
        for rotation_feature_collection in output_rotation_feature_collections]


def _all_filter(*args, **kwargs):
    """Predicate filter than always returns True (and accepts any number of parameters)."""
    return True


def _is_in_plate_pair_sequence(plate_pairs):
    """Returns a predicate filter function that returns True if a moving/fixed plate pair is in 'plate_pairs'."""
    plate_pairs = list(plate_pairs)
    def filter_is_in_plate_pair_sequence(fixed_plate_id, moving_plate_id, rotation_sequence):
        return (moving_plate_id, fixed_plate_id) in plate_pairs
    return filter_is_in_plate_pair_sequence


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
            
            if threshold_rotation_accuracy_degrees <= 0 or threshold_rotation_accuracy_degrees > 90:
                parser.error('rotation threshold must be in the range (0, 90]')
            if threshold_time_interval <= 0:
                parser.error('threshold time interval must be positive')
            
            setattr(namespace, self.dest, (threshold_rotation_accuracy_degrees, threshold_time_interval))
    
    
    def main():
    
        __description__ = \
    """Calculate stage rotations between consecutive finite rotations in plate pairs.
    
    The results are written back to the input rotation files unless an output filename prefix is provided.
    
    Note that only the rotation features satisfying the 'plate_pairs' option (if specified) are written back.
    So if an output rotation file is empty then it means all of its rotation features were filtered out.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -p 701 70 -o stage_ -- rotations.rot
     """

        # The command-line parser.
        parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
        
        # Action to parse a list of plate pairs.
        class PlatePairsAction(argparse.Action):
            def __call__(self, parser, namespace, values, option_string=None):
                # Should be an even number of numbers.
                if len(values) % 2 != 0:
                    parser.error('the plate pairs must be supplied as *pairs* of numbers')
                
                try:
                    # Convert strings to integers.
                    integer_values = map(int, values)
                    list_of_plate_pair_tuples = zip(integer_values[::2], integer_values[1::2])
                except ValueError:
                    raise argparse.ArgumentTypeError("encountered a plate id that is not an integer")

                setattr(namespace, self.dest, list_of_plate_pair_tuples)
        
        parser.add_argument('-p', '--plate_pairs', nargs='+', action=PlatePairsAction,
                metavar='moving_plate_id fixed_plate_id',
                help='One or more moving/fixed plate pairs to limit operation to. '
                     'If not specified then defaults to all plate pairs.')
        
        parser.add_argument('-o', '--output_filename_prefix', type=str,
                metavar='output_filename_prefix',
                help='Optional output filename prefix. If one is provided then an output rotation file '
                    'is created for each input rotation file by prefixing the input filenames. '
                    'If no filename prefix is provided then the input files are overwritten.')
        
        parser.add_argument('input_rotation_filenames', type=str, nargs='+',
                metavar='input_rotation_filename',
                help='One or more rotation files to perform an operation on.')
        
        # Parse command-line options.
        args = parser.parse_args()
        
        # Read the input rotation feature collections.
        input_rotation_feature_collections = [pygplates.FeatureCollection(input_rotation_filename)
                for input_rotation_filename in args.input_rotation_filenames]
        
        # Operation.
        output_rotation_feature_collections = extract_plate_pair_stage_rotations(
                input_rotation_feature_collections,
                args.plate_pairs)
        
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
