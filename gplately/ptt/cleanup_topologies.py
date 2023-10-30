
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
# Clean up topologies, including:                                                                #
#   * Removing any regular features not referenced by topologies.                                #
#   * Restricting the time periods of referenced features to match the referencing topologies.   #
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
# Need ability to query topological sections.
PYGPLATES_VERSION_REQUIRED = pygplates.Version(21)


def remove_features_not_referenced_by_topologies(
        feature_collections,
        restrict_referenced_feature_time_periods=False,
        removed_features_collections=None):
    # Docstring in numpydoc format...
    """Remove any regular features not referenced by topological features.
    
    The results are returned as a list of pygplates.FeatureCollection (one per input feature collection).
    
    The input feature collections should contain the topological and regular features that make up the topological model.
    Ensure you at least specify all topological features in the topological model, otherwise regular features
    referenced by the missing topologies will get removed.
    
    Parameters
    ----------
    feature_collections : sequence of (str, or sequence of pygplates.Feature, or pygplates.FeatureCollection, or pygplates.Feature)
        A sequence of feature collections containing both the topological features and the regular features they reference.
        Each collection in the sequence can be a filename, or a sequence (eg, list of tuple) or features, or
        a feature collection, or even a single feature.
    restrict_referenced_feature_time_periods: bool
        Whether to restrict the time periods of features referenced by topologies such that they are
        limited by the time periods of the referencing topologies.
    removed_features_collections: list
        If an empty list is provided then it will be filled with one feature collection for each input feature collection.
        And each of these will contain any removed features.
    
    Returns
    -------
    list of pygplates.FeatureCollection
        The (potentially) modified feature collections.
        Returned list is same length as ``feature_collections``.
    """
    
    # Convert each feature collection into a list of features so we can more easily remove features
    # and insert features at arbitrary locations within each feature collection
    # (for example removing an unreferenced regular feature).
    feature_collections = [list(pygplates.FeatureCollection(feature_collection))
        for feature_collection in feature_collections]
    
    # Set of feature IDs of all topological polygons and networks.
    # Note that we only keep topological lines if they are referenced by a topological polygon or network.
    feature_ids_of_topological_polygons_and_networks = set()
    # Set of feature IDs referenced by all topological polygons and networks.
    feature_ids_referenced_by_topological_polygons_and_networks = set()
    # Set of feature IDs referenced by topological lines that are referenced by topological polygons and networks.
    feature_ids_referenced_by_topological_lines_referenced_by_topological_polygons_and_networks = set()
    
    # The features referenced by topological polygons and networks.
    topological_polygon_and_network_references = []
    # The features referenced by topological lines.
    topological_line_references = dict()
    
    # All features mapped by feature ID.
    all_features = dict()
    
    # Find all topological features and their references to regular features.
    topological_reference_visitor = _TopologicalReferenceVisitor()
    for feature_collection in feature_collections:
        for feature in feature_collection:
            feature_id = feature.get_feature_id()
            
            # Enable any feature to be looked up using its feature ID.
            all_features[feature_id] = feature
            
            # See if the current feature has a topological geometry and (if so) find the features it references.
            topology_type, topological_references = topological_reference_visitor.visit_feature(feature)
            
            # For topological lines, we'll just keep track of their references for later since we don't yet know
            # which topological lines (if any) will in turn be referenced by topological polygons and  networks.
            if topology_type == pygplates.GpmlTopologicalLine:
                topological_line_references[feature_id] = (feature, topological_references)
            # For topological polygons and networks, we know they get included straight away.
            elif (topology_type == pygplates.GpmlTopologicalPolygon or
                topology_type == pygplates.GpmlTopologicalNetwork):
                topological_polygon_and_network_references.append((feature, topological_references))
                # Add the current topological polygon or network.
                feature_ids_of_topological_polygons_and_networks.add(feature_id)
            # else it's not a topological feature, so ignore it.
    
    # Referenced time periods of referenced feature (non-topological and topological lines).
    time_periods_of_referenced_non_topological_features = dict()
    time_periods_of_referenced_topological_line_features = dict()
    
    # Set of feature IDs of topological lines that are referenced by topological polygons and networks.
    feature_ids_of_topological_lines_referenced_by_topological_polygons_and_networks = set()
    
    # Start with begin time as -inf (so that it can only get larger with max() function) and
    # end time as +inf (so that it can only get smaller with min() function).
    uninitialised_time_period = float('-inf'), float('inf')  # begin_time, end_time
    
    # Find features directly referenced by topological polygons and networks.
    for topological_polygon_and_network_feature, topological_polygon_and_network_reference in topological_polygon_and_network_references:
        topological_polygon_and_network_begin_time, topological_polygon_and_network_end_time = topological_polygon_and_network_feature.get_valid_time()
        for (reference_begin_time, reference_end_time), feature_ids_referenced_in_time_period in topological_polygon_and_network_reference:
            # Restrict the current referenced time period to that of the topological polygon (or network) valid time period.
            if reference_begin_time > topological_polygon_and_network_begin_time:
                reference_begin_time = topological_polygon_and_network_begin_time
            if reference_end_time < topological_polygon_and_network_end_time:
                reference_end_time = topological_polygon_and_network_end_time
            
            # If the referenced time period and the topological polygon/network valid time overlap.
            if reference_begin_time >= reference_end_time:
                # Add the features referenced by the current topological polygon or network.
                feature_ids_referenced_by_topological_polygons_and_networks.update(feature_ids_referenced_in_time_period)
            
                # Iterate over referenced feature IDs for the currently reference time period.
                for referenced_feature_id in feature_ids_referenced_in_time_period:
                    # Point to time periods of either topological line features or non-topological features.
                    if referenced_feature_id in topological_line_references:
                        time_periods_of_referenced_features = time_periods_of_referenced_topological_line_features
                        feature_ids_of_topological_lines_referenced_by_topological_polygons_and_networks.add(referenced_feature_id)
                    else:
                        time_periods_of_referenced_features = time_periods_of_referenced_non_topological_features
                    
                    # The current referenced feature needs to include the time period referenced by the topological polygon (or network).
                    referenced_feature_max_begin_time, referenced_feature_min_end_time = (
                        time_periods_of_referenced_features.get(referenced_feature_id, uninitialised_time_period))
                    time_periods_of_referenced_features[referenced_feature_id] = (
                        max(referenced_feature_max_begin_time, reference_begin_time),
                        min(referenced_feature_min_end_time, reference_end_time))
    
    if restrict_referenced_feature_time_periods:
        # Restrict the valid time periods of all referenced *topological line* features such that they are limited by the time periods of the referencing topologies.
        for referenced_feature_id, (referenced_feature_max_begin_time, referenced_feature_min_end_time) in iteritems(time_periods_of_referenced_topological_line_features):
            referenced_feature = all_features.get(referenced_feature_id)
            if referenced_feature:  # Referenced feature might not actually exist.
                begin_time, end_time = referenced_feature.get_valid_time()
                if begin_time > referenced_feature_max_begin_time:
                    begin_time = referenced_feature_max_begin_time
                if end_time < referenced_feature_min_end_time:
                    end_time = referenced_feature_min_end_time
                
                # If the referenced time period and the topological line feature valid time overlap.
                if begin_time >= end_time:
                    referenced_feature.set_valid_time(begin_time, end_time)
    
    # Iterate over topological lines referenced by topological polygons and networks.
    for referenced_topological_line_feature_id in feature_ids_of_topological_lines_referenced_by_topological_polygons_and_networks:
        topological_line_feature, topological_line_reference = topological_line_references[referenced_topological_line_feature_id]
        topological_line_begin_time, topological_line_end_time = topological_line_feature.get_valid_time()
        for (reference_begin_time, reference_end_time), feature_ids_referenced_in_time_period in topological_line_reference:
            # Restrict the current referenced time period to that of the topological line valid time period.
            if reference_begin_time > topological_line_begin_time:
                reference_begin_time = topological_line_begin_time
            if reference_end_time < topological_line_end_time:
                reference_end_time = topological_line_end_time
            
            # If the referenced time period and the topological line valid time overlap.
            if reference_begin_time >= reference_end_time:
                # Add the features referenced by the current topological line (which is referenced by a topological polygon or network).
                feature_ids_referenced_by_topological_lines_referenced_by_topological_polygons_and_networks.update(feature_ids_referenced_in_time_period)
                
                # Iterate over referenced feature IDs for the currently referenced time period.
                for referenced_feature_id in feature_ids_referenced_in_time_period:
                    # The current referenced feature needs to include the time period referenced by the topological line.
                    referenced_feature_max_begin_time, referenced_feature_min_end_time = (
                        time_periods_of_referenced_non_topological_features.get(referenced_feature_id, uninitialised_time_period))
                    time_periods_of_referenced_non_topological_features[referenced_feature_id] = (
                        max(referenced_feature_max_begin_time, reference_begin_time),
                        min(referenced_feature_min_end_time, reference_end_time))
    
    if restrict_referenced_feature_time_periods:
        # Restrict the valid time periods of all referenced *non-topological* features such that they are limited by the time periods of the referencing topologies.
        for referenced_feature_id, (referenced_feature_max_begin_time, referenced_feature_min_end_time) in iteritems(time_periods_of_referenced_non_topological_features):
            referenced_feature = all_features.get(referenced_feature_id)
            if referenced_feature:  # Referenced feature might not actually exist.
                begin_time, end_time = referenced_feature.get_valid_time()
                if begin_time > referenced_feature_max_begin_time:
                    begin_time = referenced_feature_max_begin_time
                if end_time < referenced_feature_min_end_time:
                    end_time = referenced_feature_min_end_time
                
                # If the referenced time period and the non-topological feature valid time overlap.
                if begin_time >= end_time:
                    referenced_feature.set_valid_time(begin_time, end_time)
    
    # Only keep features that are topological polygons and networks, and any features referenced directly
    # or indirectly by them. For example, a topological polygon might reference a topological line which
    # in turn references regular features. In this case the topological line and the features it references
    # must all be kept.
    feature_ids_to_keep = (feature_ids_of_topological_polygons_and_networks |
                           feature_ids_referenced_by_topological_polygons_and_networks |
                           feature_ids_referenced_by_topological_lines_referenced_by_topological_polygons_and_networks)
    for feature_collection in feature_collections:
        # Create an extra feature collection (containing removed features) for each feature collection.
        if removed_features_collections is not None:
            removed_features_collection = pygplates.FeatureCollection()
            removed_features_collections.append(removed_features_collection)
        
        feature_index = 0
        while feature_index < len(feature_collection):
            feature = feature_collection[feature_index]
            feature_id = feature.get_feature_id()
            if feature_id not in feature_ids_to_keep:
                del feature_collection[feature_index]
                feature_index -= 1
                
                # Keep track of the removed feature if requested.
                if removed_features_collections is not None:
                    removed_features_collection.add(feature)
            
            feature_index += 1
    
    # Return our (potentially) modified feature collections as a list of pygplates.FeatureCollection.
    return [pygplates.FeatureCollection(feature_collection)
        for feature_collection in feature_collections]


# Private helper class (has '_' prefix) to find topology-related GpmlPropertyDelegate's.
class _TopologicalReferenceVisitor(pygplates.PropertyValueVisitor):
    ALL_TIME = float('inf'), float('-inf')  # begin_time, end_time
    
    def __init__(self):
        super(_TopologicalReferenceVisitor, self).__init__()
    
    def visit_feature(self, feature):
        self.topology_type = None
        self.references = []
        self.current_time_period = self.ALL_TIME
        
        # Visit all properties in the feature to find a topological line, polygon or network.
        for property in feature:
            # Get the top-level property value (containing all times) not just a specific time.
            property_value = property.get_time_dependent_value()
            # Visit the property value.
            property_value.accept_visitor(self)
            # If we visited a topological line, polygon or network then we're finished with the current feature.
            if self.topology_type:
                break
        
        return self.topology_type, self.references
    
    def visit_gpml_constant_value(self, gpml_constant_value):
        # Visit the GpmlConstantValue's nested property value.
        gpml_constant_value.get_value().accept_visitor(self)
    
    def visit_gpml_piecewise_aggregation(self, gpml_piecewise_aggregation):
        # Only need to visit if contains a topological line, polygon or network.
        value_type = gpml_piecewise_aggregation.get_value_type()
        if (value_type == pygplates.GpmlTopologicalLine or
            value_type == pygplates.GpmlTopologicalPolygon or
            value_type == pygplates.GpmlTopologicalNetwork):
            
            # NOTE: If there's only *one* time window then we ignore its time period.
            #
            # We do this for the same reason that GPlates does this (this comment from the GPlates source code)...
            #
            # This is because GPML files created with old versions of GPlates set the time period,
            # of the sole time window, to match that of the 'feature's time period (in the topology
            # build/edit tools) - newer versions set it to *all* time (distant past/future) - in fact
            # newer versions just use a GpmlConstantValue instead of GpmlPiecewiseAggregation because
            # the topology tools cannot yet create time-dependent topology (section) lists.
            # With old versions if the user expanded the 'feature's time period *after* building/editing
            # the topology then the *un-adjusted* time window time period will be incorrect and hence
            # we need to ignore it here.
            # Those old versions were around 4 years ago (prior to GPlates 1.3) - so we really shouldn't
            # be seeing any old topologies.
            # Actually I can see there are some currently in the sample data for GPlates 2.0.
            # So as a compromise we'll ignore the reconstruction time if there's only one time window
            # (a single time window shouldn't really have any time constraints on it anyway)
            # and respect the reconstruction time if there's more than one time window
            # (since multiple time windows need non-overlapping time constraints).
            # This is especially true now that pyGPlates will soon be able to generate time-dependent
            # topologies (where the reconstruction time will need to be respected otherwise multiple
            # networks from different time periods will get created instead of just one of them).
            if len(gpml_piecewise_aggregation) == 1:
                # Assume the sole time window covers *all* time (the default).
                gpml_piecewise_aggregation[0].get_value().accept_visitor(self)
            else:
                # Visit the property value in each time window.
                for gpml_time_window in gpml_piecewise_aggregation:
                    # Restrict the time period while we're visiting the time window.
                    self.current_time_period = gpml_time_window.get_begin_time(), gpml_time_window.get_end_time()
                    gpml_time_window.get_value().accept_visitor(self)
                    self.current_time_period = self.ALL_TIME
    
    def visit_gpml_topological_line(self, gpml_topological_line):
        referenced_feature_ids = set()
        # Topological line sections are topological sections (which contain a property delegate).
        for section in gpml_topological_line.get_sections():
            referenced_feature_ids.add(section.get_property_delegate().get_feature_id())
        
        self.topology_type = pygplates.GpmlTopologicalLine
        self.references.append((self.current_time_period, referenced_feature_ids))
    
    def visit_gpml_topological_polygon(self, gpml_topological_polygon):
        referenced_feature_ids = set()
        # Topological polygon exterior sections are topological sections (which contain a property delegate).
        for exterior_section in gpml_topological_polygon.get_exterior_sections():
            referenced_feature_ids.add(exterior_section.get_property_delegate().get_feature_id())
        
        self.topology_type = pygplates.GpmlTopologicalPolygon
        self.references.append((self.current_time_period, referenced_feature_ids))
    
    def visit_gpml_topological_network(self, gpml_topological_network):
        referenced_feature_ids = set()
        # Topological network boundary sections are topological sections (which contain a property delegate).
        for boundary_section in gpml_topological_network.get_boundary_sections():
            referenced_feature_ids.add(boundary_section.get_property_delegate().get_feature_id())
        # Topological network interiors are already property delegates.
        for interior in gpml_topological_network.get_interiors():
            referenced_feature_ids.add(interior.get_feature_id())
        
        self.topology_type = pygplates.GpmlTopologicalNetwork
        self.references.append((self.current_time_period, referenced_feature_ids))


if __name__ == '__main__':
    
    import os.path
    
    
    # Check the imported pygplates version.
    if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
        print('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
                os.path.basename(__file__), pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED),
            file=sys.stderr)
        sys.exit(1)
    
    
    import argparse
    
    
    def main():
    
        __description__ = \
    """Remove any regular features not referenced by topological features.
    
    The input files should contain the topological and regular features that make up the topological model.
    Ensure you at least specify all topological features in the topological model, otherwise regular features
    referenced by the missing topologies will get removed.
    
    The results are written back to the input files unless an output filename prefix is provided.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -o cleanup_topologies_ -- topologies.gpml
     """

        # The command-line parser.
        parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
        
        parser.add_argument('-o', '--output_filename_prefix', type=str,
                metavar='output_filename_prefix',
                help='Optional output filename prefix. If one is provided then an output file '
                     'is created for each input file by prefixing the input filenames. '
                     'If no filename prefix is provided then the input files are overwritten.')
        
        parser.add_argument('-d', '--removed_features_filename_prefix', type=str,
                metavar='removed_features_filename_prefix',
                help='Option to save removed features in new files with specified filename prefix. '
                     'If specified then a file is created for each input file (that has features removed) '
                     'by prefixing the input filenames. If no filename prefix is provided then the '
                     'removed features are not saved.')
        
        parser.add_argument('-p', '--restricted_referenced_time_periods', action="store_true",
                help='If specified then restrict the time periods of features referenced by topologies such that they are '
                     'limited by the time periods of the referencing topologies (default is no restriction).')
        
        parser.add_argument('input_filenames', type=str, nargs='+',
                metavar='input_filename',
                help='One or more files containing topological features and features referenced by them.')
        
        # Parse command-line options.
        args = parser.parse_args()
        
        # Read the input feature collections.
        input_feature_collections = [pygplates.FeatureCollection(input_filename)
                for input_filename in args.input_filenames]
        
        # If we're saving the removed features then provide a list for those feature collections.
        removed_features_collections = [] if args.removed_features_filename_prefix else None
        
        # Remove features not referenced by topologies.
        output_feature_collections = remove_features_not_referenced_by_topologies(
            input_feature_collections,
            args.restricted_referenced_time_periods,
            removed_features_collections)
        
        # Write the modified feature collections to disk.
        for feature_collection_index in range(len(output_feature_collections)):
            # Each output filename is the input filename with an optional prefix prepended.
            input_filename = args.input_filenames[feature_collection_index]
            
            if args.output_filename_prefix:
                dir, file_basename = os.path.split(input_filename)
                output_filename = os.path.join(dir, '{0}{1}'.format(args.output_filename_prefix, file_basename))
            else:
                output_filename = input_filename
            output_feature_collections[feature_collection_index].write(output_filename)
            
            # Write the removed features to disk.
            if args.removed_features_filename_prefix:
                dir, file_basename = os.path.split(input_filename)
                removed_features_filename = os.path.join(dir, '{0}{1}'.format(args.removed_features_filename_prefix, file_basename))
                removed_features_collections[feature_collection_index].write(removed_features_filename)
        
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
