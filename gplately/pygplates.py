"""
A light wrapping of some 
[`pyGPlates`](https://www.gplates.org/docs/pygplates/index.html) 
classes to keep track of filenames.

Each object listed here will have a `self.filenames` attribute.
"""

import pygplates as _pygplates
from pygplates import *
import warnings as _warnings
_warnings.simplefilter('always', ImportWarning)

def _is_string(value):
    # convert sets to list
    if type(value) is set:
        value = list(value)

    # check for strings inside a list
    if type(value) is list:
        bl = []
        for val in value:
            bl.append( type(val) is str )
        return all(bl)
    
    # if no list, check if string
    else:
        return type(value) is str


class RotationModel(_pygplates.RotationModel):
    """A class that wraps the 
    [`pyGPlates.RotationModel` class](https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel#pygplates.RotationModel). 
    This queries a finite rotation of a moving plate relative to any other plate, 
    optionally between two instants in geological time. 

    See [Plate reconstruction hierarchy](https://www.gplates.org/docs/pygplates/pygplates_foundations.html#pygplates-foundations-plate-reconstruction-hierarchy).

    This class provides an easy way to query rotations in any of the four 
    combinations of total/stage and equivalent/relative rotations using 
    [`get_rotation()`](https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel#pygplates.RotationModel.get_rotation). 

    Reconstruction trees can also be created at any instant 
    of geological time and these are cached internally depending on a 
    user-specified cache size parameter pass to `gplately.pygplates.RotationModel.__init__()`. 
    The reconstruction_tree_cache_size parameter of those methods controls the 
    size of an internal least-recently-used cache of reconstruction trees 
    (evicts least recently requested reconstruction tree when a new 
    reconstruction time is requested that does not currently exist in the cache). 
    This enables reconstruction trees associated with different reconstruction 
    times to be re-used instead of re-creating them, provided they have not been 
    evicted from the cache. This benefit also applies when querying rotations with 
    [`get_rotation()`](https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel#pygplates.RotationModel.get_rotation) 
    since it, in turn, requests reconstruction trees.


    This wrapping of `pygplates.RotationModel` contains all 
    [`pygplates.RotationModel` functionality](https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel#pygplates.RotationModel),
    and in addition tracks the names of files from which the rotation feature(s) are read
    using the `gplately.pygplates.RotationModel.filenames` attribute.  


    """
    def __init__(self, rotation_features):
        """**A RotationModel object can be constructed in three ways.**

        ---
        **1. Create from rotation feature collection(s) and/or rotation filename(s)**
        --------------------------------------------------------------------------
        
        Parameters
        ----------
        rotation_features : instance of `pygplates.FeatureCollection` or `str` or instance of `pygplates.Feature` or sequence of `pygplates.Feature` or sequence of any combination of those four types
            A rotation feature collection, or rotation filename, or rotation feature, 
            or sequence of rotation features, or a sequence (eg, `list` or `tuple`) of 
            any combination of those four types.

        reconstruction_tree_cache_size : int, default 150
            Number of reconstruction trees to cache internally. Defaults to 150.

        extend_total_reconstruction_poles_to_distant_past : bool, default False
            Extend each moving plate sequence back infinitely far into the distant 
            past such that reconstructed geometries will not snap back to their 
            present day positions when the reconstruction time is older than 
            the oldest times specified in the rotation features (defaults to False).

        default_anchor_plate_id : int, default 0
            The default anchored plate id to use when 
            [`get_rotation()`](https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel#pygplates.RotationModel.get_rotation)  
            and [`get_reconstruction_tree()`](https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel#pygplates.RotationModel.get_reconstruction_tree)
            are called without specifying their `anchor_plate_id` parameter. Defaults to 0.

        Raises
        ------
        OpenFileForReadingError 
            If any file is not readable (when filenames specified)

        FileFormatNotSupportedError 
            If any file format (identified by the filename extensions) 
            does not support reading (when filenames specified)


        Note that `rotation_features` can be a rotation `FeatureCollection` or a rotation filename or a rotation Feature or a sequence of rotation features, or a sequence (eg, `list` or `tuple`) of any combination of those four types.

        If any rotation filenames are specified then this method uses `FeatureCollection` internally to read the rotation files.


        Example
        -------
        Load a rotation file and some rotation adjustments (as a collection of rotation features) into a rotation model:

            rotation_adjustments = pygplates.FeatureCollection()
            ...
            rotation_model = pygplates.RotationModel(['rotations.rot', rotation_adjustments])


        ---
        **2. Create from an existing rotation model but adapt it with a potentially different cache size and/or default anchor plate ID**
        ---------------------------------------------------------------------------------------------------------------------------------

        Parameters
        ----------
        rotation_model : instance of `pygplates.RotationModel`
            An existing rotation model.

        reconstruction_tree_cache_size : int, default 2
            Number of reconstruction trees to cache internally. 
            Defaults to 2 - this is much lower than the usual default 
            cache size since the existing rotation model likely 
            already has a sizeable cache anyway - and if you are 
            leaving this at its default value then you are presumably 
            only interested in changing the default anchor plate ID 
            (not increasing the cache size).

        default_anchor_plate_id : int, defaults to the default anchor plate of `rotation_model`
            The default anchored plate id to use when 
            [`get_rotation()`](https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel#pygplates.RotationModel.get_rotation)
            and [`get_reconstruction_tree()`](https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel#pygplates.RotationModel.get_reconstruction_tree)
            are called without specifying their `anchor_plate_id` parameter. 
            Defaults to the default anchor plate of `rotation_model`.


        This is useful if you want to use an existing rotation model but with a 
        larger cache size or a different default anchor plate ID:

        Example
        -------
        The below example changes the default anchor plate ID:

            rotation_model = pygplates.RotationModel(rotation_files)
            ...
            rotation_model_anchor_1 = pygplates.RotationModel(rotation_model, default_anchor_plate_id=1)

        ---
        **3. Return an existing rotation model as a convenience**
        -------------------------------------------------------

        This is useful when defining your own function that accepts 
        rotation features or a rotation model. It avoids the hassle 
        of having to explicitly test for each source type:

            def my_function(rotation_features_or_model):
            # The appropriate constructor (__init__) overload is chosen depending on argument type.
            rotation_model = pygplates.RotationModel(rotation_features_or_model)
            ...

        Parameters
        ----------
        rotation_model : instance of `pygplates.RotationModel` or `gplately.pygplates.RotationModel`
            An existing rotation model.

        ---
        """
        super(RotationModel, self).__init__(rotation_features)
        self.filenames = []

        # update filename list
        if _is_string(rotation_features) and type(rotation_features) is list:
            self.filenames = rotation_features
        elif _is_string(rotation_features) and type(rotation_features) is str:
            self.filenames = [rotation_features]
        elif rotation_features is None:
            self.filenames = []
        elif isinstance(rotation_features, RotationModel):
            self.filenames = rotation_features.filenames
        elif hasattr(rotation_features, "filenames"):
            self.filenames = rotation_features.filenames
        else:
            msg = "\nRotationModel: No filename associated with {} in __init__".format(type(rotation_features))
            msg += "\n ensure pygplates is imported from gplately. Run,"
            msg += "\n from gplately import pygplates"
            _warnings.warn(msg, ImportWarning)
            self.filenames = []

class Feature(_pygplates.Feature):
    """A class that wraps the `pyGPlates.Feature` class. This contains tools to query and set 
    geological or plate-tectonic feature properties defined by the 
    [GPlates Geological Information Model (GPGIM)](https://www.gplates.org/docs/gpgim/).
    A feature consists of a collection of `properties`, a `feature type` and a `feature id`.

    The following operations for iterating over the properties in a feature are supported:

    Iterating over `feature` `properties`
    -------------------------------------
    * `len(f)`: Number of properties in feature `f`

    * `for p in f`: Iterates over the properties `p` in feature `f`


    This wrapping of `pygplates.Feature` contains all `pygplates.Feature` functionality,
    and in addition tracks the names of files from which the feature(s) are read 
    using the `gplately.pygplates.Feature.filenames` attribute.  


    Creating `feature`s
    -------------------

    The following methods provide convenient ways to create features:

    * [`create_reconstructable_feature()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.create_reconstructable_feature)

    * [`create_topological_feature()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.create_topological_feature)

    * [`create_tectonic_section()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.create_tectonic_section)

    * [`create_flowline()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.create_flowline)

    * [`create_motion_path()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.create_motion_path)

    * [`create_total_reconstruction_sequence()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.create_total_reconstruction_sequence)

    The following methods return the feature type and feature id:

    * [`get_feature_type()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_feature_type)

    * [`get_feature_id()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_feature_id)


    Working with `feature` properties
    ---------------------------------

    The following methods provide generic support for adding, removing, setting and getting properties:

    * [`add()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.add)

    * [`remove()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.remove)

    * [`set()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set)

    * [`get()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get)

    * [`get_value()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_value)


    Setting and getting `feature` geometries
    ----------------------------------------

    The following methods provide a convenient way to set and get feature geometry:

    * [`set_geometry()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_geometry)

    * [`get_geometry()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_geometry)

    * [`get_geometries()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_geometries)

    * [`get_all_geometries()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_all_geometries)


    Setting and getting feature topological geometry
    ------------------------------------------------

    The following methods provide a convenient way to set and get feature topological geometry (which can be a topological line, polygon or network):

    * [`set_topological_geometry()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_topological_geometry)

    * [`get_topological_geometry()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_topological_geometry)

    * [`get_topological_geometries()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_topological_geometries)

    * [`get_all_topological_geometries()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_all_topological_geometries)


    Setting and getting attributes imported from a shapefile
    --------------------------------------------------------

    The following methods provide a convenient way to set and get attributes imported from a Shapefile:

    * [`set_shapefile_attribute()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_shapefile_attribute)

    * [`set_shapefile_attributes()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_shapefile_attributes)

    * [`get_shapefile_attribute()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_shapefile_attribute)

    * [`get_shapefile_attributes()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_shapefile_attributes)


    Setting and getting enumeration properties
    ------------------------------------------

    The following methods provide a convenient way to set and get enumeration properties:

    * [`set_enumeration()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_enumeration)

    * [`get_enumeration()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_enumerationß)



    Setting and getting string, floating-point, integer and boolean properties
    --------------------------------------------------------------------------

    The following methods provide a convenient way to set and get string, floating-point, integer and boolean properties:

    * [`set_string()`](set_string())

    * [`get_string()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_string)

    * [`set_double()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_double)

    * [`get_double()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_double)

    * [`set_integer()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_integer)

    * [`get_integer()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_integer)

    * [`set_boolean()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_boolean)

    * [`get_boolean()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_boolean)


    Setting and getting common `feature` `properties`
    -------------------------------------------------

    The following methods provide a convenient way to set and get some of the properties that are common to many feature types:

    * [`set_name()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_name)

    * [`get_name()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_name)

    * [`set_description()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_description)

    * [`get_description()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_description)

    * [`set_valid_time()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_valid_time)

    * [`get_valid_time()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_valid_time)

    * [`is_valid_at_time()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.is_valid_at_time)

    * [`set_reconstruction_plate_id()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_reconstruction_plate_id)

    * [`get_reconstruction_plate_id()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_reconstruction_plate_id)

    * [`set_conjugate_plate_id()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_reconstruction_plate_id)

    * [`get_conjugate_plate_id()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_reconstruction_plate_id)

    * [`set_left_plate()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_left_plate)

    * [`get_left_plate()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_left_plate)

    * [`set_right_plate()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_right_plate)

    * [`get_right_plate()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_right_plate)

    * [`set_relative_plate()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_relative_plate)

    * [`get_relative_plate()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_relative_plate)

    * [`set_times()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_times)

    * [`get_times()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_times)

    * [`set_reconstruction_method()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_reconstruction_method)

    * [`get_reconstruction_method()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_reconstruction_method)

    * [`set_geometry_import_time()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_geometry_import_time)

    * [`get_geometry_import_time()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_geometry_import_time)

    * [`set_total_reconstruction_pole()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.set_total_reconstruction_pole)

    * [`get_total_reconstruction_pole()`](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.get_total_reconstruction_pole)

    For other properties the generic `set()`, `get()` and `get_value()` methods will still need to be used.

    A feature can be deep copied using `clone()`.

    """
    def __init__(self, feature):
        """

        Parameters
        ----------
        feature_type : instance of `pygplates.FeatureType`
            The type of feature. See 
            [here](https://www.gplates.org/docs/pygplates/generated/pygplates.featuretype#pygplates.FeatureType) 
            for a list of pygplates feature types. 

        feature_id : instance of `pygplates.FeatureId`
            The [feature identifier](https://www.gplates.org/docs/pygplates/generated/pygplates.featureid#pygplates.FeatureID).

        verify_information_model : instance of `VerifyInformationModel.yes` or `VerifyInformationModel.no`
            Specify whether to check `feature_type` with the information model (default) or not.

        Raises
        ------
        ImportWarning
            If neither a `str`, `list` of `str`, `gplately.pygplates.Feature` or `None` is passed, no
            `Feature` filenames will be collected, and the user will be alerted of this.

        InformationModelError 
            if `verify_information_model` is `VerifyInformationModel.yes` and `feature_type` is not a recognised feature type.

        """
        super(Feature, self).__init__(feature)
        self.filenames = []

        # update filename list
        if _is_string(feature) and type(feature) is list:
            self.filenames = feature
        elif _is_string(feature) and type(feature) is str:
            self.filenames = [feature]
        elif feature is None:
            self.filenames = []
        elif isinstance(feature, Feature):
            self.filenames = feature.filenames
        elif hasattr(feature, "filenames"):
            self.filenames = feature.filenames
        else:
            msg = "\nFeature: No filename associated with {} in __init__".format(type(feature))
            msg += "\n ensure pygplates is imported from gplately. Run,"
            msg += "\n from gplately import pygplates"
            _warnings.warn(msg, ImportWarning)
            self.filenames = []


    def add(self, feature):
        """Adds a property (or properties) to this feature. See original docs 
        [here](https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature.add).

        Parameters
        ----------
        property_name : instance of `pygplates.PropertyName`
            The name of the property (or properties) to add.

        property_value : instance of `pygplates.PropertyValue` or sequence (eg, `list` or `tuple`) of `pygplates.PropertyValue`
            The value (or values) of the property (or properties) to add.

        verify_information_model : instance of `VerifyInformationModel.yes` or `VerifyInformationModel.no`
            Specify whether to check `feature_type` with the information model (default) or not.


        Returns
        -------
        property_added : `Property` or `list` of `Property` depending on whether `property_value` is a `PropertyValue` or sequence of `PropertyValue`
            The property (or properties) added to the feature.
        """
        super().add(feature)
        if isinstance(feature, Feature):
            self.filenames.extend(feature.filenames)
        elif _is_string(feature):
            self.filenames.extend(feature)
        elif hasattr(feature, "filenames"):
            self.filenames.extend(feature.filenames)
        else:
            msg = "\nFeature: No filename associated with {} in add".format(type(feature))
            msg += "\n ensure pygplates is imported from gplately. Run,"
            msg += "\n from gplately import pygplates"
            _warnings.warn(msg, ImportWarning)

    def clone(self):
        """Create a duplicate of this `Feature` instance.

        This creates a new `Feature` instance with cloned versions of this feature’s `properties`. 
        The cloned feature is created with its own unique `pygplates.FeatureId`.

        Returns
        -------
        Feature : instance of `gplately.pygplates.Feature`
            The cloned `Feature` instance.

        """
        feat = super().clone()
        feat.filenames = self.filenames

class FeatureCollection(_pygplates.FeatureCollection):
    """A class that wraps the 
    [`pyGPlates.FeatureCollection`](https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection) 
    class. This aggregates a set of features into a collection. 
    This is traditionally so that a group of  features can be loaded, saved or 
    unloaded in a single operation.

    This wrapping of `pygplates.FeatureCollection` contains all 
    [`pygplates.FeatureCollection` functionality](https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection),
    and in addition tracks the names of files from which the feature 
    collection(s) are read using the 
    `gplately.pygplates.FeatureCollection.filenames` attribute.  

    Examples
    --------

    For example, to read coastline features from a file:

        coastline_feature_collection = pygplates.FeatureCollection('coastlines.gpml')

    And to write coastline features to a file:

        coastline_feature_collection = pygplates.FeatureCollection(coastline_features)
        coastline_feature_collection.write('coastlines.gpml')

    To create a new feature collection from a sequence of `features`:

        feature_collection = pygplates.FeatureCollection([feature1, feature2])

        # ...is the equivalent of...

        feature_collection = pygplates.FeatureCollection()
        feature_collection.add(feature1)
        feature_collection.add(feature2)

    The following feature collection file formats are currently supported:

    ---

    |         **File Format**        | **Filename  Extension** | **Supports  Read** | **Supports Write** |
    |:------------------------------:|:-----------------------:|:------------------:|:------------------:|
    |     GPlates Markup Language    |         ‘.gpml’         |         Yes        |         Yes        |
    |         Compressed GPML        |  ‘.gpmlz’ or ‘.gpml.gz’ |         Yes        |         Yes        |
    |          PLATES4 line          |     ‘.dat’ or ‘.pla’    |         Yes        |         Yes        |
    |        PLATES4 rotation        |          ‘.rot’         |         Yes        |         Yes        |
    |        GPlates rotation        |         ‘.grot’         |         Yes        |         Yes        |
    |         ESRI Shapefile         |          ‘.shp’         |         Yes        |         Yes        |
    |             GeoJSON            |  ‘.geojson’ or ‘.json’  |         Yes        |         Yes        |
    |           GeoPackage           |         ‘.gpkg’         |         Yes        |         Yes        |
    |             OGR GMT            |          ‘.gmt’         |         Yes        |         Yes        |
    |             GMT xy             |          ‘.xy’          |         No         |         Yes        |
    | GMAP Virtual Geomagnetic Poles |          ‘.vgp’         |         Yes        |         No         |

    ---

    In the future, support will be added to enable users to implement and register readers/writers for other file formats (or their own non-standard file formats).

    
    Operations for accessing features
    ---------------------------------

    The following operations for accessing the features are supported:


    * `len(fc)` : Number of features in feature collection `fc`.

    * `for f in fc` : Iterates over the features `f` in feature collection `fc`.

    * `fc[i]` : The feature of fc at index `i`. 


    For example:

        num_features = len(feature_collection)
        features_in_collection = [feature for feature in feature_collection]
        # assert(num_features == len(features_in_collection))

    """
    def __init__(self, features=None):
        """

        Parameters
        ----------
        features : instance of `Feature` or `str` or a sequence (eg, `list` or `tuple`) of `Feature`
            An optional filename, or sequence of features, or a single feature


        Raises
        ------
        OpenFileForReadingError 
            If file is not readable (if filename specified).

        FileFormatNotSupportedError 
            If file format (identified by the filename extension) does not support reading (when filename specified).

        """
        super(FeatureCollection, self).__init__(features)
        self.filenames = []

        # update filename list
        if _is_string(features) and type(features) is list:
            self.filenames = features
        elif _is_string(features) and type(features) is str:
            self.filenames = [features]
        elif features is None:
            self.filenames = []
        elif isinstance(features, FeatureCollection):
            self.filenames = features.filenames
        elif hasattr(features, "filenames"):
            self.filenames = features.filenames
        else:
            msg = "\nFeatureCollection: No filename associated with {} in __init__".format(type(features))
            msg += "\n ensure pygplates is imported from gplately. Run,"
            msg += "\n from gplately import pygplates"
            _warnings.warn(msg, ImportWarning)
            self.filenames = []

    def add(self, features):
        """Adds one or more features to this collection.

        Parameters
        ----------
        feature : instance of `Feature` or sequence (eg, `list` or `tuple`) of `Feature`
            One or more features to add.

        A feature collection is an unordered collection of features so there is no concept 
        of where a feature is inserted in the sequence of features.

            feature_collection.add(feature)
            feature_collection.add([feature1, feature2])
        """
        super().add(features)

        # update filename list
        if isinstance(features, FeatureCollection):
            self.filenames.extend(features.filenames)
        elif _is_string(features):
            self.filenames.extend(features)
        elif hasattr(features, "filenames"):
            self.filenames.extend(features.filenames)
        else:
            msg = "\nFeatureCollection: No filename associated with {} in add".format(type(features))
            msg += "\n ensure pygplates is imported from gplately. Run,"
            msg += "\n from gplately import pygplates"
            _warnings.warn(msg, ImportWarning)

    def clone(self):
        """Create a duplicate of this feature collection instance.

        This creates a new `FeatureCollection` instance with cloned versions of 
        this collection’s features. And the cloned features (in the cloned 
        collection) are each created with a unique `FeatureId`.

        Returns
        -------
        feature_collection : instance of `gplately.pygplates.FeatureCollection`
            The cloned `FeatureCollection`. 

        """
        fc = super().clone()
        fc.filenames = self.filenames
        return fc