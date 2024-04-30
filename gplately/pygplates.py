"""
A light wrapping of some 
[`pyGPlates`](https://www.gplates.org/docs/pygplates/index.html) 
classes to keep track of filenames.

Each object listed here will have a `self.filenames` attribute.
"""

import warnings as _warnings
from copy import copy
from typing import List, Union

import pygplates as _pygplates
from pygplates import *

_warnings.simplefilter("always", ImportWarning)


def _is_string(value):
    # convert sets to list
    if type(value) is set:
        value = list(value)

    # check for strings inside a list
    if type(value) is list:
        bl = []
        for val in value:
            bl.append(type(val) is str)
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

    def __init__(self, rotation_features, default_anchor_plate_id=0):
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
        super(RotationModel, self).__init__(
            rotation_features, default_anchor_plate_id=default_anchor_plate_id
        )
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
            msg = "\nRotationModel: No filename associated with {} in __init__".format(
                type(rotation_features)
            )
            msg += "\n ensure pygplates is imported from gplately. Run,"
            msg += "\n from gplately import pygplates"
            _warnings.warn(msg, ImportWarning)
            self.filenames = []


class Feature(_pygplates.Feature):
    """A class that wraps the `pyGPlates.Feature` class. This contains tools to query and set
    geological or plate-tectonic feature properties defined by the
    [GPlates Geological Information Model (GPGIM)](https://www.gplates.org/docs/gpgim/).
    A feature consists of a collection of `properties`, a `feature type` and a `feature id`.

    See the link below for inherited methods

    https://www.gplates.org/docs/pygplates/generated/pygplates.feature

    """

    #
    # this class seems unfinished. need to implement properly in the future.
    #

    def __init__(
        self,
        feature_type: _pygplates.FeatureType = _pygplates.FeatureType.gpml_unclassified_feature,
        feature_id: str = None,
        verify_information_model=_pygplates.VerifyInformationModel.yes,
        *,
        filenames: Union[str, List[str]] = [],
        feature: "Feature" = None,
    ):
        """
        Notes
        -----
        The signature of this constructor has been changed since gplately 1.3.0 to be compatible with pygplates.
        THe 'filenames' and 'feature' parameters must be given as keyword argument.

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

        filenames: `str` or `list` of `str`
            The filenames being associated with this feature.

        feature: instance of `gplately.pygplates.Feature`
            The other "Feature" object

        Raises
        ------
        ImportWarning
            If neither a `str` nor `list` of `str` is passed, no
            `Feature` filenames will be collected, and the user will be alerted of this.

        InformationModelError
            if `verify_information_model` is `VerifyInformationModel.yes` and `feature_type` is not a recognised feature type.

        """

        # bugfix: gplately.pygplates.Feature is not compatible with pygplates.Feature
        # see https://github.com/GPlates/gplately/issues/150
        # this gplately.pygplates.Feature class seems not completed yet. for example, the clone() method returns nothing. It looks unfinished.
        # Why is a feature associated with multiple file names?
        super().__init__(feature_type, feature_id, verify_information_model)
        self.filenames = []

        # try the best to detect backward compatibility issue
        if not isinstance(feature_type, _pygplates.FeatureType):
            raise Exception(
                "The __init__() signature has been changed. The first positional argument(besides self) is 'feature_type' now. "
                + "Check the online documentation https://gplates.github.io/gplately/pygplates.html"
            )

        # update filename list
        if isinstance(filenames, list) and all(
            isinstance(filename, str) for filename in filenames
        ):
            self.filenames = filenames
        elif isinstance(filenames, str):
            self.filenames = [filenames]
        else:
            msg = (
                f"\nFeature: No filename associated with {type(filenames)} in __init__"
                + "\n ensure pygplates is imported from gplately. Run,"
                + "\n from gplately import pygplates"
            )
            _warnings.warn(msg, ImportWarning)
            self.filenames = []

        if feature:
            self.filenames = feature.filenames
            # TODO: also need to copy everything else in the other feature into this feature

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
            msg = "\nFeature: No filename associated with {} in add".format(
                type(feature)
            )
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
    """A thin wrap of [`pyGPlates.FeatureCollection`](https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection).
    The derived FeatureCollection class contains a "filenames" attribute to track the names of files from which the feature collection is loaded.

    See the doc of base class at https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection.

    """

    def __init__(self, features=None, *, filenames: List[str] = []):
        """The constructor is compatible with the base class pygplates.FeatureCollection.

        Parameters
        ----------
        features : instance of `Feature` or `str` or a sequence (eg, `list` or `tuple`) of `Feature`
            An optional filename, or sequence of features, or a single feature

        filenames: a list of file path strings (must be given as keyword argument)
            file paths

        Raises
        ------
        OpenFileForReadingError
            If file is not readable (if filename specified).

        FileFormatNotSupportedError
            If file format (identified by the filename extension) does not support reading (when filename specified).

        """
        super().__init__(features)
        self._filenames = filenames

    @classmethod
    def from_file_list(cls, filenames: List[str] = []):
        """class method to load a feature collection from multiple files."""

        if not (
            isinstance(filenames, list)
            and all(isinstance(filename, str) for filename in filenames)
        ):
            raise Exception(
                f"The 'filenames' parameter must be a list of file path strings."
            )

        fc = _pygplates.FeatureCollection()
        for filename in filenames:
            fc.add(_pygplates.FeatureCollection(filename))

        ret = cls([f for f in fc], filenames=filenames)
        return ret

    @property
    def filenames(self):
        return self._filenames

    @filenames.setter
    def filenames(self, filenames):
        self._filenames = filenames

    @filenames.deleter
    def filenames(self):
        del self._filenames

    def add(self, feature=[], *, filename: str = None):
        """load one more file into the feature collection

        Parameters
        ----------

        feature : (Feature or sequence (eg, list or tuple) of Feature)
            one or more features to add

        filename : str
            file path string(must be passed as keyword argument)

        """
        if isinstance(feature, str):
            raise Exception(
                "The filename argument must be passed as keyword argument. The 'feature' argument must be Feature or sequence of Features."
            )

        if feature:
            super().add(feature)
        if filename:
            super().add(_pygplates.FeatureCollection(filename))
            self.filenames.append(filename)

    def clone(self):
        """Create a duplicate of this feature collection instance.

        Returns
        -------
        feature_collection : instance of `gplately.pygplates.FeatureCollection`
            The cloned `FeatureCollection`.

        """
        return FeatureCollection([f for f in self], filenames=copy(self.filenames))
