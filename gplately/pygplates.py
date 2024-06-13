"""
A light wrapping of some [`pyGPlates`](https://www.gplates.org/docs/pygplates/index.html) 
classes to keep track of filenames. Each object listed here will have a `self.filenames` attribute.
"""

import warnings as _warnings
from copy import copy
from typing import List, Union

import pygplates as _pygplates
from pygplates import *

_warnings.simplefilter("always", ImportWarning)


class RotationModel(_pygplates.RotationModel):
    """A class that wraps the [`pyGPlates.RotationModel` class](https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel#pygplates.RotationModel)."""

    def __init__(
        self,
        rotation_features,
        reconstruction_tree_cache_size=150,
        extend_total_reconstruction_poles_to_distant_past=False,
        default_anchor_plate_id=0,
    ):
        # N.B. if rotation_features is a RotationModel, the constructor
        # will not accept 'extend_total_reconstruction_poles_to_distant_past'
        # as an argument
        if isinstance(rotation_features, _pygplates.RotationModel):
            super(RotationModel, self).__init__(
                rotation_features,
                reconstruction_tree_cache_size=reconstruction_tree_cache_size,
                default_anchor_plate_id=default_anchor_plate_id,
            )
        else:
            super(RotationModel, self).__init__(
                rotation_features,
                reconstruction_tree_cache_size=reconstruction_tree_cache_size,
                extend_total_reconstruction_poles_to_distant_past=extend_total_reconstruction_poles_to_distant_past,
                default_anchor_plate_id=default_anchor_plate_id,
            )

        if isinstance(rotation_features, str):
            self._filenames = [rotation_features]
        elif (
            hasattr(rotation_features, "__iter__")
            and all(isinstance(f, str) for f in rotation_features)
        ):
            self._filenames = list(rotation_features)
        else:
            self._filenames = []

    @property
    def filenames(self):
        return self._filenames

    @filenames.setter
    def filenames(self, filenames):
        self._filenames = filenames

    @filenames.deleter
    def filenames(self):
        del self._filenames


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
