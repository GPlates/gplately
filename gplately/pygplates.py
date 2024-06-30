#
#    Copyright (C) 2024 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

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
        elif hasattr(rotation_features, "__iter__") and all(
            isinstance(f, str) for f in rotation_features
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
    """A class that wraps the
    [`pyGPlates.FeatureCollection`](https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection)
    class. This aggregates a set of features into a collection.
    This is traditionally so that a group of  features can be loaded, saved or
    unloaded in a single operation.

    This wrapping of `pygplates.FeatureCollection` contains all
    [`pygplates.FeatureCollection` functionality](https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection),
    and in addition tracks the names of files from which the feature
    collection(s) are read using the
    'gplately.pygplates.FeatureCollection.filenames' attribute.

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
            msg = "\nFeatureCollection: No filename associated with {} in __init__".format(
                type(features)
            )
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
            msg = "\nFeatureCollection: No filename associated with {} in add".format(
                type(features)
            )
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
