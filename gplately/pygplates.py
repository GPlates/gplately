"""
A light wrapping of some pyGPlates classes to keep track of filenames

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

    def __init__(self, rotation_features):
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

    def __init__(self, feature):
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
        feat = super().clone()
        feat.filenames = self.filenames

class FeatureCollection(_pygplates.FeatureCollection):

    def __init__(self, features=None):
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
        fc = super().clone()
        fc.filenames = self.filenames
        return fc