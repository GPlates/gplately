"""
A light wrapping of some pyGPlates classes to keep track of filenames

Each object listed here will have a `self.filenames` attribute.
"""

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


import pygplates as _pygplates

class RotationModel(_pygplates.RotationModel):

    def __init__(self, rotation_features):
        super(RotationModel, self).__init__(rotation_features)
        self.filenames = rotation_features

class Feature(_pygplates.Feature):

    def __init__(self, feature):
        super(Feature, self).__init__(feature)
        self.filenames = [feature]

    def add(self, feature):
        super().add(feature)
        if isinstance(feature, Feature):
            self.filenames.extend(feature.filenames)
        elif _is_string(feature):
            self.filenames.extend(feature)
        else:
            raise ValueError("Unsupported type: ", type(feature))

    def clone(self):
        feat = super().clone()
        feat.filenames = self.filenames

class FeatureCollection(_pygplates.FeatureCollection):

    def __init__(self, features=None):
        super(FeatureCollection, self).__init__(features)

        # update filename list
        if _is_string(features) and type(features) is list:
            self.filenames = features
        elif _is_string(features) and type(features) is str:
            self.filenames = [features]
        elif features is None:
            self.filenames = []
        elif isinstance(features, FeatureCollection):
            self.filenames.extend(features.filenames)
        else:
            raise ValueError("Unsupported type:", type(features))

    def add(self, features):
        super().add(features)

        # update filename list
        if isinstance(features, FeatureCollection):
            self.filenames.extend(features.filenames)
        elif _is_string(features):
            self.filenames.extend(features)
        else:
            raise ValueError("Unsupported type: ", type(features))

    def clone(self):
        fc = super().clone()
        fc.filenames = self.filenames
        return fc