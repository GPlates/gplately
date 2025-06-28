#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
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

"""This sub-module contains functions for manipulating GPML (`.gpml`, `.gpmlz`) files,
as well as `pygplates.Feature` and `pygplates.FeatureCollection` objects.
"""
import os

import pygplates

__all__ = [
    "create_feature_dict",
    "extract_feature",
    "get_topological_references",
    "is_topological",
]


def extract_feature(id, features):
    r"""Return the feature with the given feature ID.

    Searches through `features` and returns the feature with a feature ID
    matching `id`. The order of `features` is not preserved, so multiple
    features sharing the same feature ID should be avoided.

    Parameters
    ----------
    id : pygplates.FeatureId or str
        The feature ID (or string representation) to search for.
    features : valid argument for pygplates.FeaturesFunctionArgument
        `features` may be a single `pygplates.Feature`, a
        `pygplates.FeatureCollection`, a `str` filename,
        or a (potentially nested) sequence of any combination of the above.

    Returns
    -------
    pygplates.Feature or None
        The matching feature if found, otherwise None.

    Raises
    ------
    TypeError
        If `id` is not one of {`pygplates.FeatureId`, `str`},
        or if `features` is of an invalid type.
    OSError (including FileNotFoundError, IsADirectoryError)
        If `features` points to a non-existent or unrecognised file.
    """
    features = _parse_features_function_arguments(features)

    if isinstance(id, pygplates.FeatureId):
        id_type = "fid"
    elif isinstance(id, str):
        id_type = "str"
    else:
        raise TypeError("Invalid feature ID: `{}`".format(id))

    for feature in features:
        feature_id = feature.get_feature_id()
        if id_type == "fid":
            cond = id == feature_id
        else:
            cond = id == feature_id.get_string()
        if cond:
            return feature

    # If the requested feature is not found, either return None
    # or maybe raise a ValueError like list.index() does
    # e.g. [1, 2, 3].index(4)
    # raise ValueError("Feature ID `{}` is not in `features`".format(id))
    return None


def create_feature_dict(features, id_type=pygplates.FeatureId):
    r"""Create a dictionary mapping feature IDs to features.

    Feature IDs can be either `str` or `pygplates.FeatureId`,
    according to `id_type`.

    Parameters
    ----------
    features : valid argument for pygplates.FeaturesFunctionArgument
        `features` may be a single `pygplates.Feature`, a
        `pygplates.FeatureCollection`, a `str` filename,
        or a (potentially nested) sequence of any combination of the above.
    id_type : {`pygplates.FeatureId`, `str`}, optional
        By default, dictionary keys will be of type `pygplates.FeatureId`;
        pass `id_type=str` to use string representations instead.

    Returns
    -------
    dict
        Dictionary for looking up features by their feature ID.
        N.B. the result will be an empty dictonary  if `features` is empty.

    Raises
    ------
    TypeError
        If `id_type` is not one of {`pygplates.FeatureId`, `str`},
        or if `features` is of an invalid type.
    OSError (including FileNotFoundError, IsADirectoryError)
        If `features` points to a non-existent or unrecognised file.
    """
    if id_type not in {str, pygplates.FeatureId}:
        raise TypeError(
            "Invalid `key_type` value: `"
            + str(id_type)
            + "` (must be one of `pygplates.FeatureId`, `str`)"
        )

    features = _parse_features_function_arguments(features)
    if id_type is str:
        return {i.get_feature_id().get_string(): i for i in features}
    else:
        return {i.get_feature_id(): i for i in features}


def get_topological_references(features, id_type=pygplates.FeatureId):
    r"""Create a dictionary mapping topological feature IDs to
    referenced features.

    The resulting dictionary maps each topological
    `pygplates.FeatureId` in `features` to the set of
    `pygplates.FeatureId` referenced by its topological
    geometry or geometries.

    Parameters
    ----------
    features : valid argument for pygplates.FeaturesFunctionArgument
        `features` may be a single `pygplates.Feature`, a
        `pygplates.FeatureCollection`, a `str` filename,
        or a (potentially nested) sequence of any combination of the above.
    id_type : {`pygplates.FeatureId`, `str`}, optional
        By default, feature IDs will be of type `pygplates.FeatureId`;
        pass `id_type=str` to use string representations instead.

    Returns
    -------
    dict
        N.B. Dictionary keys include all topological features in `features`,
        and only topological features.

    Raises
    ------
    TypeError
        If `id_type` is not one of {`pygplates.FeatureId`, `str`},
        or if `features` is of an invalid type.
    OSError (including FileNotFoundError, IsADirectoryError)
        If `features` points to a non-existent or unrecognised file.
    """
    features = _parse_features_function_arguments(features)
    features_dict = create_feature_dict(features)
    results = {}
    for feature in features:
        if not is_topological(feature):
            continue
        references = _get_topological_references(feature, features_dict)
        if len(references) == 0:
            continue
        id = feature.get_feature_id()
        if id_type is str:
            id = id.get_string()
            references = set(i.get_string() for i in references)
        results[id] = references
    return results


def is_topological(feature):
    r"""Determine whether a feature contains a topological geometry.

    Parameters
    ----------
    feature : pygplates.Feature

    Returns
    -------
    bool
        True if `feature` contains a topological geometry, else False.

    Raises
    ------
    TypeError
        If `feature` is a type other than `pygplates.Feature`
        (more precisely, if its type does not implement
        `get_all_topological_geometries`).
    """
    try:
        return len(feature.get_all_topological_geometries()) > 0
    except AttributeError as e:
        raise TypeError(
            "`is_topological` not implemented for `{}`".format(type(feature))
        ) from e


def _get_topological_references(feature, features_dict):
    if not isinstance(feature, pygplates.Feature):
        raise TypeError("Invalid feature type: `{}`".format(type(feature)))
    topologies = feature.get_all_topological_geometries()
    referenced_ids = set()
    for topology in topologies:
        if isinstance(topology, pygplates.GpmlTopologicalLine):
            sections = topology.get_sections()
        else:
            sections = topology.get_boundary_sections()
        for section in sections:
            feature_id = section.get_property_delegate().get_feature_id()
            referenced_ids.add(feature_id)
    referenced_features = set()
    for id in referenced_ids:
        if id in features_dict:
            referenced_features.add(features_dict[id])
    for referenced_feature in referenced_features:
        referenced_ids = referenced_ids.union(
            _get_topological_references(referenced_feature, features_dict)
        )
    return referenced_ids


def _parse_features_function_arguments(features):
    r"""Load features using `pygplates.FeaturesFunctionArgument`.

    This function also tries to translate some of the exceptions
    raised by `pygplates.FeaturesFunctionArgument` into regular
    Python exceptions.
    """
    try:
        features = pygplates.FeatureCollection(
            pygplates.FeaturesFunctionArgument(features).get_features()
        )
    except pygplates.FileFormatNotSupportedError as e:
        print("Invalid filename: `{}`".format(features))
        if not os.path.exists(features):
            raise FileNotFoundError("File does not exist") from e
        if os.path.isdir(features):
            raise IsADirectoryError("File is a directory") from e
        raise OSError("Unrecognised file format") from e
    except pygplates.OpenFileForReadingError as e:
        raise FileNotFoundError(
            "Could not find input file(s): `{}`".format(features)
        ) from e
    except Exception as e:
        if str(type(e)) == "<class 'Boost.Python.ArgumentError'>":
            # This is the easiest way of catching Boost.Python.ArgumentError,
            # since it cannot be directly imported into Python
            raise TypeError("Invalid argument type: `{}`".format(type(features))) from e
        else:
            raise e
    return features


def _load_FeatureCollection(features_or_files):
    """Return a pygplates.FeatureCollection containing features loaded from one or more files or pygplates.FeatureCollection(s)."""
    if features_or_files is None:
        return None

    if not pygplates.FeaturesFunctionArgument.contains_features(features_or_files):
        raise TypeError(
            "Expected a feature collection, or filename, or feature, or sequence of features, "
            "or a sequence (eg, list or tuple) of any combination of those four types."
        )

    return pygplates.FeatureCollection(
        pygplates.FeaturesFunctionArgument(features_or_files).get_features()
    )
