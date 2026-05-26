#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
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
import abc
import pygplates  # type: ignore
from typing import List, Union


class FeatureTransformer(metaclass=abc.ABCMeta):
    """Abstract base class for feature transformers."""

    @classmethod
    def __subclasshook__(cls, subclass):
        return (
            hasattr(subclass, "transform")
            and callable(subclass.transform)
            or NotImplemented
        )

    @abc.abstractmethod
    def transform(self, feature: pygplates.Feature) -> pygplates.Feature:  # type: ignore
        """Transform a single feature.

        This abstract method must be implemented by all subclasses.
        Implementations may modify the feature in place and/or return
        a new feature object.

        Parameters
        ----------
        feature : pygplates.Feature
            The feature to transform.

        Returns
        -------
        pygplates.Feature
            The transformed feature (may be the same object or a new instance).
        """

        raise NotImplementedError


def transform_feature_collection(
    feature_collection: pygplates.FeatureCollection, transformers: List[FeatureTransformer]  # type: ignore
):
    """Transform a feature collection with a list of transformers.

    Apply a sequence of transformers to each feature in the collection.
    Transformers are applied sequentially and update features in place.

    Parameters
    ----------
    feature_collection : pygplates.FeatureCollection
        The feature collection to be transformed.
    transformers : List[FeatureTransformer]
        List of FeatureTransformer instances to be applied to each feature
        in the collection, in order.

    Returns
    -------
    pygplates.FeatureCollection
        The transformed feature collection (same object, modified in place).
    """
    transformed_feature_collection = pygplates.FeatureCollection()  # type: ignore
    for feature in feature_collection:
        for transformer in transformers:
            transformed_feature = transformer.transform(feature)
            transformed_feature_collection.add(transformed_feature)
    # If the transformers modify features in place,
    # the features in the input feature collection is already updated.
    # The transformed_feature_collection is just a new collection containing the same feature objects.
    # If transformers return new feature objects(not modifying in place),
    # we return a new feature collection containing the transformed new features.
    return transformed_feature_collection


class SetReconstructionPlateIDTransformer(FeatureTransformer):
    """Set the reconstruction plate ID of a feature to a specified value.

    This transformer updates the reconstruction plate ID property of features,
    which is essential for plate kinematic reconstructions.

    Parameters
    ----------
    plate_id : int
        The reconstruction plate ID to assign to features.
        Typically a positive integer (0 is valid but represents an unspecified plate).

    Examples
    --------
    Create a transformer that sets plate ID to 801:

    >>> transformer = SetReconstructionPlateIDTransformer(plate_id=801)

    See Also
    --------
    SetValidTimeTransformer : Transformer for setting valid time intervals
    """

    def __init__(self, plate_id: int):
        self._plate_id = plate_id

    def transform(self, feature: pygplates.Feature) -> pygplates.Feature:  # type: ignore
        feature.set_reconstruction_plate_id(self._plate_id)  # type: ignore
        return feature


class SetValidTimeTransformer(FeatureTransformer):
    """Set the valid time of a feature to specified begin and end times.

    This transformer allows partial updates to feature valid time intervals.
    If only one time bound is specified, the existing bound is preserved.

    Parameters
    ----------
    begin_time : float or None, optional
        The begin time (valid from) for the feature in Ma (millions of years ago).
        If None, the existing begin time is preserved. Default is None.
    end_time : float or None, optional
        The end time (valid until) for the feature in Ma (millions of years ago).
        If None, the existing end time is preserved. Default is None.

    Examples
    --------
    Set both begin and end times:

    >>> transformer = SetValidTimeTransformer(begin_time=100.0, end_time=50.0)

    Set only the begin time, preserving the existing end time:

    >>> transformer = SetValidTimeTransformer(begin_time=100.0)

    See Also
    --------
    SetReconstructionPlateIDTransformer : Transformer for setting reconstruction plate ID
    """

    def __init__(self, begin_time: Union[float, None] = None, end_time: Union[float, None] = None):  # type: ignore
        self._begin_time = begin_time
        self._end_time = end_time

    def transform(self, feature: pygplates.Feature) -> pygplates.Feature:  # type: ignore
        if self._begin_time is None and self._end_time is None:
            return feature
        begin_time = float("inf")
        end_time = float("-inf")
        if self._begin_time is None:
            valid_time = feature.get_valid_time(None)  # type: ignore
            if valid_time:
                begin_time = valid_time[0]
        else:
            begin_time = self._begin_time
        if self._end_time is None:
            valid_time = feature.get_valid_time(None)  # type: ignore
            if valid_time:
                end_time = valid_time[1]
        else:
            end_time = self._end_time
        feature.set_valid_time(begin_time, end_time)  # type: ignore
        return feature
