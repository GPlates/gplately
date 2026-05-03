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
    @classmethod
    def __subclasshook__(cls, subclass):
        return (
            hasattr(subclass, "transform")
            and callable(subclass.transform)
            or NotImplemented
        )

    @abc.abstractmethod
    def transform(self, feature: pygplates.Feature) -> pygplates.Feature:  # type: ignore
        """This abstract method must be implemented in subclass.

        :param feature: pygplates.Feature

        :returns: new pygplates.Feature after transformation
        """

        raise NotImplementedError


def transform_feature_collection(
    feature_collection: pygplates.FeatureCollection, transformers: List[FeatureTransformer]  # type: ignore
):
    """Transform a feature collection with a list of transformers."""
    for feature in feature_collection:
        for transformer in transformers:
            feature = transformer.transform(feature)
    return feature_collection


class SetReconstructionPlateIDTransformer(FeatureTransformer):
    """A feature transformer that sets the reconstruction plate ID of a feature to a specified value."""

    def __init__(self, plate_id: int):
        self._plate_id = plate_id

    def transform(self, feature: pygplates.Feature) -> pygplates.Feature:  # type: ignore
        feature.set_reconstruction_plate_id(self._plate_id)  # type: ignore
        return feature


class SetValidTimeTransformer(FeatureTransformer):
    """A feature transformer that sets the valid time of a feature to a specified value."""

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
