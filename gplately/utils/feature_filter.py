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
import logging
import re
from typing import List

import pygplates  # type: ignore

logger = logging.getLogger("gplately")


class FeatureFilter(metaclass=abc.ABCMeta):
    @classmethod
    def __subclasshook__(cls, subclass):
        return (
            hasattr(subclass, "should_keep")
            and callable(subclass.should_keep)
            or NotImplemented
        )

    @abc.abstractmethod
    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        """This abstract method must be implemented in subclass.

        :param feature: pygplates.Feature

        :returns: true if the feature should be kept; false otherwise
        """

        raise NotImplementedError


class FeatureNameFilter(FeatureFilter):
    """filter features by name

    for example:
        FeatureNameFilter(['Africa', 'Asia']) -- keep features who name contains 'Africa' or 'Asia'
        FeatureNameFilter(['Africa', 'Asia'], exclude=True) -- keep features who name does not contain 'Africa' or 'Asia'
        FeatureNameFilter(['Africa', 'Asia'], exact_match=True) -- keep features who name is 'Africa' or 'Asia'
        FeatureNameFilter(['Africa', 'Asia'], exclude=True, exact_match=True) -- keep features who name is not 'Africa' or 'Asia'
        FeatureNameFilter(['Africa', 'Asia'], exclude=True, exact_match=True, case_sensitive=True) -- keep features who name is not 'Africa' or 'Asia' (case sensitive)
    """

    def __init__(
        self, names: List[str], exact_match=False, case_sensitive=False, exclude=False
    ):
        self.names = names
        self.exact_match = exact_match
        self.case_sensitive = case_sensitive
        self.exclude = exclude

    def check_name(self, name_1: str, name_2: str) -> bool:
        """check if two names are the same or name_2 contains name_1"""
        if not self.case_sensitive:
            name_1_tmp = name_1.lower()
            name_2_tmp = name_2.lower()
        else:
            name_1_tmp = name_1
            name_2_tmp = name_2
        if self.exact_match:
            return name_1_tmp == name_2_tmp
        else:
            return name_1_tmp in name_2_tmp

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        if self.exclude:
            for name in self.names:
                if self.check_name(name, feature.get_name()):  # type: ignore
                    return False
            return True
        else:
            for name in self.names:
                if self.check_name(name, feature.get_name()):  # type: ignore
                    return True
            return False


class PlateIDFilter(FeatureFilter):
    """filter features by plate ID

    for example:
        PlateIDFilter([101,201,301]) -- keep features whose plate id is 101 or 201 or 301
        PlateIDFilter([101,201,301], exclude=True) -- keep features whose plate id is not 101 nor 201 nor 301

    """

    def __init__(self, pids: List[int], exclude=False):
        self.pids = pids
        self.exclude = exclude

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        if not self.exclude and feature.get_reconstruction_plate_id() in self.pids:  # type: ignore
            return True
        if self.exclude and feature.get_reconstruction_plate_id() not in self.pids:  # type: ignore
            return True
        return False


class BirthAgeFilter(FeatureFilter):
    """filter features by the time of appearance

    for example:
        BirthAgeFilter(500) -- keep features which appreared before 500 Ma
        BirthAgeFilter(500, keep_older=False) --  keep features which appreared after 500 Ma

    :param age: the age criterion
    :param keep_older: if True, return True when the feature appeared before the age criterion. If False, otherwise.

    """

    def __init__(self, age: float, keep_older=True):
        self.age = age
        self.keep_older = keep_older

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        valid_time = feature.get_valid_time(None)
        if valid_time:
            begin_time, _ = valid_time
            if self.keep_older:
                if (
                    begin_time == pygplates.GeoTimeInstant.create_distant_past()
                    or begin_time > self.age
                ):
                    return True
            else:
                if (
                    begin_time != pygplates.GeoTimeInstant.create_distant_past()
                    and begin_time < self.age
                ):
                    return True
        else:
            logger.warning(
                f"Feature {feature.get_feature_type()} {feature.get_name()} does not have valid time, and will be excluded by BirthAgeFilter."  # type: ignore
            )
        return False


class EndTimeFilter(FeatureFilter):
    """filter features by the time of disappearance

    for example:
        EndTimeFilter(500) -- keep features which disappeared before 500 Ma
        EndTimeFilter(500, disappear_before=False) --  keep features which disappeared after 500 Ma, including features that have not disappeared yet (end time is distant future)

    :param age: the age criterion
    :param disappear_before: if True, return True when the feature disappeared before the age criterion. If False, otherwise.

    """

    def __init__(self, age: float, disappear_before=True):
        self.age = age
        self.disappear_before = disappear_before

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        valid_time = feature.get_valid_time(None)
        if valid_time:
            _, end_time = valid_time
            if self.disappear_before:
                if (
                    end_time != pygplates.GeoTimeInstant.create_distant_future()
                    and end_time > self.age
                ):
                    return True
            else:
                if (
                    end_time == pygplates.GeoTimeInstant.create_distant_future()
                    or end_time < self.age
                ):
                    return True
        else:
            logger.warning(
                f"Feature {feature.get_feature_type()} {feature.get_name()} does not have valid time, and will be excluded by EndTimeFilter."  # type: ignore
            )
        return False


class FeatureTypeFilter(FeatureFilter):
    """filter features by the feature type(regular expression)

    examples:
        - gplately filter input_file output_file -t gpml:Basin
        - gplately filter input_file output_file  -t "gpml:IslandArc|gpml:Basin"

    :param feature_type_re: the regular expression to match the feature type

    """

    def __init__(self, feature_type_re: str):
        self._feature_type_re = feature_type_re

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        feature_type = feature.get_feature_type()
        if re.fullmatch(self._feature_type_re, feature_type.to_qualified_string()):
            logger.debug(
                f"feature type match: {self._feature_type_re} {feature_type.to_qualified_string()}"
            )
            return True
        else:
            logger.debug(
                f"feature type not match: {self._feature_type_re} {feature_type.to_qualified_string()}"
            )
            return False


def filter_feature_collection(
    feature_collection: pygplates.FeatureCollection, filters: List[FeatureFilter]  # type: ignore
):
    """Filter feature collection by various criteria."""
    new_feature_collection = pygplates.FeatureCollection()  # type: ignore
    for feature in feature_collection:
        keep_flag = True
        for filter in filters:
            if not filter.should_keep(feature):
                keep_flag = False
                break
        if keep_flag:
            new_feature_collection.add(feature)
    return new_feature_collection
