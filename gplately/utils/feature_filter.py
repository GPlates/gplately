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
from typing import List, Tuple, Union

import pygplates
import shapely  # type: ignore

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
                    begin_time == pygplates.GeoTimeInstant.create_distant_past()  # type: ignore
                    or begin_time > self.age
                ):
                    return True
            else:
                if (
                    begin_time != pygplates.GeoTimeInstant.create_distant_past()  # type: ignore
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
        EndTimeFilter(100) -- keep features which disappeared before 100 Ma
        EndTimeFilter(100, disappear_before=False) --  keep features which disappeared after 100 Ma, including features that have not disappeared yet (end time is distant future)

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
                    end_time != pygplates.GeoTimeInstant.create_distant_future()  # type: ignore
                    and end_time > self.age
                ):
                    return True
            else:
                if (
                    end_time == pygplates.GeoTimeInstant.create_distant_future()  # type: ignore
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


class PropertyExistsFilter(FeatureFilter):
    """filter features by the existence of a property

    Depending on the value of not_exists, this filter can be used to either keep features that have a specific property, or keep features that do not have that property. For example, if we want to keep features that do not have the property "gpml:subductionPolarity", we can use PropertyExistsFilter("gpml:subductionPolarity", not_exists=True).

    :param property_name: the name of the property to check (case-insensitive), such as gpml:subductoinPolarity
    :param not_exists: if True, keep features that do not have the property; if False, keep features that have the property.

    """

    def __init__(self, property_name: str, not_exists=False):
        self._property_name = property_name
        self._not_exists = not_exists

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        for property in feature:
            if self._not_exists:
                if (
                    property.get_name().to_qualified_string().lower()
                    == self._property_name.lower()
                ):
                    return False
            else:
                if (
                    property.get_name().to_qualified_string().lower()
                    == self._property_name.lower()
                ):
                    return True
        if self._not_exists:
            return True
        else:
            return False


class PropertyValueFilter(FeatureFilter):
    """filter features by the value of a property

    Depending on the value of not_match, this filter can be used to either keep features that have the specific property and its value matches the specified value, or keep features that either do not have that property or its value does not match the specified value.

    For example, if we want to keep features whose "gpml:subductionPolarity" property value is not "Unknown", we can use PropertyValueFilter("gpml:subductionPolarity", "Unknown", not_match=True).

    :param property_name: the name of the property to check (case-insensitive), such as gpml:subductoinPolarity
    :param property_value: the value of the property to check
    :param not_match: if True, keep features whose property value does not match the specified value, also keep features without the specified property; if False, keep features that have at least one property which has the specified value.

    """

    def __init__(
        self,
        property_name: str,
        property_value,
        not_match=False,
    ):
        self._property_name = property_name
        self._property_value = property_value
        self._not_match = not_match

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        for property in feature:
            if (
                property.get_name().to_qualified_string().lower()
                == self._property_name.lower()
            ):
                if self._not_match:
                    if property.get_value() == self._property_value:
                        return False
                else:
                    if property.get_value() == self._property_value:
                        return True
        if self._not_match:
            return True
        else:
            return False


class RegionOfInterestFilter(FeatureFilter):
    """filter features by whether they are in a region of interest

    for example:
        RegionOfInterestFilter(polygon) -- keep features that are in the polygon
        RegionOfInterestFilter(polygon, exclude=True) -- keep features that are not in the polygon

    """

    def __init__(
        self,
        region_of_interest: Union[Tuple[float, float, float, float], pygplates.PolygonOnSphere],  # type: ignore
        exclude=False,
    ):
        if isinstance(region_of_interest, pygplates.PolygonOnSphere):  # type: ignore
            self._region_of_interest = region_of_interest
        elif isinstance(region_of_interest, tuple) and len(region_of_interest) == 4:
            (left, right, bottom, top) = region_of_interest
            self._region_of_interest = pygplates.PolygonOnSphere((lat, lon) for lon, lat in [(left, bottom), (left, top), (right, top), (right, bottom)])  # type: ignore
        else:
            raise ValueError(
                "region_of_interest should be either a tuple of (left, right, bottom, top) or a pygplates.PolygonOnSphere"
            )

        self._exclude = exclude

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        # Implementation for checking if feature is in the region of interest
        roi_feature = pygplates.Feature()  # type: ignore
        roi_feature.set_geometry(self._region_of_interest)  # type: ignore
        plate_partitioner = pygplates.PlatePartitioner(pygplates.FeatureCollection([roi_feature]), pygplates.RotationModel([]))  # type: ignore
        inside_geometries = []
        outside_geometries = []
        plate_partitioner.partition_geometry(feature.get_geometries(), inside_geometries, outside_geometries)  # type: ignore
        if self._exclude:
            return len(inside_geometries) == 0
        else:
            return len(outside_geometries) == 0


def filter_feature_collection(
    feature_collection: pygplates.FeatureCollection, filters: List[FeatureFilter]  # type: ignore
):
    """Filter a feature collection using a list of filters.

    :param feature_collection: the input feature collection to be filtered
    :param filters: a list of filters to apply. A feature will be kept if it passes all the filters in the list.
    :returns: a new feature collection after filtering
    """
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
