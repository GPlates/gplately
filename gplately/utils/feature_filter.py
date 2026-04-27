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

import pygplates  # type: ignore

logger = logging.getLogger("gplately")


class FeatureFilter(metaclass=abc.ABCMeta):
    _filtrate_feature_collection = None
    _residue_feature_collection = None

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

    @property
    def filtrate_feature_collection(self) -> Union[pygplates.FeatureCollection, None]:  # type: ignore
        """the feature collection of features that passed the filter."""
        return self._filtrate_feature_collection

    @property
    def residue_feature_collection(self) -> Union[pygplates.FeatureCollection, None]:  # type: ignore
        """the feature collection of features that did not pass the filter."""
        return self._residue_feature_collection


class FeatureNameFilter(FeatureFilter):
    """filter features by name, keep features with name matching any of the specified strings by default

    for example:
        FeatureNameFilter(['Africa', 'Asia']) -- keep features who name contains 'Africa' or 'Asia'
        FeatureNameFilter(['Africa', 'Asia'], reverse=True) -- keep features who name does not contain 'Africa' or 'Asia'
        FeatureNameFilter(['Africa', 'Asia'], exact_match=True) -- keep features who name is 'Africa' or 'Asia'
        FeatureNameFilter(['Africa', 'Asia'], reverse=True, exact_match=True) -- keep features who name is not 'Africa' or 'Asia'
        FeatureNameFilter(['Africa', 'Asia'], reverse=True, exact_match=True, case_sensitive=True) -- keep features who name is not 'Africa' or 'Asia' (case sensitive)
    """

    def __init__(
        self, names: List[str], exact_match=False, case_sensitive=False, reverse=False
    ):
        """filter features by name, keep features with name matching any of the specified strings by default

        :param names: a list of strings to match the feature name
        :param exact_match: if True, the feature name must be exactly the same as one of the strings in names;
            if False, the feature name only needs to contain one of the strings in names
        :param case_sensitive: if True, the name matching is case sensitive; if False, the name matching is case insensitive
        :param reverse: if False, feature with name matching one of the strings in the parameter "names" will pass the filter,
            if True, feature with name not matching any of the strings in the parameter "names" will pass the filter.
        """
        self._names = names
        self._exact_match = exact_match
        self._case_sensitive = case_sensitive
        self._reverse = reverse
        self._filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        self._residue_feature_collection = pygplates.FeatureCollection()  # type: ignore

    def _check_name(self, name_1: str, name_2: str) -> bool:
        """check if two names are the same or name_2 contains name_1"""
        if not self._case_sensitive:
            name_1_tmp = name_1.lower()
            name_2_tmp = name_2.lower()
        else:
            name_1_tmp = name_1
            name_2_tmp = name_2
        if self._exact_match:
            return name_1_tmp == name_2_tmp
        else:
            return name_1_tmp in name_2_tmp

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        if self._reverse:
            for name in self._names:
                if self._check_name(name, feature.get_name()):  # type: ignore
                    self._residue_feature_collection.add(feature)  # type: ignore
                    return False
            self._filtrate_feature_collection.add(feature)  # type: ignore
            return True
        else:
            for name in self._names:
                if self._check_name(name, feature.get_name()):  # type: ignore
                    self._filtrate_feature_collection.add(feature)  # type: ignore
                    return True
            self._residue_feature_collection.add(feature)  # type: ignore
            return False


class PlateIDFilter(FeatureFilter):
    """filter features by plate ID, keep features with plate ID in a specified list by default

    for example:
        PlateIDFilter([101,201,301]) -- keep features whose plate id is 101 or 201 or 301
        PlateIDFilter([101,201,301], reverse=True) -- keep features whose plate id is not 101 nor 201 nor 301

    """

    def __init__(self, pids: List[int], reverse=False):
        """filter features by plate ID, keep features with plate ID in a specified list by default

        :param pids: a list of plate IDs to match
        :param reverse: if False, feature with plate ID matching one of the strings in the parameter "pids" will pass the filter,
            if True, feature with plate ID not matching any of the strings in the parameter "pids" will pass the filter.
        """
        self._pids = pids
        self._reverse = reverse
        self._filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        self._residue_feature_collection = pygplates.FeatureCollection()  # type: ignore

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        if not self._reverse:
            if feature.get_reconstruction_plate_id() in self._pids:  # type: ignore
                self._filtrate_feature_collection.add(feature)
                return True
            else:
                self._residue_feature_collection.add(feature)  # type: ignore
                return False
        else:
            if feature.get_reconstruction_plate_id() not in self._pids:  # type: ignore
                self._filtrate_feature_collection.add(feature)  # type: ignore
                return True
            else:
                self._residue_feature_collection.add(feature)  # type: ignore
                return False


class BirthAgeFilter(FeatureFilter):
    """filter features by the time of appearance, keep older features by default

    for example:
        BirthAgeFilter(500) -- keep features which appreared before 500 Ma
        BirthAgeFilter(500, reverse=False) --  keep features which appreared after 500 Ma
    """

    def __init__(self, age: float, reverse=False):
        """filter features by the time of appearance, keep older features by default

        :param age: the age criterion
        :param reverse: if False, the feature older than the age criterion will pass the filter.
            If True, the feature younger than the age criterion will pass the filter.
        """

        self._age = age
        self._reverse = reverse
        self._filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        self._residue_feature_collection = pygplates.FeatureCollection()  # type: ignore

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        valid_time = feature.get_valid_time(None)
        if valid_time:
            begin_time, _ = valid_time
            if not self._reverse:  # keep older features
                if (
                    begin_time == pygplates.GeoTimeInstant.create_distant_past()  # type: ignore
                    or begin_time >= self._age
                ):
                    self._filtrate_feature_collection.add(feature)  # type: ignore
                    return True
            else:  # keep younger features
                if (
                    begin_time != pygplates.GeoTimeInstant.create_distant_past()  # type: ignore
                    and begin_time <= self._age
                ):
                    self._filtrate_feature_collection.add(feature)  # type: ignore
                    return True
        else:
            logger.warning(
                f"Feature {feature.get_feature_type()} {feature.get_name()} does not have 'valid time', and will not pass BirthAgeFilter."  # type: ignore
            )
        self._residue_feature_collection.add(feature)  # type: ignore
        return False


class EndTimeFilter(FeatureFilter):
    """filter features by the time of disappearance, keep features that disappeared before a certain time by default

    for example:
        EndTimeFilter(100) -- keep features which disappeared before 100 Ma
        EndTimeFilter(100, reverse=False) --  keep features which disappeared after 100 Ma, including features that have not disappeared yet (end time is distant future)
    """

    def __init__(self, age: float, reverse=False):
        """filter features by the time of disappearance, keep features that disappeared before a certain time by default

        :param age: the age criterion
        :param reverse: if False, feature disappeared before the age criterion will pass the filter.
            If True, feature disappeared after the age criterion will pass the filter, including features that have not disappeared yet (end time is distant future).
        """
        self._age = age
        self._reverse = reverse
        self._filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        self._residue_feature_collection = pygplates.FeatureCollection()  # type: ignore

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        valid_time = feature.get_valid_time(None)
        if valid_time:
            _, end_time = valid_time
            if not self._reverse:  # keep features disappeared before the age criterion
                if (
                    end_time != pygplates.GeoTimeInstant.create_distant_future()  # type: ignore
                    and end_time >= self._age
                ):
                    self._filtrate_feature_collection.add(feature)  # type: ignore
                    return True
            else:  # keep features disappeared after the age criterion, including features that have not disappeared yet (end time is distant future)
                if (
                    end_time == pygplates.GeoTimeInstant.create_distant_future()  # type: ignore
                    or end_time <= self._age
                ):
                    self._filtrate_feature_collection.add(feature)  # type: ignore
                    return True
        else:
            logger.warning(
                f"Feature {feature.get_feature_type()} {feature.get_name()} does not have 'valid time', and will not pass EndTimeFilter."  # type: ignore
            )
        self._residue_feature_collection.add(feature)  # type: ignore
        return False


class FeatureTypeFilter(FeatureFilter):
    """filter features by the feature type, keep features with feature type matching a specified regular expression by default"""

    def __init__(self, feature_type_re: str, reverse=False):
        self._feature_type_re = feature_type_re
        self._reverse = reverse
        self._filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        self._residue_feature_collection = pygplates.FeatureCollection()  # type: ignore

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        feature_type = feature.get_feature_type()
        if (
            not self._reverse
        ):  # keep features with feature type matching the regular expression
            if re.fullmatch(self._feature_type_re, feature_type.to_qualified_string()):
                self._filtrate_feature_collection.add(feature)  # type: ignore
                return True
            else:
                self._residue_feature_collection.add(feature)  # type: ignore
                return False
        else:  # keep features with feature type not matching the regular expression
            if not re.fullmatch(
                self._feature_type_re, feature_type.to_qualified_string()
            ):
                self._filtrate_feature_collection.add(feature)  # type: ignore
                return True
            else:
                self._residue_feature_collection.add(feature)  # type: ignore
                return False


class PropertyExistsFilter(FeatureFilter):
    """filter features by the existence of a property, keep features that have a specific property by default

    Depending on the value of reverse, this filter can be used to either keep features that have a specific property,
    or keep features that do not have that property.

    For example, if we want to keep features that do not have the property "gpml:subductionPolarity",
    we can use PropertyExistsFilter("gpml:subductionPolarity", reverse=True).
    """

    def __init__(self, property_name: str, reverse=False):
        """filter features by the existence of a property, keep features that have a specific property by default

        :param property_name: the name of the property to check (case-insensitive), such as gpml:subductoinPolarity
        :param reverse: if False, features that have the property will pass the filter;
            if True, features that do not have the property will pass the filter.
        """

        self._property_name = property_name
        self._reverse = reverse
        self._filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        self._residue_feature_collection = pygplates.FeatureCollection()  # type: ignore

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        for property in feature:
            if (
                property.get_name().to_qualified_string().lower()
                == self._property_name.lower()
            ):
                if not self._reverse:
                    self._filtrate_feature_collection.add(feature)  # type: ignore
                    return True
                else:
                    self._residue_feature_collection.add(feature)  # type: ignore
                    return False

        # if we go through all properties and do not find the property we are looking for
        if self._reverse:
            self._filtrate_feature_collection.add(feature)  # type: ignore
            return True
        else:
            self._residue_feature_collection.add(feature)  # type: ignore
            return False


class PropertyValueFilter(FeatureFilter):
    """filter features by the value of a property, keep features that have a specific property with a specific value by default

    Depending on the value of not_match, this filter can be used to either keep features that have the specific property and
    its value matches the specified value, or keep features that either do not have that property or its value does not match
    the specified value.

    For example, if we want to keep features whose "gpml:subductionPolarity" property value is not "Unknown",
    we can use PropertyValueFilter("gpml:subductionPolarity", "Unknown", reverse=True).
    """

    def __init__(
        self,
        property_name: str,
        property_value,
        reverse=False,
    ):
        """filter features by the value of a property, keep features that have a specific property with a specific value by default

        :param property_name: the name of the property to check (case-insensitive), such as gpml:subductoinPolarity
        :param property_value: the value of the property to check
        :param reverse: if False, features which has a property and its value match the specified value will pass the filter;
            if True, features which either do not have the property or have the property but its value does not match the specified value will pass the filter.
        """
        self._property_name = property_name
        self._property_value = property_value
        self._reverse = reverse
        self._filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        self._residue_feature_collection = pygplates.FeatureCollection()  # type: ignore

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        for property in feature:
            if (
                property.get_name().to_qualified_string().lower()
                == self._property_name.lower()
                and property.get_value() == self._property_value
            ):
                if self._reverse:
                    self._residue_feature_collection.add(feature)  # type: ignore
                    return False
                else:
                    self._filtrate_feature_collection.add(feature)  # type: ignore
                    return True

        # if we go through all properties and do not find the property with the specified value
        if self._reverse:
            self._filtrate_feature_collection.add(feature)  # type: ignore
            return True
        else:
            self._residue_feature_collection.add(feature)  # type: ignore
            return False


class RegionOfInterestFilter(FeatureFilter):
    """filter features by whether they are inside a region of interest, which can be defined by a bounding box or a polygon.
    By default, features that are inside the region of interest will pass the filter.

    for example:
        RegionOfInterestFilter(polygon) -- keep features that are in the polygon
        RegionOfInterestFilter(polygon, reverse=True) -- keep features that are not in the polygon
    """

    def __init__(
        self,
        region_of_interest: Union[Tuple[float, float, float, float], pygplates.PolygonOnSphere],  # type: ignore
        reverse=False,
    ):
        """filter features by whether they are inside a region of interest, which can be defined by a bounding box or a polygon.
        By default, features that are inside the region of interest will pass the filter.

        :param region_of_interest: the region of interest, can be either a tuple of (left, right, bottom, top) or a pygplates.PolygonOnSphere
        :param reverse: if False, features that are within the region of interest will pass the filter;
            if True, features that are not inside the region of interest will pass the filter.
        """
        if isinstance(region_of_interest, pygplates.PolygonOnSphere):  # type: ignore
            self._region_of_interest = region_of_interest
        elif isinstance(region_of_interest, tuple) and len(region_of_interest) == 4:
            (left, right, bottom, top) = region_of_interest
            self._region_of_interest = pygplates.PolygonOnSphere((lat, lon) for lon, lat in [(left, bottom), (left, top), (right, top), (right, bottom)])  # type: ignore
        else:
            raise ValueError(
                "region_of_interest should be either a tuple of (left, right, bottom, top) or a pygplates.PolygonOnSphere"
            )

        self._reverse = reverse
        self._filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        self._residue_feature_collection = pygplates.FeatureCollection()  # type: ignore

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        # Implementation for checking if feature is in the region of interest
        roi_feature = pygplates.Feature()  # type: ignore
        roi_feature.set_geometry(self._region_of_interest)  # type: ignore
        plate_partitioner = pygplates.PlatePartitioner(pygplates.FeatureCollection([roi_feature]), pygplates.RotationModel([]))  # type: ignore
        inside_geometries = []
        outside_geometries = []
        plate_partitioner.partition_geometry(feature.get_geometries(), inside_geometries, outside_geometries)  # type: ignore

        # TODO:
        # we may want to add an option to specify the minimum area of the inside/outside geometry to be considered as inside/outside the region of interest,
        # because for some features, such as large polygons, they may have a small portion of their geometry inside the region of interest,
        # and we may want to consider them as not inside the region of interest.
        if not self._reverse:
            if len(inside_geometries) == 0:
                self._residue_feature_collection.add(feature)  # type: ignore
                return False
            else:
                self._filtrate_feature_collection.add(feature)  # type: ignore
                return True
        else:
            if len(outside_geometries) == 0:
                self._residue_feature_collection.add(feature)  # type: ignore
                return False
            else:
                self._filtrate_feature_collection.add(feature)  # type: ignore
                return True


def filter_feature_collection(
    feature_collection: pygplates.FeatureCollection, filters: List[FeatureFilter]  # type: ignore
):
    """Filter a feature collection using a list of filters.

    :param feature_collection: the input feature collection to be filtered
    :param filters: a list of filters to apply. A feature will be kept if it passes all the filters in the list.
    :returns: the filtered feature collection
    """
    for filter in filters:
        filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        for feature in feature_collection:
            if filter.should_keep(feature):
                filtrate_feature_collection.add(feature)
        feature_collection = filtrate_feature_collection

    return feature_collection
