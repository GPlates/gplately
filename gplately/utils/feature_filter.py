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
from typing import List, Optional, Tuple, Union

import pygplates  # type: ignore

logger = logging.getLogger("gplately")


class FeatureFilter(metaclass=abc.ABCMeta):
    def __init__(self):
        self._filtrate_feature_collection: List[Optional[pygplates.Feature]] = []
        self._residue_feature_collection: List[Optional[pygplates.Feature]] = []

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

        Parameters
        ----------
        feature : pygplates.Feature

        Returns
        -------
        bool
            True if the feature should be kept; False otherwise.
        """

        raise NotImplementedError

    @property
    def filtrate_feature_collection(self):
        return pygplates.FeatureCollection(self._filtrate_feature_collection)

    @property
    def residue_feature_collection(self):
        return pygplates.FeatureCollection(self._residue_feature_collection)

    @property
    def filtrate_features_as_list(self):
        return self._filtrate_feature_collection

    @property
    def residue_features_as_list(self):
        return self._residue_feature_collection


def filter_feature_collection(
    feature_collection: pygplates.FeatureCollection, filters: List[FeatureFilter]  # type: ignore
):
    """Filter a feature collection using a list of filters.

    Parameters
    ----------
    feature_collection : pygplates.FeatureCollection
        The input feature collection to be filtered.
    filters : list of FeatureFilter
        A list of filters to apply. A feature will be kept if it passes all the filters in the list.

    Returns
    -------
    pygplates.FeatureCollection
        The filtered feature collection.
    """
    for filter in filters:
        filtrate_feature_collection = pygplates.FeatureCollection()  # type: ignore
        for feature in feature_collection:
            if filter.should_keep(feature):
                filtrate_feature_collection.add(feature)

        if filter.filtrate_feature_collection:
            feature_collection = filter.filtrate_feature_collection
        else:
            feature_collection = filtrate_feature_collection

    return feature_collection


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
        """Constructor for the feature name filter.

        Parameters
        ----------
        names : list of str
            A list of strings to match the feature name.
        exact_match : bool, optional
            If True, the feature name must be exactly the same as one of the strings in
            ``names``; if False, the feature name only needs to contain one of the strings.
            Default is False.
        case_sensitive : bool, optional
            If True, name matching is case-sensitive; if False, it is case-insensitive.
            Default is False.
        reverse : bool, optional
            If False, features with a name matching one of the strings in ``names`` will
            pass the filter; if True, features with a name not matching any string in
            ``names`` will pass the filter. Default is False.
        """
        super().__init__()
        self._names = names
        self._exact_match = exact_match
        self._case_sensitive = case_sensitive
        self._reverse = reverse

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
                    self._residue_feature_collection.append(feature)
                    return False
            self._filtrate_feature_collection.append(feature)
            return True
        else:
            for name in self._names:
                if self._check_name(name, feature.get_name()):  # type: ignore
                    self._filtrate_feature_collection.append(feature)
                    return True
            self._residue_feature_collection.append(feature)
            return False


class PlateIDFilter(FeatureFilter):
    """filter features by plate ID, keep features with plate ID in a specified list by default

    for example:
        PlateIDFilter([101,201,301]) -- keep features whose plate id is 101 or 201 or 301
        PlateIDFilter([101,201,301], reverse=True) -- keep features whose plate id is not 101 nor 201 nor 301

    """

    def __init__(self, pids: List[int], reverse=False):
        """Constructor for the plate ID filter.

        Parameters
        ----------
        pids : list of int
            A list of plate IDs to match.
        reverse : bool, optional
            If False, features with a plate ID matching one of the values in ``pids`` will
            pass the filter; if True, features with a plate ID not matching any value in
            ``pids`` will pass the filter. Default is False.
        """
        super().__init__()
        self._pids = pids
        self._reverse = reverse

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        if not self._reverse:
            if feature.get_reconstruction_plate_id() in self._pids:  # type: ignore
                self._filtrate_feature_collection.append(feature)
                return True
            else:
                self._residue_feature_collection.append(feature)
                return False
        else:
            if feature.get_reconstruction_plate_id() not in self._pids:  # type: ignore
                self._filtrate_feature_collection.append(feature)
                return True
            else:
                self._residue_feature_collection.append(feature)
                return False


class BirthAgeFilter(FeatureFilter):
    """filter features by the time of appearance, keep older features by default

    for example:
        BirthAgeFilter(500) -- keep features which appreared before 500 Ma
        BirthAgeFilter(500, reverse=False) --  keep features which appreared after 500 Ma
    """

    def __init__(self, age: float, reverse=False):
        """Constructor for the birth age filter.

        Parameters
        ----------
        age : float
            The age criterion.
        reverse : bool, optional
            If False, features older than the age criterion will pass the filter;
            if True, features younger than the age criterion will pass the filter.
            Default is False.
        """

        super().__init__()
        self._age = age
        self._reverse = reverse

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        valid_time = feature.get_valid_time(None)
        if valid_time:
            begin_time, _ = valid_time
            if not self._reverse:  # keep older features
                if (
                    begin_time == pygplates.GeoTimeInstant.create_distant_past()  # type: ignore
                    or begin_time >= self._age
                ):
                    self._filtrate_feature_collection.append(feature)
                    return True
            else:  # keep younger features
                if (
                    begin_time != pygplates.GeoTimeInstant.create_distant_past()  # type: ignore
                    and begin_time <= self._age
                ):
                    self._filtrate_feature_collection.append(feature)
                    return True
        else:
            logger.warning(
                f"Feature {feature.get_feature_type()} {feature.get_name()} does not have 'valid time', and will not pass BirthAgeFilter."  # type: ignore
            )
        self._residue_feature_collection.append(feature)
        return False


class EndTimeFilter(FeatureFilter):
    """filter features by the time of disappearance, keep features that disappeared before a certain time by default

    for example:
        EndTimeFilter(100) -- keep features which disappeared before 100 Ma
        EndTimeFilter(100, reverse=False) --  keep features which disappeared after 100 Ma, including features that have not disappeared yet (end time is distant future)
    """

    def __init__(self, age: float, reverse=False):
        """Constructor for the end time filter.

        Parameters
        ----------
        age : float
            The age criterion.
        reverse : bool, optional
            If False, features that disappeared before the age criterion will pass the filter;
            if True, features that disappeared after the age criterion will pass the filter,
            including features that have not disappeared yet (end time is distant future).
            Default is False.
        """
        super().__init__()
        self._age = age
        self._reverse = reverse

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        valid_time = feature.get_valid_time(None)
        if valid_time:
            _, end_time = valid_time
            if not self._reverse:  # keep features disappeared before the age criterion
                if (
                    end_time != pygplates.GeoTimeInstant.create_distant_future()  # type: ignore
                    and end_time >= self._age
                ):
                    self._filtrate_feature_collection.append(feature)
                    return True
            else:  # keep features disappeared after the age criterion, including features that have not disappeared yet (end time is distant future)
                if (
                    end_time == pygplates.GeoTimeInstant.create_distant_future()  # type: ignore
                    or end_time <= self._age
                ):
                    self._filtrate_feature_collection.append(feature)
                    return True
        else:
            logger.warning(
                f"Feature {feature.get_feature_type()} {feature.get_name()} does not have 'valid time', and will not pass EndTimeFilter."  # type: ignore
            )
        self._residue_feature_collection.append(feature)
        return False


class FeatureTypeFilter(FeatureFilter):
    """filter features by the feature type, keep features with feature type matching a specified regular expression by default"""

    def __init__(self, feature_type_re: str, reverse=False):
        super().__init__()
        self._feature_type_re = feature_type_re
        self._reverse = reverse

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        feature_type = feature.get_feature_type()
        if (
            not self._reverse
        ):  # keep features with feature type matching the regular expression
            if re.fullmatch(self._feature_type_re, feature_type.to_qualified_string()):
                self._filtrate_feature_collection.append(feature)
                return True
            else:
                self._residue_feature_collection.append(feature)
                return False
        else:  # keep features with feature type not matching the regular expression
            if not re.fullmatch(
                self._feature_type_re, feature_type.to_qualified_string()
            ):
                self._filtrate_feature_collection.append(feature)
                return True
            else:
                self._residue_feature_collection.append(feature)
                return False


class PropertyExistsFilter(FeatureFilter):
    """filter features by the existence of a property, keep features that have a specific property by default

    Depending on the value of reverse, this filter can be used to either keep features that have a specific property,
    or keep features that do not have that property.

    For example, if we want to keep features that do not have the property "gpml:subductionPolarity",
    we can use PropertyExistsFilter("gpml:subductionPolarity", reverse=True).
    """

    def __init__(self, property_name: str, reverse=False):
        """Constructor for the property existence filter.

        Parameters
        ----------
        property_name : str
            The name of the property to check (case-insensitive), e.g. ``gpml:subductionPolarity``.
        reverse : bool, optional
            If False, features that have the property will pass the filter;
            if True, features that do not have the property will pass the filter.
            Default is False.
        """

        super().__init__()
        self._property_name = property_name
        self._reverse = reverse

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        for property in feature:
            if (
                property.get_name().to_qualified_string().lower()
                == self._property_name.lower()
            ):
                if not self._reverse:
                    self._filtrate_feature_collection.append(feature)
                    return True
                else:
                    self._residue_feature_collection.append(feature)
                    return False

        # if we go through all properties and do not find the property we are looking for
        if self._reverse:
            self._filtrate_feature_collection.append(feature)
            return True
        else:
            self._residue_feature_collection.append(feature)
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
        """Constructor for the property value filter.

        Parameters
        ----------
        property_name : str
            The name of the property to check (case-insensitive), e.g. ``gpml:subductionPolarity``.
        property_value :
            The value of the property to check.
        reverse : bool, optional
            If False, features that have the property and whose value matches ``property_value``
            will pass the filter; if True, features that either do not have the property or have
            the property but its value does not match will pass the filter. Default is False.
        """
        super().__init__()
        self._property_name = property_name
        self._property_value = property_value
        self._reverse = reverse

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        for property in feature:
            if (
                property.get_name().to_qualified_string().lower()
                == self._property_name.lower()
                and property.get_value() == self._property_value
            ):
                if self._reverse:
                    self._residue_feature_collection.append(feature)
                    return False
                else:
                    self._filtrate_feature_collection.append(feature)
                    return True

        # if we go through all properties and do not find the property with the specified value
        if self._reverse:
            self._filtrate_feature_collection.append(feature)
            return True
        else:
            self._residue_feature_collection.append(feature)
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
        """Constructor for the region of interest filter.

        Parameters
        ----------
        region_of_interest : tuple of (float, float, float, float) or pygplates.PolygonOnSphere
            The region of interest. Either a ``(left, right, bottom, top)`` bounding box tuple
            or a ``pygplates.PolygonOnSphere``.
        reverse : bool, optional
            If False, features within the region of interest will pass the filter;
            if True, features outside the region of interest will pass the filter.
            Default is False.
        """
        super().__init__()
        if isinstance(region_of_interest, pygplates.PolygonOnSphere):  # type: ignore
            self._region_of_interest = region_of_interest
        elif isinstance(region_of_interest, tuple) and len(region_of_interest) == 4:
            left, right, bottom, top = region_of_interest
            self._region_of_interest = pygplates.PolygonOnSphere((lat, lon) for lon, lat in [(left, bottom), (left, top), (right, top), (right, bottom)])  # type: ignore
        else:
            raise ValueError(
                "region_of_interest should be either a tuple of (left, right, bottom, top) or a pygplates.PolygonOnSphere"
            )

        self._reverse = reverse

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        # Implementation for checking if feature is in the region of interest
        roi_feature = pygplates.Feature()  # type: ignore
        roi_feature.set_geometry(self._region_of_interest)  # type: ignore
        plate_partitioner = pygplates.PlatePartitioner(pygplates.FeatureCollection([roi_feature]), pygplates.RotationModel([]))  # type: ignore
        inside_geometries = []
        outside_geometries = []
        plate_partitioner.partition_geometry(feature.get_geometries(), inside_geometries, outside_geometries)  # type: ignore

        # TODO:
        # we may want to add an option to specify the minimum area of the inside/outside geometry
        # to be considered as inside/outside the region of interest,
        # because for some features, such as large polygons, they may have a portion of their geometry
        # inside the region of interest, and we may want to consider them as inside the region of interest.
        # and the same for ploylines, they may have a portion of their geometry inside the region of interest,
        #  and we may want to consider them as inside the region of interest.
        if not self._reverse:
            if len(inside_geometries) == 0:
                self._residue_feature_collection.append(feature)
                return False
            else:
                self._filtrate_feature_collection.append(feature)
                return True
        else:
            if len(outside_geometries) == 0:
                self._residue_feature_collection.append(feature)
                return False
            else:
                self._filtrate_feature_collection.append(feature)
                return True


class TopologicalFeaturesWithDuplicateSectionsFilter(FeatureFilter):
    """find features whose topological geometries have duplicated section features.

    No parameter is needed to create this filter.
    For the given feature, if any of its topological geometries has duplicated section features,
    the feature will pass the filter.
    """

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        geoms = feature.get_all_topological_geometries()
        for geom in geoms:
            _seen_feature_ids = set()
            if not isinstance(geom, pygplates.GpmlTopologicalLine):  # type: ignore
                for section in geom.get_boundary_sections():
                    feature_id = section.get_property_delegate().get_feature_id()
                    if feature_id in _seen_feature_ids:
                        logger.debug(
                            f"Duplicate feature found: {feature_id} in feature {feature.get_feature_id()}"
                        )
                        return True
                    else:
                        _seen_feature_ids.add(feature_id)
            else:
                for section in geom.get_sections():
                    feature_id = section.get_property_delegate().get_feature_id()
                    if feature_id in _seen_feature_ids:
                        logger.debug(
                            f"Duplicate feature found: {feature_id} in feature {feature.get_feature_id()}"
                        )
                        return True
                    else:
                        _seen_feature_ids.add(feature_id)

        return False


class TopologicalSectionFeaturesFilter(FeatureFilter):
    """find section features for given topological features

    Given a collection of topological features, this filter will find section features
    that are used in the topological geometries of these topological features
    from another feature collection containing the real geometries.
    """

    def __init__(self, topological_feature_collection=pygplates.FeatureCollection()):  # type: ignore
        """Construct the filter with a collection of topological features.

        The filtering will be applied to another collection of features containing the real
        geometries.

        Parameters
        ----------
        topological_feature_collection : pygplates.FeatureCollection, optional
            A collection of topological features. The filter will find section features that
            are referenced by the topological geometries of these features.
        """
        super().__init__()
        self._topological_feature_collection = topological_feature_collection
        self._the_ids_of_section_features = set()
        for feature in self._topological_feature_collection:
            for geom in feature.get_all_topological_geometries():
                if not isinstance(geom, pygplates.GpmlTopologicalLine):  # type: ignore
                    for section in geom.get_boundary_sections():
                        self._the_ids_of_section_features.add(
                            section.get_property_delegate()
                            .get_feature_id()
                            .get_string()
                        )
                else:
                    for section in geom.get_sections():
                        self._the_ids_of_section_features.add(
                            section.get_property_delegate()
                            .get_feature_id()
                            .get_string()
                        )

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        if feature.get_feature_id().get_string() in self._the_ids_of_section_features:
            return True
        else:
            return False


class FeatureIDFilter(FeatureFilter):
    """filter features by their feature IDs, keep features with feature ID in a specified list by default

    for example:
        FeatureIDFilter(['id1', 'id2', 'id3']) -- keep features whose feature ID is 'id1' or 'id2' or 'id3'
        FeatureIDFilter(['id1', 'id2', 'id3'], reverse=True) -- keep features whose feature ID is not 'id1' nor 'id2' nor 'id3'

    """

    def __init__(self, fids: List[str], reverse=False):
        """Constructor for the feature ID filter.

        Parameters
        ----------
        fids : list of str
            A list of feature IDs to match. Each feature ID should be a universally unique string.
        reverse : bool, optional
            If False, features with a feature ID matching one of the strings in ``fids`` will
            pass the filter; if True, features with a feature ID not matching any string in
            ``fids`` will pass the filter. Default is False.
        """
        super().__init__()
        self._fids = fids
        assert len(self._fids) == len(
            set(self._fids)
        ), "The feature IDs in the parameter 'fids' should be unique."
        self._reverse = reverse
        if not self._reverse:
            self._filtrate_feature_collection = [None] * len(self._fids)
        else:
            self._residue_feature_collection = [None] * len(self._fids)

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        try:
            i = self._fids.index(feature.get_feature_id().get_string())
        except ValueError:
            i = -1
        if not self._reverse:
            if i >= 0:  # feature ID matches one of the strings in the parameter "fids"
                if self._filtrate_feature_collection[i] is None:
                    self._filtrate_feature_collection[i] = feature
                else:
                    logger.warning(
                        f"Duplicate feature ID found: {feature.get_feature_id().get_string()}. Only the first one will be kept in the filtrate feature collection."
                    )
                return True
            else:
                self._residue_feature_collection.append(feature)  # type: ignore
                return False
        else:
            if i < 0:
                self._filtrate_feature_collection.append(feature)
                return True
            else:
                if self._residue_feature_collection[i] is None:
                    self._residue_feature_collection[i] = feature  # type: ignore
                else:
                    logger.warning(
                        f"Duplicate feature ID found: {feature.get_feature_id().get_string()}. Only the first one will be kept in the residue feature collection."
                    )
                return False


class ValidTimeFilter(FeatureFilter):
    """filter features by their valid time, keep features with valid time within a specified time range."""

    def __init__(
        self, begin_time: float = float("inf"), end_time: float = float("-inf")
    ):
        """Constructor for the valid time filter.

        Parameters
        ----------
        begin_time : float, optional
            The begin time of the time range (inclusive). Default is ``inf``.
        end_time : float, optional
            The end time of the time range (inclusive). Default is ``-inf``.
        """
        super().__init__()
        self._begin_time = begin_time
        self._end_time = end_time

    def should_keep(self, feature: pygplates.Feature) -> bool:  # type: ignore
        valid_time = feature.get_valid_time(None)
        if valid_time:
            begin_time, end_time = valid_time

            if begin_time <= self._begin_time or end_time >= self._end_time:
                self._filtrate_feature_collection.append(feature)
                return True

        self._residue_feature_collection.append(feature)
        return False
