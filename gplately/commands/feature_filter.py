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

import abc
import argparse
import logging
import re
from typing import List

import pygplates

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
    def should_keep(self, feature: pygplates.Feature) -> bool:
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

    def should_keep(self, feature: pygplates.Feature) -> bool:
        if self.exclude:
            for name in self.names:
                if self.check_name(name, feature.get_name()):
                    return False
            return True
        else:
            for name in self.names:
                if self.check_name(name, feature.get_name()):
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

    def should_keep(self, feature: pygplates.Feature) -> bool:
        if not self.exclude and feature.get_reconstruction_plate_id() in self.pids:
            return True
        if self.exclude and feature.get_reconstruction_plate_id() not in self.pids:
            return True
        return False


class BirthAgeFilter(FeatureFilter):
    """filter features by the time of appearance

    for example:
        BirthAgeFilter(500) -- keep features whose time of apprearance are bigger than 500
        BirthAgeFilter(500, keep_older=False) --  keep features whose time of apprearance are smaller than 500

    :param age: the age criterion
    :param keep_older: if True, return True when the feature's birth age is older than the age criterion. If False, otherwise.

    """

    def __init__(self, age: float, keep_older=True):
        self.age = age
        self.keep_older = keep_older

    def should_keep(self, feature: pygplates.Feature) -> bool:
        valid_time = feature.get_valid_time(None)
        if valid_time:
            if self.keep_older and valid_time[0] > self.age:
                return True
            if not self.keep_older and valid_time[0] < self.age:
                return True
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

    def should_keep(self, feature: pygplates.Feature) -> bool:
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
    feature_collection: pygplates.FeatureCollection, filters: List[FeatureFilter]
):
    """Filter feature collection by various criteria."""
    new_feature_collection = pygplates.FeatureCollection()
    for feature in feature_collection:
        keep_flag = True
        for filter in filters:
            if not filter.should_keep(feature):
                keep_flag = False
                break
        if keep_flag:
            new_feature_collection.add(feature)
    return new_feature_collection


help_str = "Filter feature collection by various criteria."

__description__ = f"""{help_str}

Examples: 
    - `gplately filter input_file output_file -n Africa "North America"`
        (get features whose name contains "Africa" or "North America")

    - `gplately filter input_file output_file -p 701 714 715 101`
        (get features whose plate ID is one of 701 714 715 101)
    
    - `gplately filter input_file output_file --min-birth-age 500`
        (get features whose birth age is older than 500Myr)
    
    - `gplately filter input_file output_file --max-birth-age 500`
        (get features whose birth age is younger than 500Myr)
    
    - `gplately filter input_file output_file -n Africa "North America" -p 701 714 715 101 --min-birth-age 500`
        (get features whose name conains "Africa" or "North America" and plate ID is one of 701 714 715 101 and birth age is older than 500Myr)
    
    - `gplately filter input_file output_file -t gpml:Basin`
        (get all gpml:Basin features)
    
    - `gplately filter input_file output_file -t "gpml:IslandArc|gpml:Basin"`
        (get all gpml:Basin and gpml:IslandArc features)

    See https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_feature_filter.sh for more examples. 
"""


def add_parser(subparser):
    """add feature filter command line argument parser"""
    filter_cmd = subparser.add_parser(
        "filter",
        help=help_str,
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # feature filter command arguments
    filter_cmd.set_defaults(func=run_filter_feature_collection)
    filter_cmd.add_argument("filter_input_file", type=str, help="the input file")
    filter_cmd.add_argument(
        "filter_output_file",
        type=str,
        help="the output file into which the filtered feature collection will be saved",
    )

    # filter by feature name
    name_group = filter_cmd.add_mutually_exclusive_group()
    name_group.add_argument(
        "-n",
        "--names",
        type=str,
        dest="names",
        nargs="+",
        metavar="feature_name",
        help="features whose name contains the `feature_name` will be kept",
    )
    name_group.add_argument(
        "--exclude-names",
        type=str,
        dest="exclude_names",
        nargs="+",
        metavar="feature_name",
        help="features whose name does not contain the `feature_name` will be kept",
    )

    # filter by plate ID
    pid_group = filter_cmd.add_mutually_exclusive_group()
    pid_group.add_argument(
        "-p",
        "--pids",
        type=int,
        dest="pids",
        nargs="+",
        metavar="pid",
        help="features whose PID are in the `pids` will be kept",
    )
    pid_group.add_argument(
        "--exclude-pids",
        type=int,
        dest="exclude_pids",
        nargs="+",
        metavar="pid",
        help="features whose PID are not in the `exclude_pids` will be kept",
    )

    # filter by birth age
    birth_age_group = filter_cmd.add_mutually_exclusive_group()
    birth_age_group.add_argument(
        "-a",
        "--min-birth-age",
        type=float,
        dest="min_birth_age",
        metavar="min_birth_age",
        help="the features whose birth age is older than the `min_birth_age` will be kept",
    )
    birth_age_group.add_argument(
        "--max-birth-age",
        type=float,
        dest="max_birth_age",
        metavar="max_birth_age",
        help="the features whose birth age is younger than the `min_birth_age` will be kept",
    )

    filter_cmd.add_argument(
        "--case-sensitive",
        dest="case_sensitive",
        action="store_true",
        help="flag to indicate if the `filter by feature name` is case sensitive(not apply to feature type)",
    )
    filter_cmd.add_argument(
        "--exact-match",
        dest="exact_match",
        action="store_true",
        help="flag to indicate if the `filter by feature name` need to be exact match(not apply to feature type)",
    )

    # filter by feature type
    filter_cmd.add_argument(
        "-t",
        "--feature-type",
        type=str,
        dest="feature_type_re",
        metavar="feature_type",
        help="the feature type regular expression; features whose type matches the regular expression will be kept",
    )


def run_filter_feature_collection(args):
    """Filter the input feature collection according to command line arguments."""
    input_feature_collection = pygplates.FeatureCollection(args.filter_input_file)

    filters = []
    if args.names:
        filters.append(
            FeatureNameFilter(
                args.names,
                exact_match=args.exact_match,
                case_sensitive=args.case_sensitive,
            )
        )
    elif args.exclude_names:
        filters.append(
            FeatureNameFilter(
                args.exclude_names,
                exclude=True,
                exact_match=args.exact_match,
                case_sensitive=args.case_sensitive,
            )
        )

    if args.pids:
        filters.append(PlateIDFilter(args.pids))
    elif args.exclude_pids:
        filters.append(PlateIDFilter(args.exclude_pids, exclude=True))

    # print(args.max_birth_age)
    if args.max_birth_age is not None:
        filters.append(BirthAgeFilter(args.max_birth_age, keep_older=False))
    elif args.min_birth_age is not None:
        filters.append(BirthAgeFilter(args.min_birth_age))

    if args.feature_type_re is not None:
        filters.append(FeatureTypeFilter(args.feature_type_re))

    new_fc = filter_feature_collection(
        input_feature_collection,
        filters,
    )

    if len(new_fc) == 0:
        logger.warn(
            "No feature is left after the filtering. The output feature collection will be empty."
        )

    new_fc.write(args.filter_output_file)
    logger.info(
        f"Done! The filtered feature collection has been saved to {args.filter_output_file}."
    )
