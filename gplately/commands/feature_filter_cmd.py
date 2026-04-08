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


import argparse
import logging
import pygplates
from..utils.feature_filter import (
    FeatureNameFilter,
    PlateIDFilter,
    BirthAgeFilter,
    FeatureTypeFilter,
    filter_feature_collection,
)

logger = logging.getLogger("gplately")

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
    input_feature_collection = pygplates.FeatureCollection(args.filter_input_file) # type: ignore

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
        logger.warning(
            "No feature is left after the filtering. The output feature collection will be empty."
        )

    new_fc.write(args.filter_output_file)
    logger.info(
        f"Done! The filtered feature collection has been saved to {args.filter_output_file}."
    )
