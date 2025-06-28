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
import re

import pygplates

logger = logging.getLogger("gplately")

help_str = "Reset the feature type for the selected features. "

__description__ = f"""{help_str}

Example usage: 
    - `gplately reset_feature_type -s gpml:ClosedContinentalBoundary -t gpml:UnclassifiedFeature input_file output_file`
        (change all gpml:ClosedContinentalBoundary to gpml:UnclassifiedFeature)
        
    - `gplately reset_feature_type -s "gpml:ContinentalFragment|gpml:Coastline" -t gpml:UnclassifiedFeature input_file output_file`
        (change all gpml:ContinentalFragment and gpml:Coastline to gpml:UnclassifiedFeature)
        
    - `gplately reset_feature_type -s ".*" -t gpml:UnclassifiedFeature input_file output_file` 
        (change all feature types to gpml:UnclassifiedFeature)     

    See https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_reset_feature_type.sh for more examples. 
"""


def add_parser(subparser):
    """add `reset_feature_type` command line argument parser"""
    reset_feature_type_cmd = subparser.add_parser(
        "reset_feature_type",
        help=help_str,
        add_help=True,
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    reset_feature_type_cmd.set_defaults(func=reset_feature_type)
    reset_feature_type_cmd.add_argument("input_file", type=str, help="the input file")
    reset_feature_type_cmd.add_argument(
        "output_file",
        type=str,
        help="the output file into which the new features will be saved",
    )

    reset_feature_type_cmd.add_argument(
        "-s",
        "--select-feature-type",
        type=str,
        dest="select_feature_type_re",
        metavar="feature_type_re",
        help="the regular expression to select features by featuer type",
    )

    reset_feature_type_cmd.add_argument(
        "-t",
        "--set-feature-type",
        type=str,
        dest="set_feature_type",
        metavar="feature_type",
        help="the feature type to be set to, such as gpml:UnclassifiedFeature, gpml:Coastline, etc. See https://www.gplates.org/docs/gpgim/#FeatureClassList",
    )

    reset_feature_type_cmd.add_argument(
        "--keep-feature-id",
        dest="keep_feature_id",
        help="flag to indicate if we use the same feature IDs after resetting the feature type",
        action="store_true",
    )

    reset_feature_type_cmd.add_argument(
        "--verify-information-model",
        dest="verify_information_model",
        help="flag to indicate if we need to keep GPGIM integrity",
        action="store_true",
    )


def reset_feature_type(args):
    try:
        pygplates.Feature(
            pygplates.FeatureType.create_from_qualified_string(args.set_feature_type),
            verify_information_model=pygplates.VerifyInformationModel.yes,  # in case the user provided a wrong feature type
        )
    except pygplates.InformationModelError:
        logger.error(f"Invalid feature type: {args.set_feature_type}")
        logger.error(
            f"You have asked to reset features to an invalid feature type. It cannot be done."
        )
        return

    if args.verify_information_model:
        verify_information_model = pygplates.VerifyInformationModel.yes
    else:
        verify_information_model = pygplates.VerifyInformationModel.no

    new_fc = pygplates.FeatureCollection()
    for f in pygplates.FeatureCollection(args.input_file):
        feature_type = f.get_feature_type()
        if re.fullmatch(
            args.select_feature_type_re, feature_type.to_qualified_string()
        ):
            if args.keep_feature_id:
                new_f = pygplates.Feature(
                    pygplates.FeatureType.create_from_qualified_string(
                        args.set_feature_type
                    ),
                    f.get_feature_id(),
                )
            else:
                new_f = pygplates.Feature(
                    pygplates.FeatureType.create_from_qualified_string(
                        args.set_feature_type
                    )
                )
            for p in f:
                try:
                    new_f.add(
                        p.get_name(),
                        p.get_time_dependent_value(),
                        verify_information_model,
                    )
                except pygplates.InformationModelError:
                    logger.error(
                        f"The feature type cannot be changed to {args.set_feature_type}. "
                        + "Some properties in the features are not allowed in the new feature type."
                    )
                    return
            new_fc.add(new_f)
        else:
            new_fc.add(f)

    new_fc.write(args.output_file)
    logger.info(
        f"Reset feature type successfully. The new features have been saved in {args.output_file}."
    )
