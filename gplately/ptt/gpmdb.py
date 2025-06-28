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

"""
This sub-module retrieves paleomagnetic data from http://www.gpmdb.net,
then create and save GPlates VGP features in a .gpmlz file.
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pygplates
import requests
from plate_model_manager import PlateModelManager

# you need pygplates and pandas to run this script
# using gplately env is recommended
# `micromamba activate gplately`

DEFAULT_GPMDB_SERVER_URL = "http://www.gpmdb.net"
DATA_CACHE_DIR = "data-cache"
QUERY_DATA_URL = "get_query_data/"
QUERY_DATA_FILENAME = "query-data.json"
PMAG_RESULT_URL = "get_PMAGRESULT_data/?fmt=json"
PMAG_RESULT_FILENAME = "pmag-result.json"


def create_vgp_feature(
    RESULTNO,
    PLAT,
    PLONG,
    INC,
    DEC,
    DP,
    DM,
    LOMAGAGE,
    HIMAGAGE,
    LOWAGE,
    HIGHAGE,
    sample_site_position,
    plate_id,
):
    """function to create VGP feature"""
    # Paper reference or URL/DOI.
    PAPER_URL = "https://doi.org/10.1016/j.earscirev.2022.104258"

    # Create feature
    vgp_feature = pygplates.Feature(pygplates.FeatureType.gpml_virtual_geomagnetic_pole)
    vgp_feature.set_name(str(RESULTNO))
    vgp_feature.set_description(PAPER_URL)
    vgp_feature.set(
        pygplates.PropertyName.gpml_average_sample_site_position,
        pygplates.GmlPoint(sample_site_position),
    )
    vgp_feature.set(
        pygplates.PropertyName.gpml_pole_position,
        pygplates.GmlPoint(pygplates.PointOnSphere(PLAT, PLONG)),
    )
    vgp_feature.set(
        pygplates.PropertyName.gpml_pole_a95, pygplates.XsDouble(math.sqrt(DP + DM))
    )
    vgp_feature.set(pygplates.PropertyName.gpml_pole_dp, pygplates.XsDouble(DP))
    vgp_feature.set(pygplates.PropertyName.gpml_pole_dm, pygplates.XsDouble(DM))
    vgp_feature.set(
        pygplates.PropertyName.gpml_average_inclination, pygplates.XsDouble(INC)
    )
    vgp_feature.set(
        pygplates.PropertyName.gpml_average_declination, pygplates.XsDouble(DEC)
    )
    vgp_feature.set(
        pygplates.PropertyName.gpml_average_age,
        pygplates.XsDouble((LOMAGAGE + HIMAGAGE) / 2.0),
    )
    vgp_feature.set_valid_time(HIGHAGE, LOWAGE)

    vgp_feature.set_reconstruction_plate_id(plate_id)

    return vgp_feature


def assign_plate_id(df, static_polygon_file, rotation_file):
    """assign plate ids for sites"""
    pids = []
    sites = []
    plate_partitioner = pygplates.PlatePartitioner(
        static_polygon_file,
        rotation_file,
    )
    for index, row in df.iterrows():
        # print(index)
        sample_site_position = pygplates.PointOnSphere(row["SLAT"], row["SLONG"])
        reconstructed_static_polygon = plate_partitioner.partition_point(
            sample_site_position
        )
        if reconstructed_static_polygon:
            plate_id = (
                reconstructed_static_polygon.get_feature().get_reconstruction_plate_id()
            )
        else:
            plate_id = 0
        sites.append(sample_site_position)
        pids.append(plate_id)
    return sites, pids


class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"error: {message}\n")
        self.print_help()
        sys.exit(1)


def add_arguments(parser: argparse.ArgumentParser):
    """add command line argument parser"""
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.description = __description__

    parser.set_defaults(func=main)

    parser.add_argument(
        "-m", "--model", type=str, dest="model_name", default="Muller2022"
    )
    parser.add_argument("-o", "--outfile", type=str, dest="outfile")
    parser.add_argument(
        "--use-cached-data",
        action="store_true",
        help="use cached data for debugging purpose",
    )
    parser.add_argument(
        "--gpmdb-server-url",
        type=str,
        dest="gpmdb_server_url",
        default=DEFAULT_GPMDB_SERVER_URL,
        help="the GPMDB server URL ",
    )


__description__ = """Retrieve paleomagnetic data from https://www.gpmdb.net, create GPlates-compatible VGP features and save the VGP features in a .gpmlz file.

    The two URLs being used are 
        - https://www.gpmdb.net/get_query_data/
        - https://www.gpmdb.net/get_PMAGRESULT_data/?fmt=json

    This command will create two files.
        - the .gpmlz file -- contains the GPlates-compatible VGP features
        - the pmag-result.csv -- contains the raw paleomagnetic data

    Usage example: gplately gpmdb -m zahirovic2022 -o test.gpmlz

    The default reconstruction model being used to assign plate IDs is "Muller2022". User can choose to specify the model with "-m/--model". 
    The avaliable model names can be found with command `pmm ls` (plate-model-manager https://pypi.org/project/plate-model-manager/; use gplately conda env).

    User can specify the output .gpmlz file name with "-o/--outfile". By default, the output file name will be "vgp_features_{model name}.gmplz".

    The `gplately gpmdb` will use model "Muller2022" and create the output file "vgp_features_Muller2022.gmplz". 

    The file name for the raw paleomagnetic data is always "pmag-result.csv".

    """


def main(args):
    # get query data
    if not args.use_cached_data or not os.path.isfile(
        f"{DATA_CACHE_DIR}/{QUERY_DATA_FILENAME}"
    ):
        try:
            response = requests.get(
                f"{args.gpmdb_server_url}/{QUERY_DATA_URL}", verify=False
            )
            query_data = response.json()
        except (
            requests.exceptions.JSONDecodeError,
            requests.exceptions.ConnectionError,
        ):
            print(
                f"FATAL: The {args.gpmdb_server_url}/{QUERY_DATA_URL} did not return valid data. Check and make sure the website is up and running!"
            )
            sys.exit(1)
        Path(DATA_CACHE_DIR).mkdir(parents=True, exist_ok=True)
        with open(f"{DATA_CACHE_DIR}/{QUERY_DATA_FILENAME}", "w+") as outfile:
            outfile.write(json.dumps(query_data))
    else:
        with open(f"{DATA_CACHE_DIR}/{QUERY_DATA_FILENAME}", "r") as infile:
            query_data = json.load(infile)

    columns = []
    for c in query_data["columns"]:
        if isinstance(c, list):
            columns += c
        elif isinstance(c, str):
            columns.append(c)
        else:
            raise Exception(f"Invalid comlumn type: {type(c)}")

    df_query = pd.DataFrame(np.array(query_data["data"]), columns=columns)
    df_query = df_query.sort_values(by=["RESULTNO"], ignore_index=True)

    # get pmag-result data
    if not args.use_cached_data or not os.path.isfile(
        f"{DATA_CACHE_DIR}/{PMAG_RESULT_FILENAME}"
    ):
        try:
            response = requests.get(
                f"{args.gpmdb_server_url}/{PMAG_RESULT_URL}", verify=False
            )
            pmagresult_data = response.json()
        except (
            requests.exceptions.JSONDecodeError,
            requests.exceptions.ConnectionError,
        ):
            print(
                f"FATAL: The {args.gpmdb_server_url}/{PMAG_RESULT_URL} did not return valid data. Check and make sure the website is up and running!"
            )
            sys.exit(1)
        with open(f"{DATA_CACHE_DIR}/{PMAG_RESULT_FILENAME}", "w+") as outfile:
            outfile.write(json.dumps(pmagresult_data))
    else:
        with open(f"{DATA_CACHE_DIR}/{PMAG_RESULT_FILENAME}", "r") as infile:
            pmagresult_data = json.load(infile)

    columns = []
    for c in pmagresult_data["columns"]:
        if isinstance(c, list):
            columns += c
        elif isinstance(c, str):
            columns.append(c)
        else:
            raise Exception(f"Invalid comlumn type: {type(c)}")

    df_pmagresult = pd.DataFrame(np.array(pmagresult_data["data"]), columns=columns)
    df_pmagresult = df_pmagresult.sort_values(by=["RESULTNO"], ignore_index=True)

    df_pmagresult["LOWAGE"] = df_query["LOWAGE"]
    df_pmagresult["HIGHAGE"] = df_query["HIGHAGE"]

    df = df_pmagresult[
        [
            "RESULTNO",
            "SLAT",
            "SLONG",
            "PLAT",
            "PLONG",
            "INC",
            "DEC",
            "DP",
            "DM",
            "LOMAGAGE",
            "HIMAGAGE",
            "LOWAGE",
            "HIGHAGE",
        ]
    ]

    pm_manager = PlateModelManager()
    model = pm_manager.get_model(args.model_name, data_dir=".")

    sites, pids = assign_plate_id(
        df,
        static_polygon_file=model.get_static_polygons(),
        rotation_file=model.get_rotation_model(),
    )

    count = 0
    vgp_features = []
    for index, row in df.iterrows():
        # print(index)
        if row["DP"] is not None and row["DM"] is not None and row["INC"] is not None:
            vgp_feature = create_vgp_feature(
                RESULTNO=row["RESULTNO"],
                PLAT=row["PLAT"],
                PLONG=row["PLONG"],
                INC=row["INC"],
                DEC=row["DEC"],
                DP=row["DP"],
                DM=row["DM"],
                LOMAGAGE=row["LOMAGAGE"],
                HIMAGAGE=row["HIMAGAGE"],
                LOWAGE=row["LOWAGE"],
                HIGHAGE=row["HIGHAGE"],
                sample_site_position=sites[index],
                plate_id=pids[index],
            )

            # Add VGP feature to list.
            vgp_features.append(vgp_feature)
        else:
            count += 1
            # print(f"ignore row: {row}")

    print(f"{count} rows have been ignored due to None values.")

    # Save VGP features to file.
    if not args.outfile:
        outfile_name = f"vgp_features_{args.model_name}.gpmlz"
    else:
        outfile_name = args.outfile

    pygplates.FeatureCollection(vgp_features).write(outfile_name)

    df_pmagresult = df_pmagresult.drop(["ROCKUNITNO"], axis=1)
    df_pmagresult.to_csv("pmag-result.csv", index=False)

    print(
        f"The files {outfile_name} and pmag-result.csv have been created successfully."
    )


if __name__ == "__main__":
    # The command-line parser.
    parser = argparse.ArgumentParser(
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # add arguments
    add_arguments(parser)

    # Parse command-line options.
    args = parser.parse_args()

    # call main function
    main(args)
