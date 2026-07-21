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

"""
This sub-module retrieves paleomagnetic data from http://www.gpmdb.net,
then create and save GPlates VGP features in a .gpmlz file.
"""

import argparse
import json
import logging
import math
import os
import sys
from pathlib import Path
from typing import cast

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false

import numpy as np
import pandas as pd
import pygplates
import requests
from plate_model_manager import PlateModelManager

# you need pygplates and pandas to run this script
# using gplately env is recommended
# `micromamba activate gplately`

DEFAULT_GPMDB_SERVER_URL = "https://www.gpmdb.net/api/search/"
DATA_CACHE_DIR = "data-cache"
QUERY_DATA_FILENAME = "gpmdb.json"


logger = logging.getLogger("gplately")


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
    vgp_feature.set_name(str(RESULTNO))  # type: ignore
    vgp_feature.set_description(PAPER_URL)  # type: ignore
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
    vgp_feature.set_valid_time(HIGHAGE, LOWAGE)  # type: ignore

    vgp_feature.set_reconstruction_plate_id(plate_id)  # type: ignore

    return vgp_feature


def assign_plate_id(df, static_polygon_file, rotation_file):
    """assign plate ids for sites"""
    pids = []
    sites = []
    plate_partitioner = pygplates.PlatePartitioner(
        static_polygon_file,
        rotation_file,
    )
    for _, row in df.iterrows():
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
        "-m",
        "--model",
        type=str,
        dest="model_name",
        default="Zahirovic2022",
        help="the plate reconstruction model name to assign plate IDs for the VGP features. Default is 'Zahirovic2022'. The available model names can be found with command `pmm ls` (plate-model-manager https://pypi.org/project/plate-model-manager/).",
    )
    parser.add_argument("-o", "--outfile", type=str, dest="outfile")
    parser.add_argument(
        "--use-cached-data",
        action="store_true",
        help="use cached data to create VGP features. By default, the script will retrieve fresh data from the server. If you want to use cached data, add this command line argument. The cached data is stored in the 'data-cache' directory and the file name is 'gpmdb.json'.",
    )
    parser.add_argument(
        "--gpmdb-server-url",
        type=str,
        dest="gpmdb_server_url",
        default=DEFAULT_GPMDB_SERVER_URL,
        help="the GPMDB server URL to retrieve data from. Default is https://www.gpmdb.net/api/search/.",
    )


__description__ = """Retrieve paleomagnetic data from https://www.gpmdb.net, create GPlates-compatible VGP features, and save them to a .gpmlz file.

    Data source URL:
        - https://www.gpmdb.net/api/search/

    This command creates two files:
        - vgp_features_<model_name>.gpmlz (or the file set with -o/--outfile), containing GPlates-compatible VGP features
        - data-cache/gpmdb.json, containing the raw paleomagnetic data

    Usage example:
        gplately gpmdb -m Zahirovic2022 -o test.gpmlz

    Notes:
        - The default reconstruction model used to assign plate IDs is "Zahirovic2022".
        - You can specify a different model with -m/--model.
        - Available model names can be listed with "pmm ls"
          (plate-model-manager: https://pypi.org/project/plate-model-manager/).
        - You can specify the output .gpmlz filename with -o/--outfile.
        - If -o/--outfile is not provided, the output filename is
          "vgp_features_<model_name>.gpmlz".
        - The cached raw data filename is always "gpmdb.json".

    If https://www.gpmdb.net is not accessible, you can use backup data, for example:
        gplately gpmdb --gpmdb-server-url https://repo.gplates.org/webdav/gplately/gpmdb.json
    """


def main(args):
    # get query data
    if not args.use_cached_data or not os.path.isfile(
        f"{DATA_CACHE_DIR}/{QUERY_DATA_FILENAME}"
    ):
        try:
            response = requests.get(f"{args.gpmdb_server_url}", timeout=30)
            query_data = response.json()
        except (
            requests.exceptions.JSONDecodeError,
            requests.exceptions.ConnectionError,
        ):
            logger.error(
                f"FATAL: The {args.gpmdb_server_url} did not return valid data. Check and make sure the website is up and running!"
            )
            sys.exit(1)
        Path(DATA_CACHE_DIR).mkdir(parents=True, exist_ok=True)
        with open(f"{DATA_CACHE_DIR}/{QUERY_DATA_FILENAME}", "w+") as outfile:
            logger.info(
                f"Saving retrieved data to {DATA_CACHE_DIR}/{QUERY_DATA_FILENAME}."
            )
            outfile.write(json.dumps(query_data))
    else:
        logger.warning(
            "Using cached data. To retrieve fresh data from the server, remove the --use-cached-data command line argument."
        )
        with open(f"{DATA_CACHE_DIR}/{QUERY_DATA_FILENAME}", "r") as infile:
            query_data = json.load(infile)

    columns = query_data["columns"]
    data_array = []
    for row in query_data["data"]:
        if len(row) != len(columns):
            logger.warning(
                f"Invalid data row: {row}. The number of values in the row does not match the number of columns."
            )
            continue
        row_array = [row[key] for key in row]
        data_array.append(row_array)

    df_query = pd.DataFrame(np.array(data_array), columns=columns)
    df_query = df_query.sort_values(by=["RESULTNO"], ignore_index=True)

    # Get a new DataFrame with only the critical columns for creating VGP
    # features. Rows with missing values in any critical column are dropped.
    critical_columns = [
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

    df = cast(pd.DataFrame, df_query[critical_columns]).copy()

    # Ensure numeric columns are numeric so invalid/null-like values become NaN
    # and are dropped.
    for col in critical_columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=critical_columns).reset_index(drop=True)

    logger.warning(
        f"{len(df_query) - len(df)} rows have been dropped due to missing critical values. Only {len(df)} rows are used to create VGP features."
    )

    pm_manager = PlateModelManager()
    model = pm_manager.get_model(args.model_name, data_dir=".")

    sites, pids = assign_plate_id(
        df,
        static_polygon_file=model.get_static_polygons(),
        rotation_file=model.get_rotation_model(),
    )

    count = 0
    vgp_features = []
    for i, (_, row) in enumerate(df.iterrows()):
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
                sample_site_position=sites[i],
                plate_id=pids[i],
            )

            # Add VGP feature to list.
            vgp_features.append(vgp_feature)
        else:
            count += 1
            # print(f"ignore row: {row}")

    # Save VGP features to file.
    if not args.outfile:
        outfile_name = f"vgp_features_{args.model_name}.gpmlz"
    else:
        outfile_name = args.outfile

    pygplates.FeatureCollection(vgp_features).write(outfile_name)

    print(f"The file {outfile_name} has been created successfully.")


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
