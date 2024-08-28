#
#    Copyright (C) 2024 The University of Sydney, Australia
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

import logging
import os

logger = logging.getLogger("gplately")

disable_dev_warning = (
    "DISABLE_GPLATELY_DEV_WARNING" in os.environ
    and os.environ["DISABLE_GPLATELY_DEV_WARNING"].lower() == "true"
)


def print_dev_warning(version: str):
    if not disable_dev_warning:
        print()
        print(
            "##################################################################################################"
        )
        print(
            f"""
            WARNING: 👈👈 👀👀 
            You are using a DEV version 👉({version})👈 GPlately.     
            Some functionalities in the DEV version have not been tested thoroughly, 
            and may break your code or produce wrong results due to 
            its unstable nature(DEV in progress). Proceed With Caution!!!
            You might also need to install the DEV version plate_model_manager 
            from https://github.com/michaelchin/plate-model-manager.

            🚫 To disable this warning, set USING_DEV_VERSION to False in __init__.py or
            set DISABLE_GPLATELY_DEV_WARNING environment variable to true.
            👉`export DISABLE_GPLATELY_DEV_WARNING=true`👈
            """
        )
        print(
            "##################################################################################################"
        )
        print()


def print_using_source_code_warning(version: str):
    if not disable_dev_warning:
        if os.path.isdir(
            f"{os.path.dirname(os.path.realpath(__file__))}/../.git"
        ) or not os.path.isfile(
            f"{os.path.dirname(os.path.realpath(__file__))}/../../../../bin/gplately"
        ):
            logger.warn(
                f"The location of GPlately currently in use is {os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}. "
                + f"It seems that you are using GPlately source code directly or installed editable package with `pip install -e .`, "
                + f"the version number(❗{version}❗) 👉may not be accurate👈 in these cases."
            )
