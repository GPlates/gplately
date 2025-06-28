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

import requests
from plate_model_manager import PlateModelManager

from ..exceptions import UnableToGetModelList

logger = logging.getLogger("gplately")

help_str = "Show a list of available reconstruction models."

__description__ = f"""{help_str}

Example usage: 
    - gplately list
    - gplately list -m Merdith2021
"""


def add_parser(subparser):
    """add 'list model' command line argument parser"""
    list_cmd = subparser.add_parser(
        "list",
        help=help_str,
        add_help=True,
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # feature filter command arguments
    list_cmd.set_defaults(func=run_list_models)
    list_cmd.add_argument("-m", "--model", type=str, dest="model", nargs=1)


def run_list_models(args):
    if args.model:
        print()
        print("Layers:")
        for l in get_layer_names(args.model[0]):
            print(f"    {l}")
        print()
        print(f"Model name: {args.model[0]}")
        print(f"Model URL: {get_model_url(args.model[0])}")
        print()
    else:
        print()
        print("Models:")
        for n in get_model_names():
            print(f"    {n}")
        print()


def get_model_names():
    """return a list of model names from the servers.
    the function will try a list of urls one by one. if the first url is not working, try the second one.
    if the second one is still not working, try the thrid one. Repeat until the end of the list.
    """
    model_list_urls = [
        "https://repo.gplates.org/webdav/pmm/config/gplately_model_list.json",
        "https://www.earthbyte.org/webdav/pmm/config/gplately_model_list.json",
        "https://portal.gplates.org/static/pmm/config/gplately_model_list.json",
    ]
    names = []
    for url in model_list_urls:
        try:
            response = requests.get(url, timeout=(5, 5))
            if response.status_code == 200:
                models = response.json()
                mnames = PlateModelManager().get_available_model_names()
                for model in models:
                    if model in mnames:
                        names.append(model)
                return names
            else:
                continue
        except:
            continue
    logger.error(
        "Unable to get a list of model names from the servers. Check the network connection. See the server list below."
    )
    logger.info(model_list_urls)
    raise UnableToGetModelList


def get_layer_names(model: str):
    """Given model name, return the layer names in the model."""
    m = PlateModelManager().get_model(model)
    if m:
        return m.get_avail_layers()
    else:
        return []


def get_model_url(model: str):
    """Given model name, return the URL to the model files."""
    m = PlateModelManager().get_model(model)
    if m:
        cfg = m.get_cfg()
    else:
        cfg = {}
    if "URL" in cfg:
        return cfg["URL"]
    else:
        return ""
