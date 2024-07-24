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

import abc
import argparse
from typing import List

import pygplates

__description__ = """Show a list of available reconstruction models.

Example usage: 
    - gplately list
    - gplately list -m Merdith2021
"""


def add_parser(subparser):
    """add 'list model' command line argument parser"""
    list_cmd = subparser.add_parser(
        "list",
        help=__description__,
        add_help=True,
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # feature filter command arguments
    list_cmd.set_defaults(func=run_list_models)
    list_cmd.add_argument("-m", "--model", type=str, dest="model", nargs=1)


def run_list_models(args):
    print("run_list_models")
    if args.model:
        print(args.model)
