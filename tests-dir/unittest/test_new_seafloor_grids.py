#!/usr/bin/env python3

import os

import pygplates as _pygplates
from common import *

import gplately

print(gplately.__file__)


from gplately.seafloor_grids import make_seafloor_grids


def main():
    model_name = "merdith2021"
    make_seafloor_grids(model_name, initial_time=410, youngest_time=400)


if __name__ == "__main__":
    main()
