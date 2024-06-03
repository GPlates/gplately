#!/usr/bin/env python3

import os

import pygplates as _pygplates
from common import *

import gplately

print(gplately.__file__)


from gplately.seafloor_grids import make_seafloor_grids


def main():
    model_name = "merdith2021"
    times = list(range(410, 400, -1))
    make_seafloor_grids(model_name, times)


if __name__ == "__main__":
    main()
