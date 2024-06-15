#!/usr/bin/env python3

import datetime
import os
import time

import pygplates as _pygplates
from common import *

import gplately

print(gplately.__file__)


from gplately.seafloor_grids import make_seafloor_grids


def main():
    model_name = "merdith2021"

    start = time.time()

    make_seafloor_grids(model_name, initial_time=410, youngest_time=400)

    end = time.time()
    hours_minutes_seconds = str(datetime.timedelta(seconds=end - start)).split(":")
    print(
        f"Completed creating age grids in {hours_minutes_seconds[0]} Hours, {hours_minutes_seconds[1]} Minutes, {hours_minutes_seconds[2].split('.')[0]} Seconds "
    )


if __name__ == "__main__":
    main()
