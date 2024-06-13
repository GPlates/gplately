#!/usr/bin/env python3

import os

import pygplates as _pygplates
from common import *

import gplately

print(gplately.__file__)


def main():
    fc = gplately.pygplates.FeatureCollection()
    feature = gplately.pygplates.Feature()
    feature.set_name("feature-1")
    feature.set_geometry(_pygplates.PointOnSphere(0, 0))
    fc.add(feature)
    fc.write("fc-1.gpml")

    fc = gplately.pygplates.FeatureCollection()
    feature = gplately.pygplates.Feature()
    feature.set_name("feature-2")
    feature.set_geometry(_pygplates.PointOnSphere(10, 10))
    fc.add(feature)
    fc.write("fc-2.gpml")

    fc = gplately.pygplates.FeatureCollection.from_file_list(
        filenames=["fc-1.gpml", "fc-2.gpml"]
    )

    fc = gplately.pygplates.FeatureCollection()
    feature = gplately.pygplates.Feature()
    feature.set_name("feature-3")
    feature.set_geometry(_pygplates.PointOnSphere(-10, -10))
    fc.add(feature)
    fc.write("fc-3.gpml")

    fc = gplately.pygplates.FeatureCollection.from_file_list(
        filenames=["fc-1.gpml", "fc-2.gpml"]
    )

    fc.add(filename="fc-3.gpml")
    fc.write("fc-4.gpml")

    print(fc.filenames)

    fc_t = fc.clone()
    fc_t.write("fc-5.gpml")

    os.remove("fc-1.gpml")
    os.remove("fc-2.gpml")
    os.remove("fc-3.gpml")
    os.remove("fc-4.gpml")
    os.remove("fc-5.gpml")


if __name__ == "__main__":
    main()
