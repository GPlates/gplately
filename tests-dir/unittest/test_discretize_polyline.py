#!/usr/bin/env python3

import math
import os
import sys

import pygplates
from common import MODEL_REPO_DIR, save_fig

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import gplately
from gplately.lib import discretize_polyline

print(gplately.__file__)


def degree_to_radian(point):
    return [math.radians(point[0]), math.radians(point[1])]


def radian_to_degree(point):
    return [math.degrees(point[0]), math.degrees(point[1])]


def main(input_file, distance, output_file, strict_mode=True):
    """discretize the polylines according to the given distance

    :param input_file: the input polyline file
    :param distance: the distance(in radians) between two adjacent points
    :param output_file: save the points into a file
    :returns: no return. the result will be saved to output_file
    """
    cob_fc = pygplates.FeatureCollection(input_file)
    new_fc = pygplates.FeatureCollection()
    for f in cob_fc:
        f_name = f.get_name()
        # print(f_name)
        for geom in f.get_all_geometries():
            assert isinstance(geom, pygplates.PolylineOnSphere)
            # convert lat, lon from degrees to radians
            line_points = [degree_to_radian(p) for p in geom.to_lat_lon_list()]  # type: ignore

            # print(len(line_points))
            assert len(line_points) > 1

            new_line_points = discretize_polyline(
                line_points, distance, strict=strict_mode
            )

            # print(len(new_line_points))

            if len(new_line_points) > 1:
                nf = pygplates.Feature()
                nf.set_geometry(
                    pygplates.MultiPointOnSphere(
                        [radian_to_degree(p) for p in new_line_points]
                    )
                )
                nf.set_name(f_name)  # type: ignore
                new_fc.add(nf)
    new_fc.write(output_file)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        km = int(sys.argv[1])
    else:
        km = 200

    print(km)

    radians = km / pygplates.Earth.equatorial_radius_in_kms

    print(radians)

    # when use degree for the distance
    # degrees=2
    # radians=math.radians(degrees)

    out_file = "./output/test_discretize_polyline_results.gpmlz"
    in_file = f"{os.path.dirname(os.path.realpath(__file__))}/../data/test_discretize_polyline.gpmlz"
    main(in_file, radians, out_file)
    print(f"The points have been saved to {out_file} ")

    out_file = "./output/test_discretize_polyline_results_not_strict.gpmlz"
    main(in_file, radians, out_file, strict_mode=False)
    print(f"The points have been saved to {out_file} ")
