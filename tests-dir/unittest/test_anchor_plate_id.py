#!/usr/bin/env python3

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys

sys.path.insert(0, "../")
import gplately


def main():
    gdownload = gplately.download.DataServer("Clennett2020")
    (
        C2020_rotation_file,
        C2020_topology_features,
        C2020_static_polygons,
    ) = gdownload.get_plate_reconstruction_files()

    C2020_501 = gplately.PlateReconstruction(
        C2020_rotation_file,
        C2020_topology_features,
        C2020_static_polygons,
        default_anchor_plate_id=501,
    )

    C2020_101 = gplately.PlateReconstruction(
        C2020_rotation_file,
        C2020_topology_features,
        C2020_static_polygons,
        default_anchor_plate_id=101,
    )

    gplot501 = gplately.PlotTopologies(C2020_501, time=130)
    gplot101 = gplately.PlotTopologies(C2020_101, time=130)

    fig, ax = plt.subplots(
        figsize=(10, 5), subplot_kw={"projection": ccrs.Robinson()}, dpi=120
    )

    gplot501.plot_all_topologies(ax=ax, color="red")
    gplot101.plot_all_topologies(ax=ax, color="blue")

    plt.show()


if __name__ == "__main__":
    main()
