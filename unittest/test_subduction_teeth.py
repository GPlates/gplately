#!/usr/bin/env python3

import sys

import cartopy
import matplotlib.pyplot as plt

sys.path.insert(0, "../")
import gplately


def main():
    gdownload = gplately.download.DataServer("Muller2019")
    (
        C2020_rotation_file,
        C2020_topology_features,
        C2020_static_polygons,
    ) = gdownload.get_plate_reconstruction_files()

    fig, ax = plt.subplots(
        subplot_kw={"projection": cartopy.crs.PlateCarree()}, figsize=(12, 8)
    )
    ax.set_extent([50, 105, -10, 40], crs=cartopy.crs.PlateCarree())
    plate_model = gplately.PlateReconstruction(
        C2020_rotation_file, C2020_topology_features, C2020_static_polygons
    )
    plot_plates = gplately.PlotTopologies(plate_model, time=100)
    plot_plates.time = 100
    plot_plates.plot_ridges_and_transforms(ax, color="r")
    plot_plates.plot_trenches(ax, color="b")
    plot_plates.plot_faults(ax, color="k")
    plot_plates.plot_subduction_teeth(ax, color="green")

    plt.show()


if __name__ == "__main__":
    main()
