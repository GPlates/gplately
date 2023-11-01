#!/usr/bin/env python3

import sys

import cartopy
import matplotlib.pyplot as plt
from plate_model_manager import PlateModelManager

sys.path.insert(0, "../..")
import gplately

# test plotting subduction teeth


def main():
    pm_manager = PlateModelManager()
    model = pm_manager.get_model("Muller2019")
    model.set_data_dir("plate-model-repo")

    plt.figure(figsize=(12, 8))

    ax1 = plt.subplot(211, projection=cartopy.crs.PlateCarree())
    ax2 = plt.subplot(212, projection=cartopy.crs.Robinson())

    plate_model = gplately.PlateReconstruction(
        model.get_rotation_model(),
        topology_features=model.get_layer("Topologies"),
        static_polygons=model.get_layer("StaticPolygons"),
    )
    plot_plates = gplately.PlotTopologies(plate_model, time=100)
    plot_plates.time = 100

    ax1.set_extent([50, 105, -10, 40], crs=cartopy.crs.PlateCarree())
    plot_plates.plot_ridges_and_transforms(ax1, color="r")
    plot_plates.plot_trenches(ax1, color="b")
    plot_plates.plot_faults(ax1, color="k")
    plot_plates.plot_subduction_teeth(ax1, color="green")

    ax2.set_global()
    plot_plates.plot_ridges_and_transforms(ax2, color="r")
    plot_plates.plot_trenches(ax2, color="b")
    plot_plates.plot_faults(ax2, color="k")
    plot_plates.plot_subduction_teeth(ax2, color="green")

    plt.show()


if __name__ == "__main__":
    main()
