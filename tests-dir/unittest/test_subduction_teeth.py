#!/usr/bin/env python3

import os
import sys

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import cartopy
import matplotlib.pyplot as plt
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModelManager

import gplately

print(gplately.__file__)

# test plotting subduction teeth


def main(show=True):
    pm_manager = PlateModelManager()
    model = pm_manager.get_model("Muller2019", data_dir=MODEL_REPO_DIR)

    plt.figure(figsize=(12, 8))

    ax1 = plt.subplot(211, projection=cartopy.crs.PlateCarree())  # type: ignore
    ax2 = plt.subplot(212, projection=cartopy.crs.Robinson())  # type: ignore
    assert model
    plate_model = gplately.PlateReconstruction(
        model.get_rotation_model(),
        topology_features=model.get_layer("Topologies"),
        static_polygons=model.get_layer("StaticPolygons"),
    )
    plot_plates = gplately.PlotTopologies(plate_model, time=100)
    plot_plates.time = 100

    ax1.set_extent([50, 105, -10, 40], crs=cartopy.crs.PlateCarree())  # type: ignore
    plot_plates.plot_ridges(ax1, color="r")
    plot_plates.plot_trenches(ax1, color="b")
    plot_plates.plot_faults(ax1, color="k")
    plot_plates.plot_subduction_teeth(ax1, color="green")

    ax2.set_global()  # type: ignore
    plot_plates.plot_ridges(ax2, color="r")
    plot_plates.plot_trenches(ax2, color="b")
    plot_plates.plot_faults(ax2, color="k")
    plot_plates.plot_subduction_teeth(ax2, color="green")

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
