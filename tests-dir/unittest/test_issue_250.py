#!/usr/bin/env python3
import os
import sys

os.environ["GPLATELY_DEBUG"] = "1"
os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModelManager

import gplately

print(gplately.__file__)

# https://github.com/GPlates/gplately/issues/250


def main(show=True):
    # Call GPlately's PlateModelManager object and request data from the Müller et al. 2019 study
    pm_manager = PlateModelManager()
    muller2019_model = pm_manager.get_model("Muller2019", data_dir=MODEL_REPO_DIR)
    if not muller2019_model:
        raise Exception("Failed to get reconstruction model!")
    rotation_model = muller2019_model.get_rotation_model()
    topology_features = muller2019_model.get_topologies()
    static_polygons = muller2019_model.get_static_polygons()

    model = gplately.PlateReconstruction(
        rotation_model, topology_features, static_polygons
    )

    # Obtain features for the PlotTopologies object with PlateModelManager
    coastlines = muller2019_model.get_layer("Coastlines")
    continents = muller2019_model.get_layer("ContinentalPolygons")
    COBs = muller2019_model.get_layer("COBs")

    # Call the PlotTopologies object
    gplot = gplately.plot.PlotTopologies(
        model, coastlines=coastlines, continents=continents, COBs=COBs
    )

    gplot.time = 100  # Ma

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111, projection=ccrs.Mollweide(190))
    ax.set_global()  # type: ignore
    gplot.plot_ridges_and_transforms(ax, color="red")
    gplot.plot_ridges(ax, color="black")
    gplot.plot_transforms(ax, color="black")
    gplot.plot_misc_transforms(ax, color="black")

    print(gplot.ridges)
    print(gplot.transforms)
    print(gplot.get_ridges_and_transforms())
    print(gplot.get_misc_transforms())
    print(gplot.get_ridges())
    print(gplot.get_transforms())

    plt.title(f"{gplot.time} Ma")

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
