#!/usr/bin/env python3
# import matplotlib

# matplotlib.use("QtAgg")

import os
import sys

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import (
    PlateModel,
    PlateModelManager,
)

import gplately
from gplately import PlateReconstruction, PlotTopologies

print(gplately.__file__)

MODEL_NAME = "Muller2025"


def main(show=True):
    try:
        model = PlateModelManager().get_model(MODEL_NAME, data_dir=MODEL_REPO_DIR)
    except:
        model = PlateModel(MODEL_NAME, data_dir=MODEL_REPO_DIR, readonly=True)

    if not model:
        return

    age = 55

    test_model = PlateReconstruction(
        model.get_rotation_model(),
        topology_features=model.get_layer("Topologies"),
        static_polygons=model.get_layer("StaticPolygons"),
        plate_model=model,
    )
    gplot = PlotTopologies(
        test_model,
        coastlines=model.get_layer("Coastlines"),
        COBs=model.get_layer("COBs", return_none_if_not_exist=True),
        continents=model.get_layer("ContinentalPolygons"),
        time=age,
    )

    fig = plt.figure(figsize=(10, 5), dpi=96)
    ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=180))
    ax.set_global()

    gplot.plot_misc_transforms(ax=ax, color="red", linewidth=0.5)
    gplot.misc_transforms
    gplot.get_misc_transforms()
    gplot.get_transforms()

    plt.title(f"{age} Ma")

    if show:
        # LOOK HERE! 👀👀 👇👇
        # If the figure did not show up, you need to set your matplotlib plotting backend properly.
        # On Windows, you may install PyQt and do
        # import matplotlib
        # matplotlib.use('QtAgg')

        # if you are interested in finding what backends available on your computer and what is your current backend, do the following
        # import matplotlib.rcsetup as rcsetup
        # print(rcsetup.all_backends) # get all available backends
        # import matplotlib
        # matplotlib.get_backend() # your current backend
        #
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
