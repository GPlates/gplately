#!/usr/bin/env python3
# import matplotlib

# matplotlib.use("QtAgg")

import os
import sys

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModel, PlateModelManager

import gplately
from gplately import PlateReconstruction, PlotTopologies

print(gplately.__file__)

# test the plot function with the new PlateModel class

# MODEL_NAME = "Clennett2020"
MODEL_NAME = "Muller2019"
# MODEL_NAME = "merdith2021"


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

    # age = 100
    # gplot.time = age

    fig = plt.figure(figsize=(10, 5), dpi=96)
    ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=180))

    all_flag = 1
    plot_flag = {
        "continent_ocean_boundaries": 0,
        "coastlines": 0,
        "trenches": 0,
        "subduction_teeth": 0,
        "ridges": 0,
        "all_topologies": 0,
        "all_topological_sections": 0,
        "plate_polygon_by_id": 1,
        "unclassified_features": 0,
        "slab_edges": 0,
        "passive_continental_boundaries": 0,
        "extended_continental_crusts": 0,
        "continental_crusts": 0,
        "sutures": 0,
        "orogenic_belts": 0,
        "transitional_crusts": 0,
        "terrane_boundaries": 0,
        "inferred_paleo_boundaries": 0,
        "fracture_zones": 0,
        "faults": 0,
        "continental_rifts": 0,
        "misc_boundaries": 0,
        "transforms": 0,
        "continents": 0,
        "topological_plate_boundaries": 0,
    }

    for key in plot_flag:
        if key == "plate_polygon_by_id":
            continue
        if key == "continents":
            if plot_flag["continents"]:
                gplot.plot_continents(ax, color="grey", facecolor="0.8")
            continue
        if key == "coastlines":
            if plot_flag["coastlines"]:
                gplot.plot_coastlines(ax, edgecolor="blue", facecolor="0.5")
            continue
        if all_flag or plot_flag[key]:
            getattr(gplot, f"plot_{key}")(
                ax, color=list(np.random.choice(range(256), size=3) / 256)
            )

    ax.set_global()  # type: ignore

    if gplot.topologies:
        ids = set([f.get_reconstruction_plate_id() for f in gplot.topologies])
        for id in ids:
            if all_flag or plot_flag["plate_polygon_by_id"]:
                color = list(np.random.choice(range(256), size=3) / 256)
                gplot.plot_plate_polygon_by_id(
                    ax,
                    id,
                    facecolor=color + [0.5],
                    edgecolor=color,
                )

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
