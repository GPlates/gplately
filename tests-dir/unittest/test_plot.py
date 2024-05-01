#!/usr/bin/env python3

import os
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from plate_model_manager import PlateModel, PlateModelManager

if "GPLATELY_DEBUG" in os.environ and os.environ["GPLATELY_DEBUG"].lower() == "true":
    sys.path.insert(0, f"{os.path.dirname(os.path.realpath(__file__))}/../..")

from common import MODEL_REPO_DIR, save_fig

import gplately
from gplately import PlateReconstruction, PlotTopologies

print(gplately.__file__)

# test the plot function with the new PlateModel class

# MODEL_NAME = "Clennett2020"
MODEL_NAME = "Muller2019"


def main(show=True):
    try:
        model = PlateModelManager().get_model(MODEL_NAME, data_dir=MODEL_REPO_DIR)
    except:
        model = PlateModel(MODEL_NAME, data_dir=MODEL_REPO_DIR, readonly=True)

    age = 55

    test_model = PlateReconstruction(
        model.get_rotation_model(),
        topology_features=model.get_layer("Topologies"),
        static_polygons=model.get_layer("StaticPolygons"),
    )
    gplot = PlotTopologies(
        test_model,
        coastlines=model.get_layer("Coastlines"),
        COBs=model.get_layer("COBs"),
        time=age,
    )

    fig = plt.figure(figsize=(10, 5), dpi=96)
    ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=180))

    all_flag = 1
    plot_flag = {
        "continent_ocean_boundaries": 0,
        "coastlines": 0,
        "ridges_and_transforms": 0,
        "trenches": 0,
        "subduction_teeth": 0,
        "ridges": 1,
        "all_topologies": 0,
        "all_topological_sections": 0,
        "plate_polygon_by_id": 0,
        "unclassified_features": 0,
        "misc_transforms": 0,
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
    }

    for key in plot_flag:
        if key == "plate_polygon_by_id":
            continue
        if all_flag or plot_flag[key]:
            getattr(gplot, f"plot_{key}")(
                ax, color=list(np.random.choice(range(256), size=3) / 256)
            )

    ax.set_global()

    ids = set([f.get_reconstruction_plate_id() for f in gplot.topologies])
    for id in ids:
        if all_flag or plot_flag["plate_polygon_by_id"]:
            gplot.plot_plate_polygon_by_id(
                ax, id, facecolor="None", edgecolor="lightgreen"
            )
    plt.title(f"{age} Ma")

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
