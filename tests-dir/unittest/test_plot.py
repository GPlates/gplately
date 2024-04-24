#!/usr/bin/env python3

import os
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from plate_model_manager import PlateModelManager

if "GPLATELY_DEBUG" in os.environ and os.environ["GPLATELY_DEBUG"].lower() == "true":
    sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from common import MODEL_REPO_DIR, save_fig

import gplately
from gplately import PlateReconstruction, PlotTopologies

print(gplately.__file__)

# test the plot function with the new PlateModel class

# MODEL_NAME = "Clennett2020"
MODEL_NAME = "Muller2019"


def main(show=True):
    pm_manager = PlateModelManager()

    age = 55
    model = pm_manager.get_model(MODEL_NAME, data_dir=MODEL_REPO_DIR)

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

    plot_flag = {
        "continent_ocean_boundaries": 0,
        "coastlines": 0,
        "ridges_and_transforms": 0,
        "trenches": 0,
        "subduction_teeth": 0,
        "ridges": 1,
        "all_topologies": 0,
        "all_topological_sections": 0,
        "plot_plate_polygon_by_id": 0,
    }

    if plot_flag["continent_ocean_boundaries"]:
        gplot.plot_continent_ocean_boundaries(ax, color="cornflowerblue")
    if plot_flag["coastlines"]:
        gplot.plot_coastlines(ax, color="black")
    if plot_flag["ridges_and_transforms"]:
        gplot.plot_ridges_and_transforms(ax, color="red")
    if plot_flag["trenches"]:
        gplot.plot_trenches(ax, color="orange")
    if plot_flag["subduction_teeth"]:
        gplot.plot_subduction_teeth(ax, color="orange")
    if plot_flag["ridges"]:
        gplot.plot_ridges(ax, color="green")
    if plot_flag["all_topologies"]:
        gplot.plot_all_topologies(ax, color="red")
    if plot_flag["all_topological_sections"]:
        gplot.plot_all_topological_sections(ax, color="red")
    ax.set_global()

    ids = set([f.get_reconstruction_plate_id() for f in gplot.topologies])
    for id in ids:
        if plot_flag["plot_plate_polygon_by_id"]:
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
