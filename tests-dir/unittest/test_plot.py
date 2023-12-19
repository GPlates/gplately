#!/usr/bin/env python3

import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from plate_model_manager import PlateModelManager

sys.path.insert(0, "../..")
from common import save_fig

from gplately import PlateReconstruction, PlotTopologies

# test the plot function with the new PlateModel class

# MODEL_NAME = "Clennett2020"
MODEL_NAME = "Muller2019"


def main(show=True):
    pm_manager = PlateModelManager()

    age = 55
    model = pm_manager.get_model(MODEL_NAME)
    model.set_data_dir("plate-model-repo")

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

    gplot.plot_continent_ocean_boundaries(ax, color="cornflowerblue")
    gplot.plot_coastlines(ax, color="black")
    gplot.plot_ridges_and_transforms(ax, color="red")
    gplot.plot_trenches(ax, color="orange")
    gplot.plot_subduction_teeth(ax, color="orange")
    ax.set_global()

    ids = set([f.get_reconstruction_plate_id() for f in gplot.topologies])
    for id in ids:
        gplot.plot_plate_id(ax, id, facecolor="None", edgecolor="lightgreen")
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
