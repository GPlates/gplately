#!/usr/bin/env python3

import random
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModel, PlateModelManager

import gplately
from gplately import PlateReconstruction, PlotTopologies

print(gplately.__file__)

MODELS = [
    "Alfonso2024",
    "Muller2019",
    "Clennett2020",
    "Muller2022",
    "Muller2016",
    "Merdith2021",
    "seton2012",
]


def plot_model(ax, model_name):
    print(f"plotting {model_name}")
    try:
        model = PlateModelManager().get_model(model_name, data_dir=MODEL_REPO_DIR)
    except:
        model = PlateModel(model_name, data_dir=MODEL_REPO_DIR, readonly=True)

    age = random.randint(0, 140)

    test_model = PlateReconstruction(
        model.get_rotation_model(),
        topology_features=model.get_layer("Topologies"),
        static_polygons=model.get_layer("StaticPolygons"),
        plate_model_name=model_name,
    )
    gplot = PlotTopologies(
        test_model,
        coastlines=model.get_layer("Coastlines"),
        COBs=model.get_layer("COBs", return_none_if_not_exist=True),
        continents=model.get_layer(
            "ContinentalPolygons", return_none_if_not_exist=True
        ),
        time=age,
    )

    try:
        gplot.plot_continent_ocean_boundaries(
            ax, color=list(np.random.choice(range(256), size=3) / 256)
        )
    except:
        pass

    gplot.plot_coastlines(ax, color=list(np.random.choice(range(256), size=3) / 256))
    gplot.plot_all_topologies(
        ax, color=list(np.random.choice(range(256), size=3) / 256)
    )
    gplot.plot_continents(ax, color="grey", facecolor="0.8")

    ax.set_global()

    ax.set_title(f"{model_name}({age}Ma)", fontsize=8)


def main(show=True):
    n_column = 4
    n_row = len(MODELS) // n_column + 1
    # print(n_row)
    fig, axs = plt.subplots(
        n_row,
        n_column,
        subplot_kw={"projection": ccrs.Robinson(central_longitude=180)},
        figsize=(n_column * 3, n_row * 2),
    )

    for i in range(len(MODELS)):
        plot_model(axs.flatten()[i], MODELS[i])

    i = len(MODELS)
    while i < len(axs.flatten()):
        axs.flatten()[i].remove()
        i += 1

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
