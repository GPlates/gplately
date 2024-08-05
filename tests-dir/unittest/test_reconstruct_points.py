#!/usr/bin/env python3
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModelManager

import gplately

print(gplately.__file__)


def main(show=True):
    pm_manager = PlateModelManager()
    muller2019_model = pm_manager.get_model("Muller2019", data_dir=MODEL_REPO_DIR)
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

    pt_lons = np.array([140.0, 150.0, 160.0])
    pt_lats = np.array([-30.0, -40.0, -50.0])

    gpts = gplately.Points(model, pt_lons, pt_lats)

    rlons = np.empty((21, pt_lons.size))
    rlats = np.empty((21, pt_lons.size))

    for time in range(0, 21):
        rlons[time], rlats[time] = gpts.reconstruct(time, return_array=True)

    gplot.time = 0  # present day

    fig = plt.figure(figsize=(6, 8))
    ax1 = fig.add_subplot(111, projection=ccrs.Mercator(190))
    ax1.set_extent([130, 180, -60, -10])

    gplot.plot_coastlines(ax1, color="0.8")

    for i in range(0, len(pt_lons)):
        ax1.plot(rlons[:, i], rlats[:, i], "o", transform=ccrs.PlateCarree())

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
