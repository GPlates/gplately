#!/usr/bin/env python3
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from common import save_fig
from plate_model_manager import PlateModelManager

sys.path.insert(0, "../..")
import gplately


def main(show=True):
    # Call GPlately's PlateModelManager object and request data from the Müller et al. 2019 study
    pm_manager = PlateModelManager()
    muller2019_model = pm_manager.get_model("Muller2019", data_dir="plate-model-repo")
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

    gplot.time = 10  # Ma

    # Download all Muller et al. 2019 netCDF age grids with PlateModelManager. This is returned as a Raster object.
    agegrid = gplately.Raster(
        data=muller2019_model.get_raster("AgeGrids", int(gplot.time))
    )

    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(111, projection=ccrs.Mollweide(190))

    gplot.plot_continents(ax1, facecolor="0.8")
    gplot.plot_coastlines(ax1, color="0.5")
    gplot.plot_ridges_and_transforms(ax1, color="red")
    gplot.plot_trenches(ax1, color="k")
    gplot.plot_subduction_teeth(ax1, color="k")
    im = gplot.plot_grid(ax1, agegrid.data, cmap="YlGnBu", vmin=0, vmax=200)
    gplot.plot_plate_motion_vectors(
        ax1, spacingX=10, spacingY=10, normalise=True, zorder=10, alpha=0.5
    )

    fig.colorbar(im, orientation="horizontal", shrink=0.4, pad=0.05, label="Age (Ma)")
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
