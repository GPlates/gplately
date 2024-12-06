#!/usr/bin/env python3

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from plate_model_manager import PlateModelManager

import gplately

print(gplately.__file__)

# This file is used to debug github Issue #214.
# https://github.com/GPlates/gplately/issues/214
# Once the problem has been fixed, the exception should go away and the plot should be successful.


def main():
    # Call GPlately's PlateModelManager object and request data from the Müller et al. 2016 study
    pm_manager = PlateModelManager()
    muller2016_model = pm_manager.get_model("Muller2016", data_dir=".")
    rotation_model = muller2016_model.get_rotation_model()
    topology_features = muller2016_model.get_topologies()
    static_polygons = muller2016_model.get_static_polygons()

    model = gplately.PlateReconstruction(
        rotation_model, topology_features, static_polygons
    )

    # Obtain features for the PlotTopologies object with PlateModelManager
    coastlines = muller2016_model.get_layer("Coastlines")
    continents = None
    COBs = muller2016_model.get_layer("COBs")

    # Call the PlotTopologies object
    gplot = gplately.plot.PlotTopologies(
        model, coastlines=coastlines, continents=continents, COBs=COBs
    )

    gplot.time = 10  # Ma

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111, projection=ccrs.Mollweide(190))

    ax.set_global()
    ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=2,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )

    gplot.plot_continent_ocean_boundaries(ax, color="red")

    plt.title(f"{gplot.time} Ma")

    plt.show()


if __name__ == "__main__":
    main()
