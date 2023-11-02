#!/usr/bin/env python3


import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from plate_model_manager import PlateModelManager

from gplately import PlateReconstruction, PlotTopologies

# This is a simple example of how to use the Plate Model Manager with GPlately


def main():
    # Here is how to use PlateModelManager to create PlateReconstruction and PlotTopologies objects
    pm_manager = PlateModelManager()
    model = pm_manager.get_model("Muller2019")
    model.set_data_dir("plate-model-repo")

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

    # Now do some plotting
    fig = plt.figure(figsize=(12, 6), dpi=72)
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
    plt.show()


if __name__ == "__main__":
    main()
