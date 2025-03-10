#!/usr/bin/env python3


import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from plate_model_manager import PlateModelManager

from gplately import PlateReconstruction, PlotTopologies

# this example demonstrates how to use the plate_model_manager with GPlately.


def main():
    # use PlateModelManager to create PlateReconstruction and PlotTopologies objects
    model = PlateModelManager().get_model("Muller2019")
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

    fig = plt.figure(figsize=(12, 6), dpi=72)
    ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=180))
    ax.set_global()

    # now use PlotTopologies object to plot some model data
    gplot.plot_continent_ocean_boundaries(ax, color="cornflowerblue")
    gplot.plot_coastlines(ax, color="black")
    gplot.plot_ridges(ax, color="red")
    gplot.plot_trenches(ax, color="orange")
    gplot.plot_subduction_teeth(ax, color="orange")

    plt.title(f"{age} Ma")

    # save the map as a .png file
    output_file = f"working_with_pmm.png"
    plt.gcf().savefig(output_file, dpi=120, bbox_inches="tight")  # transparent=True)
    print(f"Done! The {output_file} has been saved.")
    plt.close(plt.gcf())


if __name__ == "__main__":
    main()
