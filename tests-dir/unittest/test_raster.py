#!/usr/bin/env python3
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModelManager

sys.path.insert(0, "../..")
import gplately


def main(show=True):
    pm_manager = PlateModelManager()
    plate_model = pm_manager.get_model("Muller2019", data_dir="plate-model-repo")

    rotation_model = plate_model.get_rotation_model()
    topology_features = plate_model.get_topologies()
    static_polygons = plate_model.get_static_polygons()

    coastlines = plate_model.get_layer("Coastlines")
    continents = plate_model.get_layer("ContinentalPolygons")
    COBs = plate_model.get_layer("COBs")

    model = gplately.PlateReconstruction(
        rotation_model, topology_features, static_polygons
    )
    gplot = gplately.PlotTopologies(model, coastlines, continents, COBs)

    age_grid_raster = gplately.Raster(
        data=plate_model.get_raster("AgeGrids", 100),
        plate_reconstruction=model,
        extent=[-180, 180, -90, 90],
    )

    xx, yy = np.meshgrid(np.linspace(-180, 180, 180), np.linspace(-90, 90, 90))

    values = age_grid_raster.query(lons=xx.flatten(), lats=yy.flatten())
    print(values[~np.isnan(values)])

    fig = plt.figure(figsize=(10, 5), dpi=96)
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    ax.scatter(
        xx.flatten(),
        yy.flatten(),
        c=values,
        marker="s",
        s=1,
        transform=ccrs.PlateCarree(),
        cmap="YlGnBu",
        vmax=200,
        vmin=0,
    )

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
