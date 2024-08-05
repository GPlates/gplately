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
    # Call GPlately's PlateModelManager object and request data from the Müller et al. 2019 study
    pm_manager = PlateModelManager()
    muller2019_model = pm_manager.get_model("Muller2019", data_dir=MODEL_REPO_DIR)
    rotation_model = muller2019_model.get_rotation_model()
    topology_features = muller2019_model.get_topologies()
    static_polygons = muller2019_model.get_static_polygons()

    model = gplately.PlateReconstruction(
        rotation_model, topology_features, static_polygons
    )

    time = 50  # Ma
    tessellation_threshold_radians = np.radians(10)
    subduction_data = model.tessellate_subduction_zones(
        time,
        tessellation_threshold_radians=tessellation_threshold_radians,
        ignore_warnings=True,
    )
    ridge_data = model.tessellate_mid_ocean_ridges(
        time,
        tessellation_threshold_radians=tessellation_threshold_radians,
        ignore_warnings=True,
    )

    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(111, projection=ccrs.Mollweide())
    ax1.set_global()

    ax1.scatter(
        subduction_data[:, 0],  # longitude
        subduction_data[:, 1],  # latitude
        color="blue",
        s=5,
        transform=ccrs.PlateCarree(),
    )
    ax1.scatter(
        ridge_data[:, 0],  # longitude
        ridge_data[:, 1],  # latitude
        color="red",
        s=5,
        transform=ccrs.PlateCarree(),
    )

    plt.title(f"{time} Ma")

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
