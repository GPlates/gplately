#!/usr/bin/env python3

import os
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModelManager

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import gplately

print(gplately.__file__)

MODEL_NAME = "Clennett2020"
# MODEL_NAME = "Muller2019"

# test reconstruction with different anchor PIDs.


def main(show=True):
    pm_manger = PlateModelManager()
    model = pm_manger.get_model(MODEL_NAME, data_dir=MODEL_REPO_DIR)
    if not model:
        raise Exception(f"Unable to get model {MODEL_NAME}!!!")

    C2020_rotation_file = model.get_rotation_model()
    C2020_topology_features = model.get_topologies()
    C2020_static_polygons = model.get_static_polygons()

    time = 80
    anchor_plate_id_1 = 501
    anchor_plate_id_2 = 101

    C2020_1 = gplately.PlateReconstruction(
        C2020_rotation_file,
        C2020_topology_features,
        C2020_static_polygons,
        anchor_plate_id=anchor_plate_id_1,
    )

    C2020_2 = gplately.PlateReconstruction(
        C2020_rotation_file,
        C2020_topology_features,
        C2020_static_polygons,
        anchor_plate_id=anchor_plate_id_2,
    )

    C2020_coastlines = model.get_coastlines()

    # Anchor plate should default to the anchor plate of 'gplately.PlateReconstruction'.
    gplot_1 = gplately.PlotTopologies(C2020_1, coastlines=C2020_coastlines, time=time)
    gplot_2 = gplately.PlotTopologies(C2020_2, coastlines=C2020_coastlines, time=time)

    C2020_present_day_agegrid_1 = gplately.Raster(
        data=model.get_raster("AgeGrids", time=0), plate_reconstruction=C2020_1
    )
    C2020_present_day_agegrid_2 = gplately.Raster(
        data=model.get_raster("AgeGrids", time=0), plate_reconstruction=C2020_2
    )

    C2020_reconstructed_agegrid_1 = C2020_present_day_agegrid_1.reconstruct(
        time, threads=4
    )
    C2020_reconstructed_agegrid_2 = C2020_present_day_agegrid_2.reconstruct(
        time, threads=4
    )

    # Show two plots, one for each anchor plate ID.
    fig, axs = plt.subplots(
        2,
        1,
        figsize=(10, 10),
        subplot_kw={"projection": ccrs.Robinson()},
        dpi=120,
    )
    ax_1 = axs[0]
    ax_2 = axs[1]

    # First anchor plate plot.
    gplot_1.plot_all_topologies(ax=ax_1, color="red")
    gplot_1.plot_coastlines(ax=ax_1, color="0.5")
    gplot_1.plot_grid(
        ax_1, C2020_reconstructed_agegrid_1, cmap="YlGnBu", vmin=0, vmax=200
    )
    ax_1.set_global()
    ax_1.set_title(f"Anchor plate {anchor_plate_id_1}", fontsize=8)

    # Second anchor plate plot.
    gplot_2.plot_all_topologies(ax=ax_2, color="blue")
    gplot_2.plot_coastlines(ax=ax_2, color="0.5")
    gplot_2.plot_grid(
        ax_2, C2020_reconstructed_agegrid_2, cmap="YlGnBu", vmin=0, vmax=200
    )
    ax_2.set_global()
    ax_2.set_title(f"Anchor plate {anchor_plate_id_2}", fontsize=8)

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
