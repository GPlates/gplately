#!/usr/bin/env python3

import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from plate_model_manager import PlateModelManager

sys.path.insert(0, "../")
from common import MODEL_REPO_DIR, save_fig

import gplately

MODEL_NAME = "Clennett2020"
# MODEL_NAME = "Muller2019"


def main(show=True):
    pm_manger = PlateModelManager()
    model = pm_manger.get_model(MODEL_NAME, data_dir=MODEL_REPO_DIR)
    if not model:
        raise Exception(f"Unable to get model {MODEL_NAME}!!!")

    C2020_rotation_file = model.get_rotation_model()
    C2020_topology_features = model.get_topologies()
    C2020_static_polygons = model.get_static_polygons()

    C2020_501 = gplately.PlateReconstruction(
        C2020_rotation_file,
        C2020_topology_features,
        C2020_static_polygons,
        anchor_plate_id=501,
    )

    C2020_101 = gplately.PlateReconstruction(
        C2020_rotation_file,
        C2020_topology_features,
        C2020_static_polygons,
        anchor_plate_id=101,
    )

    gplot501 = gplately.PlotTopologies(C2020_501, time=130)
    gplot101 = gplately.PlotTopologies(C2020_101, time=130)

    fig, ax = plt.subplots(
        figsize=(10, 5), subplot_kw={"projection": ccrs.Robinson()}, dpi=120
    )

    gplot501.plot_all_topologies(ax=ax, color="red")
    gplot101.plot_all_topologies(ax=ax, color="blue")

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
