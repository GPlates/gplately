#!/usr/bin/env python3

import math
import sys
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from common import MODEL_REPO_DIR, OUTPUT_DIR, save_fig
from plate_model_manager import PlateModel, PlateModelManager

import gplately
from gplately import ptt

print(gplately.__file__)

model_name = "muller2019"
angles = [45, 50, 55, 60, 65, 70, 75]
age = 0

try:
    model = PlateModelManager().get_model(model_name, data_dir=MODEL_REPO_DIR)
except:
    model = PlateModel(model_name, data_dir=MODEL_REPO_DIR, readonly=True)


def plot_ridges_and_transforms(ax, angle):
    print(f"plotting angle:{angle}")

    (
        topologies,
        ridge_transforms,
        ridges,
        transforms,
        trenches,
        trench_left,
        trench_right,
        other,
    ) = ptt.resolve_topologies.resolve_topologies_into_features(
        model.get_rotation_model(),
        model.get_layer("Topologies"),
        age,
        transform_segment_deviation_in_radians=math.radians(angle),
    )

    for t in topologies:
        for geom in t.get_all_geometries():
            if not geom:
                continue
            lats, lons = zip(*geom.to_lat_lon_list())
            ax.plot(
                lons,
                lats,
                color="lightgrey",
                transform=ccrs.Geodetic(),
            )

    for rt in ridge_transforms:
        for geom in rt.get_all_geometries():
            if not geom:
                continue
            lats, lons = zip(*geom.to_lat_lon_list())
            ax.plot(
                lons,
                lats,
                color="black",
                transform=ccrs.Geodetic(),
            )

    for transform in transforms:
        for geom in transform.get_all_geometries():
            if not geom:
                continue
            lats, lons = zip(*geom.to_lat_lon_list())
            ax.plot(
                lons,
                lats,
                color="blue",
                transform=ccrs.Geodetic(),
            )

    for ridge in ridges:
        for geom in ridge.get_all_geometries():
            if not geom:
                continue
            lats, lons = zip(*geom.to_lat_lon_list())
            ax.plot(
                lons,
                lats,
                color="red",
                transform=ccrs.Geodetic(),
            )

    ax.set_global()

    ax.set_title(f"{model_name}({age}Ma) angle:{angle}", fontsize=16)


def main(show=True):
    n_column = 3
    n_row = len(angles) // n_column + 1
    # print(n_row)
    fig, axs = plt.subplots(
        n_row,
        n_column,
        subplot_kw={"projection": ccrs.Robinson(central_longitude=180)},
        figsize=(n_column * 6, n_row * 3),
        constrained_layout=True,
    )
    # plt.tight_layout()

    for i in range(len(angles)):
        plot_ridges_and_transforms(axs.flatten()[i], angles[i])

    i = len(angles)
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


def plot_and_save_separately():
    for angle in angles:
        fig = plt.figure(figsize=(24, 12), dpi=96)
        ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=180))
        plt.tight_layout()
        plot_ridges_and_transforms(ax, angle)
        Path(f"{OUTPUT_DIR}/ridges_and_transforms").mkdir(parents=True, exist_ok=True)
        output_file = (
            f"{OUTPUT_DIR}/ridges_and_transforms/ridges_and_transforms_{angle}.png"
        )
        fig.savefig(output_file, dpi=120, bbox_inches="tight")  # transparent=True)
        print(f"Done! The {output_file} has been saved.")
        plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save_all":
        plot_and_save_separately()
    elif len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
