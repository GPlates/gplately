#!/usr/bin/env python

import multiprocessing
from functools import partial
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import moviepy.editor as mpy
from gmt_cpt import get_cm_from_gmt_cpt
from netCDF4 import Dataset
from plate_model_manager import PlateModel, PlateModelManager

import gplately
from gplately import PlateReconstruction, PlotTopologies

age_grid_filepath = "../../../test-age-grid-output/SEAFLOOR_AGE/merdith2021_SEAFLOOR_AGE_grid_{:0.2f}Ma.nc"
times = range(400, 411, 1)

image_output_dir = "./images"
image_filepath = image_output_dir + "/agegrid-{:0.2f}-Ma.png"


def plot_map(model, time):
    test_model = PlateReconstruction(
        model.get_rotation_model(),
        topology_features=model.get_layer("Topologies"),
        static_polygons=model.get_layer("StaticPolygons"),
    )
    gplot = PlotTopologies(
        test_model,
        coastlines=model.get_layer("Coastlines"),
        continents=model.get_layer("ContinentalPolygons"),
        time=time,
    )

    fig = plt.figure(figsize=(12, 8), dpi=300)
    ax = plt.axes(projection=ccrs.PlateCarree(), frameon=True)
    img = Dataset(age_grid_filepath.format(time))  # load grid
    cb = ax.imshow(
        img.variables["z"],
        origin="lower",
        transform=ccrs.PlateCarree(),
        extent=[-180, 180, -90, 90],
        cmap=get_cm_from_gmt_cpt("agegrid.cpt"),
        vmax=350,
        vmin=0,
    )
    fig.colorbar(cb, shrink=0.5, label="Age(Ma)", orientation="horizontal", pad=0.05)
    plt.title(f"{time} Ma")
    gplot.time = time
    gplot.plot_continents(ax, facecolor="0.8", color="none")
    gplot.plot_coastlines(ax, facecolor="0.5", color="none")
    gplot.plot_all_topologies(ax, color="blue")
    # plt.show()

    Path(image_output_dir).mkdir(parents=True, exist_ok=True)
    plt.gcf().savefig(
        image_filepath.format(time), dpi=120, bbox_inches="tight"
    )  # transparent=True)
    print(f"Done! The {image_filepath.format(time)} has been saved.")
    plt.close(plt.gcf())


def main():
    try:
        model = PlateModelManager().get_model("merdith2021")
    except:
        model = PlateModel("merdith2021", readonly=True)

    try:
        num_cpus = multiprocessing.cpu_count()
    except NotImplementedError:
        num_cpus = 1

    if num_cpus > 1:
        with multiprocessing.Pool(num_cpus) as pool:
            pool.map(
                partial(plot_map, model),
                times,
            )
    else:
        for time in times:
            plot_map(model, time)

    frame_list = [image_filepath.format(time) for time in times]
    frame_list.reverse()
    # print(frame_list)
    clip = mpy.ImageSequenceClip(frame_list, fps=6)

    clip.write_videofile(
        f"agegrids.mp4",
        codec="libx264",
        # audio_codec='aac',
        ffmpeg_params=["-s", "1140x720", "-pix_fmt", "yuv420p"],
    )  # give image size here(the numbers must divide by 2)
    print("video has been created!")


if __name__ == "__main__":
    main()
