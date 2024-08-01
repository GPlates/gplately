#!/usr/bin/env python3

import math
import multiprocessing
from functools import partial
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import moviepy.editor as mpy
from plate_model_manager import PlateModel, PlateModelManager

import gplately
from gplately import ptt

print(gplately.__file__)

MODEL_REPO_DIR = "../unittest/plate-model-repo"
OUTPUT_DIR = "./output-debug-ridges-and-transforms"

model_name = "merdith2021"
angles = [65, 70]
# angles = [45, 50, 55, 60, 65, 70, 75]
ages = range(1001)

try:
    model = PlateModelManager().get_model(model_name, data_dir=MODEL_REPO_DIR)
except:
    model = PlateModel(model_name, data_dir=MODEL_REPO_DIR, readonly=True)


def plot_and_save_ridges_and_transforms(age, angle, output_path):
    fig = plt.figure(figsize=(24, 12), dpi=96)
    ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=180))
    # plt.tight_layout()
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
    output_file = f"{output_path}/{age}.png"
    fig.savefig(output_file, dpi=120)  # transparent=True)
    print(f"Done! The {output_file} has been saved.")
    plt.close(fig)


def main():
    for angle in angles:
        print(f"plotting angle:{angle}")
        output_path = f"{OUTPUT_DIR}/angle-{angle}"
        Path(output_path).mkdir(parents=True, exist_ok=True)
        try:
            num_cpus = multiprocessing.cpu_count()
        except NotImplementedError:
            num_cpus = 1

        if num_cpus > 1:
            with multiprocessing.Pool(num_cpus) as pool:
                pool.map(
                    partial(
                        plot_and_save_ridges_and_transforms,
                        angle=angle,
                        output_path=output_path,
                    ),
                    ages,
                )
        else:
            for age in ages:
                plot_and_save_ridges_and_transforms(
                    angle=angle, age=age, output_path=output_path
                )

        frame_list = [f"{output_path}/{age}.png" for age in ages]
        frame_list.reverse()
        # print(frame_list)
        clip = mpy.ImageSequenceClip(frame_list, fps=6)

        vfn = f"{OUTPUT_DIR}/ridges_and_transforms_angle-{angle}.mp4"
        clip.write_videofile(
            vfn,
            codec="libx264",
            # audio_codec='aac',
            ffmpeg_params=["-s", "2880, 1440", "-pix_fmt", "yuv420p"],
        )  # give image size here(the numbers must divide by 2)


if __name__ == "__main__":
    main()
