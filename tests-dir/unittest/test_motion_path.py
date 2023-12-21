#!/usr/bin/env python3
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plate_model_manager import PlateModelManager

sys.path.insert(0, "../..")
from common import MODEL_REPO_DIR, save_fig

import gplately


def main(show=True):
    pm_manager = PlateModelManager()
    muller2019_model = pm_manager.get_model("Matthews2016", data_dir=MODEL_REPO_DIR)

    rotation_model = muller2019_model.get_rotation_model()
    topology_features = muller2019_model.get_topologies()
    static_polygons = muller2019_model.get_static_polygons()

    coastlines = muller2019_model.get_layer("Coastlines")

    # Call the PlateReconstruction object to create a plate motion model
    model = gplately.reconstruction.PlateReconstruction(
        rotation_model, topology_features, static_polygons
    )

    # Call the PlotTopologies object
    time = 0
    gplot = gplately.PlotTopologies(
        model, coastlines=coastlines, continents=None, COBs=None, time=time
    )

    # *****************************************do the first plot***************************

    lats = np.array([19])
    lons = np.array([-155])

    # Create the time array for the motion path - must be float
    start_reconstruction_time = 0.0
    time_step = 2.0
    max_reconstruction_time = 100.0
    time_array = np.arange(
        start_reconstruction_time, max_reconstruction_time + time_step, time_step
    )

    # Creating the motion path and obtaining rate plot arrays using the gplately Points alternative
    gpts = gplately.reconstruction.Points(model, lons, lats, time=0.0)
    lon, lat, _, _ = gpts.motion_path(
        time_array, anchor_plate_id=2, return_rate_of_motion=True
    )

    fig = plt.figure(figsize=(12, 6))
    gs = GridSpec(2, 2, figure=fig)
    ax = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_title("Motion path of the Hawaiian-Emperor Chain")

    # --- Limit map extent
    lon_min = 130.0
    lon_max = 220
    lat_min = 15.0
    lat_max = 60.0
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.gridlines(
        draw_labels=True, xlocs=np.arange(-180, 180, 20), ylocs=np.arange(-90, 90, 10)
    )

    # --- Make sure motion path longitudes are wrapped correctly to the dateline
    lon360 = gplately.tools.correct_longitudes_for_dateline(lon)

    # --- plot motion path
    ax.scatter(
        lon360,
        lat,
        100,
        marker=".",
        c=time_array,
        cmap=plt.cm.inferno,
        edgecolor="k",
        transform=ccrs.PlateCarree(),
        vmin=time_array[0],
        vmax=time_array[-1],
        zorder=4,
    )

    gplot.plot_coastlines(ax, color="DarkGrey")

    # ***************************** do the second plot **********************************

    lats = np.array([30, 30, 30])
    lons = np.array([73, 78, 83])

    # Plate ID attributes
    # motion of the Indian craton (501) with respect to the Northern European craton (301)
    plate_ID = 501
    anchor_plate_ID = 301

    # Create the time array for the motion path
    start_reconstruction_time = 0.0
    time_step = 5.0
    max_reconstruction_time = 105.0
    time_array = np.arange(
        start_reconstruction_time, max_reconstruction_time + 1, time_step
    )

    # Get the latitudes and longitudes of all points along the motion path
    lon, lat, step_times, step_rates_of_motion = model.create_motion_path(
        lons, lats, time_array, plate_ID, anchor_plate_ID, return_rate_of_motion=True
    )

    ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree(central_longitude=180))
    ax2.set_title("Motion path of the Indian craton to the North European Craton")

    gplot.plot_coastlines(ax2, color="DarkGrey")

    # --- Limit map extent
    lon_min = 0.0
    lon_max = 100
    lat_min = -40.0
    lat_max = 40.0
    ax2.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    ax2.gridlines(
        draw_labels=True, xlocs=np.arange(-180, 180, 20), ylocs=np.arange(-90, 90, 10)
    )
    # --- Make sure negtive motion path longitudes are wrapped correctly to the dateline
    lons = gplately.tools.correct_longitudes_for_dateline(lons)

    # --- plot motion path
    for i in np.arange(0, len(lons)):
        ax2.scatter(
            lon[:, i],
            lat[:, i],
            200,
            marker="*",
            c=time_array,
            cmap=plt.cm.inferno,
            edgecolor="C{}".format(i),
            linewidth=0.5,
            transform=ccrs.PlateCarree(),
            vmin=time_array[0],
            vmax=time_array[-1],
            zorder=4,
        )

    # ********************************do the third plot**********************************

    ax3 = fig.add_subplot(gs[1, :])
    ax3.set_title("Rate of motion between seed points and the North European Craton")
    ax3.plot(step_times, step_rates_of_motion)
    ax3.set_xlabel("Time (Ma)", fontsize=12)
    ax3.set_ylabel("Rate of motion (km/Myr)", fontsize=12)
    ax3.invert_xaxis()
    ax3.set_xlim([max_reconstruction_time, 0])
    ax3.grid(alpha=0.3)

    plt.subplots_adjust(wspace=0.25)

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
