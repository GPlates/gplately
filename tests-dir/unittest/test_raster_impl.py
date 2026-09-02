import os

from gplately.plot import get_topo_cmap

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModelManager
from matplotlib import image
from plate_model_manager import PresentDayRasterManager

import gplately
from gplately.auxiliary import get_gplot, get_plate_reconstruction
from gplately.utils.log_utils import turn_on_debug_logging

print(gplately.__file__)

turn_on_debug_logging()

cpt_file = "agegrid.cpt"
if not os.path.isfile(cpt_file):
    import urllib.request

    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/GPlates/gplately/refs/heads/master/tests-dir/unittest/create-age-grids-video/agegrid.cpt",
        cpt_file,
    )
from gplately.plot.gmt_cpt import get_cmap_from_gmt_cpt

model_name = "Muller2025"
plate_model = PlateModelManager().get_model(model_name, data_dir=MODEL_REPO_DIR)
np.set_printoptions(suppress=True)
assert plate_model
age_grid_raster_0_ma = gplately.Raster(
    data=plate_model.get_raster("AgeGrids", 0),
    plate_reconstruction=get_plate_reconstruction(
        model_name, model_repo_dir=MODEL_REPO_DIR
    ),
)


def test_raster_reconstruct_lon_0_360():
    age_grid_raster_0_360 = age_grid_raster_0_ma.to_longitude_positive_360()
    age_grid_raster_0_360.save_to_netcdf4("./output/age_grid_raster_0_360.nc")

    age_grid_raster_180 = age_grid_raster_0_360.to_longitude_signed_180()
    age_grid_raster_180.save_to_netcdf4("./output/age_grid_raster_180.nc")

    age_grid_raster_0_360_50ma = age_grid_raster_0_360.reconstruct(50)
    age_grid_raster_0_360_50ma.save_to_netcdf4("./output/age_grid_raster_0_360_50ma.nc")

    age_grid_raster_180_50ma = age_grid_raster_180.reconstruct(50)
    age_grid_raster_180_50ma.save_to_netcdf4("./output/age_grid_raster_180_50ma.nc")
    print(
        "Saved age_grid_raster_0_360.nc, age_grid_raster_180.nc, age_grid_raster_0_360_50ma.nc, age_grid_raster_180_50ma.nc"
    )


def test_raster_query():
    pmm_model = PlateModelManager().get_model("Muller2025", data_dir=MODEL_REPO_DIR)

    assert pmm_model

    age_grid_raster = gplately.Raster(data=pmm_model.get_raster("AgeGrids", 50))

    lons = np.linspace(-180, 180, 91)
    lats = np.linspace(-90, 90, 46)
    values = age_grid_raster.query(
        lons=lons, lats=lats, region_of_interest=8.8, pointwise=False
    )

    # print(values[~np.isnan(values)])

    fig = plt.figure(figsize=(8, 4), dpi=96)
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    xx, yy = np.meshgrid(lons, lats)
    # plot the values being returned from raster query
    ax.scatter(
        xx,
        yy,
        c=values.reshape(xx.shape),
        marker="o",
        s=2,
        transform=ccrs.PlateCarree(),
        cmap=get_cmap_from_gmt_cpt(cpt_file),
        vmax=200,
        vmin=0,
    )
    plt.show()


def test_raster_clip_by_extent():
    pmm_model = PlateModelManager().get_model("Muller2025", data_dir=MODEL_REPO_DIR)

    assert pmm_model

    age_grid_raster = gplately.Raster(data=pmm_model.get_raster("AgeGrids", 100))

    clipped_raster = age_grid_raster.clip_by_extent([-50, 50, -80, 40])

    # plot the clipped raster
    fig = plt.figure(figsize=(8, 4), dpi=96)
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    clipped_raster.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=get_cmap_from_gmt_cpt(cpt_file),
        vmax=200,
        vmin=0,
    )

    plt.show()


def test_raster_clip_by_polygons():
    pmm_model = PlateModelManager().get_model("Zahirovic2022", data_dir=MODEL_REPO_DIR)
    from gplately.utils.feature_filter import (
        FeatureNameFilter,
        PolygonAreaFilter,
        filter_feature_collection,
    )

    # find the polygon feature with area greater than 7 million and named "Australia",
    # which is the mainland of Australia, and use the Australia mainland polygon to clip a raster.
    features = filter_feature_collection(
        gplately.load_feature_collection(pmm_model.get_coastlines()),  # type: ignore
        [
            FeatureNameFilter(
                ["Australia"],
                exact_match=True,
                case_sensitive=True,
            ),
            PolygonAreaFilter(7e6, reverse=False),
        ],
    )
    assert (
        len(features) == 1
    ), f"Expected to find exactly one feature for the mainland of Australia, but found {len(features)} features."

    australia_mainland_geometry = features[0].get_geometry()
    etopo_downscaled = gplately.Raster(
        data=image.imread(PresentDayRasterManager().get_raster("ETOPO1_tif")),
        resample=(0.5, 0.5),
        origin="upper",
    )
    topo15_downscaled = gplately.Raster(
        data=PresentDayRasterManager().get_raster("topography"),
        resample=(0.5, 0.5),
    )

    # now the austrlia is in a global map
    raster_australia = etopo_downscaled.clip_by_polygons(
        [australia_mainland_geometry], fill_value="lightblue"
    )
    raster_australia_15 = topo15_downscaled.clip_by_polygons(
        [australia_mainland_geometry]
    )

    lat_lon = australia_mainland_geometry.to_lat_lon_array()
    min_lat = np.min(lat_lon[:, 0])
    max_lat = np.max(lat_lon[:, 0])
    min_lon = np.min(lat_lon[:, 1])
    max_lon = np.max(lat_lon[:, 1])
    australia_extent = (
        int(np.floor(min_lon)) - 1,
        int(np.ceil(max_lon)) + 1,
        int(np.floor(min_lat)) - 1,
        int(np.ceil(max_lat)) + 1,
    )

    clipped_australia = raster_australia.clip_by_extent(australia_extent)
    clipped_australia_15 = raster_australia_15.clip_by_extent(australia_extent)

    fig, axs = plt.subplots(
        2,
        3,
        figsize=(10, 8),
        gridspec_kw={"width_ratios": [2, 2, 1]},
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    fig.tight_layout()
    ax_1 = axs[0, 0]
    ax_2 = axs[0, 1]
    ax_3 = axs[0, 2]
    ax_4 = axs[1, 0]
    ax_5 = axs[1, 1]
    ax_6 = axs[1, 2]

    # plot the original raster
    ax_1.set_global()
    etopo_downscaled.plot(ax=ax_1)
    ax_1.set_title("Tif Raster")
    gl = ax_1.gridlines(
        draw_labels=True,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )
    gl.right_labels = False
    gl.top_labels = False

    # plot the clipped raster
    ax_2.set_global()
    raster_australia.plot(
        ax=ax_2,
    )
    gl = ax_2.gridlines(
        draw_labels=False,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )

    # plot the clipped raster with tighter extent
    ax_3.set_extent(australia_extent)
    clipped_australia.plot(ax=ax_3)

    gl = ax_3.gridlines(
        draw_labels=True,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )
    gl.left_labels = False
    gl.top_labels = False

    # plot the original raster
    ax_4.set_global()
    topo15_downscaled.plot(ax=ax_4, cmap=get_topo_cmap(), vmin=-10927, vmax=8726)
    ax_4.set_title("netCDF Raster")
    gl = ax_4.gridlines(
        draw_labels=True,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )
    gl.right_labels = False
    gl.top_labels = False

    # plot the clipped raster
    ax_5.set_global()
    raster_australia_15.plot(ax=ax_5, cmap=get_topo_cmap(), vmin=-10927, vmax=8726)

    gl = ax_5.gridlines(
        draw_labels=False,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )

    # plot the clipped raster with tighter extent
    ax_6.set_extent(australia_extent)
    clipped_australia_15.plot(
        ax=ax_6,
        cmap=get_topo_cmap(),
        vmin=-10927,
        vmax=8726,
    )
    gl = ax_6.gridlines(
        draw_labels=True,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )
    gl.left_labels = False
    gl.top_labels = False

    plt.show()


def test_upper_origin_raster_reconstruction():
    # flip the raster to upper origin for reconstruction testing
    age_grid_raster_0_ma.data = np.flipud(age_grid_raster_0_ma.data)
    age_grid_raster_0_ma.lats = np.flip(age_grid_raster_0_ma.lats)

    age_grid_raster_50_ma = age_grid_raster_0_ma.reconstruct(50)
    age_grid_raster_50_ma.save_to_netcdf4(
        "./output/age_grid_raster_upper_origin_50_ma.nc"
    )
    print("Saved age_grid_raster_upper_origin_50_ma.nc")


def test_raster_reconstruction(
    use_spatial_tree: bool = False, use_old_implementation: bool = False
):
    anchor_pid = 0
    from gplately.auxiliary import get_gplot

    gplot = get_gplot(
        "alfonso2024",
        time=0,
        default_anchor_plate_id=anchor_pid,
        model_repo_dir="plate-model-repo",
    )
    model = gplot.plate_reconstruction

    etopo = gplately.Raster(
        data=image.imread(PresentDayRasterManager().get_raster("ETOPO1_tif")),
        plate_reconstruction=model,
    )
    # etopo.lats = etopo.lats[::-1]
    etopo.data = np.flipud(etopo.data)

    etopo_downscaled = etopo.resample(0.5, 0.5)
    assert etopo_downscaled

    fig = plt.figure(figsize=(10, 5), dpi=96)
    # plot 0Ma
    ax_1 = fig.add_subplot(221, projection=ccrs.PlateCarree())
    etopo_downscaled.plot(ax_1)
    gplot.plot_continents(ax=ax_1, facecolor="none", edgecolor="black", linewidth=0.5)
    ax_1.set_title(f"0 Ma")

    # reconstruct to 50 Ma
    r_50 = etopo_downscaled.reconstruct(
        50,
        fill_value="white",
        anchor_plate_id=anchor_pid,
        use_spatial_tree=use_spatial_tree,
        use_old_implementation=use_old_implementation,
    )

    # plot 50 Ma
    ax_2 = fig.add_subplot(222, projection=ccrs.PlateCarree())
    r_50.plot(ax_2)
    gplot.time = 50
    gplot.plot_continents(ax=ax_2, facecolor="none", edgecolor="black", linewidth=0.5)
    ax_2.set_title(f"50 Ma - APID:{anchor_pid}")

    # reconstruct to 200 Ma with PID 0
    r_200 = etopo_downscaled.reconstruct(
        200,
        fill_value="white",
        anchor_plate_id=anchor_pid,
        use_spatial_tree=use_spatial_tree,
        use_old_implementation=use_old_implementation,
    )
    # plot 200 Ma
    ax_3 = fig.add_subplot(223, projection=ccrs.PlateCarree())
    r_200.plot(ax_3)
    gplot.time = 200
    gplot.plot_continents(ax=ax_3, facecolor="none", edgecolor="black", linewidth=0.5)
    ax_3.set_title(f"200 Ma - APID:{anchor_pid}")

    # reconstruct to 200 Ma with PID 701701
    anchor_pid = 701701
    r_200_701701 = etopo_downscaled.reconstruct(
        200,
        fill_value="white",
        anchor_plate_id=anchor_pid,
        use_spatial_tree=use_spatial_tree,
        use_old_implementation=use_old_implementation,
    )

    gplot = get_gplot(
        "alfonso2024",
        time=200,
        default_anchor_plate_id=anchor_pid,
        model_repo_dir="plate-model-repo",
    )

    # plot 200 Ma
    ax_4 = fig.add_subplot(224, projection=ccrs.PlateCarree())
    gplot.time = 200
    r_200_701701.plot(ax_4)
    gplot.plot_continents(ax=ax_4, facecolor="none", edgecolor="black", linewidth=0.5)
    ax_4.set_title(f"200 Ma - APID:{anchor_pid}")


def test_create_raster_from_points():
    import pygmt  # pyright: ignore[reportMissingImports]

    relief = pygmt.datasets.load_earth_relief(
        resolution="10m", region=[-180, 180, -80, 80]
    )

    # randomly subsample the regular grid into "scattered" points
    rng = np.random.default_rng(0)
    ny, nx = relief.shape
    n_points = 5000
    iy = rng.integers(0, ny, n_points)
    ix = rng.integers(0, nx, n_points)

    lon = relief.lon.values[ix]
    lat = relief.lat.values[iy]
    val = relief.values[iy, ix]

    # Create a Raster from the points
    raster = gplately.Raster.from_points(
        lon=lon,
        lat=lat,
        values=val,
        spacing="0.5d",  # dataset covers only a small region, so use fine spacing
    )
    # quick visual check
    fig = pygmt.Figure()
    raster.plot(ax_or_fig=fig, use_gmt=True)
    fig.plot(
        x=lon,
        y=lat,
        style="c0.05c",
        fill="black",
        label="Sample points",
    )
    fig.colorbar(frame='af+l"Earth Relief (m)"')
    fig.legend()
    fig.show()


def main(argv):
    pass


if __name__ == "__main__":
    # main(sys.argv)
    test_raster_reconstruct_lon_0_360()
