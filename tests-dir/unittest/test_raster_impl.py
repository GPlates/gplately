import os

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
from gplately.mapping.gmt_cpt import get_cmap_from_gmt_cpt

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


def main(argv):
    pass


if __name__ == "__main__":
    # main(sys.argv)
    test_raster_reconstruct_lon_0_360()
