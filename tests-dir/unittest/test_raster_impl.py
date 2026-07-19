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
from gplately.auxiliary import get_plate_reconstruction
from gplately.utils.log_utils import turn_on_debug_logging

print(gplately.__file__)

# turn_on_debug_logging()

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
age_grid_raster_0 = gplately.Raster(
    data=plate_model.get_raster("AgeGrids", 0),
    plate_reconstruction=get_plate_reconstruction(
        model_name, model_repo_dir=MODEL_REPO_DIR
    ),
)


def test_raster_longitude_conversion():
    np.set_printoptions(suppress=True)
    assert plate_model
    age_grid_raster = gplately.Raster(
        data=plate_model.get_raster("AgeGrids", 0),
        plate_reconstruction=get_plate_reconstruction(
            model_name, model_repo_dir=MODEL_REPO_DIR
        ),
    )
    # print(age_grid_raster.shape)
    # print(age_grid_raster.lons)
    # print(age_grid_raster.lats)

    age_grid_raster_0_360 = age_grid_raster.to_longitude_positive_360()
    # print(age_grid_raster_0_360.lons)
    # print(age_grid_raster_0_360.shape)
    # print(age_grid_raster_0_360.lats)
    age_grid_raster_0_360.save_to_netcdf4("./output/age_grid_raster_0_360.nc")

    age_grid_raster_180 = age_grid_raster_0_360.to_longitude_signed_180()
    # print(age_grid_raster_180.lons)
    # print(age_grid_raster_180.shape)
    # print(age_grid_raster_180.lats)
    age_grid_raster_180.save_to_netcdf4("./output/age_grid_raster_180.nc")

    age_grid_raster_0_360_50ma = age_grid_raster_0_360.reconstruct(50)
    # print(age_grid_raster_0_360_50ma.lons)
    # print(age_grid_raster_0_360_50ma.shape)
    # print(age_grid_raster_0_360_50ma.lats)
    age_grid_raster_0_360_50ma.save_to_netcdf4("./output/age_grid_raster_0_360_50ma.nc")

    age_grid_raster_180_50ma = age_grid_raster_180.reconstruct(50)
    # print(age_grid_raster_180_50ma.lons)
    # print(age_grid_raster_180_50ma.shape)
    # print(age_grid_raster_180_50ma.lats)
    age_grid_raster_180_50ma.save_to_netcdf4("./output/age_grid_raster_180_50ma.nc")


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
    age_grid_raster_0.data = np.flipud(age_grid_raster_0.data)
    age_grid_raster_0.lats = np.flip(age_grid_raster_0.lats)

    age_grid_raster_50_ma = age_grid_raster_0.reconstruct(50)
    age_grid_raster_50_ma.save_to_netcdf4(
        "./output/age_grid_raster_upper_origin_50_ma.nc"
    )


def test_raster_reconstruction():
    anchor_pid = 0
    model = gplately.auxiliary.get_plate_reconstruction(
        "Muller2019", model_repo_dir="plate-model-repo"
    )

    etopo = gplately.Raster(
        data=image.imread(PresentDayRasterManager().get_raster("ETOPO1_tif")),
        plate_reconstruction=model,
    )
    etopo.lats = etopo.lats[::-1]

    etopo_downscaled = etopo.resample(0.5, 0.5)
    assert etopo_downscaled

    fig = plt.figure(figsize=(10, 5), dpi=96)
    # plot 0Ma
    ax_1 = fig.add_subplot(221, projection=ccrs.PlateCarree())
    etopo_downscaled.imshow(ax_1, interpolation="none")
    ax_1.set_title(f"0 Ma")

    # reconstruct to 50 Ma
    r_50 = etopo_downscaled.reconstruct(
        50, fill_value="white", anchor_plate_id=anchor_pid
    )

    # plot 50 Ma
    ax_2 = fig.add_subplot(222, projection=ccrs.PlateCarree())
    r_50.imshow(ax_2, interpolation="none")  # type: ignore
    ax_2.set_title(f"50 Ma - APID:{anchor_pid}")

    # reconstruct to 100 Ma
    r_100 = etopo_downscaled.reconstruct(
        100, fill_value="white", anchor_plate_id=anchor_pid
    )

    # plot 100 Ma
    ax_3 = fig.add_subplot(223, projection=ccrs.PlateCarree())
    r_100.imshow(ax_3, interpolation="none")  # type: ignore
    ax_3.set_title(f"100 Ma - APID:{anchor_pid}")

    # reconstruct to 200 Ma
    r_200 = etopo_downscaled.reconstruct(
        200, fill_value="white", anchor_plate_id=anchor_pid
    )

    # plot 200 Ma
    ax_4 = fig.add_subplot(224, projection=ccrs.PlateCarree())
    r_200.imshow(ax_4, interpolation="none")  # type: ignore
    ax_4.set_title(f"200 Ma - APID:{anchor_pid}")


def main(argv):
    pass


if __name__ == "__main__":
    # main(sys.argv)
    test_raster_longitude_conversion()
