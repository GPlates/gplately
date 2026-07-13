#!/usr/bin/env python3
import os
import sys

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModelManager

import gplately
from gplately.utils.longitude_convert import (
    to_longitude_positive_360,
    to_longitude_signed_180,
    upwrap_antimeridian_wraparound,
)
from gplately.auxiliary import get_plate_reconstruction

print(gplately.__file__)

# gplately.turn_on_debug_logging()


def test_to_longitude_positive_360():

    # Global, no duplicate
    lons = np.array([-180, -90, 0, 90])
    grid = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
    g, lo = to_longitude_positive_360(grid, lons)
    assert np.allclose(lo, [0, 90, 180, 270])
    assert np.allclose(g, [[3, 4, 1, 2], [7, 8, 5, 6]])

    # Global, duplicate closing column
    lons = np.array([-180, -90, 0, 90, 180])  # global
    grid = np.array([[1, 2, 3, 4, 1], [5, 6, 7, 8, 5]])
    g, lo = to_longitude_positive_360(grid, lons)
    assert np.allclose(lo, [0, 90, 180, 270, 360])
    assert np.allclose(g, [[3, 4, 1, 2, 3], [7, 8, 5, 6, 7]])

    # regional, on the edge of 0-360, but not crossing the antimeridian
    lons = np.array([-160, -100, -50, 0])
    grid = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
    g, lo = to_longitude_positive_360(grid, lons)
    assert np.allclose(lo, [200, 260, 310, 360])
    assert np.allclose(g, [[1, 2, 3, 4], [5, 6, 7, 8]])

    try:
        # Genuine regional crossing -> raise ValueError
        lons = np.array([-10, 0, 10])  # regional, crosses antimeridian in 0-360
        grid = np.array([[1, 2, 3]])
        to_longitude_positive_360(grid, lons)
    except ValueError as e:
        print(e)

    # Regional, no crossing
    lons = np.array([180, 190, 210, 360])
    grid = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
    g, lo = to_longitude_signed_180(grid, lons)
    assert np.allclose(lo, [-180, -170, -150, 0])
    assert np.allclose(g, [[1, 2, 3, 4], [5, 6, 7, 8]])

    # Global, 0-360 input, no duplicate
    lons = np.array([0, 90, 180, 270])
    grid = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
    g, lo = to_longitude_signed_180(grid, lons)
    assert np.allclose(lo, [-180, -90, 0, 90])
    assert np.allclose(g, [[3, 4, 1, 2], [7, 8, 5, 6]])

    # Global, 0-360 input, duplicate closing column
    lons = np.array([0, 90, 180, 270, 360])
    grid = np.array([[1, 2, 3, 4, 1], [5, 6, 7, 8, 5]])
    g, lo = to_longitude_signed_180(grid, lons)
    assert np.allclose(lo, [-180, -90, 0, 90, 180])
    assert np.allclose(g, [[3, 4, 1, 2, 3], [7, 8, 5, 6, 7]])

    # Regional crossing -> should error
    try:
        lons = np.array([170, 175, 180, 185, 190])  # i.e. straddles antimeridian
        grid = np.array([[1, 2, 3, 4, 5]])
        to_longitude_signed_180(grid, lons)
        # ValueError: Converting this regional grid...
    except ValueError as e:
        print(e)

    # Test upwrap_antimeridian function
    lons = [160, 170, 180, -170, -160]
    unwrapped_lons = upwrap_antimeridian_wraparound(lons)
    assert np.allclose(unwrapped_lons, [160, 170, 180, 190, 200])

    lons = np.array([340, 350, 360, 10, 20])
    unwrapped_lons = upwrap_antimeridian_wraparound(lons)
    assert np.allclose(unwrapped_lons, [-20, -10, 0, 10, 20])

    print(
        "All tests passed for to_longitude_positive_360, to_longitude_signed_180, and upwrap_antimeridian."
    )


def test_raster_longitude_conversion():
    model_name = "Muller2025"
    plate_model = PlateModelManager().get_model(model_name, data_dir=MODEL_REPO_DIR)

    assert plate_model
    np.set_printoptions(suppress=True)

    age_grid_raster = gplately.Raster(
        data=plate_model.get_raster("AgeGrids", 0),
        plate_reconstruction=get_plate_reconstruction(
            model_name, model_repo_dir=MODEL_REPO_DIR
        ),
    )
    print(age_grid_raster.shape)
    print(age_grid_raster.lons)
    print(age_grid_raster.lats)

    age_grid_raster_0_360 = age_grid_raster.to_longitude_positive_360()
    print(age_grid_raster_0_360.lons)
    print(age_grid_raster_0_360.shape)
    print(age_grid_raster_0_360.lats)
    age_grid_raster_0_360.save_to_netcdf4("./output/age_grid_raster_0_360.nc")

    age_grid_raster_180 = age_grid_raster_0_360.to_longitude_signed_180()
    print(age_grid_raster_180.lons)
    print(age_grid_raster_180.shape)
    print(age_grid_raster_180.lats)
    age_grid_raster_180.save_to_netcdf4("./output/age_grid_raster_180.nc")

    age_grid_raster_0_360_50ma = age_grid_raster_0_360.reconstruct(50)
    print(age_grid_raster_0_360_50ma.lons)
    print(age_grid_raster_0_360_50ma.shape)
    print(age_grid_raster_0_360_50ma.lats)
    age_grid_raster_0_360_50ma.save_to_netcdf4("./output/age_grid_raster_0_360_50ma.nc")

    age_grid_raster_180_50ma = age_grid_raster_180.reconstruct(50)
    print(age_grid_raster_180_50ma.lons)
    print(age_grid_raster_180_50ma.shape)
    print(age_grid_raster_180_50ma.lats)
    age_grid_raster_180_50ma.save_to_netcdf4("./output/age_grid_raster_180_50ma.nc")


def main(argv):
    if len(argv) == 2 and argv[1] == "save":
        show = False
    else:
        show = True

    pm_manager = PlateModelManager()
    plate_model = pm_manager.get_model("Muller2019", data_dir=MODEL_REPO_DIR)

    assert plate_model

    rotation_model = plate_model.get_rotation_model()
    topology_features = plate_model.get_topologies()
    static_polygons = plate_model.get_static_polygons()

    model = gplately.PlateReconstruction(
        rotation_model, topology_features, static_polygons
    )

    age_grid_raster = gplately.Raster(
        data=plate_model.get_raster("AgeGrids", 100),
        plate_reconstruction=model,
        extent=(-180, 180, -90, 90),
    )

    xx, yy = np.meshgrid(np.linspace(-180, 180, 180), np.linspace(-90, 90, 90))

    values = age_grid_raster.query(
        lons=xx.flatten(), lats=yy.flatten(), region_of_interest=8.8
    )

    # print(values[~np.isnan(values)])

    fig = plt.figure(figsize=(10, 5), dpi=96)
    ax_1 = fig.add_subplot(121, projection=ccrs.PlateCarree())

    # plot the values being returned from raster query
    ax_1.scatter(
        xx.flatten(),
        yy.flatten(),
        c=values,
        marker="s",
        s=1,
        transform=ccrs.PlateCarree(),
        cmap="YlGnBu",
        vmax=200,
        vmin=0,
    )

    clipped_raster = age_grid_raster.clip_by_extent([-50, 50, -80, 40])

    # plot the clipped raster
    ax_2 = fig.add_subplot(122, projection=ccrs.PlateCarree())
    clipped_raster.plot(
        ax=ax_2,
        transform=ccrs.PlateCarree(),
        cmap="YlGnBu",
        vmax=200,
        vmin=0,
    )

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    # main(sys.argv)
    # test_to_longitude_positive_360()
    test_raster_longitude_conversion()
