import copy
import os
import warnings

import numpy as np
import pygplates
import pytest
from gplately.utils.longitude_convert import (
    to_longitude_positive_360,
    to_longitude_signed_180,
    upwrap_antimeridian_wraparound,
)
from conftest import (
    gplately_merdith_raster,
    gplately_merdith_static_geometries,
    gplately_plate_reconstruction_object,
)
from conftest import gplately_raster_object as graster
from conftest import logger, muller_2019_model, pt_lat, pt_lon

import gplately

logger.info(__name__)
# ========================================= <gplately.Raster> =========================================

""" 
A series of automated tests that ensure GPlately's <DataServer> object collects the 0 and 100 Ma
age grids from Müller et al. (2019) to initialise the <Raster> object. The following 
methods in the object are tested:

    - __init__ 
    - interpolate
    - resample
    - resize
    - fill_NaNs
    - reconstruct
        For the test, the present-day Müller et al. 2019 age grid is reconstructed to its
        configuration at 50 Ma.
"""


def test_longitude_conversions():
    """
    Test the longitude conversion functions in gplately.utils.longitude_convert module.
    The functions tested are:
        - to_longitude_positive_360
        - to_longitude_signed_180
        - upwrap_antimeridian_wraparound
    The tests cover various scenarios, including global and regional grids, with and without duplicate columns, and crossing the antimeridian.
    """
    # Global, no duplicate
    lons = np.array([-180, -90, 0, 90])
    grid = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
    g, lo = to_longitude_positive_360(grid, lons)
    assert np.allclose(lo, [0, 90, 180, 270])
    assert np.allclose(g, [[3, 4, 1, 2], [7, 8, 5, 6]])

    # Global, duplicate closing column
    lons = np.array([-180, -90, 0, 90, 180])
    grid = np.array([[1, 2, 3, 4, 1], [5, 6, 7, 8, 5]])
    g, lo = to_longitude_positive_360(grid, lons)
    assert np.allclose(lo, [0, 90, 180, 270, 360])
    assert np.allclose(g, [[3, 4, 1, 2, 3], [7, 8, 5, 6, 7]])

    # Regional, on the edge of 0-360, but not crossing the antimeridian
    lons = np.array([-160, -100, -50, 0])
    grid = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
    g, lo = to_longitude_positive_360(grid, lons)
    assert np.allclose(lo, [200, 260, 310, 360])
    assert np.allclose(g, [[1, 2, 3, 4], [5, 6, 7, 8]])

    # Genuine regional crossing -> raise ValueError
    lons = np.array([-10, 0, 10])
    grid = np.array([[1, 2, 3]])
    with pytest.raises(ValueError):
        to_longitude_positive_360(grid, lons)

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
    lons = np.array([170, 175, 180, 185, 190])
    grid = np.array([[1, 2, 3, 4, 5]])
    with pytest.raises(ValueError):
        to_longitude_signed_180(grid, lons)

    # upwrap_antimeridian_wraparound
    lons = [160, 170, 180, -170, -160]
    unwrapped_lons = upwrap_antimeridian_wraparound(lons)
    assert np.allclose(unwrapped_lons, [160, 170, 180, 190, 200])

    lons = np.array([340, 350, 360, 10, 20])
    unwrapped_lons = upwrap_antimeridian_wraparound(lons)
    assert np.allclose(unwrapped_lons, [-20, -10, 0, 10, 20])


# TEST LINEAR POINT DATA INTERPOLATION (WITH PT. COORDS FROM CONFTEST)
def test_point_interpolation(graster):
    interpolated_points = graster.interpolate(
        pt_lon,
        pt_lat,
        method="linear",
        return_indices=False,
    )
    assert interpolated_points.any(), "Unable to interpolate points"


def test_bilinear_interpolation():
    array = np.array([[0, 1], [1, 2]], dtype=float)
    graster = gplately.Raster(data=array, extent=(0, 1, 0, 1))

    ilon = ilat = 2.0 / 3
    result = 1.0 + 1.0 / 3
    interp = graster.interpolate(ilon, ilat, method="linear")
    assert np.isclose(interp, result), "Linear interpolation in x direction failed"

    # get interpolation coordinates
    interp, (ci, cj) = graster.interpolate(
        ilon, ilat, method="linear", return_indices=True
    )
    assert ci == 1 and cj == 1, "Indices of interpolation are incorrect"


def test_nearest_neighbour_interpolation():
    array = np.array([[0, 1], [1, 2]], dtype=float)
    graster = gplately.Raster(data=array, extent=(0, 1, 0, 1))

    ilon = ilat = 2.0 / 3
    result = 2
    interp = graster.interpolate(ilon, ilat, method="nearest")
    assert np.isclose(interp, result), "Linear interpolation in x direction failed"

    # get interpolation coordinates
    interp, (ci, cj) = graster.interpolate(
        ilon, ilat, method="nearest", return_indices=True
    )
    assert ci == 1 and cj == 1, "Indices of interpolation are incorrect"


def test_clip_by_polygons_rgb_uses_channel_aware_default_fill_value():
    data = np.full((3, 3, 3), 255, dtype=np.uint8)
    raster = gplately.Raster(data=data, extent=(0.0, 2.0, 0.0, 2.0))
    polygon = pygplates.PolygonOnSphere(
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]
    )

    clipped = raster.clip_by_polygons([polygon])

    assert np.array_equal(clipped.data[0, 0], [255, 255, 255])
    assert np.array_equal(clipped.data[-1, -1], [0, 0, 0])


def test_clip_by_polygons_preserves_upper_origin():
    data = np.full((3, 3), 255, dtype=np.uint8)
    raster = gplately.Raster(data=data, extent=(0.0, 2.0, 0.0, 2.0), origin="upper")
    polygon = pygplates.PolygonOnSphere(
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]
    )

    clipped = raster.clip_by_polygons([polygon], fill_value=0)

    assert clipped.origin == "upper"
    assert clipped.data[-1, 0] == 255
    assert clipped.data[0, -1] == 0


def test_rotate_reference_frames_ignores_masked_fill_values():
    fill_value = np.float64(9.969209968386869e36)
    data = np.linspace(0, 100, 100).reshape(10, 10).astype(float)
    data[3:6, 3:6] = fill_value
    masked_data = np.ma.masked_values(data, fill_value)

    raster = gplately.Raster(data=masked_data, extent=(-180, 180, -90, 90))
    empty_rotation_model = pygplates.FeatureCollection()
    rotated = raster.rotate_reference_frames(
        grid_spacing_degrees=0.5,
        reconstruction_time=0,
        from_rotation_features_or_model=empty_rotation_model,
        to_rotation_features_or_model=empty_rotation_model,
    )
    rotated_data = np.asarray(rotated.data)

    assert not np.any(np.abs(rotated_data) > 1e20)


def test_rotate_reference_frames_defaults_to_input_resolution():
    data = np.arange(24, dtype=float).reshape(4, 6)
    raster = gplately.Raster(data=data, extent=(10.0, 40.0, -20.0, 20.0))
    empty_rotation_model = pygplates.FeatureCollection()

    rotated = raster.rotate_reference_frames(
        reconstruction_time=0,
        from_rotation_features_or_model=empty_rotation_model,
        to_rotation_features_or_model=empty_rotation_model,
    )

    assert rotated.data.shape == raster.data.shape
    assert np.allclose(rotated.extent, raster.extent)


# TEST AGE GRID RESIZING (AT RESOLUTIONS OF RES_X = 1000, RES_Y = 400)
def test_resizing(graster):
    resized_agegrid = graster.resize(1000, 400, return_array=True)
    assert np.shape(resized_agegrid) == (400, 1000), "Unable to rezise"


# TEST THAT RESAMPLED AXES LAND EXACTLY ON THE REQUESTED EXTENT
@pytest.mark.parametrize("spacing", [0.05, 0.1, 0.2, 0.25, 0.5, 1.0])
def test_resample_axis_endpoints_are_exact(spacing):
    """Resampled axes must not drift past the poles or the dateline.

    Building these axes with ``np.arange`` accumulates floating-point error, so a
    0.2-degree global latitude axis used to end at 90.00000000000256 and trip the
    pole-clipping warning in ``sample_grid``.
    """
    raster = gplately.Raster(
        data=np.zeros((181, 361)), extent=(-180.0, 180.0, -90.0, 90.0)
    )

    with warnings.catch_warnings():
        warnings.filterwarnings("error", message="Invalid values encountered in lat")
        resampled = raster.resample(spacing, spacing)

    assert resampled.lats[0] == -90.0 and resampled.lats[-1] == 90.0
    assert resampled.lons[0] == -180.0 and resampled.lons[-1] == 180.0


# TEST FILLING NaNs IN AGE GRIDS
def test_fill_NaNs(graster):
    no_NaNs = graster.fill_NaNs(return_array=True)
    assert not np.isnan(no_NaNs).all(), "Unable to fill NaNs"


def test_reconstruct(graster):
    graster_downsampled = graster.resample(1.0, 1.0)
    reconstructed_raster = graster_downsampled.reconstruct(50)
    assert np.shape(reconstructed_raster) == np.shape(
        graster_downsampled
    ), "The shape of the reconstructed raster does not match the original raster."

    # test inplace reconstruction with anchor_plate_id=701
    ra = graster_downsampled.copy()
    ra.reconstruct(100, inplace=True, anchor_plate_id=701)
    assert ra.plate_reconstruction is not None
    assert ra.plate_reconstruction.rotation_model
    assert ra.plate_reconstruction.rotation_model.get_default_anchor_plate_id() == 701

    # test return a new Raster obj with anchor_plate_id=701
    ret = graster_downsampled.reconstruct(150, return_array=False, anchor_plate_id=701)
    assert isinstance(ret, gplately.Raster)
    assert ret.plate_reconstruction is not None
    assert ret.plate_reconstruction.rotation_model
    assert ret.plate_reconstruction.rotation_model.get_default_anchor_plate_id() == 701
    assert graster_downsampled.plate_reconstruction is not None
    assert graster_downsampled.plate_reconstruction.rotation_model
    assert (
        graster_downsampled.plate_reconstruction.rotation_model.get_default_anchor_plate_id()
        != 701
    )


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a large file from Internet and takes long time to finish. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_reverse_reconstruct(
    gplately_merdith_raster,
    gplately_merdith_static_geometries,
):
    continents = gplately_merdith_static_geometries[1]
    original_data = np.array(gplately_merdith_raster.data)

    gplately_merdith_raster.reconstruct(
        50,
        partitioning_features=continents,
        inplace=True,
    )
    gplately_merdith_raster.reconstruct(
        0,
        partitioning_features=continents,
        inplace=True,
    )
    diff = gplately_merdith_raster.data - original_data
    # RMS error after reconstructing and reverse reconstructing
    # should be fairly small if reconstruction is working well
    rmse = np.sqrt(np.nanmean(diff**2))
    assert rmse < 250.0  # make sure RMSE is within a reasonable limit


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase is computationally heavy and can be slow. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
def test_reverse_reconstruct_1(
    graster,
    gplately_muller_static_geometries,
):
    continents = gplately_muller_static_geometries[1]
    original_data = np.array(graster.data)

    graster.reconstruct(
        50,
        partitioning_features=continents,
        inplace=True,
    )
    graster.reconstruct(
        0,
        partitioning_features=continents,
        inplace=True,
    )
    diff = graster.data - original_data
    # RMS error after reconstructing and reverse reconstructing
    # should be fairly small if reconstruction is working well
    rmse = np.sqrt(np.nanmean(diff**2))
    logger.info(f"test_reverse_reconstruct_1: the rmse is {rmse}.")
    assert rmse < 250.0  # make sure RMSE is within a reasonable limit


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
def test_copy(graster):
    new = graster.copy()
    assert np.allclose(new.data, graster.data, equal_nan=True)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
def test_new(graster):
    new = gplately.Raster(graster)
    assert np.allclose(new.data, graster.data, equal_nan=True)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
def test_array(graster):
    arr = np.array(graster)
    assert np.allclose(arr, graster.data, equal_nan=True)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
def test_add(graster):
    data = graster.data
    other = gplately.Raster(data=np.ones_like(graster.data))

    assert np.allclose((graster + 1).data, data + 1, equal_nan=True)
    assert np.allclose((1 + graster).data, 1 + data, equal_nan=True)
    assert isinstance(graster + 1, gplately.Raster)
    assert isinstance(1 + graster, gplately.Raster)
    assert np.allclose(graster + other, data + 1, equal_nan=True)
    assert np.allclose(other + graster, 1 + data, equal_nan=True)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
def test_sub(graster):
    data = graster.data
    other = gplately.Raster(data=np.ones_like(graster.data))

    assert np.allclose((graster - 1).data, data - 1, equal_nan=True)
    assert np.allclose((1 - graster).data, 1 - data, equal_nan=True)
    assert isinstance(graster - 1, gplately.Raster)
    assert isinstance(1 - graster, gplately.Raster)
    assert np.allclose(graster - other, data - 1, equal_nan=True)
    assert np.allclose(other - graster, 1 - data, equal_nan=True)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
def test_mul(graster):
    data = graster.data
    other = gplately.Raster(data=np.full_like(graster.data, 2))

    assert np.allclose((graster * 2).data, data * 2, equal_nan=True)
    assert np.allclose((2 * graster).data, 2 * data, equal_nan=True)
    assert isinstance(graster * 2, gplately.Raster)
    assert isinstance(2 * graster, gplately.Raster)
    assert np.allclose(graster * other, data * 2, equal_nan=True)
    assert np.allclose(other * graster, 2 * data, equal_nan=True)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
def test_truediv(graster):
    data = graster.data
    other = gplately.Raster(data=np.full_like(graster.data, 2))

    assert np.allclose((graster / 2).data, data / 2, equal_nan=True)
    assert np.allclose((2 / graster).data, 2 / data, equal_nan=True)
    assert isinstance(graster / 2, gplately.Raster)
    assert isinstance(2 / graster, gplately.Raster)
    assert np.allclose(graster / other, data / 2, equal_nan=True)
    assert np.allclose(other / graster, 2 / data, equal_nan=True)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
def test_floordiv(graster):
    data = graster.data
    other = gplately.Raster(data=np.full_like(graster.data, 2))

    assert np.allclose((graster // 2).data, data // 2, equal_nan=True)
    assert np.allclose((2 // graster).data, 2 // data, equal_nan=True)
    assert isinstance(graster // 2, gplately.Raster)
    assert isinstance(2 // graster, gplately.Raster)
    assert np.allclose(graster // other, data // 2, equal_nan=True)
    assert np.allclose(other // graster, 2 // data, equal_nan=True)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
@pytest.mark.filterwarnings("ignore:invalid value", "ignore:divide by zero")
def test_mod(graster):
    data = graster.data.astype(np.int_)
    graster.data = data

    other = gplately.Raster(data=np.full_like(graster.data, 2))

    assert np.allclose((graster % 2).data, data % 2, equal_nan=True)
    assert np.allclose((2 % graster).data, 2 % data, equal_nan=True)
    assert isinstance(graster % 2, gplately.Raster)
    assert isinstance(2 % graster, gplately.Raster)
    assert np.allclose(graster % other, data % 2, equal_nan=True)
    assert np.allclose(other % graster, 2 % data, equal_nan=True)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="No need to run this testcase everytime. Set GPLATELY_TEST_LEVEL>=1 to activate it.",
)
@pytest.mark.filterwarnings("ignore:invalid value", "ignore:divide by zero")
def test_pow(graster):
    graster = graster.copy()
    data = graster.data
    data = data + np.nanmin(data)
    data = data / np.nanmax(data)
    graster.data = data

    other = gplately.Raster(data=np.full_like(graster.data, 2))

    assert np.allclose((graster**2).data, data**2, equal_nan=True)
    assert np.allclose((2**graster), 2**data, equal_nan=True)
    assert isinstance(graster**2, gplately.Raster)
    assert isinstance(2**graster, gplately.Raster)
    assert np.allclose(graster**other, data**2, equal_nan=True)
    assert np.allclose(other**graster, 2**data, equal_nan=True)
