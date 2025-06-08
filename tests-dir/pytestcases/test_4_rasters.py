import os

import numpy as np
import pytest
from conftest import gplately_merdith_raster, gplately_merdith_static_geometries
from conftest import gplately_raster_object as graster
from conftest import logger, pt_lat, pt_lon

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


# TEST AGE GRID RESIZING (AT RESOLUTIONS OF RES_X = 1000, RES_Y = 400)
def test_resizing(graster):
    resized_agegrid = graster.resize(1000, 400, return_array=True)
    assert np.shape(resized_agegrid) == (400, 1000), "Unable to rezise"


# TEST FILLING NaNs IN AGE GRIDS
def test_fill_NaNs(graster):
    no_NaNs = graster.fill_NaNs(return_array=True)
    assert not np.isnan(no_NaNs).all(), "Unable to fill NaNs"


def test_reconstruct(graster):
    reconstructed_raster = graster.reconstruct(50)
    assert np.shape(reconstructed_raster) == np.shape(
        graster
    ), "Unable to reconstruct age grid"


@pytest.mark.skipif(
    int(os.getenv("TEST_LEVEL", 0)) < 1, reason="requires TEST_LEVEL higher than 0"
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


def test_copy(graster):
    new = graster.copy()
    assert np.allclose(new.data, graster.data, equal_nan=True)


def test_new(graster):
    new = gplately.Raster(graster)
    assert np.allclose(new.data, graster.data, equal_nan=True)


def test_array(graster):
    arr = np.array(graster)
    assert np.allclose(arr, graster.data, equal_nan=True)


def test_add(graster):
    data = graster.data
    other = gplately.Raster(data=np.ones_like(graster.data))

    assert np.allclose((graster + 1).data, data + 1, equal_nan=True)
    assert np.allclose((1 + graster).data, 1 + data, equal_nan=True)
    assert isinstance(graster + 1, gplately.Raster)
    assert isinstance(1 + graster, gplately.Raster)
    assert np.allclose(graster + other, data + 1, equal_nan=True)
    assert np.allclose(other + graster, 1 + data, equal_nan=True)


def test_sub(graster):
    data = graster.data
    other = gplately.Raster(data=np.ones_like(graster.data))

    assert np.allclose((graster - 1).data, data - 1, equal_nan=True)
    assert np.allclose((1 - graster).data, 1 - data, equal_nan=True)
    assert isinstance(graster - 1, gplately.Raster)
    assert isinstance(1 - graster, gplately.Raster)
    assert np.allclose(graster - other, data - 1, equal_nan=True)
    assert np.allclose(other - graster, 1 - data, equal_nan=True)


def test_mul(graster):
    data = graster.data
    other = gplately.Raster(data=np.full_like(graster.data, 2))

    assert np.allclose((graster * 2).data, data * 2, equal_nan=True)
    assert np.allclose((2 * graster).data, 2 * data, equal_nan=True)
    assert isinstance(graster * 2, gplately.Raster)
    assert isinstance(2 * graster, gplately.Raster)
    assert np.allclose(graster * other, data * 2, equal_nan=True)
    assert np.allclose(other * graster, 2 * data, equal_nan=True)


def test_truediv(graster):
    data = graster.data
    other = gplately.Raster(data=np.full_like(graster.data, 2))

    assert np.allclose((graster / 2).data, data / 2, equal_nan=True)
    assert np.allclose((2 / graster).data, 2 / data, equal_nan=True)
    assert isinstance(graster / 2, gplately.Raster)
    assert isinstance(2 / graster, gplately.Raster)
    assert np.allclose(graster / other, data / 2, equal_nan=True)
    assert np.allclose(other / graster, 2 / data, equal_nan=True)


def test_floordiv(graster):
    data = graster.data
    other = gplately.Raster(data=np.full_like(graster.data, 2))

    assert np.allclose((graster // 2).data, data // 2, equal_nan=True)
    assert np.allclose((2 // graster).data, 2 // data, equal_nan=True)
    assert isinstance(graster // 2, gplately.Raster)
    assert isinstance(2 // graster, gplately.Raster)
    assert np.allclose(graster // other, data // 2, equal_nan=True)
    assert np.allclose(other // graster, 2 // data, equal_nan=True)


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
