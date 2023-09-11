import numpy as np
import gplately
from conftest import (
    reconstruction_times,
    pt_lon,
    pt_lat,
    gplately_raster_object as graster,
    gplately_merdith_raster,
    gplately_merdith_static_geometries,
)

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
        method='linear',
        return_indices=False,
    )
    assert interpolated_points.any(), "Unable to interpolate points"


def test_bilinear_interpolation():
    array = np.array([[0,1],[1,2]], dtype=float)
    graster = gplately.Raster(data=array, extent=[0,1,0,1])

    ilon = ilat = 2.0/3
    result = 1.0 + 1.0/3
    interp = graster.interpolate(ilon, ilat, method="linear")
    assert np.isclose(interp, result), "Linear interpolation in x direction failed"

    # get interpolation coordinates
    interp, (ci,cj) = graster.interpolate(ilon, ilat, method="linear", return_indices=True)
    assert ci == 1 and cj == 1, "Indices of interpolation are incorrect"

def test_nearest_neighbour_interpolation():
    array = np.array([[0,1],[1,2]], dtype=float)
    graster = gplately.Raster(data=array, extent=[0,1,0,1])

    ilon = ilat = 2.0/3
    result = 2
    interp = graster.interpolate(ilon, ilat, method="nearest")
    assert np.isclose(interp, result), "Linear interpolation in x direction failed"

    # get interpolation coordinates
    interp, (ci,cj) = graster.interpolate(ilon, ilat, method="nearest", return_indices=True)
    assert ci == 1 and cj == 1, "Indices of interpolation are incorrect"    


# TEST AGE GRID RESIZING (AT RESOLUTIONS OF RES_X = 1000, RES_Y = 400)
def test_resizing(graster):
    resized_agegrid = graster.resize(1000, 400, return_array=True)
    assert np.shape(resized_agegrid)==(400,1000), "Unable to rezise"


# TEST FILLING NaNs IN AGE GRIDS
def test_fill_NaNs(graster):
    no_NaNs = graster.fill_NaNs(return_array=True)
    assert not np.isnan(no_NaNs).all(), "Unable to fill NaNs"

def test_reconstruct(graster):
    reconstructed_raster = graster.reconstruct(50)
    assert (
        np.shape(reconstructed_raster) == np.shape(graster)
    ), "Unable to reconstruct age grid"


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
    rmse = np.sqrt(np.nanmean(diff ** 2))
    assert rmse < 250.0  # make sure RMSE is within a reasonable limit
