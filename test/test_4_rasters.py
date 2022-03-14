import pytest
import gplately
import numpy as np
from conftest import reconstruction_times, pt_lon, pt_lat
from conftest import gplately_raster_object as graster

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
    interpolated_points = graster.interpolate(pt_lon, pt_lat, method='linear', return_indices=False, return_distances=False)
    assert interpolated_points.any(), "Unable to interpolate points"


# TEST AGE GRID RESIZING (AT RESOLUTIONS OF RES_X = 1000, RES_Y = 400)
def test_resizing(graster):
    resized_agegrid = graster.resize(1000, 400, overwrite=False)
    assert np.shape(resized_agegrid)==(400,1000), "Unable to rezise"


# TEST FILLING NaNs IN AGE GRIDS
def test_fill_NaNs(graster):
    no_NaNs = graster.fill_NaNs(overwrite=False)
    assert not np.isnan(no_NaNs).all(), "Unable to fill NaNs"

def test_reconstruct(graster):
    reconstructed_raster = graster.reconstruct(50)
    assert np.shape(reconstructed_raster), "Unable to reconstruct age grid"