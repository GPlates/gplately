import pytest
import gplately
import numpy as np
import numpy.ma as ma
from conftest import reconstruction_times, pt_lon, pt_lat
from conftest import gplately_plate_reconstruction_object as model
from conftest import gplately_plot_topologies_object as gplot
from conftest import gplately_data_server_object as gdownload

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

# CALL THE RASTER OBJECT
@pytest.mark.parametrize("time", reconstruction_times)
def test_gplately_raster_object(time, model, gdownload):
    masked_age_grid = gdownload.get_age_grid(time)
    graster = gplately.Raster(model, array=masked_age_grid, extent=[-180,180,-90,90])
    assert graster, "Unable to create a <gplately.Raster> object with Müller et al. (2019) at {} Ma.".format(time)


# TEST LINEAR POINT DATA INTERPOLATION (WITH PT. COORDS FROM CONFTEST)
@pytest.mark.parametrize("time", reconstruction_times)
def test_point_interpolation(time, gdownload):
    masked_age_grid = gdownload.get_age_grid(time)
    graster = gplately.Raster(model, array=masked_age_grid, extent=[-180,180,-90,90])
    interpolated_points = graster.interpolate(
        lons=pt_lon, 
        lats=pt_lat, 
        method='linear', 
        return_indices=False, 
        return_distances=False
    )
    assert interpolated_points.any(), "Unable to interpolate points on a Müller et al. (2019) {} Ma age grid.".format(time)


# TEST AGE GRID RESAMPLING (AT X_SPACING = Y_SPACING = 2)
#@pytest.mark.parametrize("time", reconstruction_times)
#def test_resampling(time, graster):
    #resampled_agegrid = graster.resample(2,2, overwrite=False)
    #assert resampled_agegrid.any(), "Unable to resample a Müller et al. (2019) {} Ma age grid to another grid spacing.".format(time)


# TEST AGE GRID RESIZING (AT RESOLUTIONS OF RES_X = 1000, RES_Y = 400)
@pytest.mark.parametrize("time", reconstruction_times)
def test_resizing(time, gdownload):
    masked_age_grid = gdownload.get_age_grid(time)
    graster = gplately.Raster(model, array=masked_age_grid, extent=[-180,180,-90,90])
    resized_agegrid = graster.resize(1000, 400, overwrite=False)
    assert ma.shape(resized_agegrid)==(400,1000), "Unable to rezise a Müller et al. (2019) {} Ma age grid.".format(time)


# TEST FILLING NaNs IN AGE GRIDS
@pytest.mark.parametrize("time", reconstruction_times)
def test_fill_NaNs(time, gdownload):
    masked_age_grid = gdownload.get_age_grid(time)
    graster = gplately.Raster(model, array=masked_age_grid, extent=[-180,180,-90,90])
    no_NaNs = graster.fill_NaNs(overwrite=False)
    assert not np.isnan(no_NaNs).all(), "Unable to fill NaNs for the Müller et al. (2019) {} Ma age grid.".format(time)


