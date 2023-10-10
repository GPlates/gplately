import pytest
import gplately
import os
import numpy as np
import pandas as pd
import statistics
from conftest import gplately_seafloorgrid_object as seafloorgrid
from conftest import gridding_times

# ========================================= <gplately.SeafloorGrid> =========================================

""" 
A series of automated tests that ensure GPlately's <SeafloorGrid> object can be initialised and all gridding
routines can be performed. The routines must also return sensible grids - for example, spreading rate grids 
must not have an excessive amount of NaN entries (which may mean incorrect masking has occurred).
The following methods in the object are tested:

    - __init__ 
    - reconstruct_by_topologies
    - lon_lat_z_to_netcdf

These tests will be conducted using the Muller et al. 2019 model, with the gridding time range of 250-249 Ma (one
time step).
"""

zval_names = ["SEAFLOOR_AGE", "SPREADING_RATE"]

# CALL THE SEAFLOORGRID OBJECT
def test_gplately_SeafloorGrid_object(
    seafloorgrid
):
    #assert gplot, "No <gplately.PlotTopologies> object made with {}.".format(model)
    assert seafloorgrid, "Unable to create a <gplately.SeafloorGrid> object with \
    Müller et al. (2019) at a max time of {} Ma, and a min time of {} Ma.".format(gridding_times[1], gridding_times[0])

# =======================================================================================================================

# TEST RECONSTRUCTION BY TOPOLOGIES BY TESTING THE CONTENTS OF OUTPUT FILES.
# Gridding input npz files are ok if they have no NaN entries in the spreading rate column.
# Spreading rate ultimately builds the seafloor age. 
@pytest.mark.parametrize("time", gridding_times)
def test_reconstruct_by_topologies(
    time, 
    seafloorgrid
):
    _reconstruct_by_topologies(time, seafloorgrid, clean=True)


# Separate reconstruction by topologies code so that test_lat_lon_z_to_netCDF
# does not depend on the result of test_reconstruct_by_topologies
# This allows the tests to be run in parallel using pytest-xdist
def _reconstruct_by_topologies(time, seafloorgrid, clean=False):
    # This next line must run smoothly. An interruption will happen in pytest
    # if any aspect of reconstruct_by_topologies goes wrong.
    reconstruction_process = seafloorgrid.reconstruct_by_topologies()

    # Otherwise, the following lines then test the validity of written gridding inputs. 

    # Accounts for a given save directory only
    npz_gridding_input = "{:s}/{}_gridding_input_{:0.1f}Ma.npz".format(
            seafloorgrid.save_directory,
            seafloorgrid.file_collection,
            time
        )

    npz_loaded = np.load(npz_gridding_input)

    curr_data = pd.DataFrame.from_dict(

        {item: npz_loaded[item] for item in npz_loaded.files}, 
        orient='columns'
    )
    curr_data.columns = seafloorgrid.total_column_headers 

    unique_data = curr_data.drop_duplicates(

        subset=[
            "CURRENT_LONGITUDES", 
            "CURRENT_LATITUDES"
        ]
    )

    # Gridding input critical data
    age_data = np.array(
        unique_data["SEAFLOOR_AGE"].to_list()
    )
    spreading_rate_data = np.array(
        unique_data["SPREADING_RATE"].to_list()
    )

    # Ensure spreading rate is sensible at max_time; namely that 
    # it is only the initial spreading rate
    if time == np.sort(gridding_times)[-1]:
        initial_sr = np.unique(spreading_rate_data)
        assert initial_sr == seafloorgrid.initial_ocean_mean_spreading_rate, "Initial grids have not been allocated the correct initial mean ocean spreading rate."

    # If it is not the max time, but still near the start of the reconstruction tree
    # (in this test we use the first time step after max_time),
    # ensure the initial spreading rate is still the mode in the array
    elif time == np.sort(gridding_times)[-2]:
        most_common_spreading_rate = statistics.mode(spreading_rate_data)
        assert most_common_spreading_rate == seafloorgrid.initial_ocean_mean_spreading_rate, "Initial grids have not been allocated the correct initial mean ocean spreading rate."

    # Otherwise, make sure the array is not empty, and has no NaNs.
    else:
        assert np.unique(np.isnan(spreading_rate_data))[0] is False, "Some spreading rates in the {} Ma gridding input have been ascribed a NaN value.".format(time)

    if clean:
        for i in ("", "unmasked_"):
            for zval_name in zval_names:
                filename = os.path.join(
                    seafloorgrid.save_directory,
                    "{}_{}_grid_{}{}Ma.nc".format(
                        seafloorgrid.file_collection,
                        zval_name,
                        i,
                        time,
                    )
                )
                if os.path.exists(filename):
                    os.remove(filename)


# test netCDF writing
@pytest.mark.parametrize("zval_name", zval_names)
def test_lat_lon_z_to_netCDF(
    zval_name, 
    seafloorgrid
):
    time = gridding_times[0]

    # Test the creation of a masked and unmasked age grid
    _reconstruct_by_topologies(time, seafloorgrid)
    seafloorgrid.lat_lon_z_to_netCDF(
        zval_name, 
        unmasked=True
    )

    grid_output_unmasked = "{}/{}_{}_grid_unmasked_{}Ma.nc".format(
        seafloorgrid.save_directory,
        seafloorgrid.file_collection,
        zval_name,
        time
    )
    
    grid_output_dir = "{}/{}_{}_grid_{}Ma.nc".format(
        seafloorgrid.save_directory,
        seafloorgrid.file_collection,
        zval_name,
        time
    )

    age_grid_unmasked = gplately.Raster(
        data=grid_output_unmasked,
        extent=[-180,180,-90,90]
    )

    age_grid = gplately.Raster(
        data=grid_output_dir,
        extent=[-180,180,-90,90]
    )

    assert age_grid_unmasked.data.any(), "Could not produce an unmasked {} grid.".format(zval_name)
    assert age_grid.data.any(), "Could not produce a masked {} grid.".format(zval_name)







