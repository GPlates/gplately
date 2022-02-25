import pytest
import gplately
import numpy as np
from conftest import gplately_plate_reconstruction_object as model
from conftest import reconstruction_times
from conftest import gplately_data_server_object as gdownload

# ========================================= <gplately.PlateReconstruction> =====================================

""" 
A series of automated tests that ensure GPlately's <DataServer> object collects the necessary
rotation files, topology features and static polygons to initialise the <PlateReconstruction> 
object with the Müller et al. (2019) plate reconstruction model. The following methods in the 
object are tested:

    - __init__
    - tesselate_subduction_zones
    - tesselate_mid_ocean_ridges

    Using pyGPlates and Plate Tectonic Tools:
    - total_subduction_zone_length
    - total_ridge_length
    - total_continental_arc_length

    - get_point_velocities
"""

# CALLING THE PLATE RECONSTRUCTION OBJECT
def test_gplately_reconstruction_object(model):
    #assert model, "No <gplately.PlateReconstruction> object made with {}.".format(model)
    assert model, "Unable to create a <gplately.PlateReconstruction> object with Müller et al. (2019)."


# TESSELLATION OF TRENCHES AND RIDGES
"""
Contents are:
    0 - longitude of sampled point
    1 - latitude of sampled point
    2 - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
    3 - subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
    4 - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
    5 - trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
    6 - length of arc segment (in degrees) that current point is on
    7 - trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
    8 - subducting plate ID
    9 - trench plate ID
"""
@pytest.mark.parametrize(
    "time", 
    reconstruction_times
)
def test_tessellate_trenches(time, model):
    subduction_data = model.tesselate_subduction_zones(
        time, 
        ignore_warnings=False
    )
    # CONDITIONS
    assert subduction_data.any(), "There is no trench data inside Müller et al. (2019) at {} Ma.".format(time)
    assert (np.abs(subduction_data[:,0]) <= 180).all(), "Some trench lons exceed a magnitude of 180 degrees in Müller et al. (2019) at {} Ma.".format(time)
    assert (np.abs(subduction_data[:,1]) <= 90).all(), "Some trench lats exceed a magnitude of 90 degrees in Müller et al. (2019) at {} Ma.".format(time)
    assert (np.abs(subduction_data[:,2]) <= 50).all(), "Some trench convergence velocities exceed 20 cm/yr in Müller et al. (2019) at {} Ma.".format(time)
    # condition for subd convergence obliquity angle
    # condition for trench absolute velocity magnitude
    # condition for trench absolute velocity obliquity angle
    # condition for trench normal azimuth angle

@pytest.mark.parametrize(
    "time", 
    reconstruction_times
)
def test_tessellate_ridges(time, model):
    ridge_data = model.tesselate_mid_ocean_ridges(
        time, 
        ignore_warnings=False
    )
    assert ridge_data.any(), "There is no trench data inside Müller et al. (2019) at {} Ma.".format(time)


# TOTAL TRENCH AND RIDGE LENGTHS: PYGPLATES
@pytest.mark.parametrize(
    "time", 
    reconstruction_times
)
def test_pygplates_trench_length(time, model):
    total_sz_length = model.total_subduction_zone_length(
        time, 
        use_pygplates=True, 
        ignore_warnings=False
    )
    assert total_sz_length, "Could not calculate total SZ lengths via pyGPlates for Müller et al. (2019) at {} Ma.".format(time)


@pytest.mark.parametrize(
    "time", 
    reconstruction_times
)
def test_pygplates_ridge_length(time, model):
    total_ridge_length = model.total_ridge_length(
        time, 
        use_pygplates=True, 
        ignore_warnings=False
    )
    assert total_ridge_length, "Could not calculate total MOR lengths via pyGPlates for Müller et al. (2019) at {} Ma.".format(time)


# TOTAL TRENCH AND RIDGE LENGTHS: PTT
@pytest.mark.parametrize(
    "time", 
    reconstruction_times
)
def test_pygplates_trench_length(time, model):
    total_sz_length = model.total_subduction_zone_length(
        time, 
        use_ptt=True, 
        ignore_warnings=False
    )
    # Subduction zone length acceptable range adapted from Pall et al. (2018) https://doi.org/10.5194/cp-14-857-2018
    assert 50000 <= total_sz_length <= 100000, "Could not calculate total SZ lengths via PTT for Müller et al. (2019) at {} Ma.".format(time)


@pytest.mark.parametrize(
    "time", 
    reconstruction_times
)
def test_pygplates_ridge_length(time, model):
    total_ridge_length = model.total_ridge_length(
        time, 
        use_ptt=True, 
        ignore_warnings=False
    )
    assert total_ridge_length, "Could not calculate total MOR lengths via PTT for Müller et al. (2019) at {} Ma.".format(time)


# TOTAL CONTINENTAL ARC LENGTHS AT 0 AND 100 MA
@pytest.mark.parametrize(
    "time", 
    reconstruction_times
)
def test_cont_arc_length(time, model):
    gdownload = gplately.download.DataServer("Muller2019")
    continental_grid_directory = gdownload.get_age_grid(time)
    # Use 281km as the trench-arc distance 281 km; the median distance frin global analysis at the present-day (Pall et al. 2018)
    trench_arc_distance = 281
    total_cont_arc_length = model.total_continental_arc_length(
        time, 
        continental_grid_directory, 
        trench_arc_distance, 
        ignore_warnings=False
    )
    assert total_cont_arc_length, "Could not calculate total continental arc lengths for Müller et al. (2019) at {} Ma.".format(time)


# POINT VELOCITIES
@pytest.mark.parametrize(
    "time", 
    reconstruction_times
)
def test_point_velocities(time, model):
    lons = [-30]
    lats = [15]
    point_velocities = model.get_point_velocities(
        lons,
        lats, 
        time,
        delta_time=1.0
    )
    assert point_velocities.any(), "Could not calculate point data velocities for Müller et al. (2019) at {} Ma.".format(time)


def download_multifeature_gplates_sample_data(gdownload):
    johansson_2018 = gdownload.get_feature_data("Johannson2018")
    assert isinstance(johansson_2018, pygplates.FeatureCollection), "Cached GPlates 2.3 sample data not an instance of <pygplates.FeatureCollection>."


def download_multifeature_gplates_sample_data(gdownload):
    seafloor_fabric = gdownload.get_feature_data("SeafloorFabric")
    assert [isinstance(sf, pygplates.FeatureCollection) for sf in seafloor_fabric], "Cached GPlates 2.3 sample data not an instance of <pygplates.FeatureCollection>."

