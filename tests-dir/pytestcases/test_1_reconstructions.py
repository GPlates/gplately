import numpy as np
import pytest
from conftest import gplately_plate_reconstruction_object as model
from conftest import logger, muller_2019_model, reconstruction_times
from pygplates import FeatureCollection, RotationModel

import gplately

logger.info(__name__)


# ========================================= <gplately.PlateReconstruction> =====================================

""" 
A series of automated tests that ensure GPlately's <DataServer> object collects the necessary
rotation files, topology features and static polygons to initialise the <PlateReconstruction> 
object with the Müller et al. (2019) plate reconstruction model. The following methods in the 
object are tested:

    - __init__
    - tessellate_subduction_zones
    - tessellate_mid_ocean_ridges

    Using pyGPlates and Plate Tectonic Tools:
    - total_subduction_zone_length
    - total_ridge_length
    - total_continental_arc_length

    - get_point_velocities
"""

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


@pytest.mark.parametrize("time", reconstruction_times)
def test_tessellate_trenches(time, model):
    subduction_data = model.tessellate_subduction_zones(time, ignore_warnings=True)
    # CONDITIONS
    assert (
        subduction_data.any()
    ), "There is no trench data inside Muller et al. (2019) at {} Ma.".format(time)
    assert (
        np.abs(subduction_data[:, 0]) <= 180
    ).all(), "Some trench lons exceed a magnitude of 180 degrees in Muller et al. (2019) at {} Ma.".format(
        time
    )
    assert (
        np.abs(subduction_data[:, 1]) <= 90
    ).all(), "Some trench lats exceed a magnitude of 90 degrees in Muller et al. (2019) at {} Ma.".format(
        time
    )
    assert (
        np.abs(subduction_data[:, 2]) <= 50
    ).all(), "Some trench convergence velocities exceed 50 cm/yr in Muller et al. (2019) at {} Ma.".format(
        time
    )


@pytest.mark.parametrize("time", reconstruction_times)
def test_tessellate_ridges(time, model):
    ridge_data = model.tessellate_mid_ocean_ridges(time, ignore_warnings=False)
    assert (
        ridge_data.any()
    ), "There is no ridge data inside Müller et al. (2019) at {} Ma.".format(time)
    assert (
        np.abs(ridge_data[:, 0]) <= 180
    ).all(), "Some ridge lons exceed a magnitude of 180 degrees in Muller et al. (2019) at {} Ma.".format(
        time
    )
    assert (
        np.abs(ridge_data[:, 1]) <= 90
    ).all(), "Some ridge lats exceed a magnitude of 90 degrees in Muller et al. (2019) at {} Ma.".format(
        time
    )
    assert (
        np.abs(ridge_data[:, 2]) <= 50
    ).all(), "Some ridge convergence velocities exceed 50 cm/yr in Muller et al. (2019) at {} Ma.".format(
        time
    )


# TOTAL TRENCH AND RIDGE LENGTHS: PYGPLATES
@pytest.mark.parametrize("time", reconstruction_times)
def test_pygplates_trench_length(time, model):
    total_sz_length = model.total_subduction_zone_length(time, use_ptt=False)
    assert (
        50000 <= total_sz_length <= 100000
    ), "Could not calculate total SZ lengths for Muller et al. (2019) at {} Ma.".format(
        time
    )
    total_sz_length_using_ptt = model.total_subduction_zone_length(
        time, use_ptt=True, ignore_warnings=True
    )
    assert (
        50000 <= total_sz_length_using_ptt <= 100000
    ), "Could not calculate total SZ lengths for Muller et al. (2019) at {} Ma.".format(
        time
    )
    # There are some differences in subduction zone length with and without using PTT.
    diffs = {
        # PTT misses some subduction zones at 0 Ma. Eg, one is duplicated, and mislabelled as a Transform
        # (the duplicate that's attached to the subducting plate)...
        0: 4000,
        100: 200,
    }
    assert np.abs(total_sz_length - total_sz_length_using_ptt) < diffs.get(
        time,
        # Default diff when value for this time has not been computed...
        1000.0,
    ), "Could not calculate total SZ lengths for Muller et al. (2019) at {} Ma.".format(
        time
    )


@pytest.mark.parametrize("time", reconstruction_times)
def test_pygplates_ridge_length(time, model):
    total_ridge_length = model.total_ridge_length(time, use_ptt=False)
    assert (
        total_ridge_length
    ), "Could not calculate total MOR lengths for Muller et al. (2019) at {} Ma.".format(
        time
    )
    total_ridge_length_using_ptt = model.total_ridge_length(
        time, use_ptt=True, ignore_warnings=True
    )
    # There are some differences in ridge length with and without using PTT since
    # PTT uses stage rotations for spreading whereas not using PTT uses divergence velocities.
    assert (
        np.abs(total_ridge_length - total_ridge_length_using_ptt) < 3000.0
    ), "Could not calculate total MOR lengths for Muller et al. (2019) at {} Ma.".format(
        time
    )


# TOTAL CONTINENTAL ARC LENGTHS AT 0 AND 100 MA
@pytest.mark.parametrize("time", reconstruction_times)
def test_cont_arc_length(time, model, muller_2019_model):
    agegrid_path = muller_2019_model.get_raster("AgeGrids", time)
    agegrid = gplately.Raster(
        data=agegrid_path,
        plate_reconstruction=model,
        extent=(-180, 180, -90, 90),
    )

    continental_grid = np.isnan(np.array(agegrid))
    # Use 281km as the trench-arc distance 281 km; the median distance frin global analysis at the present-day (Pall et al. 2018)
    trench_arc_distance = 281
    total_cont_arc_length = model.total_continental_arc_length(
        time,
        continental_grid,
        trench_arc_distance,
        ignore_warnings=True,
    )
    approx_lengths = {
        0: 52200,
        100: 46000,
    }
    err_msg = (
        "Could not calculate total continental arc lengths for Muller et "
        + "al. (2019) at {} Ma.".format(time)
    )
    if time in approx_lengths.keys():
        # Check arc length approximately matches precomputed value
        cmp = approx_lengths[time]
        diff = np.abs(total_cont_arc_length - cmp)
        assert diff <= 1000.0, err_msg
    else:
        # Value for this time has not been computed, so just make sure no
        # errors were encountered
        assert total_cont_arc_length >= 0.0, err_msg


# POINT VELOCITIES
@pytest.mark.parametrize("time", reconstruction_times)
def test_point_velocities(time, model):
    lons = [-30]
    lats = [15]
    point_velocities = model.get_point_velocities(lons, lats, time, delta_time=1.0)
    assert (
        point_velocities.any()
    ), "Could not calculate point data velocities for Muller et al. (2019) at {} Ma.".format(
        time
    )


def test_rotation_model_copy(model):
    rot_model = model.rotation_model
    rot_model_copy = RotationModel(rot_model)
    rot1 = rot_model.get_rotation(50, 101)
    rot2 = rot_model_copy.get_rotation(50, 101)
    assert rot1 == rot2


def test_pickle_reconstruction(model, muller_2019_model):
    import pickle

    # Also create a model from actual pygplates objects (instead of filenames).
    pygplates_rotation_model = RotationModel(muller_2019_model.get_rotation_model())
    pygplates_topology_features = [
        FeatureCollection(f) for f in muller_2019_model.get_topologies()
    ]
    pygplates_static_polygons = [
        FeatureCollection(f) for f in muller_2019_model.get_static_polygons()
    ]
    pygplates_model = gplately.PlateReconstruction(
        rotation_model=pygplates_rotation_model,
        topology_features=pygplates_topology_features,
        static_polygons=pygplates_static_polygons,
    )

    # Test both the model created using PlateModelManager (ie, using filenames) and
    # the model created using actual pygplates objects (ie, not filenames).
    #
    # Pickling of the former will be faster than the latter.
    for m in (model, pygplates_model):
        pickled_model = pickle.loads(pickle.dumps(m))

        time = 100.0
        plate_id = 701

        # Test pickled rotation model, topology features and static polygons.
        # Since these are handled specially when pickling.
        assert pickled_model.rotation_model.get_rotation(
            time, plate_id
        ) == m.rotation_model.get_rotation(time, plate_id)
        assert pickled_model.topology_features and len(m.topology_features) == len(
            pickled_model.topology_features
        )
        assert pickled_model.static_polygons and len(m.static_polygons) == len(
            pickled_model.static_polygons
        )

        # The snapshots in the pickled model get rebuilt (rather than pickled).
        pickled_topological_snapshot = pickled_model.topological_snapshot(
            time, anchor_plate_id=plate_id
        )
        assert pickled_topological_snapshot.get_reconstruction_time() == time
        assert pickled_topological_snapshot.get_anchor_plate_id() == plate_id
        pickled_static_polygons_snapshot = pickled_model.static_polygons_snapshot(
            time, anchor_plate_id=plate_id
        )
        assert pickled_static_polygons_snapshot.get_reconstruction_time() == time
        assert pickled_static_polygons_snapshot.get_anchor_plate_id() == plate_id


def test_auxiliary_get_plate_reconstruction():
    from gplately import auxiliary

    m = auxiliary.get_plate_reconstruction("Muller2022")
    assert isinstance(m, gplately.PlateReconstruction)
