import math

import numpy as np
import pygplates
import pytest
from conftest import (
    anchor_plate_ids,
    gplately_plate_reconstruction_object,
    gplately_points_lonlat,
    gplately_points_object,
    logger,
    muller_2019_model,
    reconstruction_times,
)

import gplately

logger.info(__name__)
# ========================================= <gplately.Points> =====================================

""" 
A series of automated tests that ensure GPlately's <Points> object is initialised. The object must
process an array of latitudinal and longitudinal coordinates into point data that:
    - are reconstructable through geological time with respect to another geological element or 
    the absolute reference frame, 
    - are ascribed to kinetic topological plates,
    - have point velocity data (the velocities of these plates)

The following methods in the object are tested:

    - __init__
    - reconstruct
    - reconstruct_to_birth_age
    - plate_velocity
        Ensures plate velocities are <insert condition here>. 

"""


# POINTS CREATION AT NON-ZERO TIME AND NON-ZERO ANCHOR PLATE
def test_point_creation_non_zero_time_and_anchor_plate(muller_2019_model):
    # Longitude and latitude of the Hawaiian-Emperor Seamount chain seed points
    # at 50 Ma relative to anchor plate 701 (via 'model_anchored_701').
    rlons = np.array([-121.41333619, -151.95768647])
    rlats = np.array([0.46020662, 36.19123148])

    # Load the model files.
    #
    # Note: We do it once (outside the loop) to avoid converting files to a rotation model or
    #       feature collections in each loop iteration (since that is fairly slow).
    rotation_model = pygplates.RotationModel(muller_2019_model.get_rotation_model())
    topology_features = pygplates.FeatureCollection()
    for f in muller_2019_model.get_topologies():
        topology_features.add(pygplates.FeatureCollection(f))
    static_polygons = pygplates.FeatureCollection()
    for f in muller_2019_model.get_static_polygons():
        static_polygons.add(pygplates.FeatureCollection(f))

    # Try specifying point plate IDs and having them assigned.
    # Also use anchor plate 701, but either via the PlateReconstruction or the Points object (not both).
    for point_plate_ids, model_anchor_plate_id, points_anchor_plate_id in [
        ((901, 901), None, 701),
        ((901, 901), 701, None),
        (None, None, 701),
        (None, 701, None),
    ]:
        model = gplately.PlateReconstruction(
            rotation_model=rotation_model,
            topology_features=topology_features,
            static_polygons=static_polygons,
            anchor_plate_id=model_anchor_plate_id,
        )
        gpts = gplately.Points(
            model,
            lons=rlons,
            lats=rlats,
            time=50,
            plate_id=point_plate_ids,
            anchor_plate_id=points_anchor_plate_id,
        )

        # Note: The geometry (point) in each "feature" is always present-day.
        #       But the lons and lats are locations at the initial 'time' (ie, 50 Ma).
        present_day_lat_lons = np.array(
            [feature.get_geometry().to_lat_lon() for feature in gpts.features]
        )
        present_day_lons = present_day_lat_lons[:, 1]
        present_day_lats = present_day_lat_lons[:, 0]

        # Note: This data was verified using GPlates.
        assert present_day_lons == pytest.approx([-155.4696, 164.3])
        assert present_day_lats == pytest.approx([19.8202, 53.5])

        assert (gpts.plate_id == [901, 901]).all()
        # Anchor plate is 701 regardless of whether set through PlateReconstruction or Points.
        assert gpts.anchor_plate_id == 701
        # And cannot change the anchor plate of a Points object.
        with pytest.raises(AttributeError):
            gpts.anchor_plate_id = 801

        for return_array in (False, True):
            # Reconstruct to 100 Ma.
            #
            # Note: If we don't specify the anchor plate then it will default to 701.
            rlonslats = gpts.reconstruct(100, return_array=return_array)
            if return_array:
                rlons_100, rlats_100 = rlonslats
            else:
                rlons_100, rlats_100 = rlonslats.lons, rlonslats.lats
            # Note: This data was verified using GPlates.
            assert rlons_100 == pytest.approx([-99.22843546, -136.30698742])
            assert rlats_100 == pytest.approx([-28.44815591, 0.20291823])
            #
            # Test anchor plate zero (ie, a non-default anchor plate for 'gpts').
            rlonslats = gpts.reconstruct(
                100, anchor_plate_id=0, return_array=return_array
            )
            if return_array:
                rlons_100, rlats_100 = rlonslats
            else:
                rlons_100, rlats_100 = rlonslats.lons, rlonslats.lats
            # Note: This data was verified using GPlates.
            assert rlons_100 == pytest.approx([-117.22543964, -150.55606233])
            assert rlats_100 == pytest.approx([-3.62350339, 28.85941415])


# POINTS CREATION WITH REMOVAL OF UNRECONSTRUCTABLE POINTS
def test_point_creation_remove_unreconstructable_points(
    gplately_plate_reconstruction_object,
):
    # Longitude and latitude of the Hawaiian-Emperor Seamount chain seed points
    # at 50 Ma relative to anchor plate 701 (via 'model_anchored_701').
    rlons = np.array([-121.41333619, -151.95768647])
    rlats = np.array([0.46020662, 36.19123148])

    gpts = gplately.Points(
        gplately_plate_reconstruction_object,
        lons=rlons,
        lats=rlats,
        time=50,
        plate_id=901,
        age=40,
    )
    # No points removed by default.
    assert gpts.size == 2

    gpts = gplately.Points(
        gplately_plate_reconstruction_object,
        lons=rlons,
        lats=rlats,
        time=50,
        plate_id=901,
        age=40,
        remove_unreconstructable_points=True,
    )
    # Both points removed (since their age is less than 'time').
    assert gpts.size == 0

    gpts = gplately.Points(
        gplately_plate_reconstruction_object,
        lons=rlons,
        lats=rlats,
        time=50,
        plate_id=901,
        age=np.array([40, 100]),
        remove_unreconstructable_points=True,
    )
    # First point removed (since its age is less than 'time').
    assert gpts.size == 1

    unreconstructable_point_indices = []
    gpts = gplately.Points(
        gplately_plate_reconstruction_object,
        lons=rlons,
        lats=rlats,
        time=50,
        plate_id=901,
        age=np.array([40, 100]),
        remove_unreconstructable_points=unreconstructable_point_indices,
    )
    # First point removed (since its age is less than 'time').
    assert gpts.size == 1
    # And index of first point returned in our list.
    assert unreconstructable_point_indices == [0]

    unreconstructable_point_indices = []
    gpts = gplately.Points(
        gplately_plate_reconstruction_object,
        lons=rlons,
        lats=rlats,
        time=50,
        plate_id=None,  # assign plate IDs
        age=None,  # assign ages
        remove_unreconstructable_points=unreconstructable_point_indices,
    )
    # No points removed since both auto-assigned ages (83 and 100 Ma) were greater than (or equal to) 'time' (50 Ma).
    #
    # The fact itself, that both points were auto-assigned, means the associated reconstructed static polygons
    # themselves existed at 'time', which means the polygon ages must have been greater than (or equal to) 'time'.
    assert gpts.size == 2
    assert not unreconstructable_point_indices

    unreconstructable_point_indices = []
    gpts = gplately.Points(
        gplately_plate_reconstruction_object,
        lons=rlons,
        lats=rlats,
        time=50,
        plate_id=None,  # assign plate IDs
        age=np.array([100, 40]),
        remove_unreconstructable_points=unreconstructable_point_indices,
    )
    # Second point removed.
    # Although its plate ID is successfully assigned, its age is explicitly provided but is less than 'time'.
    assert gpts.size == 1
    # And index of second point returned in our list.
    assert unreconstructable_point_indices == [1]


# POINTS RECONSTRUCTION
@pytest.mark.parametrize("time", reconstruction_times)
@pytest.mark.parametrize("anchor_plate_id", anchor_plate_ids)
def test_point_reconstruction(
    time,
    anchor_plate_id,
    gplately_plate_reconstruction_object,
    gplately_points_lonlat,
):
    # Dict of (time, anchor plate ID) to (rlons, rlats).
    #
    # Note: This data was verified using GPlates.
    approx_rlons_rlats_dict = {
        (0, 0): ([-155.4696, 164.3], [19.8202, 53.5]),
        (0, 701): ([-155.4696, 164.3], [19.8202, 53.5]),
        (100, 0): ([-117.22543964, -150.55606233], [-3.62350339, 28.85941415]),
        (100, 701): ([-99.22843546, -136.30698742], [-28.44815591, 0.20291823]),
    }

    # Create points that:
    # 1. exist for all time,
    # 2. exist back to 50 Ma,
    # 3. have ages assigned using the model's static polygons.
    for age in (np.inf, 50, None):

        lons, lats = gplately_points_lonlat
        gpts = gplately.Points(
            gplately_plate_reconstruction_object,
            lons=lons,
            lats=lats,
            age=age,
        )

        approx_rlons, approx_rlats = approx_rlons_rlats_dict[time, anchor_plate_id]
        approx_point_indices = (0, 1)

        if age is None:
            # Ages were assigned. The first point assigned 83 Ma. The second assigned 100 Ma.
            # So, at 100 Ma, only the second point will be in the output.
            if time == 100:
                # Remove first point.
                approx_rlons = np.delete(approx_rlons, 0)
                approx_rlats = np.delete(approx_rlats, 0)
                approx_point_indices = approx_point_indices[1:]
        elif age == 50:
            # Both points appeared at 50 Ma.
            # So, at 100 Ma, neither point will be in the output.
            if time == 100:
                approx_rlons = np.delete(approx_rlons, (0, 1))
                approx_rlats = np.delete(approx_rlats, (0, 1))
                approx_point_indices = ()

        for return_array in (False, True):
            for return_point_indices in (False, True):
                returned_object = gpts.reconstruct(
                    time,
                    anchor_plate_id=anchor_plate_id,
                    return_array=return_array,
                    return_point_indices=return_point_indices,
                )
                if return_array:
                    if return_point_indices:
                        rlons, rlats, point_indices = returned_object
                    else:
                        rlons, rlats = returned_object
                else:
                    if return_point_indices:
                        rpoints, point_indices = returned_object
                        rlons, rlats = rpoints.lons, rpoints.lats
                    else:
                        rpoints = returned_object
                        rlons, rlats = rpoints.lons, rpoints.lats
                assert rlons == pytest.approx(approx_rlons)
                assert rlats == pytest.approx(approx_rlats)
                if return_point_indices:
                    assert (point_indices == approx_point_indices).all()
                # If 'return_array' is False then we have a Points object.
                # Try reconstructing that, in turn, to 'time'.
                # We should get the same reconstructed results.
                if not return_array:
                    rlons, rlats = rpoints.reconstruct(
                        time, anchor_plate_id=anchor_plate_id, return_array=True
                    )
                    assert rlons == pytest.approx(approx_rlons)
                    assert rlats == pytest.approx(approx_rlats)
                    rpoints_recon, rpoint_indices = rpoints.reconstruct(
                        time,
                        anchor_plate_id=anchor_plate_id,
                        return_array=False,
                        return_point_indices=True,
                    )
                    # Note: When 'age=50' and 'time=100' we have no reconstructed points.
                    if age == 50 and time == 100:
                        assert rpoints.size == 0
                    assert rpoints_recon.lons == pytest.approx(approx_rlons)
                    assert rpoints_recon.lats == pytest.approx(approx_rlats)
                    # All points in 'rpoints' should get reconstructed to 'time' (because 'rpoints.time=time').
                    assert (rpoint_indices == np.arange(rpoints.size, dtype=int)).all()


# POINTS RECONSTRUCTION TO BIRTH AGE
@pytest.mark.parametrize("anchor_plate_id", anchor_plate_ids)
def test_point_reconstruction_to_birth_age(
    anchor_plate_id, gplately_plate_reconstruction_object, gplately_points_lonlat
):
    # Dict of anchor plate ID to (rlons, rlats).
    #
    # Note: This data was verified using GPlates.
    approx_rlons_rlats_dict = {
        0: ([-130.74465555, -150.5560233], [9.64809664, 28.85941415]),
        701: ([-121.41333619, -136.30698742], [0.46020662, 0.20291823]),
    }
    approx_rlons, approx_rlats = approx_rlons_rlats_dict[anchor_plate_id]

    lons, lats = gplately_points_lonlat

    # Points exist for all time.
    gpts = gplately.Points(
        gplately_plate_reconstruction_object,
        lons=lons,
        lats=lats,
        age=np.inf,
    )

    rlons, rlats = gpts.reconstruct_to_birth_age(
        ages=[50, 100], anchor_plate_id=anchor_plate_id
    )
    assert rlons == pytest.approx(approx_rlons)
    assert rlats == pytest.approx(approx_rlats)

    # Points assigned ages using static polygons.
    # First point age should be 83 Ma.
    # Second point age should be 100 Ma.
    gpts = gplately.Points(
        gplately_plate_reconstruction_object,
        lons=lons,
        lats=lats,
        age=None,
    )

    # First point exists at 50 Ma.
    # Second point exists at 100 Ma.
    rlons, rlats, point_indices = gpts.reconstruct_to_birth_age(
        ages=[50, 100], anchor_plate_id=anchor_plate_id, return_point_indices=True
    )
    assert rlons == pytest.approx(approx_rlons)
    assert rlats == pytest.approx(approx_rlats)
    assert (point_indices == [0, 1]).all()

    # Change reconstruct age of first point to 100 Ma (from 50).
    # Now it should not get reconstructed (it's appearance age is 83 Ma).
    rlons, rlats, point_indices = gpts.reconstruct_to_birth_age(
        ages=[100, 100], anchor_plate_id=anchor_plate_id, return_point_indices=True
    )
    assert rlons == pytest.approx(approx_rlons[1:])
    assert rlats == pytest.approx(approx_rlats[1:])
    assert (point_indices == [1]).all()

    # Points appear at 90 Ma.
    gpts = gplately.Points(
        gplately_plate_reconstruction_object,
        lons=lons,
        lats=lats,
        age=90,
    )

    # Second point should not get reconstructed (90 < 100).
    rlons, rlats, point_indices = gpts.reconstruct_to_birth_age(
        ages=[50, 100], anchor_plate_id=anchor_plate_id, return_point_indices=True
    )
    assert rlons == pytest.approx(approx_rlons[:1])
    assert rlats == pytest.approx(approx_rlats[:1])
    assert (point_indices == [0]).all()

    # Neither point should get reconstructed (90 < 100).
    rlons, rlats, point_indices = gpts.reconstruct_to_birth_age(
        ages=[100, 100], anchor_plate_id=anchor_plate_id, return_point_indices=True
    )
    assert rlons.size == 0
    assert rlats.size == 0
    assert point_indices.size == 0


# TESTING PLATE VELOCITY CALCULATIONS
@pytest.mark.parametrize("time", reconstruction_times)
def test_plate_velocity(
    time, gplately_plate_reconstruction_object, gplately_points_lonlat
):
    # Dict of time to (vlons, vlats, rlons, rlats).
    #
    # Note: This data was verified using GPlates.
    approx_vlons_vlats_rlons_rlats_dict = {
        0: (
            [-6.45695260, -5.57032570],
            [2.70292492, 1.92700367],
            [-155.4696, 164.3],
            [19.8202, 53.5],
        ),
        100: (
            [-3.31739410, -3.06703484],
            [-1.82746963, -2.29404105],
            [-117.2254, -150.5560],
            [-3.6235, 28.8594],
        ),
    }

    # Create points that:
    # 1. exist for all time,
    # 2. exist back to 50 Ma,
    # 3. have ages assigned using the model's static polygons.
    for age in (np.inf, 50, None):

        lons, lats = gplately_points_lonlat
        gpts = gplately.Points(
            gplately_plate_reconstruction_object,
            lons=lons,
            lats=lats,
            age=age,
        )

        approx_vlons, approx_vlats, approx_rlons, approx_rlats = (
            approx_vlons_vlats_rlons_rlats_dict[time]
        )
        approx_point_indices = (0, 1)

        if age is None:
            # Ages were assigned. The first point assigned 83 Ma. The second assigned 100 Ma.
            # So, at 100 Ma, only the second point will be in the output.
            if time == 100:
                # Remove first point.
                approx_vlons = np.delete(approx_vlons, 0)
                approx_vlats = np.delete(approx_vlats, 0)
                approx_rlons = np.delete(approx_rlons, 0)
                approx_rlats = np.delete(approx_rlats, 0)
                approx_point_indices = approx_point_indices[1:]
        elif age == 50:
            # Both points appeared at 50 Ma.
            # So, at 100 Ma, neither point will be in the output.
            if time == 100:
                approx_vlons = np.delete(approx_vlons, (0, 1))
                approx_vlats = np.delete(approx_vlats, (0, 1))
                approx_rlons = np.delete(approx_rlons, (0, 1))
                approx_rlats = np.delete(approx_rlats, (0, 1))
                approx_point_indices = ()

        vlons, vlats = gpts.plate_velocity(time)
        assert vlons == pytest.approx(approx_vlons)
        assert vlats == pytest.approx(approx_vlats)

        vlons, vlats, rlons, rlats = gpts.plate_velocity(
            time, return_reconstructed_points=True
        )
        assert vlons == pytest.approx(approx_vlons)
        assert vlats == pytest.approx(approx_vlats)
        assert rlons == pytest.approx(approx_rlons)
        assert rlats == pytest.approx(approx_rlats)

        vlons, vlats, point_indices = gpts.plate_velocity(
            time, return_point_indices=True
        )
        assert vlons == pytest.approx(approx_vlons)
        assert vlats == pytest.approx(approx_vlats)
        assert (point_indices == approx_point_indices).all()


def test_point_attributes(gplately_points_object):
    attr = np.arange(0, gplately_points_object.size)
    gplately_points_object.add_attributes(FROMAGE=attr, TOAGE=attr)


def test_pickle_Points(gplately_points_object):
    import pickle

    gpts_dump = pickle.dumps(gplately_points_object)
    gpts_load = pickle.loads(gpts_dump)

    attr = np.arange(0, gplately_points_object.size)
    gplately_points_object.add_attributes(FROMAGE=attr, TOAGE=attr)
    gpts_dump = pickle.dumps(gplately_points_object)
    gpts_load = pickle.loads(gpts_dump)


def test_change_ancbor_plate(gplately_points_object):
    gplately_points_object.rotate_reference_frames(
        50, from_rotation_reference_plate=0, to_rotation_reference_plate=101
    )


def test_reconstruct_points_func():
    data_1 = {
        "points": "95, 54, 142, -33",
        "time": 140,
        "model": "Muller2019",
    }
    data_2 = {"lons": "95, -117.26, 142", "lats": "54, 32.7, -33", "time": 140}
    data_3 = {
        "lons": [95, -117.26, 142],
        "lats": [54, 32.7, -33],
        "times": [140, 100, 50],
        "model": "Muller2019",
    }
    data_4 = {
        "lats": "50, 10, 50",
        "lons": "-100, 160, 100",
        "time": "100",
        "model": "PALEOMAP",
    }
    ret = gplately.reconstruct_points(
        lons=[95, -117.26, 142],
        lats=[54, 32.7, -33],
        model_name="Muller2019",
        times=[140, 100, 50],
    )
    logger.info(ret)
    assert len(ret) == 3
    assert ret[0]["time"] == 140
    assert ret[1]["time"] == 100
    assert ret[2]["time"] == 50

    assert math.isclose(ret[0]["lons"][0], 62.6938, abs_tol=0.001)
    assert math.isnan(ret[0]["lons"][1])
    assert math.isclose(ret[0]["lons"][2], 126.7291, abs_tol=0.001)
    assert math.isclose(ret[1]["lons"][0], 69.9307, abs_tol=0.001)
    assert math.isnan(ret[1]["lons"][1])
    assert math.isclose(ret[1]["lons"][2], 141.4215, abs_tol=0.001)
    assert math.isclose(ret[2]["lons"][0], 84.3967, abs_tol=0.001)
    assert math.isnan(ret[2]["lons"][1])
    assert math.isclose(ret[2]["lons"][2], 137.0091, abs_tol=0.001)

    assert math.isclose(ret[0]["lats"][0], 58.8486, abs_tol=0.001)
    assert math.isnan(ret[0]["lats"][1])
    assert math.isclose(ret[0]["lats"][2], -61.6615, abs_tol=0.001)
    assert math.isclose(ret[1]["lats"][0], 56.6657, abs_tol=0.001)
    assert math.isnan(ret[1]["lats"][1])
    assert math.isclose(ret[1]["lats"][2], -55.9227, abs_tol=0.001)
    assert math.isclose(ret[2]["lats"][0], 56.2327, abs_tol=0.001)
    assert math.isnan(ret[2]["lats"][1])
    assert math.isclose(ret[2]["lats"][2], -55.5198, abs_tol=0.001)

    ret = gplately.reconstruct_points(
        lons=[95, 142],
        lats=[54, -33],
        model_name="Muller2019",
        times=140,
    )

    logger.info(ret)

    ret = gplately.reconstruct_points(
        lons=[62.6938, 126.7291],
        lats=[58.8486, -61.6615],
        model_name="Muller2019",
        times=140,
        reverse=True,
    )
    assert len(ret) == 1
    assert math.isclose(ret[0]["lats"][0], 54, abs_tol=0.001)
    assert math.isclose(ret[0]["lats"][1], -33, abs_tol=0.001)
    assert math.isclose(ret[0]["lons"][0], 95, abs_tol=0.001)
    assert math.isclose(ret[0]["lons"][1], 142, abs_tol=0.001)

    logger.info(ret)

    ret = gplately.reconstruct_points(
        lons=[95, -117.26, 142],
        lats=[54, 32.7, -33],
        model_name="Muller2019",
        times=[140, 100, 50],
        ignore_valid_time=True,
    )

    logger.info(ret)

    ret = gplately.reconstruct_points(
        lons=[95, -117.26, 142],
        lats=[54, 32.7, -33],
        model_name="Muller2019",
        times=[140, 100, 50],
        anchor_plate_id=701,
    )

    logger.info(ret)

    ret = gplately.reconstruct_points(
        lons=[95, -117.26, 142],
        lats=[54, 32.7, -33],
        model_name="Muller2019",
        times=[140, 100, 50],
        pids=401,
    )

    logger.info(ret)

    ret = gplately.reconstruct_points(
        lons=[95, -117.26, 142],
        lats=[54, 32.7, -33],
        model_name="Muller2019",
        times=[140, 100, 50],
        pids=[401, 901, 801],
    )

    logger.info(ret)

    ret = gplately.reconstruct_points(
        lons=[95, -117.26, 142],
        lats=[54, 32.7, -33],
        model_name="Muller2019",
        times=[140, 100, 50],
        pids=[401, 901, 801],
        valid_time=(600, 0),
    )

    logger.info(ret)

    ret = gplately.reconstruct_points(
        lons=[95, -117.26, 142],
        lats=[54, 32.7, -33],
        model_name="Muller2019",
        times=[140, 100, 5],
        pids=[401, 901, 801],
        valid_time=[(600, 0), (10, 0), (100, 0)],
    )

    logger.info(ret)
