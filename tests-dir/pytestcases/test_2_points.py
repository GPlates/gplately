import pytest
import pygplates
import gplately
import numpy as np
from conftest import reconstruction_times, anchor_plate_ids
from conftest import gplately_points_object as gpts


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
def test_point_creation(muller_2019_model):
    # Longitude and latitude of the Hawaiian-Emperor Seamount chain seed points
    # at 50 Ma relative to anchor plate 701 (via 'model_anchored_701').
    rlons = np.array([-121.41333619, -151.95768647])
    rlats = np.array([0.46020662, 36.19123148])

    # Load the model files.
    #
    # Note: We do it once (outside the loop) to avoid converting files to a rotation model or
    #       feature collections in each loop iteration (since that is fairly slow).
    rotation_model = pygplates.RotationModel(muller_2019_model.get_rotation_model())
    topology_features = [
        pygplates.FeatureCollection(file) for file in muller_2019_model.get_topologies()
    ]
    static_polygons = [
        pygplates.FeatureCollection(file)
        for file in muller_2019_model.get_static_polygons()
    ]

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

        # Reconstruct to 100 Ma.
        #
        # Note: If we don't specify the anchor plate then it will default to 701.
        rlons_100, rlats_100 = gpts.reconstruct(100, return_array=True)
        # Note: This data was verified using GPlates.
        assert rlons_100 == pytest.approx([-99.22843546, -136.30698742])
        assert rlats_100 == pytest.approx([-28.44815591, 0.20291823])
        #
        # Test anchor plate zero (ie, a non-default anchor plate for 'gpts').
        rlons_100, rlats_100 = gpts.reconstruct(
            100, anchor_plate_id=0, return_array=True
        )
        # Note: This data was verified using GPlates.
        assert rlons_100 == pytest.approx([-117.22543964, -150.55606233])
        assert rlats_100 == pytest.approx([-3.62350339, 28.85941415])


# POINTS RECONSTRUCTION
@pytest.mark.parametrize("time", reconstruction_times)
@pytest.mark.parametrize("anchor_plate_id", anchor_plate_ids)
def test_point_reconstruction(time, anchor_plate_id, gpts):
    rlons, rlats = gpts.reconstruct(
        time, anchor_plate_id=anchor_plate_id, return_array=True
    )
    # Dict of (time, anchor plate ID) to (rlons, rlats).
    #
    # Note: This data was verified using GPlates.
    approx_rlons_rlats_dict = {
        (0, 0): ([-155.4696, 164.3], [19.8202, 53.5]),
        (0, 701): ([-155.4696, 164.3], [19.8202, 53.5]),
        (100, 0): ([-117.22543964, -150.55606233], [-3.62350339, 28.85941415]),
        (100, 701): ([-99.22843546, -136.30698742], [-28.44815591, 0.20291823]),
    }
    approx_rlons, approx_rlats = approx_rlons_rlats_dict[time, anchor_plate_id]
    assert rlons == pytest.approx(approx_rlons)
    assert rlats == pytest.approx(approx_rlats)


# POINTS RECONSTRUCTION TO BIRTH AGE
@pytest.mark.parametrize("anchor_plate_id", anchor_plate_ids)
def test_point_reconstruction_to_birth_age(anchor_plate_id, gpts):
    rlons, rlats = gpts.reconstruct_to_birth_age(
        ages=[50, 100], anchor_plate_id=anchor_plate_id
    )
    # Dict of anchor plate ID to (rlons, rlats).
    #
    # Note: This data was verified using GPlates.
    approx_rlons_rlats_dict = {
        0: ([-130.74465555, -150.5560233], [9.64809664, 28.85941415]),
        701: ([-121.41333619, -136.30698742], [0.46020662, 0.20291823]),
    }
    approx_rlons, approx_rlats = approx_rlons_rlats_dict[anchor_plate_id]
    assert rlons == pytest.approx(approx_rlons)
    assert rlats == pytest.approx(approx_rlats)


# TESTING PLATE VELOCITY CALCULATIONS
@pytest.mark.parametrize("time", reconstruction_times)
def test_plate_velocity(time, gpts):
    plate_vel = gpts.plate_velocity(time, delta_time=1)
    assert (
        plate_vel
    ), "Unable to calculate plate velocities of point data at {} Ma with Muller et al. (2019).".format(
        time
    )


def test_point_attributes(gpts):
    attr = np.arange(0, gpts.size)
    gpts.add_attributes(FROMAGE=attr, TOAGE=attr)


def test_pickle_Points(gpts):
    import pickle

    gpts_dump = pickle.dumps(gpts)
    gpts_load = pickle.loads(gpts_dump)

    attr = np.arange(0, gpts.size)
    gpts.add_attributes(FROMAGE=attr, TOAGE=attr)
    gpts_dump = pickle.dumps(gpts)
    gpts_load = pickle.loads(gpts_dump)


def test_change_ancbor_plate(gpts):
    gpts.rotate_reference_frames(
        50, from_rotation_reference_plate=0, to_rotation_reference_plate=101
    )
