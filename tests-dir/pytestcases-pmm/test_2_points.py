import pytest
import gplately
import numpy as np
from conftest import reconstruction_times
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
    - plate_velocity
        Ensures plate velocities are <insert condition here>. 

"""

# TESTING THE POINTS OBJECT
@pytest.mark.parametrize("time", reconstruction_times)
def test_point_reconstruction(time, gpts):
    rlons, rlats = gpts.reconstruct(time, return_array=True, anchor_plate_id=0)
    assert (rlons, rlats), "Unable to reconstruct point data to {} Ma with Muller et al. (2019).".format(time)

       
# TESTING PLATE VELOCITY CALCULATIONS
@pytest.mark.parametrize("time", reconstruction_times)
def test_plate_velocity(time, gpts):
    plate_vel = gpts.plate_velocity(time, delta_time=1)
    assert plate_vel, "Unable to calculate plate velocities of point data at {} Ma with Muller et al. (2019).".format(time)

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