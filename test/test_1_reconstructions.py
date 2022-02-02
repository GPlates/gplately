import pytest
import gplately

from .conftest import generic_reconstuction_object as model

def test_gplately_reconstruction_object(model):

    time = 50 #Ma
    subduction_data = model.tesselate_subduction_zones(time, ignore_warnings=True)
    ridge_data = model.tesselate_mid_ocean_ridges(time, ignore_warnings=True)

    assert ridge_data.any(), "There is no data inside ridge_data"


