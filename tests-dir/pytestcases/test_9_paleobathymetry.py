import numpy as np
import pytest
from conftest import logger

import gplately
from gplately.grids.paleobathymetry import (
    AGE_DEPTH_MODELS,
    age_to_basement_depth,
    dutkiewicz_2017_sediment_thickness,
    paleobathymetry,
    sediment_isostatic_correction,
)

logger.info(__name__)


@pytest.mark.parametrize("model", AGE_DEPTH_MODELS)
def test_age_to_basement_depth_shape_and_sign(model):
    age = np.array([np.nan, 0.0, 5.0, 50.0, 150.0])
    depth = age_to_basement_depth(age, model=model)

    assert depth.shape == age.shape
    assert np.isnan(depth[0])
    # basement depth is negative (below sea level) everywhere age is defined
    assert np.all(depth[1:] < 0)
    # older crust is deeper (more negative) than younger crust
    assert np.all(np.diff(depth[1:]) < 0)


def test_age_to_basement_depth_unknown_model():
    with pytest.raises(ValueError):
        age_to_basement_depth(np.array([10.0]), model="not-a-model")


def test_age_to_basement_depth_gdh1_reference_values():
    # Stein & Stein (1992) GDH1, checked against the published formula by hand.
    age = np.array([0.0, 20.0, 100.0])
    depth = age_to_basement_depth(age, model="gdh1")
    expected = np.array(
        [
            -2600.0,
            -(2600.0 + 365.0 * np.sqrt(20.0)),
            -(5651.0 - 2473.0 * np.exp(-0.0278 * 100.0)),
        ]
    )
    np.testing.assert_allclose(depth, expected)


def test_dutkiewicz_sediment_thickness_positive_and_masked():
    age = np.array([np.nan, 0.0, 100.0, 100.0])
    distance_km = np.array(
        [np.nan, 500.0, 500.0, 3500.0]
    )  # last value exceeds the clamp

    thickness = dutkiewicz_2017_sediment_thickness(age, distance_km)

    assert np.isnan(thickness[0])
    assert np.all(thickness[1:] > 0)
    # clamped distance (3500 -> 3000 km) must match evaluating at exactly the clamp.
    clamped = dutkiewicz_2017_sediment_thickness(np.array([100.0]), np.array([3000.0]))
    np.testing.assert_allclose(thickness[3], clamped[0])


def test_sediment_isostatic_correction_nan_propagates_and_nonnegative():
    sediment_thickness_m = np.array([np.nan, 0.0, 500.0, 5000.0])
    correction = sediment_isostatic_correction(sediment_thickness_m)

    assert np.isnan(correction[0])
    assert np.all(correction[1:] >= 0.0)


def test_paleobathymetry_combines_basement_and_sediment():
    basement_depth_m = np.array([-3000.0, -4000.0])
    sediment_thickness_m = np.array([0.0, 500.0])

    result = paleobathymetry(basement_depth_m, sediment_thickness_m)
    correction = sediment_isostatic_correction(sediment_thickness_m)

    np.testing.assert_allclose(
        result, basement_depth_m + sediment_thickness_m - correction
    )
    # zero sediment thickness means paleobathymetry equals basement depth exactly.
    assert result[0] == basement_depth_m[0]


def test_public_api_exports():
    for name in (
        "AGE_DEPTH_MODELS",
        "DUTKIEWICZ_2017_SEDIMENT_THICKNESS",
        "age_to_basement_depth",
        "dutkiewicz_2017_sediment_thickness",
        "sediment_isostatic_correction",
        "paleobathymetry",
    ):
        assert hasattr(gplately, name)
