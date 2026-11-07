"""Real CLI/API interface tests on a small synthetic topology (no model downloads).

These compare assembly, routing and serialization with the existing public APIs.
They do not establish independent scientific accuracy or modern/PTT equivalence.
"""

import hashlib
import importlib
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pygplates
import pytest

import gplately

DATA = Path(__file__).resolve().parents[1] / "data" / "subduction_cli"
ROTATIONS = DATA / "rotations.rot"
TOPOLOGIES = DATA / "topologies.gpml"
PTT = importlib.import_module("gplately.ptt.subduction_convergence")
BASE_SCALARS = [
    "ConvergenceVelocityMagnitude",
    "ConvergenceObliquityDegrees",
    "TrenchAbsoluteVelocityMagnitude",
    "TrenchAbsoluteObliquityDegrees",
    "SubductingLengthDegrees",
    "TrenchNormalAzimuthDegrees",
    "SubductingPlateId",
    "TrenchPlateId",
]
OPTIONAL_SCALARS = {
    "output_distance_to_nearest_edge_of_trench": [
        "DistanceToNearestEdgeOfTrenchDegrees"
    ],
    "output_distance_to_start_edge_of_trench": ["DistanceToStartEdgeOfTrenchDegrees"],
    "output_convergence_velocity_components": [
        "ConvergenceVelocityOrthogonal",
        "ConvergenceVelocityParallel",
    ],
    "output_trench_absolute_velocity_components": [
        "TrenchAbsoluteVelocityOrthogonal",
        "TrenchAbsoluteVelocityParallel",
    ],
    "output_subducting_absolute_velocity": [
        "SubductingAbsoluteVelocityMagnitude",
        "SubductingAbsoluteObliquityDegrees",
    ],
    "output_subducting_absolute_velocity_components": [
        "SubductingAbsoluteVelocityOrthogonal",
        "SubductingAbsoluteVelocityParallel",
    ],
    "output_trench_normal": ["TrenchNormalX", "TrenchNormalY", "TrenchNormalZ"],
}
ALL_FLAGS = dict.fromkeys(OPTIONAL_SCALARS, True)
FLAG_CASES = (
    [pytest.param({}, id="base")]
    + [pytest.param({flag: True}, id=flag) for flag in OPTIONAL_SCALARS]
    + [pytest.param(ALL_FLAGS, id="all")]
)


@pytest.fixture(scope="module")
def reconstruction():
    return gplately.PlateReconstruction(str(ROTATIONS), str(TOPOLOGIES))


def run_cli(
    folder,
    method="pygplates",
    alias="sc",
    flags=None,
    times=(0, 1),
    increment=1,
    anchor=0,
    delta=1,
    degrees=5,
    kms=None,
    extra=(),
    gpml=True,
):
    folder.mkdir(parents=True, exist_ok=True)
    arguments = [
        sys.executable,
        "-m",
        "gplately",
        alias,
        str(folder / "result"),
        "-r",
        str(ROTATIONS),
        "-m",
        str(TOPOLOGIES),
        "-t",
        str(times[0]),
        str(times[1]),
        "-i",
        str(increment),
        "-a",
        str(anchor),
        "-v",
        str(delta),
    ]
    if gpml:
        arguments += ["-g", str(folder / "coverage.gpml")]
    if method is not None:
        arguments += ["--method", method]
    if kms is not None:
        arguments += ["-k", str(kms)]
    elif degrees is not None:
        arguments += ["-d", str(degrees)]
    if flags:
        # Output order is defined by the API, independent of CLI flag order.
        arguments += ["-x"] + list(reversed(list(flags)))
    arguments += list(extra)
    environment = dict(os.environ, MPLBACKEND="Agg")
    # Use the same production source as the test process, including when
    # running these regressions against a separate untouched source tree.
    source = str(Path(gplately.__file__).resolve().parent.parent)
    environment["PYTHONPATH"] = source + os.pathsep + environment.get("PYTHONPATH", "")
    result = subprocess.run(
        arguments,
        cwd=folder,
        env=environment,
        capture_output=True,
        text=True,
        timeout=60,
    )
    (folder / "stdout.log").write_text(result.stdout)
    (folder / "stderr.log").write_text(result.stderr)
    return result


def assert_outputs(folder, expected_by_time, flags):
    coverage = list(pygplates.FeatureCollection(str(folder / "coverage.gpml")))
    assert len(coverage) == len(expected_by_time)
    assert sorted(p.name for p in folder.glob("result_*.xy")) == sorted(
        f"result_{time:.2f}.xy" for time in expected_by_time
    )
    scalar_names = BASE_SCALARS + [
        name
        for flag, names in OPTIONAL_SCALARS.items()
        if flags.get(flag)
        for name in names
    ]
    for time, expected in expected_by_time.items():
        assert expected.ndim == 2 and len(expected) > 0
        assert expected.shape[1] == 2 + len(scalar_names)
        assert np.isfinite(expected).all()
        path = folder / f"result_{time:.2f}.xy"
        actual = np.loadtxt(path, ndmin=2)
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
        for line in path.read_text().splitlines():
            fields = line.split()
            assert len(fields) == expected.shape[1]
            assert fields[8:10] == ["201", "101"]
        feature = next(
            item
            for item in coverage
            if item.get_valid_time() == (time + 0.5, time - 0.5)
        )
        geometry, scalars = feature.get_geometry(
            coverage_return=pygplates.CoverageReturn.geometry_and_scalars
        )
        np.testing.assert_allclose(
            geometry.to_lat_lon_array(),
            expected[:, [1, 0]],
            rtol=0,
            atol=1e-12,
        )
        expected_scalars = {
            pygplates.ScalarType.create_gpml(name): expected[:, index + 2]
            for index, name in enumerate(scalar_names)
        }
        assert set(scalars) == set(expected_scalars)
        for scalar, values in expected_scalars.items():
            np.testing.assert_allclose(scalars[scalar], values, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("flags", FLAG_CASES)
def test_modern_cli_schema_and_gpml(tmp_path, reconstruction, flags):
    result = run_cli(tmp_path, flags=flags)
    assert result.returncode == 0, result.stderr
    expected = {
        time: reconstruction.tessellate_subduction_zones(
            time,
            math.radians(5),
            return_geodataframe=False,
            use_ptt=False,
            anchor_plate_id=0,
            velocity_delta_time=1,
            **flags,
        )
        for time in [0, 1]
    }
    assert_outputs(tmp_path, expected, flags)


@pytest.mark.parametrize(
    "alias", ["subduction-convergence", "sc", "subduction_convergence"]
)
def test_modern_cli_aliases(tmp_path, reconstruction, alias):
    result = run_cli(tmp_path, alias=alias, times=(0, 0))
    assert result.returncode == 0, result.stderr
    expected = reconstruction.tessellate_subduction_zones(
        0, math.radians(5), use_ptt=False, return_geodataframe=False
    )
    assert_outputs(tmp_path, {0: expected}, {})


def test_modern_cli_km_anchor_velocity_and_inclusive_times(tmp_path, reconstruction):
    threshold = math.radians(4)
    result = run_cli(
        tmp_path,
        kms=threshold * pygplates.Earth.equatorial_radius_in_kms,
        anchor=101,
        delta=0.5,
        increment=2,
        times=(0, 2),
    )
    assert result.returncode == 0, result.stderr
    expected = {
        time: reconstruction.tessellate_subduction_zones(
            time,
            threshold,
            use_ptt=False,
            return_geodataframe=False,
            anchor_plate_id=101,
            velocity_delta_time=0.5,
        )
        for time in [0, 2]
    }
    assert_outputs(tmp_path, expected, {})


def test_modern_cli_empty_time_has_no_text_file(tmp_path):
    result = run_cli(tmp_path, flags=ALL_FLAGS, times=(20, 20))
    assert result.returncode == 0, result.stderr
    assert not list(tmp_path.glob("result_*.xy"))
    assert len(pygplates.FeatureCollection(str(tmp_path / "coverage.gpml"))) == 0


def test_modern_cli_preserves_default_spacing(tmp_path, reconstruction):
    result = run_cli(tmp_path, degrees=None, times=(0, 0))
    assert result.returncode == 0, result.stderr
    expected = reconstruction.tessellate_subduction_zones(
        0,
        math.radians(PTT.DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES),
        use_ptt=False,
        return_geodataframe=False,
    )
    assert_outputs(tmp_path, {0: expected}, {})


def test_modern_cli_text_without_optional_gpml(tmp_path, reconstruction):
    result = run_cli(tmp_path, times=(0, 0), gpml=False)
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "coverage.gpml").exists()
    expected = reconstruction.tessellate_subduction_zones(0, math.radians(5))
    assert len(expected) > 0
    np.testing.assert_allclose(
        np.loadtxt(tmp_path / "result_0.00.xy", ndmin=2),
        expected,
        rtol=1e-12,
        atol=1e-12,
    )


@pytest.mark.parametrize("method", [None, "ptt"], ids=["default", "explicit-ptt"])
def test_legacy_cli_matches_unchanged_public_ptt(tmp_path, method):
    result = run_cli(tmp_path, method=method, flags=ALL_FLAGS)
    assert result.returncode == 0, result.stderr
    expected = {
        time: np.asarray(
            PTT.subduction_convergence(
                pygplates.RotationModel(str(ROTATIONS)),
                list(pygplates.FeatureCollection(str(TOPOLOGIES))),
                math.radians(5),
                time,
                velocity_delta_time=1,
                anchor_plate_id=0,
                **ALL_FLAGS,
            )
        )
        for time in [0, 1]
    }
    assert_outputs(tmp_path, expected, ALL_FLAGS)


def test_both_methods_preserve_fixture_files(tmp_path):
    before = {
        path: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [ROTATIONS, TOPOLOGIES]
    }
    for method in ["ptt", "pygplates"]:
        result = run_cli(tmp_path / method, method=method, times=(0, 0))
        assert result.returncode == 0, result.stderr
    assert before == {
        path: hashlib.sha256(path.read_bytes()).hexdigest() for path in before
    }


@pytest.mark.parametrize(
    "arguments,message",
    [
        ({"method": "unsupported"}, "invalid choice"),
        ({"times": (2, 1)}, "First (young) value in time range"),
        ({"increment": 0}, "not a positive number"),
        ({"flags": {"unsupported_output": True}}, "not a valid output parameter"),
        ({"degrees": 0}, "Sampling distance must be finite and positive"),
        ({"kms": 0}, "Sampling distance must be finite and positive"),
        ({"times": (0, float("inf"))}, "Time range must be finite"),
        ({"increment": float("inf")}, "Time increment must be finite and positive"),
        ({"degrees": float("inf")}, "Sampling distance must be finite and positive"),
        ({"delta": float("nan")}, "Velocity delta time must be finite and positive"),
        (
            {"times": (1, 1), "increment": 1e-320},
            "Time increment cannot advance reconstruction time",
        ),
    ],
)
def test_cli_rejects_invalid_arguments(tmp_path, arguments, message):
    result = run_cli(tmp_path, **arguments)
    assert result.returncode != 0
    assert message in result.stderr
    assert not list(tmp_path.glob("result_*.xy"))
