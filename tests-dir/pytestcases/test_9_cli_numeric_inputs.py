# SPDX-License-Identifier: GPL-2.0-only
"""Numeric CLI guards, with small original synthetic native inputs."""

from pathlib import Path
import sys

import numpy as np
import pygplates
import pytest

from gplately.__main__ import main
from gplately.commands import _plate_model_args as model_args
from gplately.commands import subduction_convergence as command
from gplately.reconstruction import PlateReconstruction

DATA = Path(__file__).resolve().parents[1] / "data" / "cli_numeric"


def _run(
    monkeypatch,
    tmp_path,
    *extra,
    subcommand="subduction-convergence",
    gpml=False,
    model=False,
):
    args = ["gplately", subcommand]
    if model:
        args += ["--model", "synthetic-fixture"]
    else:
        args += ["-r", str(DATA / "rotations.rot"), "-m", str(DATA / "topologies.gpml")]
    args += ["-t", "0", "1", "-i", "1", "-v", "1", "-w"]
    if "-d" not in extra and "-k" not in extra:
        args += ["-d", "5"]
    args += ["-x", "output_trench_normal", "output_distance_to_nearest_edge_of_trench"]
    if gpml:
        args += ["-g", str(tmp_path / "coverage.gpml")]
    args += [*extra, "--", str(tmp_path / "convergence")]
    monkeypatch.setattr(sys, "argv", args)
    main()


def _forbid_reconstruction(*args, **kwargs):
    pytest.fail("Invalid numeric input reached reconstruction")


class _ObservedNativeLimit(Exception):
    """A deliberately interrupted native prefix, never full CLI completion."""


def _observe_native(monkeypatch, limit=3):
    original = PlateReconstruction.tessellate_subduction_zones
    calls = []

    def observed(self, time, *args, **kwargs):
        # Delegate every call to the real implementation and native pyGPlates.
        result = original(self, time, *args, **kwargs)
        calls.append((time, kwargs))
        if len(calls) >= limit:
            raise _ObservedNativeLimit("Completed the bounded real native prefix")
        return result

    monkeypatch.setattr(PlateReconstruction, "tessellate_subduction_zones", observed)
    return calls


@pytest.mark.parametrize(
    "extra",
    [
        ["-t", "nan", "1"],
        ["-t", "0", "inf"],
        ["-i", "nan"],
        ["-i", "inf"],
        ["-v", "nan"],
        ["-d", "nan"],
        ["-k", "nan"],
        ["--convergence-threshold=nan"],
        ["--convergence-threshold=-inf"],
    ],
)
def test_nonfinite_numbers_fail_before_reconstruction(monkeypatch, tmp_path, extra):
    monkeypatch.setattr(command, "PlateReconstruction", _forbid_reconstruction)
    with pytest.raises(ValueError, match="finite"):
        _run(monkeypatch, tmp_path, *extra)
    assert not list(tmp_path.iterdir())


def test_initial_stall_fails_before_reconstruction(monkeypatch, tmp_path):
    monkeypatch.setattr(command, "PlateReconstruction", _forbid_reconstruction)
    with pytest.raises(ValueError, match="advance"):
        _run(monkeypatch, tmp_path, "-t", "10", "11", "-i", "1e-20")
    assert not list(tmp_path.iterdir())


def test_later_stall_fails_before_repeating_a_native_time(monkeypatch, tmp_path):
    calls = _observe_native(monkeypatch)
    with pytest.raises(ValueError, match="advance.*earlier times"):
        _run(
            monkeypatch,
            tmp_path,
            "-t",
            "0.9999999999999999",
            "2",
            "-i",
            "1.1102230246251565e-16",
        )
    assert [time for time, _ in calls] == [0.9999999999999999, 1.0]
    # Earlier XY output is retained; the failure does not promise rollback.
    assert (tmp_path / "convergence_1.00.xy").is_file()


@pytest.mark.parametrize("young,old", [(10.0, 10.0), (1.0, 1.0000000000000002)])
def test_failed_increment_at_native_endpoint_finishes_once(
    monkeypatch, tmp_path, young, old
):
    calls = _observe_native(monkeypatch)
    _run(monkeypatch, tmp_path, "-t", str(young), str(old), "-i", "1e-20")
    assert [time for time, _ in calls] == [young]
    assert len(list(tmp_path.glob("*.xy"))) == 1


def test_native_endpoint_equivalence_does_not_stop_advancing_samples(
    monkeypatch, tmp_path
):
    calls = _observe_native(monkeypatch, limit=2)
    # Native tolerance admits both samples (and can admit more). Observe only
    # this bounded prefix; equality at the first sample must not stop the CLI.
    with pytest.raises(_ObservedNativeLimit):
        _run(
            monkeypatch,
            tmp_path,
            "-t",
            "1",
            "1.0000000000000002",
            "-i",
            "2.220446049250313e-16",
        )
    assert [time for time, _ in calls] == [1.0, 1.0000000000000002]


@pytest.mark.parametrize(
    "subcommand,use_ptt",
    [
        ("subduction-convergence", False),
        ("sc", False),
        ("subduction_convergence", True),
    ],
)
def test_native_outputs_aliases_and_current_coverage_times(
    monkeypatch, tmp_path, subcommand, use_ptt
):
    calls = _observe_native(monkeypatch)
    extra = ["--use-ptt"] if use_ptt else []
    _run(monkeypatch, tmp_path, *extra, subcommand=subcommand, gpml=True)
    assert [time for time, _ in calls] == [0.0, 1.0]
    for time in (0, 1):
        path = tmp_path / f"convergence_{time:0.2f}.xy"
        rows = np.atleast_2d(np.loadtxt(path))
        assert len(rows) > 0
        assert rows.shape[1] == 14
        assert np.isfinite(rows).all()
        for line in path.read_text().splitlines():
            assert all(value == str(int(value)) for value in line.split()[8:10])
    features = pygplates.FeatureCollection(str(tmp_path / "coverage.gpml"))
    assert {feature.get_valid_time() for feature in features} == {
        (0.5, -0.5),
        (1.5, 0.5),
    }


def test_decimal_terminal_uses_native_time_comparison(monkeypatch, tmp_path):
    calls = _observe_native(monkeypatch, limit=6)
    _run(monkeypatch, tmp_path, "-t", "0", "0.3", "-i", "0.1", gpml=True)
    assert [time for time, _ in calls] == [0.0, 0.1, 0.2, 0.30000000000000004]
    assert len(list(tmp_path.glob("*.xy"))) == 4
    assert len(pygplates.FeatureCollection(str(tmp_path / "coverage.gpml"))) == 4


def test_negative_finite_convergence_threshold_remains_supported(monkeypatch, tmp_path):
    calls = _observe_native(monkeypatch)
    _run(monkeypatch, tmp_path, "--convergence-threshold=-1")
    assert [time for time, _ in calls] == [0.0, 1.0]
    assert [kwargs["convergence_threshold_in_cm_per_yr"] for _, kwargs in calls] == [
        -1.0,
        -1.0,
    ]


def test_normal_single_point_keeps_current_coverage_interval(monkeypatch, tmp_path):
    calls = _observe_native(monkeypatch)
    _run(monkeypatch, tmp_path, "-t", "10", "10", gpml=True)
    assert [time for time, _ in calls] == [10.0]
    features = pygplates.FeatureCollection(str(tmp_path / "coverage.gpml"))
    assert [feature.get_valid_time() for feature in features] == [(10.5, 9.5)]


@pytest.mark.parametrize("option", ["-i", "-v"])
def test_existing_positive_domain_rejections(monkeypatch, tmp_path, option):
    monkeypatch.setattr(command, "PlateReconstruction", _forbid_reconstruction)
    with pytest.raises(SystemExit) as error:
        _run(monkeypatch, tmp_path, option, "0")
    assert error.value.code == 1
    assert not list(tmp_path.iterdir())


def test_nonfinite_validation_precedes_model_lookup(monkeypatch, tmp_path):
    def forbidden_lookup(*args, **kwargs):
        pytest.fail("Non-finite input reached the model resolver")

    monkeypatch.setattr(
        model_args, "_get_rotation_and_topology_files", forbidden_lookup
    )
    with pytest.raises(ValueError, match="finite"):
        _run(monkeypatch, tmp_path, "-i", "nan", model=True)
    assert not list(tmp_path.iterdir())


def test_valid_arguments_reach_model_transport_and_real_native_calculation(
    monkeypatch, tmp_path
):
    # Model transport only is doubled; no live network/cache behaviour claim.
    lookups = []

    def fixture_lookup(name, repo):
        lookups.append(name)
        return [str(DATA / "rotations.rot")], [str(DATA / "topologies.gpml")]

    monkeypatch.setattr(model_args, "_get_rotation_and_topology_files", fixture_lookup)
    calls = _observe_native(monkeypatch)
    _run(monkeypatch, tmp_path, subcommand="sc", model=True)
    assert lookups == ["synthetic-fixture"]
    assert [time for time, _ in calls] == [0.0, 1.0]
    assert len(list(tmp_path.glob("*.xy"))) == 2
