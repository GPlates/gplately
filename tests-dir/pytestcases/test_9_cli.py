import math
import sys
import warnings

import numpy as np
import pygplates
import pytest

import gplately
from gplately.__main__ import main as gplately_main
from gplately.ptt import subduction_convergence as ptt_subduction_convergence

# Use the Müller et al. (2019) plate model at 0 and 1 Ma.

_EXTRA_OUTPUT_PARAMETERS = (
    "output_trench_normal",
    "output_distance_to_nearest_edge_of_trench",
)


def _sort_rows(data):
    # Before pyGPlates 1.1, the modern method's sample order can differ between runs
    # (it comes from pyGPlates' plate boundary statistics), so compare rows regardless of order.
    return data[np.lexsort(data.T[::-1])]


def _run_gplately(monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", ["gplately", *args])
    gplately_main()


def _run_subduction_convergence(
    monkeypatch, rotation_files, topology_files, output_prefix, output_gpml, *args
):
    _run_gplately(
        monkeypatch,
        "subduction-convergence",
        "-r",
        *rotation_files,
        "-m",
        *topology_files,
        "-t",
        "0",
        "1",
        "-d",
        "0.5",
        "-x",
        *_EXTRA_OUTPUT_PARAMETERS,
        "-g",
        output_gpml,
        "-w",
        *args,
        "--",
        output_prefix,
    )


def test_subduction_convergence(monkeypatch, tmp_path, muller_2019_model):
    rotation_files = muller_2019_model.get_rotation_model()
    topology_files = muller_2019_model.get_topologies()
    output_prefix = str(tmp_path / "convergence")
    output_gpml = str(tmp_path / "convergence.gpml")

    _run_subduction_convergence(
        monkeypatch, rotation_files, topology_files, output_prefix, output_gpml
    )

    plate_reconstruction = gplately.PlateReconstruction(rotation_files, topology_files)
    for time in (0, 1):
        expected = plate_reconstruction.tessellate_subduction_zones(
            time,
            math.radians(0.5),
            output_distance_to_nearest_edge_of_trench=True,
            output_trench_normal=True,
        )
        assert len(expected) > 0

        # Extra columns are in the order of the API (not the command line).
        output_filename = f"{output_prefix}_{time:0.2f}.xy"
        np.testing.assert_allclose(
            _sort_rows(np.loadtxt(output_filename)), _sort_rows(expected)
        )

        # Plate IDs should be written as integers.
        with open(output_filename) as output_file:
            for line in output_file:
                for plate_id in line.split()[8:10]:
                    assert plate_id == str(int(plate_id))

    # One coverage feature per time.
    assert len(pygplates.FeatureCollection(output_gpml)) == 2


def test_subduction_convergence_use_ptt_matches_previous_cli(
    monkeypatch, tmp_path, muller_2019_model
):
    rotation_files = muller_2019_model.get_rotation_model()
    topology_files = muller_2019_model.get_topologies()

    output_prefix = str(tmp_path / "convergence")
    output_gpml = str(tmp_path / "convergence.gpml")
    _run_subduction_convergence(
        monkeypatch,
        rotation_files,
        topology_files,
        output_prefix,
        output_gpml,
        "--use-ptt",
    )

    # What the CLI did before '--use-ptt' was added.
    previous_output_prefix = str(tmp_path / "previous")
    previous_output_gpml = str(tmp_path / "previous.gpml")
    with warnings.catch_warnings():
        # Same as '-w'.
        warnings.simplefilter("ignore")
        ptt_subduction_convergence.subduction_convergence_over_time(
            previous_output_prefix,
            "xy",
            rotation_files,
            topology_files,
            math.radians(0.5),
            0,
            1,
            1,
            anchor_plate_id=0,
            output_gpml_filename=previous_output_gpml,
            **dict.fromkeys(_EXTRA_OUTPUT_PARAMETERS, True),
        )

    for time in (0, 1):
        # Same values in the same order (zeros are written as '0.0' instead of '0').
        np.testing.assert_array_equal(
            np.loadtxt(f"{output_prefix}_{time:0.2f}.xy"),
            np.loadtxt(f"{previous_output_prefix}_{time:0.2f}.xy"),
        )

    def _coverages(filename):
        coverages = []
        for feature in pygplates.FeatureCollection(filename):
            geometry, scalars = feature.get_geometry(
                coverage_return=pygplates.CoverageReturn.geometry_and_scalars
            )
            coverages.append(
                (
                    feature.get_valid_time(),
                    geometry.to_lat_lon_list(),
                    {
                        name.to_qualified_string(): values
                        for name, values in scalars.items()
                    },
                )
            )
        return coverages

    assert _coverages(output_gpml) == _coverages(previous_output_gpml)


@pytest.mark.parametrize(
    "subcommand", ["subduction-convergence", "sc", "subduction_convergence"]
)
def test_subduction_convergence_threshold_requires_modern_method(
    monkeypatch, tmp_path, subcommand
):
    # Fails before reading any files, so they don't need to exist.
    # Also checks each subcommand name runs the new command.
    with pytest.raises(ValueError, match="--use-ptt"):
        _run_gplately(
            monkeypatch,
            subcommand,
            "-r",
            "rotations.rot",
            "-m",
            "topologies.gpml",
            "--use-ptt",
            "--convergence-threshold",
            "0",
            "--",
            str(tmp_path / "convergence"),
        )
