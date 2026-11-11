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


def _run_ridge_spreading_rate(
    monkeypatch, rotation_files, topology_files, output_prefix, *args
):
    _run_gplately(
        monkeypatch,
        "ridge-spreading-rate",
        "-r",
        *rotation_files,
        "-m",
        *topology_files,
        "-d",
        "0.5",
        *args,
        "--",
        output_prefix,
    )


@pytest.mark.parametrize("use_ptt", [False, True])
def test_ridge_spreading_rate(monkeypatch, tmp_path, muller_2019_model, use_ptt):
    rotation_files = muller_2019_model.get_rotation_model()
    topology_files = muller_2019_model.get_topologies()
    output_prefix = str(tmp_path / "spreading_rate")
    output_gpml = str(tmp_path / "spreading_rate.gpml")

    _run_ridge_spreading_rate(
        monkeypatch,
        rotation_files,
        topology_files,
        output_prefix,
        "-t",
        "0",
        "1",
        "-g",
        output_gpml,
        *(["--use-ptt"] if use_ptt else []),
    )

    plate_reconstruction = gplately.PlateReconstruction(rotation_files, topology_files)
    coverage_features = list(pygplates.FeatureCollection(output_gpml))
    assert len(coverage_features) == 2
    for time, coverage_feature in zip((0, 1), coverage_features):
        expected = plate_reconstruction.tessellate_mid_ocean_ridges(
            time,
            math.radians(0.5),
            use_ptt=use_ptt,
            output_obliquity_and_normal_and_left_right_plates=True,
        )
        assert len(expected) > 0

        output_filename = f"{output_prefix}_{time:0.2f}.xy"
        output_data = np.loadtxt(output_filename)
        np.testing.assert_allclose(_sort_rows(output_data), _sort_rows(expected))

        # Plate IDs should be written as integers.
        with open(output_filename) as output_file:
            for line in output_file:
                for plate_id in line.split()[6:8]:
                    assert plate_id == str(int(plate_id))

        # The coverage has the same points and values as the output file.
        assert coverage_feature.get_valid_time() == (time + 0.5, time - 0.5)
        geometry, scalars = coverage_feature.get_geometry(
            coverage_return=pygplates.CoverageReturn.geometry_and_scalars
        )
        np.testing.assert_allclose(
            geometry.to_lat_lon_array(), output_data[:, [1, 0]], atol=1e-10
        )
        for column, name in enumerate(
            [
                "SpreadingVelocityMagnitude",
                "SpreadingObliquityDegrees",
                "SpreadingLengthDegrees",
                "SpreadingNormalAzimuthDegrees",
                "LeftPlateId",
                "RightPlateId",
            ],
            start=2,
        ):
            np.testing.assert_allclose(
                scalars[pygplates.ScalarType.create_gpml(name)], output_data[:, column]
            )


def test_ridge_spreading_rate_options(monkeypatch, tmp_path, muller_2019_model):
    rotation_files = muller_2019_model.get_rotation_model()
    topology_files = muller_2019_model.get_topologies()
    output_prefix = str(tmp_path / "diverging")

    # Find all diverging plate boundaries.
    _run_ridge_spreading_rate(
        monkeypatch,
        rotation_files,
        topology_files,
        output_prefix,
        "-t",
        "0",
        "0",
        "--all-plate-boundaries",
        "--include-transform-segments",
        "--divergence-threshold",
        "0",
    )

    expected = gplately.PlateReconstruction(
        rotation_files, topology_files
    ).tessellate_mid_ocean_ridges(
        0,
        math.radians(0.5),
        spreading_feature_types=None,
        transform_segment_deviation_in_radians=None,
        divergence_threshold_in_cm_per_yr=0.0,
        output_obliquity_and_normal_and_left_right_plates=True,
    )
    assert len(expected) > 0
    np.testing.assert_allclose(
        _sort_rows(np.loadtxt(f"{output_prefix}_0.00.xy")), _sort_rows(expected)
    )


@pytest.mark.parametrize("subcommand", ["ridge-spreading-rate", "rsr"])
def test_ridge_spreading_rate_threshold_requires_modern_method(
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
            "--divergence-threshold",
            "0",
            "--",
            str(tmp_path / "spreading_rate"),
        )


@pytest.mark.parametrize(
    "subcommand", ["subduction-convergence", "ridge-spreading-rate"]
)
def test_coverage_valid_times_follow_time_increment(
    monkeypatch, tmp_path, muller_2019_model, subcommand
):
    output_gpml = str(tmp_path / "coverages.gpml")
    _run_gplately(
        monkeypatch,
        subcommand,
        "-r",
        *muller_2019_model.get_rotation_model(),
        "-m",
        *muller_2019_model.get_topologies(),
        "-t",
        "0",
        "4",
        "-i",
        "2",
        "-g",
        output_gpml,
        "--",
        str(tmp_path / "output"),
    )

    # Consecutive coverages neither overlap nor leave gaps.
    assert [
        feature.get_valid_time() for feature in pygplates.FeatureCollection(output_gpml)
    ] == [(1, -1), (3, 1), (5, 3)]


@pytest.mark.parametrize(
    "subcommand, args",
    [
        # Use '--use-ptt' since its output order is the same between runs (see '_sort_rows').
        ("subduction-convergence", ["-t", "0", "0", "--use-ptt", "-w"]),
        ("ridge-spreading-rate", ["-t", "0", "0", "--use-ptt"]),
        ("resolve-topologies", ["-t", "0"]),
    ],
)
def test_model_option(monkeypatch, tmp_path, muller_2019_model, subcommand, args):
    # Output from specifying the rotation and topology files...
    files_dir = tmp_path / "files"
    files_dir.mkdir()
    _run_gplately(
        monkeypatch,
        subcommand,
        "-r",
        *muller_2019_model.get_rotation_model(),
        "-m",
        *muller_2019_model.get_topologies(),
        *args,
        "--",
        str(files_dir / "output"),
    )

    # ...should be the same as from specifying the model (already downloaded by the fixture).
    model_dir = tmp_path / "model"
    model_dir.mkdir()
    _run_gplately(
        monkeypatch,
        subcommand,
        "--model",
        "Muller2019",
        "--model-repo",
        muller_2019_model.data_dir,
        *args,
        "--",
        str(model_dir / "output"),
    )

    # Compare the output files (ignoring GPlates layer files, which can contain file paths).
    filenames = sorted(
        path.name for path in files_dir.iterdir() if not path.name.endswith(".xml")
    )
    assert filenames
    assert filenames == sorted(
        path.name for path in model_dir.iterdir() if not path.name.endswith(".xml")
    )
    for filename in filenames:
        if filename.endswith(".xy"):
            assert (files_dir / filename).read_bytes() == (
                model_dir / filename
            ).read_bytes()
        elif filename.endswith(".shp"):
            # Resolved topologies get new feature IDs each run (and their order can differ),
            # so compare the geometries.
            assert _shapefile_geometries(files_dir / filename) == _shapefile_geometries(
                model_dir / filename
            )
        else:
            # The other shapefile components (attributes contain the new feature IDs).
            assert filename.endswith((".dbf", ".shx", ".prj")), filename


def _shapefile_geometries(filename):
    return sorted(
        tuple(geometry.to_lat_lon_list())
        for feature in pygplates.FeatureCollection(str(filename))
        for geometry in feature.get_all_geometries()
    )


@pytest.mark.parametrize(
    "args",
    [
        # Both a model and files.
        ["--model", "Muller2019", "-r", "rotations.rot"],
        # Neither.
        [],
        # Only rotation files.
        ["-r", "rotations.rot"],
        # Only topology files (eg, mistakenly using '-m' for a model name).
        ["-m", "Muller2019"],
        # A model folder without a model.
        ["-r", "rotations.rot", "-m", "topologies.gpml", "--model-repo", "models"],
    ],
)
@pytest.mark.parametrize(
    "subcommand",
    ["subduction-convergence", "ridge-spreading-rate", "resolve-topologies"],
)
def test_model_option_usage_errors(monkeypatch, tmp_path, subcommand, args):
    # Fails before reading any files (or downloading the model).
    with pytest.raises(SystemExit) as exc_info:
        _run_gplately(
            monkeypatch,
            subcommand,
            *args,
            "-t",
            "0",
            *(["0"] if subcommand != "resolve-topologies" else []),
            "--",
            str(tmp_path / "output"),
        )
    assert exc_info.value.code != 0
