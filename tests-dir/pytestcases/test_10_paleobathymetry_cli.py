"""Argparse-level tests for the CLI subcommands added alongside gplately.paleobathymetry /
gplately.sediment_thickness (#449): 'paleobathymetry', 'generate-distance-grids',
'generate-sediment-grids' and 'generate-passive-margins'.

These build a parser, register one subcommand's arguments via its module's add_parser(),
parse a fixed argv list, and assert on the resulting namespace -- they never call the
command's actual func (which would need real files, a network connection, and a plate
model). That makes them fast and offline, unlike tests-dir/test-cli.sh, whose coverage of
these same subcommands needs a full end-to-end run instead. See gh issue #463.
"""

import argparse

import pytest

from gplately.commands import continent_contouring, paleobathymetry, sediment_thickness


def _build_subparser(add_parser_func):
    """Register a single subcommand via *add_parser_func* on a fresh top-level parser."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    add_parser_func(subparsers)
    return parser


# =============================================================================
# 'paleobathymetry' / 'pb'
# =============================================================================


def test_paleobathymetry_parser_defaults():
    parser = _build_subparser(paleobathymetry.add_parser)
    args = parser.parse_args(["paleobathymetry", "outdir"])

    assert args.command == "paleobathymetry"
    assert args.output_dir == "outdir"
    assert args.model_name is None
    assert args.plate_model_repo == "plate-model-repo"
    assert args.age_grid_template is None
    assert args.min_time == 0
    assert args.max_time == 100
    assert args.time_step == 1
    assert args.grid_spacing == 0.5
    assert args.anchor_plate_id is None
    assert args.proximity_filenames == []
    assert args.rotation_filenames == []
    assert args.topology_filenames == []
    assert args.max_reconstruction_time is None
    assert args.clamp_distance_km == 3000.0
    assert args.route_around_continents is False
    assert args.continent_obstacle_filenames == []
    assert args.shortest_path_grid_depth == 6
    assert args.age_depth_model == "gdh1"
    assert args.richards_table is None
    assert args.pybacktrack is False
    assert args.static_polygons is None
    assert args.present_day_age_grid is None
    assert callable(args.func)


def test_paleobathymetry_parser_all_flags():
    parser = _build_subparser(paleobathymetry.add_parser)
    args = parser.parse_args(
        [
            "paleobathymetry",
            "outdir",
            "-m",
            "muller2025",
            "-f",
            "my-repo",
            "--age-grid-template",
            "agegrids/age_{time:.0f}Ma.nc",
            "-e",
            "0",
            "-s",
            "50",
            "--time-step",
            "2",
            "-r",
            "0.25",
            "-a",
            "701",
            "--proximity-features",
            "cobs1.gpml",
            "cobs2.gpml",
            "--rotations",
            "a.rot",
            "b.rot",
            "--topologies",
            "topo.gpml",
            "--max-reconstruction-time",
            "180",
            "--clamp-distance-km",
            "2000",
            "--route-around-continents",
            "--continent-obstacles",
            "coastlines.gpml",
            "--shortest-path-grid-depth",
            "8",
            "--age-depth-model",
            "rhcw18",
            "--richards-table",
            "table.txt",
            "--pybacktrack",
            "--static-polygons",
            "statics.gpml",
            "--present-day-age-grid",
            "age_0Ma.nc",
        ]
    )

    assert args.output_dir == "outdir"
    assert args.model_name == "muller2025"
    assert args.plate_model_repo == "my-repo"
    assert args.age_grid_template == "agegrids/age_{time:.0f}Ma.nc"
    assert args.min_time == 0
    assert args.max_time == 50
    assert args.time_step == 2
    assert args.grid_spacing == 0.25
    assert args.anchor_plate_id == 701
    assert args.proximity_filenames == ["cobs1.gpml", "cobs2.gpml"]
    assert args.rotation_filenames == ["a.rot", "b.rot"]
    assert args.topology_filenames == ["topo.gpml"]
    assert args.max_reconstruction_time == 180
    assert args.clamp_distance_km == 2000
    assert args.route_around_continents is True
    assert args.continent_obstacle_filenames == ["coastlines.gpml"]
    assert args.shortest_path_grid_depth == 8
    assert args.age_depth_model == "rhcw18"
    assert args.richards_table == "table.txt"
    assert args.pybacktrack is True
    assert args.static_polygons == "statics.gpml"
    assert args.present_day_age_grid == "age_0Ma.nc"


def test_paleobathymetry_alias_pb_resolves_same_subcommand():
    parser = _build_subparser(paleobathymetry.add_parser)
    args = parser.parse_args(["pb", "outdir"])
    assert args.output_dir == "outdir"
    assert callable(args.func)


def test_paleobathymetry_rejects_unknown_age_depth_model():
    parser = _build_subparser(paleobathymetry.add_parser)
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["paleobathymetry", "outdir", "--age-depth-model", "no-such-model"]
        )


def test_paleobathymetry_requires_output_dir():
    parser = _build_subparser(paleobathymetry.add_parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["paleobathymetry"])


# =============================================================================
# 'generate-distance-grids' / 'gdg'
# =============================================================================


def test_generate_distance_grids_parser_defaults():
    parser = _build_subparser(sediment_thickness.add_parser)
    args = parser.parse_args(["generate-distance-grids", "outdir"])

    assert args.command == "generate-distance-grids"
    assert args.output_dir == "outdir"
    assert args.model_name is None
    assert args.min_time == 0
    assert args.max_time == 100
    assert args.time_step == 1
    assert args.grid_spacing == 0.5
    assert args.anchor_plate_id is None
    assert args.proximity_filenames == []
    assert args.rotation_filenames == []
    assert args.topology_filenames == []
    assert args.max_reconstruction_time is None
    assert args.clamp_distance_km == 3000.0
    assert args.route_around_continents is False
    assert args.continent_obstacle_filenames == []
    assert args.shortest_path_grid_depth == 6
    assert callable(args.func)


def test_generate_distance_grids_parser_all_flags():
    parser = _build_subparser(sediment_thickness.add_parser)
    args = parser.parse_args(
        [
            "generate-distance-grids",
            "outdir",
            "-m",
            "muller2025",
            "-f",
            "my-repo",
            "--age-grid-template",
            "agegrids/age_{time:.0f}Ma.nc",
            "-e",
            "10",
            "-s",
            "60",
            "--time-step",
            "5",
            "-r",
            "1.0",
            "-a",
            "701",
            "--proximity-features",
            "cobs.gpml",
            "--rotations",
            "a.rot",
            "--topologies",
            "topo.gpml",
            "--max-reconstruction-time",
            "170",
            "--clamp-distance-km",
            "1500",
            "--route-around-continents",
            "--continent-obstacles",
            "coastlines.gpml",
            "--shortest-path-grid-depth",
            "7",
        ]
    )

    assert args.output_dir == "outdir"
    assert args.model_name == "muller2025"
    assert args.plate_model_repo == "my-repo"
    assert args.age_grid_template == "agegrids/age_{time:.0f}Ma.nc"
    assert args.min_time == 10
    assert args.max_time == 60
    assert args.time_step == 5
    assert args.grid_spacing == 1.0
    assert args.anchor_plate_id == 701
    assert args.proximity_filenames == ["cobs.gpml"]
    assert args.rotation_filenames == ["a.rot"]
    assert args.topology_filenames == ["topo.gpml"]
    assert args.max_reconstruction_time == 170
    assert args.clamp_distance_km == 1500
    assert args.route_around_continents is True
    assert args.continent_obstacle_filenames == ["coastlines.gpml"]
    assert args.shortest_path_grid_depth == 7


def test_generate_distance_grids_alias_gdg_resolves_same_subcommand():
    parser = _build_subparser(sediment_thickness.add_parser)
    args = parser.parse_args(["gdg", "outdir"])
    assert args.output_dir == "outdir"
    assert callable(args.func)


# =============================================================================
# 'generate-sediment-grids' / 'gsg'
# =============================================================================


def test_generate_sediment_grids_parser_defaults():
    parser = _build_subparser(sediment_thickness.add_parser)
    args = parser.parse_args(
        ["generate-sediment-grids", "outdir", "--distance-grids-dir", "distances"]
    )

    assert args.command == "generate-sediment-grids"
    assert args.output_dir == "outdir"
    assert args.distance_grids_dir == "distances"
    assert args.model_name is None
    assert args.min_time == 0
    assert args.max_time == 100
    assert args.time_step == 1
    assert args.grid_spacing == 0.5
    assert args.anchor_plate_id is None
    assert callable(args.func)


def test_generate_sediment_grids_requires_distance_grids_dir():
    # --distance-grids-dir is the one flag marked required=True across all four parsers --
    # this is real argparse-enforced validation, not just a default, so it's worth pinning.
    parser = _build_subparser(sediment_thickness.add_parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["generate-sediment-grids", "outdir"])


def test_generate_sediment_grids_alias_gsg_resolves_same_subcommand():
    parser = _build_subparser(sediment_thickness.add_parser)
    args = parser.parse_args(["gsg", "outdir", "--distance-grids-dir", "distances"])
    assert args.output_dir == "outdir"
    assert args.distance_grids_dir == "distances"
    assert callable(args.func)


# =============================================================================
# 'generate-passive-margins' / 'gpm'
# =============================================================================


def test_generate_passive_margins_parser_defaults():
    parser = _build_subparser(continent_contouring.add_parser)
    args = parser.parse_args(["generate-passive-margins", "outdir"])

    assert args.command == "generate-passive-margins"
    assert args.output_dir == "outdir"
    assert args.model_name is None
    assert args.plate_model_repo == "plate-model-repo"
    assert args.rotation_filenames == []
    assert args.topology_filenames == []
    assert args.continent_filenames == []
    assert args.min_time == 0
    assert args.max_time == 100
    assert args.time_step == 1
    assert args.point_spacing == 0.25
    assert args.area_threshold_km2 == 0.0
    assert args.buffer_and_gap_km == 0.0
    assert args.exclusion_area_threshold_km2 == 800000.0
    assert args.max_distance_active_margin_km == 500.0
    assert args.anchor_plate_id is None
    assert callable(args.func)


def test_generate_passive_margins_parser_all_flags():
    parser = _build_subparser(continent_contouring.add_parser)
    args = parser.parse_args(
        [
            "generate-passive-margins",
            "outdir",
            "-m",
            "muller2025",
            "-f",
            "my-repo",
            "--rotations",
            "a.rot",
            "--topologies",
            "topo.gpml",
            "--continents",
            "continents.gpml",
            "-e",
            "0",
            "-s",
            "200",
            "--time-step",
            "10",
            "-r",
            "0.5",
            "--area-threshold-km2",
            "1000",
            "--buffer-and-gap-km",
            "50",
            "--exclusion-area-threshold-km2",
            "500000",
            "--max-distance-active-margin-km",
            "300",
            "-a",
            "701",
        ]
    )

    assert args.output_dir == "outdir"
    assert args.model_name == "muller2025"
    assert args.plate_model_repo == "my-repo"
    assert args.rotation_filenames == ["a.rot"]
    assert args.topology_filenames == ["topo.gpml"]
    assert args.continent_filenames == ["continents.gpml"]
    assert args.min_time == 0
    assert args.max_time == 200
    assert args.time_step == 10
    assert args.point_spacing == 0.5
    assert args.area_threshold_km2 == 1000
    assert args.buffer_and_gap_km == 50
    assert args.exclusion_area_threshold_km2 == 500000
    assert args.max_distance_active_margin_km == 300
    assert args.anchor_plate_id == 701


def test_generate_passive_margins_alias_gpm_resolves_same_subcommand():
    parser = _build_subparser(continent_contouring.add_parser)
    args = parser.parse_args(["gpm", "outdir"])
    assert args.output_dir == "outdir"
    assert callable(args.func)


def test_generate_passive_margins_requires_output_dir():
    parser = _build_subparser(continent_contouring.add_parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["generate-passive-margins"])
