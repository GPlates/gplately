# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

GPlately is an object-oriented Python interface to pyGPlates for plate tectonic reconstruction (points, lines, polygons, and rasters through deep geologic time, plate velocities, subduction/spreading rates, and paleomap plotting).

`pygplates` is a required native dependency and is **not** pip-installable — it must come from conda-forge or the project's Docker image. Local dev/test environments are conda/micromamba-based (see `tests-dir/test-env.yml`, `conda/`, `docker/`).

## Commands

Install in editable mode (env must already have pygplates, e.g. via conda-forge):
```
pip install -e .
pip install -e ".[dev]"   # adds black, isort, bumpver, pip-tools, pytest
```

Run the standard automated suite (what CI runs, matches `pytest.ini`'s `testpaths`):
```
python -m pytest -vv tests-dir/pytestcases
```

Run a single test file or test:
```
python -m pytest tests-dir/pytestcases/test_2_points.py
python -m pytest tests-dir/pytestcases/test_2_points.py::test_name
```

Run the full test suite, including slow/network-heavy cases gated behind `GPLATELY_TEST_LEVEL` (see `os.getenv("GPLATELY_TEST_LEVEL", 0)` checks in `tests-dir/pytestcases/test_4_rasters.py`):
```
GPLATELY_TEST_LEVEL=100 python -m pytest -vv tests-dir/pytestcases
```

Other test entry points (see `tests-dir/readme.md`):
```
./tests-dir/test-cli.sh                       # CLI smoke tests
./tests-dir/unittest/run_all.sh               # human-verified scripts, not run in CI
./tests-dir/unittest/test-seafloor-gridding.sh
./scripts/run_all_notebooks.sh                # executes example notebooks
```

Build docs:
```
./scripts/build-sphinx-doc.sh
```

CLI entry point (`gplately = "gplately.__main__:main"`), exposes many subcommands (`list_models`, `combine`, `feature_filter`, `reset_feature_type`, `seafloor_grids`, `regrid`, `rotate_grid`, `fix_crossovers`, `remove_rotations`, `cleanup_topologies`, `convert_xy_to_gplates`, `diagnose_rotations`, `resolve_topologies`, `rotation_tools`, `separate_ridge_transform_segments`, `subduction_convergence`, `gpmdb`, `get_cli_config_example`):
```
gplately <subcommand> --help
```

## CI behaviour to know about

`.github/workflows/build_and_test.yml` first diffs the PR/push against its base and only runs the (matrixed, multi-OS/multi-Python) test job if the change touches `gplately/**`. Changes confined to notebooks, docs, or `tests-dir/` won't trigger the full matrix — keep this in mind when reasoning about whether CI will exercise a change.

## Architecture

Everything lives under `gplately/`. `gplately/__init__.py` is the public API surface — check it to see what's actually exported vs. internal.

**Core user-facing classes**, one per top-level module:
- `PlateReconstruction` (`reconstruction.py`) — wraps a rotation model + topology/static-polygon features; the central object most other classes take as a constructor argument.
- `Points` (`points.py`) — reconstructs point data through time using a `PlateReconstruction`.
- `Raster` (`raster.py`) — reconstructs gridded/raster data through time.
- `PlotTopologies` (`plot/plot_topologies.py`) — builds geopandas GeoDataFrames of reconstructed geometries for a given reconstruction time.
- `DataServer` (`data_server.py`) — legacy data-fetching interface; newer code should prefer `plate_model_manager`'s `PlateModelManager`/`PresentDayRasterManager`, which `gplately` re-exports directly.
- `SeafloorGrid` / `TopologySeafloorGrid` / `IsochronSeafloorGrid` (`grids/`) — seafloor age/spreading-rate gridding through time.
- `auxiliary.py` — `get_plate_reconstruction()` / `get_gplot()` convenience factories that wire the above together from a model name.

**Plotting** (`plot/`): `PlotEngine` (`plot_engine.py`) is an abstract base class with two backends, `CartopyPlotEngine` and `PygmtPlotEngine`; `PlotTopologies` is injected with one of these to render. Add new plotting backends by subclassing `PlotEngine`.

**Gridding** (`grids/`): `oceans.py` holds `SeafloorGrid`; `_grids.py`/`_utils.py` hold shared netCDF grid I/O helpers used across the grid classes.

**Low-level reconstruction math** (`lib/`): `reconstruct.py` (point/velocity reconstruction functions used by `Points`), plus `rotation.py`, `icosahedron.py`, `isopolate.py`, `polyline.py`, `quaternions.py` — geometric/numerical primitives, not typically touched directly by end users.

**PTT — PlateTectonicTools** (`ptt/`): a vendored/ported toolkit (topology resolution, subduction convergence, ridge spreading rate, rotation-file diagnostics, crossover fixing, etc.). Each module is both importable from Python and wired up as a `gplately` CLI subcommand in `__main__.py`.

**CLI subcommands** (`commands/`): implementations for `gplately`-CLI-only subcommands (`feature_filter_cmd.py`, `list_models.py`, `regrid.py`, `reset_feature_type.py`, `rotate_grid.py`, `seafloor_grids.py`) that aren't part of `ptt/`. Each module exposes an `add_parser(subparser)` used by `__main__.py`.

**Backward compatibility** (`deprecated/`): `gplately/__init__.py` registers these modules into `sys.modules` under old import paths (`gplately.mapping`, `gplately.pygplates`, `gplately.download`, `gplately.data`, `gplately.parallel`, `gplately.oceans`) so pre-refactor import paths keep working. When renaming/moving public functionality, add a shim here rather than breaking old imports.

**Shared utilities** (`utils/`): `check_pmm.py` enforces the `plate-model-manager` version required by the installed `gplately` (via `REQUIRED_PMM_VERSION`, checked at import time in `__init__.py`); `log_utils.py` sets up logging on import (`gplately.log`, `gplately-pytest.log`, `ptt.log` files); `feature_filter.py` / `feature_transformer.py` back the `feature_filter` CLI command; `io_utils.py` has `load_feature_collection`.

## Testing layout

`tests-dir/` (named to avoid clashing with Python packaging conventions, not `tests/`):
- `pytestcases/` — the formal, CI-run suite. Files are numbered by feature area (`test_0_imports.py`, `test_1_reconstructions.py`, `test_2_points.py`, `test_3_plottopologies.py`, `test_4_rasters.py`, `test_5_seafloorgrid.py`, `test_6_geometry.py`, `test_7_models.py`, `test_8_feature_filter.py`). Shared fixtures (reconstruction/points/raster/seafloorgrid objects built from the Müller et al. 2019 model) live in `pytestcases/conftest.py`.
- `unittest/` — scripts needing human/visual verification; explicitly excluded from `pytest.ini`'s default run and not part of CI.
- `debug_utils/` — ad hoc developer scripts for reproducing bugs/perf tuning; not a formal test suite.
- `data/` — test data files.

Tests pull real plate models/rasters at run time via `plate_model_manager` (`PlateModelManager`, `PresentDayRasterManager`) rather than committed fixtures — expect network access and caching directories (`data-cache/`, `plate-model-repo/`) to appear when running tests locally.
