import argparse
import math
import os

import netCDF4
import numpy as np
import pygplates
import pytest
from conftest import logger

import gplately

# Aliased to avoid clashing with the Step 4 `paleobathymetry` function imported above --
# these are the gplately.commands modules backing the CLI subcommands of the same name.
from gplately.commands import continent_contouring as continent_contouring_cmd
from gplately.commands import paleobathymetry as paleobathymetry_cmd
from gplately.commands import sediment_thickness as sediment_thickness_cmd
from gplately.grids._utils import (
    DEFAULT_DECIMAL_PLACES_IN_TIME,
    DEFAULT_DISTANCE_GRID_DECIMAL_PLACES_IN_TIME,
    check_time_increment_covers_times,
    distance_grid_filename,
    format_time_in_filename,
    resolve_decimal_places_in_time,
    time_index_at_or_after,
    uniform_time_step,
)
from gplately.grids.continent_contouring import (
    _VALID_TIME_EPSILON,
    generate_passive_margins,
    passive_margin_polylines,
)
from gplately.grids.paleobathymetry import (
    AGE_DEPTH_MODELS,
    age_to_basement_depth,
    dutkiewicz_2017_sediment_thickness,
    paleobathymetry,
    sediment_isostatic_correction,
    simple_paleobathymetry,
)
from gplately.grids.pybacktrack_paleobathymetry import (
    merge_pybacktrack_paleobathymetry,
    ocean_age_to_depth_function,
)
from gplately.grids.sediment_thickness import (
    _check_proximity_features,
    generate_distance_grids,
    generate_input_points_grid,
    generate_sediment_thickness_grids,
)
from gplately.lib import shortest_path

try:
    import pybacktrack  # noqa: F401

    HAS_PYBACKTRACK = True
except ImportError:
    HAS_PYBACKTRACK = False

requires_pybacktrack = pytest.mark.skipif(
    not HAS_PYBACKTRACK,
    reason="pybacktrack is an optional dependency (gplately[paleobathymetry])",
)

logger.info(__name__)


# =============================================================================
# gplately.grids.paleobathymetry -- Steps 1, 3, 4 (age_to_basement_depth,
# dutkiewicz_2017_sediment_thickness, sediment_isostatic_correction, paleobathymetry)
# =============================================================================


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


def test_age_to_basement_depth_rhcw18_negative_age_is_ridge_crest_not_sea_level():
    # Tiny negative ages (floating-point/grid-interpolation noise right at a ridge) must be
    # clamped to the table's ridge-crest depth (~-2500 m), not fall through to np.interp's
    # `left` fallback taken literally as 0 m (sea level) -- see age_to_basement_depth().
    depth = age_to_basement_depth(np.array([-0.01, 0.0]), model="rhcw18")
    assert depth[0] < -2000.0
    np.testing.assert_allclose(depth[0], depth[1], atol=1.0)


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


def test_public_api_exports_paleobathymetry():
    for name in (
        "AGE_DEPTH_MODELS",
        "DUTKIEWICZ_2017_SEDIMENT_THICKNESS",
        "age_to_basement_depth",
        "dutkiewicz_2017_sediment_thickness",
        "sediment_isostatic_correction",
        "paleobathymetry",
    ):
        assert hasattr(gplately, name)


# =============================================================================
# gplately.grids.sediment_thickness -- Step 2 (generate_distance_grids) and Step 3's
# grid-level driver (generate_sediment_thickness_grids)
# =============================================================================


@pytest.fixture(scope="module")
def synthetic_age_grid_filename(tmp_path_factory):
    # A coarse, synthetic seafloor-age grid (not real data): a smooth age field with a
    # "continental" (NaN) patch, just to exercise masking/reconstruction end to end without
    # depending on a large downloaded age grid.
    lon = np.linspace(-180.0, 180.0, 37)  # 10 degree spacing
    lat = np.linspace(-90.0, 90.0, 19)
    lon_2d, lat_2d = np.meshgrid(lon, lat)

    age = np.abs(lon_2d) / 4.0  # 0-45 Ma, smoothly varying with longitude
    continent = (lat_2d > 20) & (lat_2d < 40) & (lon_2d > -20) & (lon_2d < 20)
    age[continent] = np.nan

    path = tmp_path_factory.mktemp("sediment_thickness") / "synthetic_age_0Ma.nc"
    gplately.write_netcdf_grid(str(path), age)
    return str(path)


def test_generate_input_points_grid():
    lon_1d, lat_1d, lon_flat, lat_flat = generate_input_points_grid(10.0)
    assert lon_1d.min() == -180.0
    assert lon_1d.max() == 180.0
    assert lat_1d.min() == -90.0
    assert lat_1d.max() == 90.0
    assert lon_flat.size == lat_flat.size == lon_1d.size * lat_1d.size


def test_generate_input_points_grid_invalid_spacing():
    with pytest.raises(ValueError):
        generate_input_points_grid(0.0)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_generate_distance_grids(
    gplately_muller_reconstruction_files,
    gplately_muller_static_geometries,
    synthetic_age_grid_filename,
):
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    _, _, cobs = gplately_muller_static_geometries

    age_grid_filenames_and_times = [(synthetic_age_grid_filename, 0.0)]

    distance_grids = generate_distance_grids(
        rotation_model=rotation_model,
        proximity_features=cobs,
        topological_features=topology_features,
        age_grid_filenames_and_times=age_grid_filenames_and_times,
        grid_spacing=10.0,
        time_increment=1,
        max_reconstruction_time=5,
        clamp_distance_km=3000.0,
    )

    assert set(distance_grids.keys()) == {0.0}
    lon, lat, grid = distance_grids[0.0]
    assert grid.shape == (lat.size, lon.size)

    finite = np.isfinite(grid)
    assert finite.any()
    assert np.all(grid[finite] >= 0.0)
    assert np.all(grid[finite] <= 3000.0)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_generate_distance_grids_with_continent_obstacles(
    gplately_muller_reconstruction_files,
    gplately_muller_static_geometries,
    synthetic_age_grid_filename,
):
    # Routing around continents should never produce a *shorter* path than a great-circle
    # line, on average -- some individual grid points can come out slightly shorter due to
    # the underlying grid's node-interpolation smoothing (see gplately.lib.shortest_path),
    # but the aggregate effect across many points should be unambiguously non-negative, and
    # at least some points (where a great-circle line cuts through land) should come out
    # meaningfully longer.
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    coastlines, _, cobs = gplately_muller_static_geometries

    age_grid_filenames_and_times = [(synthetic_age_grid_filename, 0.0)]
    common_kwargs = dict(
        rotation_model=rotation_model,
        proximity_features=cobs,
        topological_features=topology_features,
        age_grid_filenames_and_times=age_grid_filenames_and_times,
        grid_spacing=10.0,
        time_increment=1,
        max_reconstruction_time=5,
        clamp_distance_km=None,
    )

    great_circle_grids = generate_distance_grids(**common_kwargs)
    routed_grids = generate_distance_grids(
        continent_obstacle_features=coastlines,
        shortest_path_grid_subdivision_depth=5,
        **common_kwargs,
    )

    _, _, great_circle_grid = great_circle_grids[0.0]
    _, _, routed_grid = routed_grids[0.0]

    finite = np.isfinite(great_circle_grid) & np.isfinite(routed_grid)
    assert finite.any()
    diff = routed_grid[finite] - great_circle_grid[finite]

    assert diff.mean() >= 0.0
    assert diff.max() > 20.0  # at least one point should show real routing effect


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_generate_distance_grids_output_directory(
    tmp_path,
    gplately_muller_reconstruction_files,
    gplately_muller_static_geometries,
    synthetic_age_grid_filename,
):
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    _, _, cobs = gplately_muller_static_geometries
    age_grid_filenames_and_times = [(synthetic_age_grid_filename, 0.0)]
    output_dir = tmp_path / "distances"

    generate_distance_grids(
        rotation_model=rotation_model,
        proximity_features=cobs,
        topological_features=topology_features,
        age_grid_filenames_and_times=age_grid_filenames_and_times,
        grid_spacing=10.0,
        time_increment=1,
        max_reconstruction_time=5,
        output_directory=str(output_dir),
    )

    assert (output_dir / "mean_distance_10.0d_0.0.nc").is_file()


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_generate_sediment_thickness_grids(
    gplately_muller_reconstruction_files,
    gplately_muller_static_geometries,
    synthetic_age_grid_filename,
):
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    _, _, cobs = gplately_muller_static_geometries
    age_grid_filenames_and_times = [(synthetic_age_grid_filename, 0.0)]

    distance_grids = generate_distance_grids(
        rotation_model=rotation_model,
        proximity_features=cobs,
        topological_features=topology_features,
        age_grid_filenames_and_times=age_grid_filenames_and_times,
        grid_spacing=10.0,
        time_increment=1,
        max_reconstruction_time=5,
        clamp_distance_km=3000.0,
    )

    sediment_grids = generate_sediment_thickness_grids(
        age_grid_filenames_and_times, distance_grids
    )

    dist_lon, dist_lat, _ = distance_grids[0.0]
    lon, lat, thickness = sediment_grids[0.0]
    np.testing.assert_array_equal(lon, dist_lon)
    np.testing.assert_array_equal(lat, dist_lat)

    finite = np.isfinite(thickness)
    assert finite.any()
    assert np.all(thickness[finite] > 0.0)


def test_generate_sediment_thickness_grids_missing_time_raises():
    with pytest.raises(KeyError):
        generate_sediment_thickness_grids(
            [("does-not-matter.nc", 0.0)],
            distance_grids={},  # no entry for time 0.0
        )


# =============================================================================
# Age grids stored with descending latitude
# =============================================================================


def _write_age_grid_verbatim(path, lon, lat, age):
    """Write an age grid to netCDF with its coordinate order exactly as given.

    ``gplately.write_netcdf_grid()`` always writes latitude ascending, so it cannot
    produce the north-up (descending latitude) grids that GDAL writes -- and therefore
    that any GeoTIFF or GIS round-trip produces. These tests need both orders, so they
    write the netCDF directly.
    """
    with netCDF4.Dataset(str(path), "w") as cdf:
        cdf.createDimension("lon", len(lon))
        cdf.createDimension("lat", len(lat))
        cdf_lon = cdf.createVariable("lon", "f8", ("lon",))
        cdf_lat = cdf.createVariable("lat", "f8", ("lat",))
        cdf_z = cdf.createVariable("z", "f4", ("lat", "lon"), fill_value=np.nan)
        cdf_lon.units = "degrees_east"
        cdf_lat.units = "degrees_north"
        cdf_lon[:] = np.asarray(lon)
        cdf_lat[:] = np.asarray(lat)
        cdf_z[:, :] = np.asarray(age)
    return str(path)


@pytest.fixture(scope="module")
def latitude_ordered_age_grids(tmp_path_factory):
    """The same synthetic age field written twice: latitude ascending, then descending.

    The field increases monotonically from south to north, so reading one as though it
    were the other mirrors every row rather than leaving the result unchanged.

    The field has no NaNs on purpose. ``sample_grid()`` interpolates with
    ``scipy.ndimage.map_coordinates(order=1)``, which blends row *i* with row *i+1* even
    at an exact node; since ``0 * NaN`` is ``NaN``, every NaN spreads one row towards the
    start of the array. That direction is south in an ascending grid and north in a
    descending one, so a masked region is never flip-invariant and would mask the
    geometric property being tested here.
    """
    lon = np.linspace(-180.0, 180.0, 37)  # 10 degree spacing
    lat = np.linspace(-90.0, 90.0, 19)
    _, lat_2d = np.meshgrid(lon, lat)

    age = (lat_2d + 90.0) / 4.0  # 0-45 Ma, increasing from south to north

    directory = tmp_path_factory.mktemp("latitude_order")
    ascending = _write_age_grid_verbatim(
        directory / "age_ascending_0Ma.nc", lon, lat, age
    )
    descending = _write_age_grid_verbatim(
        directory / "age_descending_0Ma.nc", lon, lat[::-1], age[::-1, :]
    )
    return ascending, descending


def test_descending_latitude_age_grid_stays_descending(latitude_ordered_age_grids):
    """The fixture must really give us two different storage orders.

    ``read_netcdf_grid()`` re-sorts latitude only when it also has to realign longitudes
    from 0-360, which these grids do not need, so the descending file comes back
    descending. Without this the tests below could pass vacuously.
    """
    ascending, descending = latitude_ordered_age_grids

    age_up, _, lat_up = gplately.read_netcdf_grid(ascending, return_grids=True)
    age_down, _, lat_down = gplately.read_netcdf_grid(descending, return_grids=True)

    assert lat_up[0] < lat_up[-1]
    assert lat_down[0] > lat_down[-1]
    # The same field, just stored the other way up.
    np.testing.assert_array_equal(lat_down, lat_up[::-1])
    np.testing.assert_array_equal(age_down, age_up[::-1, :])


# =============================================================================
# Refusing inputs that produce plausible but wrong numbers
# =============================================================================


def _feature_with_geometry(geometry, name):
    feature = pygplates.Feature()
    feature.set_geometry(geometry)
    feature.set_name(name)
    return feature


def test_proximity_polylines_are_accepted():
    _check_proximity_features(
        [
            _feature_with_geometry(
                pygplates.PolylineOnSphere([(0, 0), (0, 10), (10, 10)]), "a margin"
            )
        ]
    )


def test_proximity_polygons_are_refused():
    """A continent-ocean boundary polygon wraps active margins as well as passive ones.

    Plate models ship exactly such polygons under the name "COBs", so this is the easy
    mistake to make, and it yields a plausible grid of wrong numbers rather than an error.
    """
    polygon = _feature_with_geometry(
        pygplates.PolygonOnSphere([(0, 0), (0, 10), (10, 10)]), "Africa COB"
    )
    with pytest.raises(ValueError, match="polygon"):
        _check_proximity_features([polygon])
    with pytest.raises(ValueError, match="Africa COB"):
        _check_proximity_features([polygon])


def test_empty_proximity_features_are_refused():
    """Nothing to measure to means every point reports half the Earth's circumference.

    Reachable through a proximity_feature_types filter that matches nothing, which produces
    a full grid of ~20,000 km rather than an error.
    """
    with pytest.raises(ValueError, match="empty"):
        _check_proximity_features([])

    # The reachable way to get here: a feature-type filter that matches nothing.
    with pytest.raises(ValueError, match="gpml:PassiveContinentalBoundary"):
        _check_proximity_features([], ["gpml:PassiveContinentalBoundary"])


def test_cli_skips_times_whose_distance_grid_was_never_written(monkeypatch, tmp_path):
    """'generate-distance-grids' writes nothing for an age grid it skipped.

    One absent file is that skip, and the remaining times should still be processed; all of
    them absent means the directory or the naming options are wrong, which is an error.
    """
    times = [0.0, 1.0]
    monkeypatch.setattr(
        sediment_thickness_cmd,
        "_resolve_age_grid_filenames_and_times",
        lambda args: ([("unused.nc", time) for time in times], None),
    )
    captured = {}
    monkeypatch.setattr(
        sediment_thickness_cmd,
        "generate_sediment_thickness_grids",
        lambda pairs, distance_grids, **kwargs: captured.update(
            pairs=pairs, distance_grids=distance_grids
        ),
    )
    monkeypatch.setattr(
        sediment_thickness_cmd,
        "read_netcdf_grid",
        lambda path, **kwargs: (np.zeros((2, 2)), np.zeros(2), np.zeros(2)),
    )

    args = argparse.Namespace(
        distance_grids_dir=str(tmp_path),
        grid_spacing=0.5,
        output_dir=str(tmp_path),
        decimal_places_in_time=None,
    )

    # Nothing written at all: that is a mistake, not a skip.
    with pytest.raises(Exception, match="No distance grids found"):
        sediment_thickness_cmd._run_generate_sediment_grids(args)

    # Only 1 Ma written: 0 Ma was skipped upstream, so carry on with what exists.
    (tmp_path / distance_grid_filename(0.5, 1.0, 1)).write_text("")
    sediment_thickness_cmd._run_generate_sediment_grids(args)
    assert [time for _, time in captured["pairs"]] == [1.0]
    assert sorted(captured["distance_grids"]) == [1.0]


def test_cli_does_not_take_proximity_features_from_the_plate_model():
    """The plate model's COBs layer must not be used as a silent fallback."""

    class PlateModelWithCOBs:
        def get_rotation_model(self):
            return ["rotations.rot"]

        def get_layer(self, name, return_none_if_not_exist=False):
            return [f"{name.lower()}.gpml"]

    args = argparse.Namespace(
        rotation_filenames=None, topology_filenames=None, proximity_filenames=[]
    )
    with pytest.raises(Exception, match="--proximity-features"):
        sediment_thickness_cmd._resolve_rotation_topology_proximity_files(
            args, PlateModelWithCOBs()
        )


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_an_all_nan_age_grid_produces_no_output(
    gplately_muller_reconstruction_files,
    gplately_muller_static_geometries,
    tmp_path,
):
    """An age grid with nothing usable in it used to write a grid of NaN.

    A file full of NaN is indistinguishable from a real result until something reads it, so
    the original workflow warns and moves on instead.
    """
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    _, _, cobs = gplately_muller_static_geometries

    empty_age_grid = str(tmp_path / "all_nan_0Ma.nc")
    gplately.write_netcdf_grid(empty_age_grid, np.full((19, 37), np.nan))

    results = generate_distance_grids(
        rotation_model=rotation_model,
        proximity_features=cobs,
        topological_features=topology_features,
        age_grid_filenames_and_times=[(empty_age_grid, 0.0)],
        grid_spacing=10.0,
        max_reconstruction_time=3,
        output_directory=str(tmp_path),
    )

    assert results == {}
    assert list(tmp_path.glob("mean_distance_*.nc")) == []


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_adjacent_time_slices_do_not_overlap(
    gplately_muller_reconstruction_files, gplately_muller_static_geometries
):
    """Each slice's features must end just before the next slice's begin.

    Without the epsilon both slices are valid at the instant they share, and GPlates draws
    two sets of margins on top of each other at every interval boundary.
    """
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    _, continental_polygons, _ = gplately_muller_static_geometries
    time_step = 1.0

    result = generate_passive_margins(
        rotation_model=rotation_model,
        continent_features=continental_polygons,
        topological_features=topology_features,
        times=[0.0, time_step],
        point_spacing_degrees=4.0,
        time_step=time_step,
    )

    valid_times = [
        feature.get_valid_time() for feature in result["continent_contour_features"]
    ]
    boundary = 0.0 + 0.5 * time_step  # shared by the 0 Ma and 1 Ma slices
    valid_at_boundary = [
        (begin, end) for begin, end in valid_times if end <= boundary <= begin
    ]
    # Exactly one slice may claim the boundary instant, never both -- and at least one
    # feature must reach it, or this would pass simply by finding no features.
    assert valid_at_boundary
    assert all(
        math.isclose(end, boundary, rel_tol=0, abs_tol=1e-9)
        for begin, end in valid_at_boundary
    ), valid_at_boundary
    assert _VALID_TIME_EPSILON > 0


# =============================================================================
# Output time step vs reconstruction time increment
# =============================================================================


def test_uniform_time_step_recovers_the_output_spacing():
    assert uniform_time_step([0.0, 10.0, 20.0]) == 10.0
    assert uniform_time_step([0.0, 0.5, 1.0]) == 0.5
    assert uniform_time_step([20.0, 0.0, 10.0]) == 10.0  # order must not matter
    assert uniform_time_step([5.0]) == 1.0  # one time: any positive increment will do


def test_uniform_time_step_rejects_uneven_times():
    """pyBacktrack takes one increment, so uneven times cannot be honoured silently."""
    with pytest.raises(ValueError, match="evenly spaced"):
        uniform_time_step([0.0, 10.0, 25.0])


def test_time_increment_must_divide_the_output_times():
    """Each age grid's walk starts at its time snapped up to a multiple of the increment.

    That snapping is what puts every age grid on one shared grid of times, so it stays --
    but a time that is not a multiple of the increment is then reconstructed from the wrong
    starting time, which has to be refused rather than silently tolerated.
    """
    check_time_increment_covers_times([0.0, 1.0, 2.0], 1.0)
    check_time_increment_covers_times([0.0, 0.5, 1.0], 0.5)
    check_time_increment_covers_times([0.0, 10.0, 20.0], 2.0)

    with pytest.raises(ValueError, match="not a multiple of"):
        check_time_increment_covers_times([0.0, 0.5, 1.0], 1.0)


@pytest.mark.parametrize("bad", [0, -1.0])
def test_time_increment_must_be_positive(bad):
    """0 used to surface as a bare ZeroDivisionError, and a negative walked backwards."""
    with pytest.raises(ValueError, match="must be positive"):
        check_time_increment_covers_times([0.0, 1.0], bad)


def test_times_from_a_fractional_step_snap_to_themselves():
    """Accumulated float error must not push a valid time up onto the next increment.

    The times a fractional step produces are not exact: 0 + 3 * 0.1 is
    0.30000000000000004, and 0.30000000000000004 / 0.1 is 3.0000000000000004, which a plain
    ceil() sends to 4 -- reconstructing a 0.3 Ma grid from 0.4 Ma. The validator and the
    snap share one rule so that every time the validator accepts maps to itself.
    """
    times = [i * 0.1 for i in range(11)]
    assert repr(times[3]) == "0.30000000000000004"  # the case this guards

    check_time_increment_covers_times(times, 0.1)
    for time in times:
        snapped = time_index_at_or_after(time, 0.1) * 0.1
        assert math.isclose(snapped, time, rel_tol=1e-9, abs_tol=1e-9)

    # A time genuinely off the grid still snaps forward, as it must.
    assert time_index_at_or_after(0.35, 0.1) == 4


@pytest.mark.parametrize(
    "missing",
    ["output_directory", "static_polygon_filename", "present_day_age_grid_filename"],
)
def test_step_5_prerequisites_checked_before_steps_1_to_4(
    synthetic_age_grid_filename, missing
):
    """These used to be checked at the end, after the whole pipeline had already run."""
    kwargs = dict(
        rotation_model="does-not-exist.rot",
        proximity_features="does-not-exist.gpml",
        topological_features="does-not-exist.gpml",
        age_grid_filenames_and_times=[(synthetic_age_grid_filename, 0.0)],
        output_directory="unused",
        pybacktrack=True,
        static_polygon_filename="does-not-exist.gpml",
        present_day_age_grid_filename=synthetic_age_grid_filename,
    )
    kwargs[missing] = None

    with pytest.raises(ValueError, match=missing):
        simple_paleobathymetry(**kwargs)


def test_distance_grids_validate_before_opening_any_file():
    """The check has to precede the rotation/topology file loading, not follow it.

    Every filename here is nonexistent, so anything other than the ValueError below means
    the validation happened too late to be useful.
    """
    with pytest.raises(ValueError, match="not a multiple of"):
        generate_distance_grids(
            rotation_model="does-not-exist.rot",
            proximity_features="does-not-exist.gpml",
            topological_features="does-not-exist.gpml",
            age_grid_filenames_and_times=[("does-not-exist.nc", 0.5)],
            time_increment=1,
        )


def test_uneven_times_rejected_before_steps_1_to_4_run(synthetic_age_grid_filename):
    """Step 5 cannot express uneven times, and saying so afterwards is no use."""
    with pytest.raises(ValueError, match="evenly spaced"):
        simple_paleobathymetry(
            rotation_model="does-not-exist.rot",
            proximity_features="does-not-exist.gpml",
            topological_features="does-not-exist.gpml",
            age_grid_filenames_and_times=[
                (synthetic_age_grid_filename, time) for time in (0.0, 10.0, 25.0)
            ],
            output_directory="unused",
            pybacktrack=True,
            static_polygon_filename="does-not-exist.gpml",
            present_day_age_grid_filename=synthetic_age_grid_filename,
        )


def test_pybacktrack_gets_the_output_step_not_the_reconstruction_increment(
    monkeypatch, synthetic_age_grid_filename, tmp_path
):
    """Step 5's increment is the spacing of the grids Steps 1-4 wrote.

    It used to be handed Step 2's reconstruction increment. With output every 10 Myr and
    the default 1 Myr reconstruction increment, pyBacktrack generated 0, 1, 2 ... 20 Ma and
    went looking for paleobathymetry_1Ma.nc -- files Step 4 never wrote.
    """
    from gplately.grids import pybacktrack_paleobathymetry as pybacktrack_module
    from gplately.grids import sediment_thickness as sediment_thickness_module

    lon = np.linspace(-180.0, 180.0, 37)
    lat = np.linspace(-90.0, 90.0, 19)
    times = [0.0, 10.0, 20.0]
    captured = {}

    monkeypatch.setattr(
        sediment_thickness_module,
        "generate_distance_grids",
        lambda *, age_grid_filenames_and_times, **kwargs: {
            time: (lon, lat, np.full((lat.size, lon.size), 500.0))
            for _, time in age_grid_filenames_and_times
        },
    )
    monkeypatch.setattr(
        sediment_thickness_module,
        "generate_sediment_thickness_grids",
        lambda age_grid_filenames_and_times, distance_grids, **kwargs: {
            time: (lon, lat, np.full((lat.size, lon.size), 100.0))
            for _, time in age_grid_filenames_and_times
        },
    )
    monkeypatch.setattr(
        pybacktrack_module,
        "merge_pybacktrack_paleobathymetry",
        lambda **kwargs: captured.update(kwargs),
    )

    simple_paleobathymetry(
        rotation_model="rotations.rot",
        proximity_features="cobs.gpml",
        topological_features="topologies.gpml",
        age_grid_filenames_and_times=[
            (synthetic_age_grid_filename, time) for time in times
        ],
        time_increment=1,
        output_directory=str(tmp_path),
        pybacktrack=True,
        static_polygon_filename="static_polygons.gpml",
        present_day_age_grid_filename=synthetic_age_grid_filename,
    )

    assert captured["time_increment"] == 10.0
    assert captured["oldest_time"] == 20.0
    assert captured["youngest_time"] == 0.0


@pytest.mark.parametrize(
    "add_parser, argv",
    [
        (paleobathymetry_cmd.add_parser, ["paleobathymetry", "outdir"]),
        (sediment_thickness_cmd.add_parser, ["generate-distance-grids", "outdir"]),
    ],
)
def test_time_step_and_time_increment_are_independent(add_parser, argv):
    """--time-step used to set both, so a coarse output step coarsened the sampling too."""
    parser = _build_subparser(add_parser)

    defaults = parser.parse_args(argv + ["--time-step", "10"])
    assert defaults.time_step == 10
    assert defaults.time_increment == 1

    both = parser.parse_args(argv + ["--time-step", "10", "--time-increment", "0.5"])
    assert both.time_step == 10
    assert both.time_increment == 0.5


def test_distance_grids_cli_passes_the_reconstruction_increment(monkeypatch, tmp_path):
    """The CLI must hand generate_distance_grids() --time-increment, not --time-step."""
    captured = {}
    monkeypatch.setattr(
        sediment_thickness_cmd,
        "_resolve_age_grid_filenames_and_times",
        lambda args: ([("unused.nc", 0.0)], None),
    )
    monkeypatch.setattr(
        sediment_thickness_cmd,
        "_resolve_rotation_topology_proximity_files",
        lambda args, plate_model: ("rot", "topo", "cobs"),
    )
    monkeypatch.setattr(
        sediment_thickness_cmd, "_resolve_distance_grid_kwargs", lambda args, model: {}
    )
    monkeypatch.setattr(
        sediment_thickness_cmd,
        "generate_distance_grids",
        lambda **kwargs: captured.update(kwargs),
    )

    parser = _build_subparser(sediment_thickness_cmd.add_parser)
    args = parser.parse_args(
        ["generate-distance-grids", str(tmp_path), "--time-step", "10"]
    )
    sediment_thickness_cmd._run_generate_distance_grids(args)

    assert captured["time_increment"] == 1


# =============================================================================
# Time in output filenames
# =============================================================================


def _flat_distance_grids(times, value_km=500.0):
    lon = np.linspace(-180.0, 180.0, 37)
    lat = np.linspace(-90.0, 90.0, 19)
    grid = np.full((lat.size, lon.size), value_km)
    return {time: (lon, lat, grid) for time in times}


def test_time_in_filename_matches_pybacktrack_format():
    """Our filenames and pyBacktrack's must be generated by the same rule.

    pyBacktrack builds its own format as ``{time:.Nf}`` from
    ``output_file_decimal_places_in_time`` (and the ``merge_...`` equivalent), and goes
    looking for the Steps 1-4 grids written here. If the two rules ever diverge, Step 5
    silently finds nothing to merge.
    """
    for places in range(4):
        for time in (0.0, 1.0, 2.5, 10.25, 137.0):
            assert format_time_in_filename(time, places) == "{:.{}f}".format(
                time, places
            )


def test_distance_grid_filename_reproduces_original_workflow_names():
    """The default must not rename what the predicting-sediment-thickness workflow wrote."""
    assert (
        distance_grid_filename(0.5, 0.0, DEFAULT_DISTANCE_GRID_DECIMAL_PLACES_IN_TIME)
        == "mean_distance_0.5d_0.0.nc"
    )
    assert (
        distance_grid_filename(0.1, 137.0, DEFAULT_DISTANCE_GRID_DECIMAL_PLACES_IN_TIME)
        == "mean_distance_0.1d_137.0.nc"
    )


def test_integer_times_keep_their_original_filenames(
    synthetic_age_grid_filename, tmp_path
):
    """Whole-number times must be named exactly as before this parameter existed."""
    times = [0.0, 1.0, 2.0]
    generate_sediment_thickness_grids(
        [(synthetic_age_grid_filename, time) for time in times],
        _flat_distance_grids(times),
        output_directory=str(tmp_path),
    )
    assert sorted(path.name for path in tmp_path.glob("*.nc")) == [
        "sediment_thickness_0Ma.nc",
        "sediment_thickness_1Ma.nc",
        "sediment_thickness_2Ma.nc",
    ]


def test_colliding_times_raise_before_doing_any_work(
    synthetic_age_grid_filename, tmp_path
):
    """A fractional time step at zero decimal places used to overwrite silently.

    0, 0.5, 1, 1.5 and 2 Ma all formatted to three distinct names, so five grids were
    computed and three written. The error has to arrive before the computation, not after.
    """
    times = [0.0, 0.5, 1.0, 1.5, 2.0]
    with pytest.raises(ValueError, match="decimal_places_in_time"):
        generate_sediment_thickness_grids(
            [(synthetic_age_grid_filename, time) for time in times],
            _flat_distance_grids(times),
            output_directory=str(tmp_path),
        )
    assert list(tmp_path.iterdir()) == []


def test_fractional_times_are_written_once_decimal_places_allow_it(
    synthetic_age_grid_filename, tmp_path
):
    times = [0.0, 0.5, 1.0, 1.5, 2.0]
    generate_sediment_thickness_grids(
        [(synthetic_age_grid_filename, time) for time in times],
        _flat_distance_grids(times),
        output_directory=str(tmp_path),
        decimal_places_in_time=1,
    )
    assert sorted(path.name for path in tmp_path.glob("*.nc")) == [
        "sediment_thickness_0.0Ma.nc",
        "sediment_thickness_0.5Ma.nc",
        "sediment_thickness_1.0Ma.nc",
        "sediment_thickness_1.5Ma.nc",
        "sediment_thickness_2.0Ma.nc",
    ]


def test_no_output_directory_means_no_filename_constraint(synthetic_age_grid_filename):
    """Times that could not be named distinctly are still fine in memory."""
    times = [0.0, 0.5, 1.0]
    results = generate_sediment_thickness_grids(
        [(synthetic_age_grid_filename, time) for time in times],
        _flat_distance_grids(times),
    )
    assert sorted(results) == times


@pytest.mark.parametrize("bad", [-1, 1.5, "1", True])
def test_decimal_places_in_time_must_be_a_non_negative_int(bad):
    with pytest.raises(ValueError, match="non-negative integer"):
        resolve_decimal_places_in_time(bad, DEFAULT_DECIMAL_PLACES_IN_TIME)


def test_decimal_places_in_time_accepts_numpy_integers():
    """Times and steps often come from numpy, so np.int64 must not be rejected."""
    assert resolve_decimal_places_in_time(np.int64(2), 0) == 2
    assert resolve_decimal_places_in_time(np.int32(0), 1) == 0


def test_age_grid_times_may_be_an_iterator(synthetic_age_grid_filename, tmp_path):
    """The filename check walks the input before the main loop, so it must be materialised.

    A generator would otherwise be exhausted by the check, leaving the main loop with
    nothing to do -- no grids written, no error raised.
    """
    times = [0.0, 1.0]
    results = generate_sediment_thickness_grids(
        ((synthetic_age_grid_filename, time) for time in times),
        _flat_distance_grids(times),
        output_directory=str(tmp_path),
    )
    assert sorted(results) == times
    assert len(list(tmp_path.glob("*.nc"))) == 2


def test_cli_rejects_colliding_times_before_reading_distance_grids(
    monkeypatch, tmp_path
):
    """The CLI reads every distance grid before handing over to the library.

    Without its own up-front check it would do all that I/O and only then hit the
    library's ValueError, which is exactly the wasted work the check exists to avoid.
    """
    times = [0.0, 0.5, 1.0]
    monkeypatch.setattr(
        sediment_thickness_cmd,
        "_resolve_age_grid_filenames_and_times",
        lambda args: ([("unused.nc", time) for time in times], None),
    )

    def fail_if_read(*args, **kwargs):
        raise AssertionError("distance grids were read before the collision was caught")

    monkeypatch.setattr(sediment_thickness_cmd, "read_netcdf_grid", fail_if_read)

    args = argparse.Namespace(
        distance_grids_dir=str(tmp_path),
        grid_spacing=0.5,
        output_dir=str(tmp_path),
        decimal_places_in_time=None,
    )
    with pytest.raises(ValueError, match="decimal_places_in_time"):
        sediment_thickness_cmd._run_generate_sediment_grids(args)


def test_decimal_places_in_time_defaults_when_none():
    assert resolve_decimal_places_in_time(None, 1) == 1
    assert resolve_decimal_places_in_time(0, 1) == 0


def test_simple_paleobathymetry_leaves_each_step_its_own_default(
    monkeypatch, synthetic_age_grid_filename, tmp_path
):
    """Step 2's grids must stay findable by the 'generate-sediment-grids' subcommand.

    The distance grids carry one decimal place of time and these grids carry none, both
    inherited from the workflows they came from. Resolving the default here and pushing it
    down would rename the distance grids to something the CLI reader does not look for, so
    an unset value has to reach each step unresolved.
    """
    from gplately.grids import sediment_thickness as sediment_thickness_module

    seen = {}
    lon = np.linspace(-180.0, 180.0, 37)
    lat = np.linspace(-90.0, 90.0, 19)

    def fake_distance_grids(*, age_grid_filenames_and_times, **kwargs):
        seen["distance"] = kwargs.get("decimal_places_in_time")
        return {
            time: (lon, lat, np.full((lat.size, lon.size), 500.0))
            for _, time in age_grid_filenames_and_times
        }

    def fake_sediment_grids(age_grid_filenames_and_times, distance_grids, **kwargs):
        seen["sediment"] = kwargs.get("decimal_places_in_time")
        return {
            time: (lon, lat, np.full((lat.size, lon.size), 100.0))
            for _, time in age_grid_filenames_and_times
        }

    monkeypatch.setattr(
        sediment_thickness_module, "generate_distance_grids", fake_distance_grids
    )
    monkeypatch.setattr(
        sediment_thickness_module,
        "generate_sediment_thickness_grids",
        fake_sediment_grids,
    )

    simple_paleobathymetry(
        rotation_model="unused.rot",
        proximity_features="unused.gpml",
        topological_features="unused.gpml",
        age_grid_filenames_and_times=[(synthetic_age_grid_filename, 0.0)],
        output_directory=str(tmp_path),
    )

    assert seen == {"distance": None, "sediment": None}
    assert (tmp_path / "paleobathymetry_0Ma.nc").is_file()


@pytest.mark.skipif(not HAS_PYBACKTRACK, reason="requires the pybacktrack package")
@pytest.mark.parametrize("decimal_places", [None, 0, 2])
def test_decimal_places_forwarded_to_both_pybacktrack_parameters(
    monkeypatch, decimal_places
):
    """Step 5 names its own output, and finds ours, using two separate pyBacktrack knobs.

    Both were hardwired to 0, so a fractional time step broke the hand-off in step with
    the writers. They must now follow the same value the writers used.
    """
    import pybacktrack

    captured = {}

    def fake_reconstruct(output_file_prefix, **kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(
        pybacktrack, "reconstruct_paleo_bathymetry_grids", fake_reconstruct
    )

    merge_pybacktrack_paleobathymetry(
        output_file_prefix="out/paleobathymetry_${time}Ma.nc",
        merge_paleobathymetry_filename_format="in/paleobathymetry_${time}Ma.nc",
        rotation_filenames="rotations.rot",
        static_polygon_filename="static_polygons.gpml",
        present_day_age_grid_filename="age_0Ma.nc",
        grid_spacing_degrees=0.5,
        oldest_time=2.0,
        decimal_places_in_time=decimal_places,
    )

    expected = (
        DEFAULT_DECIMAL_PLACES_IN_TIME if decimal_places is None else decimal_places
    )
    assert captured["output_file_decimal_places_in_time"] == expected
    assert captured["merge_paleo_bathymetry_file_decimal_places_in_time"] == expected


def test_sediment_thickness_grids_invariant_to_age_grid_latitude_order(
    latitude_ordered_age_grids,
):
    """Step 3 must read the same ages from a grid whichever way up it is stored.

    ``generate_sediment_thickness_grids()`` samples the age grid onto the distance grid's
    points. Handing ``sample_grid()`` an unsigned extent reads a descending-latitude grid
    upside-down, applying southern-hemisphere ages to northern points -- silently, and
    with entirely plausible-looking output.
    """
    ascending, descending = latitude_ordered_age_grids

    lon = np.linspace(-180.0, 180.0, 37)
    lat = np.linspace(-90.0, 90.0, 19)
    distance_grids = {0.0: (lon, lat, np.full((lat.size, lon.size), 500.0))}

    _, _, thickness_up = generate_sediment_thickness_grids(
        [(ascending, 0.0)], distance_grids
    )[0.0]
    _, _, thickness_down = generate_sediment_thickness_grids(
        [(descending, 0.0)], distance_grids
    )[0.0]

    assert np.isfinite(thickness_up).any()
    np.testing.assert_allclose(thickness_down, thickness_up, rtol=1e-9)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_simple_paleobathymetry_invariant_to_age_grid_latitude_order(
    gplately_muller_reconstruction_files,
    gplately_muller_static_geometries,
    latitude_ordered_age_grids,
):
    """The same property end to end, covering the other two age-grid sampling sites.

    ``simple_paleobathymetry()`` samples the age grid once inside
    ``generate_distance_grids()`` (to decide which points are ocean, and how old each is)
    and again for Step 1's basement depth.
    """
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    _, _, cobs = gplately_muller_static_geometries
    ascending, descending = latitude_ordered_age_grids

    def run(age_grid_filename):
        return simple_paleobathymetry(
            rotation_model=rotation_model,
            proximity_features=cobs,
            topological_features=topology_features,
            age_grid_filenames_and_times=[(age_grid_filename, 0.0)],
            grid_spacing=10.0,
            time_increment=1,
            max_reconstruction_time=5,
        )[0.0]

    _, _, depth_up = run(ascending)
    _, _, depth_down = run(descending)

    assert np.isfinite(depth_up).any()
    np.testing.assert_allclose(depth_down, depth_up, rtol=1e-9)


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_simple_paleobathymetry_sediment_thickness_kwargs_can_override_max_distance_km(
    gplately_muller_reconstruction_files,
    gplately_muller_static_geometries,
    synthetic_age_grid_filename,
):
    # simple_paleobathymetry() always passes clamp_distance_km through to
    # generate_sediment_thickness_grids() as max_distance_km; a caller-supplied
    # sediment_thickness_kwargs["max_distance_km"] (a documented override) must not collide
    # with that and raise "got multiple values for keyword argument".
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    _, _, cobs = gplately_muller_static_geometries

    result = simple_paleobathymetry(
        rotation_model=rotation_model,
        proximity_features=cobs,
        topological_features=topology_features,
        age_grid_filenames_and_times=[(synthetic_age_grid_filename, 0.0)],
        grid_spacing=10.0,
        max_reconstruction_time=5,
        clamp_distance_km=3000.0,
        sediment_thickness_kwargs={"max_distance_km": 5000.0},
    )

    assert set(result.keys()) == {0.0}


def test_public_api_exports_sediment_thickness():
    for name in (
        "generate_input_points_grid",
        "generate_distance_grids",
        "generate_sediment_thickness_grids",
        "simple_paleobathymetry",
    ):
        assert hasattr(gplately, name)


# =============================================================================
# gplately.grids.pybacktrack_paleobathymetry -- Step 5 (merge_pybacktrack_paleobathymetry)
# =============================================================================


@pytest.fixture(scope="module")
def synthetic_paleobathymetry_grid_filenames(tmp_path_factory):
    # Coarse, synthetic Steps-1-4 paleobathymetry grids (not real data) at the same spacing as
    # pybacktrack's bundled reconstruction, just to exercise the merge plumbing end to end
    # without depending on a full gplately.simple_paleobathymetry() run.
    lon = np.linspace(-180.0, 180.0, 37)  # 10 degree spacing
    lat = np.linspace(-90.0, 90.0, 19)
    lon_2d, lat_2d = np.meshgrid(lon, lat)

    out_dir = tmp_path_factory.mktemp("pybacktrack_merge_input")
    for time in (0.0, 1.0):
        depth = -2000.0 - 10.0 * time - np.abs(lat_2d)  # arbitrary smooth "ocean" depth
        ocean = np.abs(lon_2d) < 90.0  # arbitrary ocean/continent split
        depth = np.where(ocean, depth, np.nan)
        gplately.write_netcdf_grid(
            str(out_dir / "paleobathymetry_{:.0f}Ma.nc".format(time)), depth
        )
    return str(out_dir / "paleobathymetry_${time}Ma.nc")


@requires_pybacktrack
def test_merge_pybacktrack_paleobathymetry(
    tmp_path, synthetic_paleobathymetry_grid_filenames
):
    import pybacktrack.bundle_data as bundle_data

    output_dir = tmp_path / "pybacktrack_output"
    output_dir.mkdir()

    gplately.merge_pybacktrack_paleobathymetry(
        output_file_prefix=str(output_dir / "paleobathymetry_${time}Ma.nc"),
        merge_paleobathymetry_filename_format=synthetic_paleobathymetry_grid_filenames,
        rotation_filenames=bundle_data.BUNDLE_RECONSTRUCTION_ROTATION_FILENAMES,
        static_polygon_filename=bundle_data.BUNDLE_RECONSTRUCTION_STATIC_POLYGON_FILENAME,
        present_day_age_grid_filename=bundle_data.BUNDLE_AGE_GRID_FILENAME,
        grid_spacing_degrees=10.0,
        oldest_time=1,
        youngest_time=0,
        time_increment=1,
        age_depth_model="gdh1",
    )

    merged_path = output_dir / "paleobathymetry_0Ma.nc"
    assert merged_path.is_file()

    merged, lon, lat = gplately.read_netcdf_grid(str(merged_path), return_grids=True)
    finite = np.isfinite(merged)
    assert finite.any()

    # pyBacktrack's whole purpose is to fill in cells the merged-in Steps-1-4 grid left as NaN
    # (submerged continental crust / crust that has since subducted); the merged output should
    # therefore have at least as much coverage as the input.
    step4, _, _ = gplately.read_netcdf_grid(
        synthetic_paleobathymetry_grid_filenames.replace("${time}", "0"),
        return_grids=True,
    )
    assert finite.sum() >= np.isfinite(step4).sum()


@requires_pybacktrack
def test_merge_pybacktrack_paleobathymetry_unknown_age_depth_model(
    synthetic_paleobathymetry_grid_filenames,
):
    import pybacktrack.bundle_data as bundle_data

    with pytest.raises(ValueError):
        gplately.merge_pybacktrack_paleobathymetry(
            output_file_prefix="unused_${time}",
            merge_paleobathymetry_filename_format=synthetic_paleobathymetry_grid_filenames,
            rotation_filenames=bundle_data.BUNDLE_RECONSTRUCTION_ROTATION_FILENAMES,
            static_polygon_filename=bundle_data.BUNDLE_RECONSTRUCTION_STATIC_POLYGON_FILENAME,
            present_day_age_grid_filename=bundle_data.BUNDLE_AGE_GRID_FILENAME,
            grid_spacing_degrees=10.0,
            oldest_time=1,
            # "parsons_sclater" used to belong here, as the one model with no pyBacktrack
            # counterpart. pyBacktrack is now given gplately's own conversion instead of one
            # of its built-in models, so every model gplately knows works -- only a name it
            # does not know is an error.
            age_depth_model="not-a-model",
        )


# -----------------------------------------------------------------------------
# Steps 1-4 and Step 5 must compute depth with the same model
# -----------------------------------------------------------------------------


@pytest.mark.parametrize("model", AGE_DEPTH_MODELS)
def test_every_age_depth_model_can_be_handed_to_pybacktrack(model):
    """pyBacktrack takes a function of age, so none of gplately's models is out of reach."""
    assert callable(ocean_age_to_depth_function(model))


@pytest.mark.parametrize("alias", ["richards", "r18", "RHCW18", "ghd1", "stein_stein"])
def test_age_depth_model_aliases_are_accepted(alias):
    assert callable(ocean_age_to_depth_function(alias))


def test_unknown_age_depth_model_is_refused():
    with pytest.raises(ValueError, match="Unknown age-depth model"):
        ocean_age_to_depth_function("not-a-model")


def test_the_age_to_depth_callable_survives_pickling():
    """pyBacktrack pickles this when use_all_cpus is set, and a closure would not survive."""
    import pickle

    function = ocean_age_to_depth_function("gdh1")
    assert pickle.loads(pickle.dumps(function))(50.0) == function(50.0)


@pytest.mark.skipif(not HAS_PYBACKTRACK, reason="requires the pybacktrack package")
@pytest.mark.parametrize("model", AGE_DEPTH_MODELS)
def test_step_5_depths_match_steps_1_to_4_exactly(model):
    """The property the whole arrangement exists for.

    Step 5's output is merged into Steps 1-4's, so any disagreement between the two shows
    up as a step change at the merge boundary. Going through gplately's own conversion
    makes them the same computation rather than two implementations that ought to match.
    """
    import pybacktrack

    ages = np.arange(0.0, 200.5, 0.5)
    ours = -np.asarray(age_to_basement_depth(ages, model=model), dtype=float)
    theirs = np.array(
        [
            pybacktrack.convert_age_to_depth(
                float(age), ocean_age_to_depth_function(model)
            )
            for age in ages
        ]
    )
    np.testing.assert_allclose(theirs, ours, rtol=0, atol=0)


@pytest.mark.skipif(not HAS_PYBACKTRACK, reason="requires the pybacktrack package")
def test_a_substituted_richards_table_reaches_step_5_too(tmp_path):
    """Replacing the RHCW18 table changes the model, so Step 5 has to be given it as well."""
    import pybacktrack

    table = tmp_path / "table.dat"
    ages = np.arange(0.0, 401.0, 1.0)
    np.savetxt(table, np.column_stack([ages, -(2400.0 + 30.0 * np.sqrt(ages))]))

    function = ocean_age_to_depth_function("rhcw18", str(table))
    sample = np.arange(0.0, 200.5, 0.5)
    ours = -np.asarray(
        age_to_basement_depth(
            sample, model="rhcw18", richards_table_filename=str(table)
        ),
        dtype=float,
    )
    theirs = np.array(
        [pybacktrack.convert_age_to_depth(float(age), function) for age in sample]
    )
    np.testing.assert_allclose(theirs, ours, rtol=0, atol=0)

    # And it is genuinely a different relationship from the shipped table.
    shipped = -np.asarray(age_to_basement_depth(sample, model="rhcw18"), dtype=float)
    assert np.nanmax(np.abs(ours - shipped)) > 100.0


@pytest.mark.skipif(not HAS_PYBACKTRACK, reason="requires the pybacktrack package")
def test_crosby09_is_not_pybacktracks_crosby_2007():
    """Why a built-in model is not used, recorded as a measurement.

    gplately's "crosby09" is Crosby & McKenzie (2009)'s empirical piecewise fit;
    pyBacktrack's CROSBY_2007 is the plate-cooling model of Crosby's 2007 thesis. The names
    line up, the models do not, and picking the built-in by name would put a step change at
    the merge boundary.
    """
    import pybacktrack

    ages = np.arange(0.0, 200.5, 0.5)
    built_in = np.array(
        [
            pybacktrack.convert_age_to_depth(
                float(age), pybacktrack.AGE_TO_DEPTH_MODEL_CROSBY_2007
            )
            for age in ages
        ]
    )
    ours = -np.asarray(age_to_basement_depth(ages, model="crosby09"), dtype=float)

    assert np.nanmax(np.abs(ours - built_in)) > 100.0
    assert abs(ours[0] - built_in[0]) > 25.0  # they differ at the ridge crest too


def test_the_age_to_depth_callable_is_memoised():
    """Called once per ocean point per decompaction step, so repeats should be cheap."""
    function = ocean_age_to_depth_function("rhcw18")
    first = function(37.5)
    info_before = ocean_age_to_depth_function("rhcw18").func.cache_info()
    assert function(37.5) == first
    info_after = ocean_age_to_depth_function("rhcw18").func.cache_info()
    assert info_after.hits > info_before.hits


@pytest.mark.parametrize("pybacktrack_enabled", [False, True])
def test_unusable_age_depth_model_reported_before_steps_1_to_4_run(
    pybacktrack_enabled, synthetic_age_grid_filename, tmp_path
):
    """Finding out from Step 4, after Step 2 has run for hours, is too late."""
    with pytest.raises(ValueError, match="Unknown age-depth model"):
        simple_paleobathymetry(
            rotation_model="does-not-exist.rot",
            proximity_features="does-not-exist.gpml",
            topological_features="does-not-exist.gpml",
            age_grid_filenames_and_times=[(synthetic_age_grid_filename, 0.0)],
            age_depth_model="not-a-model",
            output_directory=str(tmp_path),
            pybacktrack=pybacktrack_enabled,
            static_polygon_filename="does-not-exist.gpml",
            present_day_age_grid_filename=synthetic_age_grid_filename,
        )


def test_public_api_exports_pybacktrack_paleobathymetry():
    assert hasattr(gplately, "merge_pybacktrack_paleobathymetry")


# =============================================================================
# gplately.grids.continent_contouring -- dynamically-contoured passive margins
# =============================================================================


def test_passive_margin_polylines_no_subduction_zones():
    contour = pygplates.PolylineOnSphere(
        [(0.0, 0.0), (0.0, 10.0), (0.0, 20.0), (0.0, 30.0)]
    )
    margins = passive_margin_polylines(contour, [], max_distance_radians=0.1)
    assert len(margins) == 1
    assert list(margins[0].get_points()) == list(contour.get_points())


def test_passive_margin_polylines_entirely_active():
    contour = pygplates.PolylineOnSphere([(0.0, 0.0), (0.0, 10.0), (0.0, 20.0)])
    # A subduction zone running right alongside the whole contour.
    subduction_zone = pygplates.PolylineOnSphere([(1.0, 0.0), (1.0, 10.0), (1.0, 20.0)])
    margins = passive_margin_polylines(
        contour, [subduction_zone], max_distance_radians=np.radians(2.0)
    )
    assert margins == []


def test_passive_margin_polylines_splits_around_active_segment():
    # A contour running along longitude 0, latitude 0 to 30.
    # Widely-spaced points so that non-adjacent segments stay far from the subduction zone
    # even accounting for it being close to a shared vertex of two segments.
    contour = pygplates.PolylineOnSphere(
        [(0.0, lon) for lon in (0.0, 20.0, 40.0, 60.0, 80.0)]
    )
    # A subduction zone entirely within the middle segment's longitude range (40-60), not
    # touching either of its endpoints, so only that segment is classified as active.
    subduction_zone = pygplates.PolylineOnSphere([(2.0, 45.0), (2.0, 55.0)])
    margins = passive_margin_polylines(
        contour, [subduction_zone], max_distance_radians=np.radians(3.0)
    )
    # The lon [0, 20, 40] run and the lon [60, 80] run should survive as two separate
    # polylines, with the lon [40, 60] segment (active margin) removed between them.
    assert len(margins) == 2
    margin_lons = [
        [round(lon, 6) for _, lon in (p.to_lat_lon() for p in margin.get_points())]
        for margin in margins
    ]
    assert margin_lons == [[0.0, 20.0, 40.0], [60.0, 80.0]]


def test_passive_margin_polylines_merges_across_closed_ring_seam():
    # A closed ring (points[0] == points[-1], as continent contours are) with exactly one
    # active edge, located away from the arbitrary start/end point. Without the seam fix,
    # the single remaining passive arc would be incorrectly split into two polylines at the
    # array boundary (points[0]/points[-1]) instead of returned as one.
    lons = [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0, 0.0]
    ring = pygplates.PolylineOnSphere([(0.0, lon) for lon in lons])
    # Near the midpoint of the (180, 225) edge only.
    subduction_zone = pygplates.PolylineOnSphere([(0.5, 195.0), (0.5, 210.0)])

    margins = passive_margin_polylines(
        ring, [subduction_zone], max_distance_radians=np.radians(1.1)
    )

    assert len(margins) == 1
    assert len(list(margins[0].get_points())) == len(lons) - 1


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_generate_passive_margins(
    gplately_muller_reconstruction_files, gplately_muller_static_geometries
):
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    _, continental_polygons, _ = gplately_muller_static_geometries

    result = generate_passive_margins(
        rotation_model=rotation_model,
        continent_features=continental_polygons,
        topological_features=topology_features,
        times=[0.0],
        point_spacing_degrees=2.0,
    )

    assert set(result.keys()) == {
        "continent_contour_features",
        "passive_margin_features",
        "continent_masks",
    }
    assert len(result["continent_contour_features"]) > 0
    # Passive margins are what's left of the contours after removing active (near-subduction)
    # segments; a single contour can fragment into several passive-margin polylines, so there
    # can legitimately be more passive-margin features than contour features. Just check some
    # were found (Earth has plenty of passive margins).
    assert len(result["passive_margin_features"]) > 0

    mask = result["continent_masks"][0.0]
    assert mask.dtype == bool
    # Some, but not all, of the globe should be continental crust.
    assert 0.0 < mask.mean() < 1.0


@pytest.mark.skipif(
    int(os.getenv("GPLATELY_TEST_LEVEL", 0)) < 1,
    reason="This testcase downloads a full Muller2019 plate model from the Internet. Set GPLATELY_TEST_LEVEL higher than 1 to activate it.",
)
def test_generate_passive_margins_output_directory(
    tmp_path, gplately_muller_reconstruction_files, gplately_muller_static_geometries
):
    rotation_model, topology_features, _ = gplately_muller_reconstruction_files
    _, continental_polygons, _ = gplately_muller_static_geometries

    output_dir = tmp_path / "passive_margins"
    generate_passive_margins(
        rotation_model=rotation_model,
        continent_features=continental_polygons,
        topological_features=topology_features,
        times=[0.0],
        point_spacing_degrees=2.0,
        output_directory=str(output_dir),
    )

    assert (output_dir / "continent_contour_features.gpmlz").is_file()
    assert (output_dir / "passive_margin_features.gpmlz").is_file()
    assert (output_dir / "continent_mask_0.nc").is_file()


def test_public_api_exports_continent_contouring():
    for name in ("generate_passive_margins", "passive_margin_polylines"):
        assert hasattr(gplately, name)


# =============================================================================
# gplately.lib.shortest_path -- shortest-path-around-obstacles engine backing
# generate_distance_grids()'s continent_obstacle_features option
# =============================================================================


def test_routed_distance_around_obstacle_is_longer_than_great_circle():
    grid = shortest_path.Grid(subdivision_depth=5)  # ~2.8 degree spacing

    # A square "continent" straddling the equator between longitude -10 and 10.
    obstacle = pygplates.PolygonOnSphere([(20, -10), (20, 10), (-20, 10), (-20, -10)])
    obstacle_grid = grid.create_obstacle_grid([obstacle])

    source = pygplates.PointOnSphere(0.0, -30.0)
    distance_grid = obstacle_grid.create_distance_grid([source])

    # Directly on the other side of the obstacle: a great-circle line would cut through it.
    target = pygplates.PointOnSphere(0.0, 30.0)
    routed_distance = distance_grid.shortest_distance(target)
    great_circle_distance = pygplates.GeometryOnSphere.distance(source, target)

    assert routed_distance is not None
    assert routed_distance > great_circle_distance


def test_routed_distance_matches_great_circle_when_path_is_clear():
    grid = shortest_path.Grid(subdivision_depth=5)

    obstacle = pygplates.PolygonOnSphere([(20, -10), (20, 10), (-20, 10), (-20, -10)])
    obstacle_grid = grid.create_obstacle_grid([obstacle])

    source = pygplates.PointOnSphere(0.0, -30.0)
    distance_grid = obstacle_grid.create_distance_grid([source])

    # Same side as the source, well away from the obstacle: no detour needed.
    target = pygplates.PointOnSphere(0.0, -35.0)
    routed_distance = distance_grid.shortest_distance(target)
    great_circle_distance = pygplates.GeometryOnSphere.distance(source, target)

    assert routed_distance is not None
    # Grid discretisation introduces a small amount of noise, but a clear path should stay
    # close to the great-circle distance.
    assert routed_distance == pytest.approx(
        great_circle_distance, abs=math.radians(1.0)
    )


def test_target_unreachable_when_completely_enclosed():
    grid = shortest_path.Grid(subdivision_depth=5)

    # A source inside a fully-enclosed obstacle can't reach a target outside it (and vice
    # versa) -- here we enclose the *target* instead, with the source outside.
    enclosing_obstacle = pygplates.PolygonOnSphere([(5, -5), (5, 5), (-5, 5), (-5, -5)])
    obstacle_grid = grid.create_obstacle_grid([enclosing_obstacle])

    source = pygplates.PointOnSphere(0.0, -60.0)
    # distance_threshold_radians small enough that the long way around is excluded.
    distance_grid = obstacle_grid.create_distance_grid(
        [source], distance_threshold_radians=math.radians(30.0)
    )

    target = pygplates.PointOnSphere(0.0, 60.0)
    routed_distance = distance_grid.shortest_distance(target)
    assert routed_distance is None


# =============================================================================
# gplately.commands.{paleobathymetry,sediment_thickness,continent_contouring} -- CLI
# argument-parsing for 'paleobathymetry'/'pb', 'generate-distance-grids'/'gdg',
# 'generate-sediment-grids'/'gsg' and 'generate-passive-margins'/'gpm' (#463)
#
# NOTE: these were originally a separate file (test_10_paleobathymetry_cli.py), kept apart
# on the grounds that they're a different kind of test from the rest of this file -- fast,
# offline, fixture-free argparse wiring checks, versus the numeric/behavioural library tests
# above (which often need the Muller2019 fixtures and are sometimes network-gated). Merged
# in on request to keep one file per feature area; flagging the concern here rather than
# silently dropping it, in case it's worth splitting them out again later.
# =============================================================================


def _build_subparser(add_parser_func):
    """Register a single subcommand via *add_parser_func* on a fresh top-level parser."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    add_parser_func(subparsers)
    return parser


def test_paleobathymetry_parser_defaults():
    parser = _build_subparser(paleobathymetry_cmd.add_parser)
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
    parser = _build_subparser(paleobathymetry_cmd.add_parser)
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
    parser = _build_subparser(paleobathymetry_cmd.add_parser)
    args = parser.parse_args(["pb", "outdir"])
    assert args.output_dir == "outdir"
    assert callable(args.func)


def test_paleobathymetry_rejects_unknown_age_depth_model():
    parser = _build_subparser(paleobathymetry_cmd.add_parser)
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["paleobathymetry", "outdir", "--age-depth-model", "no-such-model"]
        )


def test_paleobathymetry_requires_output_dir():
    parser = _build_subparser(paleobathymetry_cmd.add_parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["paleobathymetry"])


def test_generate_distance_grids_parser_defaults():
    parser = _build_subparser(sediment_thickness_cmd.add_parser)
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
    parser = _build_subparser(sediment_thickness_cmd.add_parser)
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
    parser = _build_subparser(sediment_thickness_cmd.add_parser)
    args = parser.parse_args(["gdg", "outdir"])
    assert args.output_dir == "outdir"
    assert callable(args.func)


def test_generate_sediment_grids_parser_defaults():
    parser = _build_subparser(sediment_thickness_cmd.add_parser)
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


@pytest.mark.parametrize(
    "add_parser, argv",
    [
        (paleobathymetry_cmd.add_parser, ["paleobathymetry", "outdir"]),
        (sediment_thickness_cmd.add_parser, ["generate-distance-grids", "outdir"]),
        (
            sediment_thickness_cmd.add_parser,
            ["generate-sediment-grids", "outdir", "--distance-grids-dir", "distances"],
        ),
        (continent_contouring_cmd.add_parser, ["generate-passive-margins", "outdir"]),
    ],
)
def test_decimal_places_in_time_flag_on_every_subcommand(add_parser, argv):
    """All four subcommands write time-stamped files, so all four need the flag."""
    parser = _build_subparser(add_parser)

    assert parser.parse_args(argv).decimal_places_in_time is None
    assert (
        parser.parse_args(
            argv + ["--decimal-places-in-time", "2"]
        ).decimal_places_in_time
        == 2
    )


def test_generate_sediment_grids_requires_distance_grids_dir():
    # --distance-grids-dir is the one flag marked required=True across all four parsers --
    # this is real argparse-enforced validation, not just a default, so it's worth pinning.
    parser = _build_subparser(sediment_thickness_cmd.add_parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["generate-sediment-grids", "outdir"])


def test_generate_sediment_grids_alias_gsg_resolves_same_subcommand():
    parser = _build_subparser(sediment_thickness_cmd.add_parser)
    args = parser.parse_args(["gsg", "outdir", "--distance-grids-dir", "distances"])
    assert args.output_dir == "outdir"
    assert args.distance_grids_dir == "distances"
    assert callable(args.func)


def test_generate_passive_margins_parser_defaults():
    parser = _build_subparser(continent_contouring_cmd.add_parser)
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
    parser = _build_subparser(continent_contouring_cmd.add_parser)
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
    parser = _build_subparser(continent_contouring_cmd.add_parser)
    args = parser.parse_args(["gpm", "outdir"])
    assert args.output_dir == "outdir"
    assert callable(args.func)


def test_generate_passive_margins_requires_output_dir():
    parser = _build_subparser(continent_contouring_cmd.add_parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["generate-passive-margins"])
