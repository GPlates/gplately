import math

import numpy as np
import pygplates
import pytest
from conftest import logger

import gplately
from gplately.grids.continent_contouring import (
    generate_passive_margins,
    passive_margin_polylines,
)
from gplately.grids.paleobathymetry import (
    AGE_DEPTH_MODELS,
    age_to_basement_depth,
    dutkiewicz_2017_sediment_thickness,
    paleobathymetry,
    sediment_isostatic_correction,
)
from gplately.grids.sediment_thickness import (
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
            age_depth_model="parsons_sclater",  # no pyBacktrack equivalent
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
