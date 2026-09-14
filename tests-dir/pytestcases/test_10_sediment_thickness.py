import numpy as np
import pytest
from conftest import logger

import gplately
from gplately.grids.sediment_thickness import (
    generate_distance_grids,
    generate_input_points_grid,
    generate_sediment_thickness_grids,
)

logger.info(__name__)


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


def test_public_api_exports():
    for name in (
        "generate_input_points_grid",
        "generate_distance_grids",
        "generate_sediment_thickness_grids",
        "simple_paleobathymetry",
    ):
        assert hasattr(gplately, name)
