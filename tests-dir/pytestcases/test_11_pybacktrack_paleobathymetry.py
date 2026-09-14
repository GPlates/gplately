import numpy as np
import pytest
from conftest import logger

import gplately

pybacktrack = pytest.importorskip(
    "pybacktrack",
    reason="pybacktrack is an optional dependency (gplately[paleobathymetry])",
)

logger.info(__name__)


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


def test_public_api_export():
    assert hasattr(gplately, "merge_pybacktrack_paleobathymetry")
