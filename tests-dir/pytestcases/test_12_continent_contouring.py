import numpy as np
import pygplates
import pytest
from conftest import logger

import gplately
from gplately.grids.continent_contouring import (
    generate_passive_margins,
    passive_margin_polylines,
)

logger.info(__name__)


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


def test_public_api_exports():
    for name in ("generate_passive_margins", "passive_margin_polylines"):
        assert hasattr(gplately, name)
