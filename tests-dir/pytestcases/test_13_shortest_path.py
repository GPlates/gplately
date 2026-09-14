import math

import pygplates
import pytest
from conftest import logger

from gplately.lib import shortest_path

logger.info(__name__)


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
