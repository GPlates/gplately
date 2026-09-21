#
#    Copyright (C) 2026 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

"""Distance-to-passive-margin and predicted sediment-thickness grids for ocean crust.

This is a port of the core engine in EarthByte's `predicting-sediment-thickness
<https://github.com/EarthByte/predicting-sediment-thickness>`__ (``ocean_basin_proximity.py``
and ``predict_sediment_thickness.py``), used by both that repository's own workflow and by
`simple_paleobathymetry <https://github.com/EarthByte/simple_paleobathymetry>`__'s Steps 2-3
(see `gplately#445 <https://github.com/GPlates/gplately/issues/445>`__ /
`gplately#444 <https://github.com/GPlates/gplately/issues/444>`__).

* :func:`generate_distance_grids` -- for each ocean point in a seafloor-age grid, reconstruct
  it backward through time (from the reconstruction time to its formation age at a mid-ocean
  ridge) and compute its *lifetime-mean* distance to the nearest "proximity" feature (typically
  passive-margin continent-ocean-boundary line segments).
* :func:`generate_sediment_thickness_grids` -- combine a seafloor-age grid with the distance
  grid above into a predicted sediment-thickness grid, via
  :func:`gplately.dutkiewicz_2017_sediment_thickness`.

`generate_distance_grids`'s `proximity_features` need not be a static COB line-segment file --
:func:`gplately.generate_passive_margins`'s dynamically-contoured
``passive_margin_features`` (`gplately#446
<https://github.com/GPlates/gplately/issues/446>`__) can be passed straight in instead, which
also captures passive margins that existed in the past but not at present day (see that
function's docstring).

`generate_distance_grids`'s `continent_obstacle_features` routes distances *around* continents
(rather than a great-circle straight line that may cut through land) -- see its docstring and
:mod:`gplately.lib.shortest_path`, a port of the original ``ocean_basin_proximity.py``'s
``shortest_path.py``.

**Not yet included here:** the "topological proximity features" mode of the original
``ocean_basin_proximity.py`` (measuring distance to resolved plate-boundary sections rather
than static/reconstructed features) -- only non-topological proximity features (e.g. COB line
segments, or :func:`gplately.generate_passive_margins`'s output) are supported.
"""

import logging
import math
import os

import numpy as np
import pygplates

from ._grids import read_netcdf_grid, sample_grid, write_netcdf_grid
from ._utils import (
    DEFAULT_DECIMAL_PLACES_IN_TIME,
    DEFAULT_DISTANCE_GRID_DECIMAL_PLACES_IN_TIME,
    check_time_increment_covers_times,
    check_times_are_distinct_in_filenames,
    time_index_at_or_after,
    distance_grid_filename,
    format_time_in_filename,
    resolve_decimal_places_in_time,
)
from .paleobathymetry import dutkiewicz_2017_sediment_thickness
from ..lib import shortest_path
from ..ptt.utils.proximity_query import find_closest_geometries_to_points

_DEFAULT_PLATE_BOUNDARY_OBSTACLE_FEATURE_TYPES = (
    pygplates.FeatureType.gpml_mid_ocean_ridge,
    pygplates.FeatureType.gpml_subduction_zone,
)

logger = logging.getLogger("gplately")

__all__ = [
    "generate_input_points_grid",
    "generate_distance_grids",
    "generate_sediment_thickness_grids",
]


def generate_input_points_grid(grid_spacing_degrees):
    """Generate a global, gridline-registered lon/lat point grid at the given spacing.

    Parameters
    ----------
    grid_spacing_degrees : float
        Spacing between points, in degrees. Must be positive.

    Returns
    -------
    lon_1d, lat_1d : numpy.ndarray
        1-D coordinate arrays, from -180 to 180 (inclusive of both endpoints when they are an
        integer multiple of the spacing away from the start) and -90 to 90 respectively.
    lon_flat, lat_flat : numpy.ndarray
        The corresponding flattened (lon, lat) coordinates of every point in the grid, in
        row-major (lat-major) order matching a ``(lat, lon)``-shaped 2-D grid.
    """
    if grid_spacing_degrees <= 0:
        raise ValueError("grid_spacing_degrees must be positive.")

    num_lon = int(math.floor(360.0 / grid_spacing_degrees)) + 1
    num_lat = int(math.floor(180.0 / grid_spacing_degrees)) + 1
    lon_1d = np.linspace(-180.0, 180.0, num_lon)
    lat_1d = np.linspace(-90.0, 90.0, num_lat)

    lon_grid, lat_grid = np.meshgrid(lon_1d, lat_1d)  # each shaped (num_lat, num_lon)
    return lon_1d, lat_1d, lon_grid.ravel(), lat_grid.ravel()


def _accumulate_mean_distance_for_age_grid(
    flat_indices,
    point_lons,
    point_lats,
    point_ages,
    age_grid_time,
    num_output_points,
    rotation_model,
    proximity_features,
    topological_features,
    time_increment,
    max_reconstruction_time,
    distance_threshold_radians,
    shortest_path_grid=None,
    continent_obstacle_features=None,
    plate_boundary_obstacle_feature_types=None,
):
    """Reconstruct one age grid's ocean points backward through time, accumulating each
    point's mean distance (km) to the nearest proximity feature over its lifetime.

    Returns a 1-D array of length `num_output_points` (mean distance, NaN where a point never
    contributed a sample -- e.g. it had no defined age).
    """
    # Absolute time (Ma) at which each output grid point's crust formed at the ridge,
    # indexed by its position in the full output grid (so it stays valid however the set of
    # "currently active" points shrinks below).
    formation_time = np.full(num_output_points, np.nan, dtype="float64")
    formation_time[flat_indices] = age_grid_time + np.asarray(
        point_ages, dtype="float64"
    )

    sum_distance_km = np.zeros(num_output_points, dtype="float64")
    num_distance = np.zeros(num_output_points, dtype="int64")

    current_points = [
        pygplates.PointOnSphere(float(lat), float(lon))
        for lon, lat in zip(point_lons, point_lats)
    ]
    current_indices = np.asarray(flat_indices, dtype="int64")

    time_index = time_index_at_or_after(age_grid_time, time_increment)
    while current_points:
        time = time_index * time_increment
        if max_reconstruction_time is not None and time > max_reconstruction_time:
            break

        # Built once per time step; used below for point-stepping, and (if routing around
        # continents) to also resolve plate-boundary obstacles at this time.
        topological_model = pygplates.TopologicalModel(
            topological_features, rotation_model
        )

        # Reconstruct the (non-topological) proximity features to the current time.
        reconstructed_feature_geometries = []
        pygplates.reconstruct(
            proximity_features, rotation_model, reconstructed_feature_geometries, time
        )
        proximity_geometries = [
            reconstructed_feature_geometry.get_reconstructed_geometry()
            for reconstructed_feature_geometry in reconstructed_feature_geometries
        ]

        if shortest_path_grid is not None:
            # Distance routed *around* continent obstacles, at this reconstruction time.
            obstacle_geometries = []
            reconstructed_obstacle_feature_geometries = []
            pygplates.reconstruct(
                continent_obstacle_features,
                rotation_model,
                reconstructed_obstacle_feature_geometries,
                time,
            )
            obstacle_geometries.extend(
                reconstructed_obstacle_feature_geometry.get_reconstructed_geometry()
                for reconstructed_obstacle_feature_geometry in reconstructed_obstacle_feature_geometries
            )
            if plate_boundary_obstacle_feature_types:
                for (
                    resolved_topological_section
                ) in topological_model.topological_snapshot(
                    time
                ).get_resolved_topological_sections():
                    if (
                        resolved_topological_section.get_feature().get_feature_type()
                        not in plate_boundary_obstacle_feature_types
                    ):
                        continue
                    for (
                        shared_sub_segment
                    ) in resolved_topological_section.get_shared_sub_segments():
                        obstacle_geometries.append(
                            shared_sub_segment.get_resolved_geometry()
                        )

            obstacle_grid = shortest_path_grid.create_obstacle_grid(obstacle_geometries)
            shortest_path_distance_grid = obstacle_grid.create_distance_grid(
                proximity_geometries, distance_threshold_radians
            )

            for point_index, point in enumerate(current_points):
                distance_radians = shortest_path_distance_grid.shortest_distance(point)
                if distance_radians is None:
                    distance_radians = math.pi
                flat_index = current_indices[point_index]
                sum_distance_km[flat_index] += (
                    distance_radians * pygplates.Earth.mean_radius_in_kms
                )
                num_distance[flat_index] += 1
        else:
            # Distance (great circle) from each currently-active point to the nearest proximity
            # geometry, at this reconstruction time.
            closest_geometries = find_closest_geometries_to_points(
                current_points,
                proximity_geometries,
                distance_threshold_radians=distance_threshold_radians,
            )
            for point_index, closest_geometry in enumerate(closest_geometries):
                distance_radians = (
                    closest_geometry[0] if closest_geometry is not None else math.pi
                )
                flat_index = current_indices[point_index]
                sum_distance_km[flat_index] += (
                    distance_radians * pygplates.Earth.mean_radius_in_kms
                )
                num_distance[flat_index] += 1

        # Step the still-active points back one time increment (younger to older); drop any
        # that the topological model deactivates, or that reach their formation time.
        reconstructed_time_span = topological_model.reconstruct_geometry(
            current_points,
            initial_time=time,
            oldest_time=time + time_increment,
            youngest_time=time,
            time_increment=time_increment,
            deactivate_points=None,
        )
        reconstructed_points = reconstructed_time_span.get_geometry_points(
            time + time_increment, return_inactive_points=True
        )

        next_points = []
        next_indices = []
        if reconstructed_points:
            for point_index, point in enumerate(reconstructed_points):
                if point is None:
                    continue
                flat_index = current_indices[point_index]
                if time + time_increment > formation_time[flat_index]:
                    continue
                next_points.append(point)
                next_indices.append(flat_index)

        current_points = next_points
        current_indices = np.asarray(next_indices, dtype="int64")
        time_index += 1

    with np.errstate(invalid="ignore"):
        mean_distance_km = np.where(
            num_distance > 0, sum_distance_km / np.maximum(num_distance, 1), np.nan
        )
    return mean_distance_km


def generate_distance_grids(
    rotation_model,
    proximity_features,
    topological_features,
    age_grid_filenames_and_times,
    grid_spacing=0.5,
    time_increment=1,
    max_reconstruction_time=None,
    anchor_plate_id=0,
    clamp_distance_km=None,
    proximity_feature_types=None,
    distance_threshold_radians=None,
    continent_obstacle_features=None,
    plate_boundary_obstacle_feature_types=_DEFAULT_PLATE_BOUNDARY_OBSTACLE_FEATURE_TYPES,
    shortest_path_grid_subdivision_depth=6,
    output_directory=None,
    *,
    decimal_places_in_time=None,
):
    """For each ocean point in each age grid, compute its lifetime-mean distance (km) to the
    nearest proximity feature (typically passive continental margins).

    Each ocean point is reconstructed backward through time -- from the reconstruction time to
    the moment its crust formed at a mid-ocean ridge (its age) -- and its distance to the
    nearest (reconstructed) proximity feature is sampled at every ``time_increment``. The
    returned/written grid is the mean of those samples.

    Parameters
    ----------
    rotation_model : pygplates.RotationModel, or any argument accepted by its constructor
        The rotation model.
    proximity_features : any argument accepted by pygplates.FeaturesFunctionArgument
        Non-topological features to measure distance to -- e.g. passive-margin
        continent-ocean-boundary line segments. Must be line/point/polygon geometries, not
        topological plate boundaries.
    topological_features : any argument accepted by pygplates.FeaturesFunctionArgument
        The topological plate boundary / network features used to reconstruct ocean points
        backward through time.
    age_grid_filenames_and_times : sequence of (str, float)
        ``(age_grid_filename, reconstruction_time)`` pairs. Each age grid gives the age (Ma) of
        the ocean crust at its reconstruction time; NaN (or masked) cells are treated as
        non-ocean and excluded.
    grid_spacing : float, default: 0.5
        Spacing (degrees) of the point grid that distances are computed/output on.
    time_increment : float, default: 1
        Time increment (Myr) used both for stepping the backward reconstruction and for
        sampling distance along each point's lifetime. This is independent of the spacing of
        the times in `age_grid_filenames_and_times`: a coarser set of output times does not
        mean coarser sampling. Every one of those times must however be a multiple of this
        increment, because each age grid's walk starts at its own time snapped up onto the
        shared grid of multiples -- a time that does not land on that grid would be
        reconstructed from the wrong starting time. Must be positive.
    max_reconstruction_time : float, optional
        Do not reconstruct any point older than this (Ma). ``None`` means points are
        reconstructed back to their formation age regardless of how old that is (limited only
        by the temporal extent of `topological_features`).
    anchor_plate_id : int, default: 0
        Anchor plate used to build `rotation_model` (ignored if `rotation_model` is already a
        :class:`pygplates.RotationModel`).
    clamp_distance_km : float, optional
        If given, mean distances above this are clamped to it (the sediment-thickness
        relationship in :func:`generate_sediment_thickness_grids` is only calibrated over a
        finite distance range; see its docstring).
    proximity_feature_types : sequence of str, optional
        If given, restrict `proximity_features` to these qualified feature type names (e.g.
        ``["gpml:PassiveContinentalBoundary"]``) before measuring distance.
    distance_threshold_radians : float, optional
        Reject/ignore proximities further than this (radians). ``None`` means no threshold.
    continent_obstacle_features : any argument accepted by pygplates.FeaturesFunctionArgument, optional
        If given, distances are routed *around* these (reconstructed) obstacle geometries
        (typically continent/coastline polygons) instead of being a great-circle straight line
        -- see :mod:`gplately.lib.shortest_path`. ``None`` (the default) uses plain great-circle
        distance, which can be misleading where land lies between a point and the nearest
        margin (e.g. across a narrow isthmus).
    plate_boundary_obstacle_feature_types : sequence of pygplates.FeatureType, optional
        Only used when `continent_obstacle_features` is given: resolved plate-boundary sections
        of these feature types are *also* treated as obstacles at each time (default: mid-ocean
        ridges and subduction zones, matching the original workflow). Pass an empty sequence to
        only use `continent_obstacle_features`.
    shortest_path_grid_subdivision_depth : int, default: 6
        Only used when `continent_obstacle_features` is given: subdivision depth of the
        :class:`gplately.lib.shortest_path.Grid` used to route around obstacles -- spacing is
        ``90 / 2^depth`` degrees (default: 1.40625 degrees, matching the original workflow).
        Higher values are more accurate but slower.
    output_directory : str, optional
        If given, write each time's grid to
        ``<output_directory>/mean_distance_<grid_spacing>d_<time>.nc``.
    decimal_places_in_time : int, optional
        Decimal places of the reconstruction time in the output filenames. Defaults to 1,
        reproducing the filenames of the workflow this was ported from. Times that are not
        distinct at this resolution would overwrite each other, so a clash raises
        `ValueError` rather than silently discarding grids -- raise this value when using a
        fractional time step. It is the same rule as pyBacktrack's
        ``output_file_decimal_places_in_time``.

    Returns
    -------
    dict
        Maps each reconstruction time (as given in `age_grid_filenames_and_times`) to a
        ``(lon, lat, grid)`` tuple: 1-D longitude/latitude coordinate arrays and a 2-D
        ``(lat, lon)`` array of lifetime-mean distance in kilometres (NaN where the age grid
        has no data).

    Raises
    ------
    ValueError
        If `time_increment` is not positive, if some time in `age_grid_filenames_and_times`
        is not a multiple of it, or if two of those times would be written to the same file.
    """
    # Everything cheap is validated up front, before a single file is opened or a single
    # point reconstructed: these are the mistakes that would otherwise surface as a wrong
    # number or a missing grid at the end of a long run.
    age_grid_filenames_and_times = list(age_grid_filenames_and_times)
    age_grid_times = [time for _, time in age_grid_filenames_and_times]
    check_time_increment_covers_times(age_grid_times, time_increment)
    decimal_places_in_time = resolve_decimal_places_in_time(
        decimal_places_in_time, DEFAULT_DISTANCE_GRID_DECIMAL_PLACES_IN_TIME
    )
    if output_directory:
        check_times_are_distinct_in_filenames(
            age_grid_times,
            decimal_places_in_time,
            "mean_distance_{:.1f}d_{{}}.nc".format(grid_spacing),
        )

    rotation_model = pygplates.RotationModel(
        rotation_model, default_anchor_plate_id=anchor_plate_id
    )

    proximity_feature_list = pygplates.FeaturesFunctionArgument(
        proximity_features
    ).get_features()
    if proximity_feature_types:
        feature_types = [
            pygplates.FeatureType.create_from_qualified_string(feature_type)
            for feature_type in proximity_feature_types
        ]
        proximity_feature_list = [
            feature
            for feature in proximity_feature_list
            if feature.get_feature_type() in feature_types
        ]

    topological_feature_list = pygplates.FeaturesFunctionArgument(
        topological_features
    ).get_features()

    shortest_path_grid = None
    continent_obstacle_feature_list = None
    if continent_obstacle_features is not None:
        shortest_path_grid = shortest_path.Grid(shortest_path_grid_subdivision_depth)
        continent_obstacle_feature_list = pygplates.FeaturesFunctionArgument(
            continent_obstacle_features
        ).get_features()
        plate_boundary_obstacle_feature_types = [
            (
                pygplates.FeatureType.create_from_qualified_string(feature_type)
                if isinstance(feature_type, str)
                else feature_type
            )
            for feature_type in (plate_boundary_obstacle_feature_types or [])
        ]

    lon_1d, lat_1d, lon_flat, lat_flat = generate_input_points_grid(grid_spacing)
    num_output_points = lon_flat.size

    if output_directory:
        os.makedirs(output_directory, exist_ok=True)

    results = {}
    for age_grid_filename, age_grid_time in age_grid_filenames_and_times:
        age_grid, grid_lon, grid_lat = read_netcdf_grid(
            age_grid_filename, return_grids=True
        )
        ages = sample_grid(
            lon_flat,
            lat_flat,
            age_grid,
            method="linear",
            # First/last coordinates, not min/max: a descending latitude axis must keep
            # its sign, since that is how sample_grid() knows the row order of the grid.
            extent=(
                float(grid_lon[0]),
                float(grid_lon[-1]),
                float(grid_lat[0]),
                float(grid_lat[-1]),
            ),
        )
        valid = np.isfinite(ages)

        if not np.any(valid):
            logger.warning(
                "All input points are outside the age grid: %s", age_grid_filename
            )
            mean_distance_km = np.full(num_output_points, np.nan, dtype="float64")
        else:
            mean_distance_km = _accumulate_mean_distance_for_age_grid(
                flat_indices=np.flatnonzero(valid),
                point_lons=lon_flat[valid],
                point_lats=lat_flat[valid],
                point_ages=ages[valid],
                age_grid_time=age_grid_time,
                num_output_points=num_output_points,
                rotation_model=rotation_model,
                proximity_features=proximity_feature_list,
                topological_features=topological_feature_list,
                time_increment=time_increment,
                max_reconstruction_time=max_reconstruction_time,
                distance_threshold_radians=distance_threshold_radians,
                shortest_path_grid=shortest_path_grid,
                continent_obstacle_features=continent_obstacle_feature_list,
                plate_boundary_obstacle_feature_types=plate_boundary_obstacle_feature_types,
            )
            if clamp_distance_km is not None:
                mean_distance_km = np.where(
                    mean_distance_km > clamp_distance_km,
                    clamp_distance_km,
                    mean_distance_km,
                )

        grid = mean_distance_km.reshape(lat_1d.size, lon_1d.size)
        results[age_grid_time] = (lon_1d, lat_1d, grid)

        if output_directory:
            output_path = os.path.join(
                output_directory,
                distance_grid_filename(
                    grid_spacing, age_grid_time, decimal_places_in_time
                ),
            )
            write_netcdf_grid(output_path, grid)

    return results


def generate_sediment_thickness_grids(
    age_grid_filenames_and_times,
    distance_grids,
    output_directory=None,
    *,
    decimal_places_in_time=None,
    **sediment_thickness_kwargs,
):
    """Predict compacted sediment thickness by combining seafloor-age and distance-to-margin grids.

    For each time, samples the age grid onto the distance grid's points and evaluates
    :func:`gplately.dutkiewicz_2017_sediment_thickness`.

    Parameters
    ----------
    age_grid_filenames_and_times : sequence of (str, float)
        ``(age_grid_filename, time)`` pairs, as passed to :func:`generate_distance_grids`.
    distance_grids : dict
        The return value of :func:`generate_distance_grids`: maps time to a
        ``(lon, lat, distance_km)`` tuple, on the grid that the output sediment-thickness grid
        will also be on.
    output_directory : str, optional
        If given, write each time's grid to
        ``<output_directory>/sediment_thickness_<time>Ma.nc``.
    decimal_places_in_time : int, optional
        Decimal places of the reconstruction time in the output filenames. Defaults to 0,
        reproducing the filenames of the workflow this was ported from. Times that are not
        distinct at this resolution would overwrite each other, so a clash raises
        `ValueError` rather than silently discarding grids -- raise this value when using a
        fractional time step. It is the same rule as pyBacktrack's
        ``output_file_decimal_places_in_time``.
    **sediment_thickness_kwargs
        Passed through to
        :func:`gplately.dutkiewicz_2017_sediment_thickness` (e.g. to override
        the default Dutkiewicz et al. 2017 constants).

    Returns
    -------
    dict
        Maps each time to a ``(lon, lat, grid)`` tuple: 1-D longitude/latitude coordinate
        arrays and a 2-D ``(lat, lon)`` array of predicted sediment thickness in metres.
    """
    # Materialised because the filename check below walks it before the main loop does.
    age_grid_filenames_and_times = list(age_grid_filenames_and_times)

    decimal_places_in_time = resolve_decimal_places_in_time(
        decimal_places_in_time, DEFAULT_DECIMAL_PLACES_IN_TIME
    )
    if output_directory:
        check_times_are_distinct_in_filenames(
            [time for _, time in age_grid_filenames_and_times],
            decimal_places_in_time,
            "sediment_thickness_{}Ma.nc",
        )
        os.makedirs(output_directory, exist_ok=True)

    results = {}
    for age_grid_filename, time in age_grid_filenames_and_times:
        if time not in distance_grids:
            raise KeyError(
                f"No distance grid for time {time} (age grid {age_grid_filename!r}); "
                "did you pass the same age_grid_filenames_and_times to generate_distance_grids?"
            )
        lon, lat, distance_km = distance_grids[time]

        age_grid, grid_lon, grid_lat = read_netcdf_grid(
            age_grid_filename, return_grids=True
        )
        lon_2d, lat_2d = np.meshgrid(lon, lat)
        ages = sample_grid(
            lon_2d,
            lat_2d,
            age_grid,
            method="linear",
            # First/last coordinates, not min/max: a descending latitude axis must keep
            # its sign, since that is how sample_grid() knows the row order of the grid.
            extent=(
                float(grid_lon[0]),
                float(grid_lon[-1]),
                float(grid_lat[0]),
                float(grid_lat[-1]),
            ),
        )

        sediment_thickness_m = dutkiewicz_2017_sediment_thickness(
            ages, distance_km, **sediment_thickness_kwargs
        )
        results[time] = (lon, lat, sediment_thickness_m)

        if output_directory:
            output_path = os.path.join(
                output_directory,
                "sediment_thickness_{}Ma.nc".format(
                    format_time_in_filename(time, decimal_places_in_time)
                ),
            )
            write_netcdf_grid(output_path, sediment_thickness_m)

    return results
