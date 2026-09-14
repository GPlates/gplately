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

"""Dynamically-contoured passive margins through time.

A port of EarthByte's `continent-contouring
<https://github.com/EarthByte/continent-contouring>`__ ``create_passive_margins.py`` (the one
use case that repository contains -- see `gplately#446
<https://github.com/GPlates/gplately/issues/446>`__), built entirely on gplately's own
continent-contouring engine (:class:`gplately.ptt.continent_contours.ContinentContouring`).

For each requested time, reconstructed continental polygons are contoured into continents
(:class:`gplately.ptt.continent_contours.ContouredContinent` -- see that engine); each contour
is then split into **passive margin** segments by removing the parts that lie close to a
subduction zone (those are **active** margins). This is also what
`simple_paleobathymetry <https://github.com/EarthByte/simple_paleobathymetry>`__ uses as its
optional dynamically-contoured alternative to a static continent-ocean-boundary (COB)
line-segment file for :func:`gplately.generate_distance_grids`'s proximity target (see that
function's docstring, and *simple_paleobathymetry*'s README, "Step 2") -- unlike a static COB
file, contouring re-derives every paleo-margin at each time step, so it also captures passive
margins that existed in the past but not at present day.

* :func:`passive_margin_polylines` -- split one continent-contour polyline into its
  passive-margin segments (the core "use case" logic).
* :func:`generate_passive_margins` -- the full driver: contour continents through time and
  split each contour, using :class:`gplately.ptt.continent_contours.ContinentContouring`.
"""

import os

import pygplates

from ._grids import write_netcdf_grid
from ..ptt.continent_contours import ContinentContouring

__all__ = ["passive_margin_polylines", "generate_passive_margins"]


def passive_margin_polylines(
    contour_polyline, subduction_zone_lines, max_distance_radians
):
    """Split a continent-contour polyline into its passive-margin segments.

    A great-circle-arc segment of the contour is an **active margin** if it lies within
    `max_distance_radians` of any subduction-zone line; everything else is a **passive
    margin**. Consecutive passive-margin arcs are joined into single polylines.

    Parameters
    ----------
    contour_polyline : pygplates.PolylineOnSphere
        One contour of a contoured continent (see
        :meth:`gplately.ptt.continent_contours.ContouredContinent.get_contours`).
    subduction_zone_lines : sequence of pygplates.GeometryOnSphere
        Resolved subduction-zone geometries at the same reconstruction time.
    max_distance_radians : float
        A contour segment within this distance (radians) of a subduction zone is active;
        everything else is passive.

    Returns
    -------
    list of pygplates.PolylineOnSphere
        The passive-margin segments of `contour_polyline` (possibly empty, if the whole
        contour is an active margin).
    """
    points = list(contour_polyline.get_points())
    if len(points) < 2:
        return []

    def _near_subduction(point0, point1):
        if not subduction_zone_lines:
            return False
        segment = pygplates.PolylineOnSphere([point0, point1])
        for subduction_zone_line in subduction_zone_lines:
            if (
                pygplates.GeometryOnSphere.distance(
                    segment, subduction_zone_line, max_distance_radians
                )
                is not None
            ):
                return True
        return False

    # Continent contours are closed rings (points[0] == points[-1]). If we split them
    # starting from points[0] as-is, a passive-margin stretch that happens to straddle that
    # arbitrary start/end point would be cut into two separate output polylines instead of
    # one. Avoid this by rotating the ring to start right after an active edge (if any),
    # so the array boundary never falls in the middle of a passive stretch.
    if points[0] == points[-1] and len(points) > 2:
        unique_points = points[:-1]
        n = len(unique_points)
        for offset in range(n):
            if _near_subduction(unique_points[offset - 1], unique_points[offset]):
                rotated = unique_points[offset:] + unique_points[:offset]
                points = rotated + [rotated[0]]
                break
        # If no active edge was found, the whole ring is one passive margin; leave `points`
        # as-is (the loop below will still return it as a single closed polyline).

    margins = []
    run = [points[0]]
    for i in range(1, len(points)):
        if _near_subduction(points[i - 1], points[i]):
            if len(run) >= 2:
                margins.append(pygplates.PolylineOnSphere(run))
            run = [points[i]]
        else:
            run.append(points[i])
    if len(run) >= 2:
        margins.append(pygplates.PolylineOnSphere(run))
    return margins


def generate_passive_margins(
    rotation_model,
    continent_features,
    topological_features,
    times,
    point_spacing_degrees=0.25,
    area_threshold_square_kms=0.0,
    buffer_and_gap_distance_kms=0.0,
    exclusion_area_threshold_square_kms=800000.0,
    separation_distance_threshold_radians=None,
    max_distance_of_subduction_from_active_margin_kms=500.0,
    anchor_plate_id=0,
    time_step=1.0,
    output_directory=None,
):
    """Dynamically contour continents through time and split each contour into passive margins.

    Parameters
    ----------
    rotation_model : pygplates.RotationModel, or any argument accepted by its constructor
        The rotation model.
    continent_features : any argument accepted by pygplates.FeaturesFunctionArgument
        The continental polygons (or cratons) to contour -- *not* COB line segments or
        boundaries; see :class:`gplately.ptt.continent_contours.ContinentContouring`.
    topological_features : any argument accepted by pygplates.FeaturesFunctionArgument
        The topological plate boundary / network features, used to resolve subduction zones
        (which mark active margins) at each time.
    times : sequence of float
        Reconstruction times (Ma) to contour.
    point_spacing_degrees : float, default: 0.25
        Grid spacing (degrees) of the point grid used to contour/aggregate the continental
        polygons into continents.
    area_threshold_square_kms : float, default: 0.0
        Contoured continents smaller than this are excluded.
    buffer_and_gap_distance_kms : float, default: 0.0
        Expand continents ocean-ward by this distance (also closes small gaps between nearby
        continental polygons) before contouring.
    exclusion_area_threshold_square_kms : float, default: 800000.0
        Enclosed interior gaps (e.g. lakes) smaller than this are dropped (i.e. filled in as
        continental crust) rather than left as holes in the contour.
    separation_distance_threshold_radians : float, optional
        Continental polygons closer than this are merged into the same continent. Defaults to
        :data:`gplately.ptt.continent_contours.DEFAULT_CONTINENT_SEPARATION_DISTANCE_THRESHOLD_RADIANS`
        (a small numerical-tolerance value) if not given.
    max_distance_of_subduction_from_active_margin_kms : float, default: 500.0
        A contour segment within this distance of a subduction zone is an active margin;
        everything else is a passive margin (see :func:`passive_margin_polylines`).
    anchor_plate_id : int, default: 0
        Anchor plate used to build `rotation_model` (ignored if `rotation_model` is already a
        :class:`pygplates.RotationModel`).
    time_step : float, default: 1.0
        Used only to set each output feature's valid-time window (``time +/- 0.5 *
        time_step``), so a feature is only "visible" at its own reconstruction time when
        loaded back into GPlates. Should match the spacing of `times`, if they're evenly
        spaced.
    output_directory : str, optional
        If given, writes ``continent_contour_features.gpmlz`` and
        ``passive_margin_features.gpmlz`` (the aggregated feature collections, across all
        `times`), plus one ``continent_mask_<time>.nc`` per time (a global grid, 1.0 where
        continental crust, 0.0 elsewhere, at `point_spacing_degrees` resolution).

    Returns
    -------
    dict
        ``{"continent_contour_features": pygplates.FeatureCollection,
        "passive_margin_features": pygplates.FeatureCollection,
        "continent_masks": {time: numpy.ndarray}}`` -- the continent-ocean-boundary contours
        (every contour around continental crust, both passive and active margins), the
        passive-margin subset of those contours, and each time's continental-crust mask.
    """
    earth_radius_km = pygplates.Earth.mean_radius_in_kms

    rotation_model = pygplates.RotationModel(
        rotation_model, default_anchor_plate_id=anchor_plate_id
    )
    continent_feature_list = pygplates.FeaturesFunctionArgument(
        continent_features
    ).get_features()
    topological_feature_list = pygplates.FeaturesFunctionArgument(
        topological_features
    ).get_features()

    contourer_kwargs = {}
    if separation_distance_threshold_radians is not None:
        contourer_kwargs["continent_separation_distance_threshold_radians"] = (
            separation_distance_threshold_radians
        )

    contourer = ContinentContouring(
        rotation_model,
        continent_feature_list,
        continent_contouring_point_spacing_degrees=point_spacing_degrees,
        continent_contouring_area_threshold_steradians=(
            area_threshold_square_kms / (earth_radius_km * earth_radius_km)
        ),
        continent_contouring_buffer_and_gap_distance_radians=(
            buffer_and_gap_distance_kms / earth_radius_km
        ),
        continent_exclusion_area_threshold_steradians=(
            exclusion_area_threshold_square_kms / (earth_radius_km * earth_radius_km)
        ),
        **contourer_kwargs,
    )
    max_distance_radians = (
        max_distance_of_subduction_from_active_margin_kms / earth_radius_km
    )

    if output_directory:
        os.makedirs(output_directory, exist_ok=True)

    contour_features = []
    passive_margin_features = []
    continent_masks = {}

    for time in times:
        time = float(time)

        subduction_zone_lines = []
        resolved_topologies = []
        shared_boundary_sections = []
        pygplates.resolve_topologies(
            topological_feature_list,
            rotation_model,
            resolved_topologies,
            time,
            shared_boundary_sections,
        )
        for shared_boundary_section in shared_boundary_sections:
            if (
                shared_boundary_section.get_feature().get_feature_type()
                == pygplates.FeatureType.gpml_subduction_zone
            ):
                for (
                    shared_sub_segment
                ) in shared_boundary_section.get_shared_sub_segments():
                    subduction_zone_lines.append(
                        shared_sub_segment.get_resolved_geometry()
                    )

        continent_mask, contoured_continents = (
            contourer.get_continent_mask_and_contoured_continents(time)
        )
        continent_masks[time] = continent_mask

        for contoured_continent in contoured_continents:
            for contour in contoured_continent.get_contours():
                contour_feature = pygplates.Feature()
                contour_feature.set_geometry(contour)
                contour_feature.set_valid_time(
                    time + 0.5 * time_step, time - 0.5 * time_step
                )
                contour_features.append(contour_feature)

                for margin in passive_margin_polylines(
                    contour, subduction_zone_lines, max_distance_radians
                ):
                    margin_feature = pygplates.Feature()
                    margin_feature.set_geometry(margin)
                    margin_feature.set_valid_time(
                        time + 0.5 * time_step, time - 0.5 * time_step
                    )
                    passive_margin_features.append(margin_feature)

        if output_directory:
            mask_path = os.path.join(
                output_directory, "continent_mask_{:.0f}.nc".format(time)
            )
            write_netcdf_grid(mask_path, continent_mask.astype("float64"))

    contour_feature_collection = pygplates.FeatureCollection(contour_features)
    passive_margin_feature_collection = pygplates.FeatureCollection(
        passive_margin_features
    )

    if output_directory:
        contour_feature_collection.write(
            os.path.join(output_directory, "continent_contour_features.gpmlz")
        )
        passive_margin_feature_collection.write(
            os.path.join(output_directory, "passive_margin_features.gpmlz")
        )

    return {
        "continent_contour_features": contour_feature_collection,
        "passive_margin_features": passive_margin_feature_collection,
        "continent_masks": continent_masks,
    }
