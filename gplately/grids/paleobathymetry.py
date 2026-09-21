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

"""Paleobathymetry of oceanic crust, reconstructed from seafloor age and sediment loading.

This is a port of the physics in EarthByte's `simple_paleobathymetry
<https://github.com/EarthByte/simple_paleobathymetry>`__ workflow (see its README for the
full derivation and references). Paleobathymetry is computed as::

    paleobathymetry = basement_depth + sediment_thickness - isostatic_correction

* :func:`age_to_basement_depth` -- seafloor age -> depth to igneous basement (thermal
  subsidence; the *simple_paleobathymetry* workflow's Step 1).
* :func:`dutkiewicz_2017_sediment_thickness` -- seafloor age + distance to the nearest
  passive continental margin -> predicted compacted sediment thickness (Dutkiewicz et al.,
  2017; Step 3).
* :func:`sediment_isostatic_correction` / :func:`paleobathymetry` -- combine basement depth
  and sediment thickness into paleobathymetry, correcting for sediment loading (Sykes, 1996;
  Step 4).
* :func:`simple_paleobathymetry` -- run Steps 1-4 end to end, calling
  :mod:`gplately.grids.sediment_thickness` (`gplately#445 <https://github.com/GPlates/gplately/issues/445>`__)
  for Step 2 (each ocean point's lifetime-mean distance to the nearest passive continental
  margin), and optionally Step 5, calling
  :mod:`gplately.grids.pybacktrack_paleobathymetry` (`gplately#447
  <https://github.com/GPlates/gplately/issues/447>`__) to merge in pyBacktrack's present-day
  paleobathymetry (requires the optional ``pybacktrack`` package).

"""

import logging
import math
import os
from importlib.resources import files

import numpy as np
import pygplates

from ._grids import read_netcdf_grid, sample_grid, write_netcdf_grid
from ._utils import (
    DEFAULT_DECIMAL_PLACES_IN_TIME,
    check_times_are_distinct_in_filenames,
    format_time_in_filename,
    resolve_decimal_places_in_time,
    uniform_time_step,
)

logger = logging.getLogger("gplately")

__all__ = [
    "AGE_DEPTH_MODELS",
    "DUTKIEWICZ_2017_SEDIMENT_THICKNESS",
    "age_to_basement_depth",
    "dutkiewicz_2017_sediment_thickness",
    "sediment_isostatic_correction",
    "paleobathymetry",
    "simple_paleobathymetry",
]


# Four published thermal-subsidence (age -> basement depth) models, keyed by their canonical
# name plus the alternative spellings accepted by `age_to_basement_depth`'s `model` argument.
AGE_DEPTH_MODELS = ("gdh1", "rhcw18", "parsons_sclater", "crosby09")

_AGE_DEPTH_MODEL_ALIASES = {
    "gdh1": "gdh1",
    "ghd1": "gdh1",
    "stein_stein": "gdh1",
    "steinstein": "gdh1",
    "rhcw18": "rhcw18",
    "rchw18": "rhcw18",
    "richards": "rhcw18",
    "richards18": "rhcw18",
    "r18": "rhcw18",
    "parsons_sclater": "parsons_sclater",
    "ps": "parsons_sclater",
    "ps_tbl": "parsons_sclater",
    "parsons": "parsons_sclater",
    "crosby09": "crosby09",
    "crosby": "crosby09",
    "crosby_2009": "crosby09",
    "cmk09": "crosby09",
}

# The Richards, Hoggard, Cowton & White (2018) age-depth relationship has no closed-form
# expression, so it ships as a lookup table (age in Ma, depth in m) that we interpolate.
_DEFAULT_RICHARDS_TABLE_FILENAME = str(
    files("gplately").joinpath("data", "RHCW18_age_depth.dat")
)
_richards_table_cache = {}


def _load_richards_table(richards_table_filename=None):
    path = str(richards_table_filename or _DEFAULT_RICHARDS_TABLE_FILENAME)
    if path not in _richards_table_cache:
        table = np.loadtxt(path)
        order = np.argsort(table[:, 0])
        _richards_table_cache[path] = (table[order, 0], table[order, 1])
    return _richards_table_cache[path]


def age_to_basement_depth(age, model="gdh1", richards_table_filename=None):
    """Convert seafloor age to depth-to-basement using a thermal subsidence model.

    Depth to basement is the seafloor depth *before* any sediment is added on top: as oceanic
    lithosphere moves away from a mid-ocean ridge it cools, the mantle lithosphere thickens and
    becomes denser, and the column subsides isostatically.

    Parameters
    ----------
    age : array_like
        Seafloor age in Ma. NaN (no age, e.g. continental crust) stays NaN in the output.
    model : str, default: "gdh1"
        One of ``"gdh1"`` (Stein & Stein, 1992), ``"rhcw18"`` (Richards, Hoggard, Cowton &
        White, 2018), ``"parsons_sclater"`` (Parsons & Sclater, 1977) or ``"crosby09"``
        (Crosby, McKenzie & Sclater, 2009); see :data:`AGE_DEPTH_MODELS`. Common alternative
        spellings are also accepted.
    richards_table_filename : str, optional
        Path to the two-column (age Ma, depth m) age-depth lookup table used by the
        ``"rhcw18"`` model. Defaults to the table shipped with GPlately (the RHCW18 preferred
        parameters from the `RHCW18_Plate_Model
        <https://github.com/freddrichards/RHCW18_Plate_Model>`__ repository).

    Returns
    -------
    numpy.ndarray
        Depth to basement in metres, negative downwards, NaN where ``age`` is NaN.

    References
    ----------
    * Stein, C.A. & Stein, S. (1992). A model for the global variation in oceanic depth and
      heat flow with lithospheric age. *Nature*, 359, 123-129.
    * Parsons, B. & Sclater, J.G. (1977). An analysis of the variation of ocean floor
      bathymetry and heat flow with age. *JGR*, 82, 803-827.
    * Crosby, A.G., McKenzie, D. & Sclater, J.G. (2006); Crosby, A.G. & McKenzie, D. (2009).
      The relationship between depth, age and gravity in the oceans; and an analysis of young
      ocean depth, gravity and global residual topography. *GJI*.
    * Richards, F.D., Hoggard, M.J., Cowton, L.R. & White, N.J. (2018). Reassessing the
      thermal structure of oceanic lithosphere with revised global inventories of basement
      depths and heat flow measurements. *JGR: Solid Earth*, 123, 9136-9161.
    """
    age = np.asarray(age, dtype="float64")
    depth = np.full(age.shape, np.nan, dtype="float64")
    valid = ~np.isnan(age)

    key = _AGE_DEPTH_MODEL_ALIASES.get(str(model).strip().lower())
    if key is None:
        raise ValueError(
            f"Unknown age-depth model {model!r}. Choose one of: "
            + ", ".join(AGE_DEPTH_MODELS)
        )

    if key == "gdh1":
        neg = valid & (age < 0)
        young = valid & (age >= 0) & (age <= 20)
        old = valid & (age > 20)
        depth[neg] = -2600.0
        depth[young] = -(2600.0 + 365.0 * np.sqrt(age[young]))
        depth[old] = -(5651.0 - 2473.0 * np.exp(-0.0278 * age[old]))

    elif key == "rhcw18":
        ages_tbl, depths_tbl = _load_richards_table(richards_table_filename)
        # Ages older than the table's range are held at its deepest value (np.interp's
        # default `right`); depths are negated to match this module's sign convention.
        # Tiny negative ages (floating-point/grid-interpolation noise right at a ridge) are
        # clamped to 0 first, so they get the table's actual ridge-crest depth rather than
        # np.interp's `left` fallback (which would otherwise be taken literally as 0 m --
        # sea level -- like the other three models, this model's zero-age depth should be
        # ~-2500 m, not 0 m).
        clamped_age = np.clip(age[valid], 0.0, None)
        depth[valid] = -1.0 * np.interp(clamped_age, ages_tbl, depths_tbl)

    elif key == "parsons_sclater":
        neg = valid & (age < 0)
        young = valid & (age >= 0) & (age < 70)
        old = valid & (age >= 70)
        depth[neg] = -2500.0
        depth[young] = -(2500.0 + 350.0 * np.sqrt(age[young]))
        depth[old] = -(6400.0 - 3200.0 * np.exp(-age[old] / 62.8))

    else:  # crosby09
        young = valid & (age >= 0) & (age <= 75)
        mid = valid & (age > 75) & (age <= 160)
        old = valid & (age > 160)
        depth[young] = -(2652.0 + 324.0 * np.sqrt(age[young]))
        depth[mid] = -(
            5028.0 + 5.26 * age[mid] - 250.0 * np.sin((age[mid] - 75.0) / 30.0)
        )
        depth[old] = -5750.0

    return depth


# Dutkiewicz, A., Muller, R.D., Wang, X., O'Callaghan, S., Cannon, J. & Wright, N.M. (2017).
# Predicting sediment thickness on vanished ocean crust since 200 Ma. Geochemistry,
# Geophysics, Geosystems, 18, 4586-4603.
#
# The fitted degree-3 polynomial relates the *logarithm* of compacted sediment thickness to
# standardised seafloor age and standardised distance to the nearest passive continental
# margin. These are the published constants (also the gplately Zahirovic2022 defaults used by
# EarthByte's predicting-sediment-thickness and simple_paleobathymetry workflows); distances
# are in **kilometres**, matching the 0-3000 km calibrated range these constants were fitted
# over (simple_paleobathymetry's README describes the same constants as metres, which appears
# to be a documentation error -- see the discussion on gplately#444).
DUTKIEWICZ_2017_SEDIMENT_THICKNESS = {
    "mean_age": 61.18406823,
    "mean_distance_km": 1835.28118479,
    "variance_age": 1934.6999014,
    "variance_distance_km2": 1207521.8995806,
    "max_age": 191.87276,
    "max_distance_km": 3000.0,
    "polynomial_coefficients": (
        5.441401190368497,
        0.46893096,
        -0.07320928,
        -0.24077496,
        -0.10840657,
        0.00381672,
        0.06831728,
        0.01179914,
        0.01158149,
        -0.39880562,
    ),
}


def dutkiewicz_2017_sediment_thickness(
    age,
    distance_to_margin_km,
    mean_age=DUTKIEWICZ_2017_SEDIMENT_THICKNESS["mean_age"],
    mean_distance_km=DUTKIEWICZ_2017_SEDIMENT_THICKNESS["mean_distance_km"],
    variance_age=DUTKIEWICZ_2017_SEDIMENT_THICKNESS["variance_age"],
    variance_distance_km2=DUTKIEWICZ_2017_SEDIMENT_THICKNESS["variance_distance_km2"],
    polynomial_coefficients=DUTKIEWICZ_2017_SEDIMENT_THICKNESS[
        "polynomial_coefficients"
    ],
    max_age=DUTKIEWICZ_2017_SEDIMENT_THICKNESS["max_age"],
    max_distance_km=DUTKIEWICZ_2017_SEDIMENT_THICKNESS["max_distance_km"],
):
    """Predict compacted sediment thickness from seafloor age and distance to a passive margin.

    Implements the Dutkiewicz et al. (2017) relationship: a degree-3 polynomial, fitted to a
    global compilation of sediment-thickness measurements, relating the logarithm of compacted
    sediment thickness to standardised seafloor age and standardised distance to the nearest
    passive continental margin (older, nearer-margin crust has thicker sediment cover).

    Parameters
    ----------
    age : array_like
        Seafloor age in Ma.
    distance_to_margin_km : array_like
        Distance to the nearest passive continental margin, in kilometres -- typically the
        *lifetime-mean* distance of each point of ocean floor to the nearest passive margin,
        as it drifted from the ridge to its current position (this module does not compute
        that distance grid; see the module docstring).
    mean_age, mean_distance_km, variance_age, variance_distance_km2, polynomial_coefficients : float or sequence
        Standardisation and polynomial-fit constants. Default to the published Dutkiewicz et
        al. (2017) values (:data:`DUTKIEWICZ_2017_SEDIMENT_THICKNESS`); override to use a
        relationship retrained on a different present-day age grid or sediment-thickness
        dataset.
    max_age, max_distance_km : float or None
        Values above these are clamped before standardising, since the fit is only meaningful
        within its calibrated range. ``None`` disables clamping for that predictor.

    Returns
    -------
    numpy.ndarray
        Predicted compacted sediment thickness in metres (positive). NaN in, NaN out.

    References
    ----------
    Dutkiewicz, A., Muller, R.D., Wang, X., O'Callaghan, S., Cannon, J. & Wright, N.M. (2017).
    Predicting sediment thickness on vanished ocean crust since 200 Ma. *Geochemistry,
    Geophysics, Geosystems*, 18, 4586-4603.
    """
    age = np.asarray(age, dtype="float64")
    distance_km = np.asarray(distance_to_margin_km, dtype="float64")

    if max_age is not None:
        age = np.where(age > max_age, max_age, age)
    if max_distance_km is not None:
        distance_km = np.where(
            distance_km > max_distance_km, max_distance_km, distance_km
        )

    a = (age - mean_age) / math.sqrt(variance_age)
    d = (distance_km - mean_distance_km) / math.sqrt(variance_distance_km2)

    c = polynomial_coefficients
    log_sediment_thickness = (
        c[0] * 1.0
        + c[1] * a
        + c[2] * d
        + c[3] * a * a
        + c[4] * a * d
        + c[5] * d * d
        + c[6] * a * a * a
        + c[7] * a * a * d
        + c[8] * a * d * d
        + c[9] * d * d * d
    )
    return np.exp(log_sediment_thickness)


def sediment_isostatic_correction(sediment_thickness_m):
    """Isostatic correction (Sykes, 1996) for the seafloor depression caused by sediment load.

    Adding a sediment pile on top of basement does not raise the seafloor by the full pile
    thickness, because the added weight pushes the crust down. This is the correction to
    subtract.

    Parameters
    ----------
    sediment_thickness_m : array_like
        Sediment thickness in metres (positive).

    Returns
    -------
    numpy.ndarray
        Isostatic correction in metres (positive; clamped to >= 0), NaN where
        ``sediment_thickness_m`` is NaN.

    References
    ----------
    Sykes, T.J.S. (1996). A correction for sediment load upon the ocean floor: uniform versus
    varying sediment density estimations. *Marine Geology*, 133, 35-49.
    """
    h_km = np.asarray(sediment_thickness_m, dtype="float64") / 1000.0
    correction_m = (0.43422 * h_km - 0.010395 * h_km * h_km) * 1000.0
    # np.clip (unlike `np.where(correction_m >= 0, ...)`) propagates NaN rather than silently
    # turning masked (non-ocean) NaN cells into 0.
    return np.clip(correction_m, 0.0, None)


def paleobathymetry(basement_depth_m, sediment_thickness_m):
    """Combine basement depth and sediment thickness into isostatically-compensated paleobathymetry.

    ``paleobathymetry = basement_depth + sediment_thickness - isostatic_correction`` (see
    :func:`sediment_isostatic_correction`). This describes oceanic crust only -- crust whose
    age is defined in the seafloor-age grid that ``basement_depth_m`` was derived from (see
    :func:`age_to_basement_depth`); it does not cover submerged continental crust or crust that
    has since subducted (see the module docstring re: Step 5 / pyBacktrack).

    Parameters
    ----------
    basement_depth_m : array_like
        Depth to igneous basement in metres, negative downwards (:func:`age_to_basement_depth`).
    sediment_thickness_m : array_like
        Sediment thickness in metres, positive (e.g. :func:`dutkiewicz_2017_sediment_thickness`).

    Returns
    -------
    numpy.ndarray
        Paleobathymetry (depth of the sediment-covered ocean floor) in metres, negative
        downwards.
    """
    basement_depth_m = np.asarray(basement_depth_m, dtype="float64")
    sediment_thickness_m = np.asarray(sediment_thickness_m, dtype="float64")
    correction_m = sediment_isostatic_correction(sediment_thickness_m)
    return basement_depth_m + sediment_thickness_m - correction_m


def simple_paleobathymetry(
    rotation_model,
    proximity_features,
    topological_features,
    age_grid_filenames_and_times,
    age_depth_model="gdh1",
    grid_spacing=0.5,
    time_increment=1,
    max_reconstruction_time=None,
    anchor_plate_id=0,
    clamp_distance_km=3000.0,
    richards_table_filename=None,
    output_directory=None,
    sediment_thickness_kwargs=None,
    continent_obstacle_features=None,
    plate_boundary_obstacle_feature_types=None,
    shortest_path_grid_subdivision_depth=6,
    pybacktrack=False,
    static_polygon_filename=None,
    present_day_age_grid_filename=None,
    pybacktrack_kwargs=None,
    *,
    decimal_places_in_time=None,
):
    """Run the full *simple_paleobathymetry* workflow (Steps 1-4, optionally 5) end to end.

    For each given time: Step 2 reconstructs ocean points backward through time to compute
    their lifetime-mean distance to the nearest passive continental margin
    (:func:`gplately.generate_distance_grids`); Step 3 predicts sediment
    thickness from that distance and seafloor age
    (:func:`gplately.generate_sediment_thickness_grids`, which calls
    :func:`dutkiewicz_2017_sediment_thickness`); Step 1 converts seafloor age to basement depth
    (:func:`age_to_basement_depth`); and Step 4 combines them into paleobathymetry
    (:func:`paleobathymetry`). If `pybacktrack` is true, Step 5 additionally merges in
    pyBacktrack's present-day paleobathymetry (see
    :func:`gplately.grids.pybacktrack_paleobathymetry.merge_pybacktrack_paleobathymetry`) to
    also cover submerged continental crust and crust that has since subducted.

    Parameters
    ----------
    rotation_model, proximity_features, topological_features, age_grid_filenames_and_times, grid_spacing, time_increment, max_reconstruction_time, anchor_plate_id, clamp_distance_km
        Passed to :func:`gplately.generate_distance_grids`; see its
        docstring. `proximity_features` should be passive-margin continent-ocean-boundary line
        segments (not polygons -- see the *simple_paleobathymetry* README's Step 2 for why).
    age_depth_model : str, default: "gdh1"
        Passed to :func:`age_to_basement_depth` as `model`.
    richards_table_filename : str, optional
        Passed to :func:`age_to_basement_depth` (only used when `age_depth_model` is
        ``"rhcw18"``). Forwarded to Step 5 as well, so both sides use the same table.
    output_directory : str, optional
        If given, intermediate distance/sediment-thickness grids and the final paleobathymetry
        grids (``paleobathymetry_<time>Ma.nc``) are all written under this directory (in
        ``Distances/``, ``SedimentThickness/`` and directly here, respectively). **Required**
        if `pybacktrack` is true, since Step 5 merges the Step 4 grids by reading them back off
        disk.
    decimal_places_in_time : int, optional
        Decimal places of the reconstruction time in the filenames this writes, and in the
        names Step 5 uses to find the Step 4 grids again. Defaults to 0 for the
        sediment-thickness and paleobathymetry grids and 1 for the intermediate distance
        grids, reproducing the filenames of the workflows each was ported from -- an
        explicit value applies to all of them. Times that are not
        distinct at this resolution would overwrite each other, so a clash raises
        `ValueError` rather than silently discarding grids -- raise this value when using a
        fractional time step. It is the same rule as pyBacktrack's
        ``output_file_decimal_places_in_time``.
    sediment_thickness_kwargs : dict, optional
        Extra keyword arguments passed to :func:`dutkiewicz_2017_sediment_thickness` (via
        :func:`gplately.generate_sediment_thickness_grids`), e.g. to
        override the default Dutkiewicz et al. (2017) constants.
    continent_obstacle_features, plate_boundary_obstacle_feature_types, shortest_path_grid_subdivision_depth
        Passed to :func:`gplately.generate_distance_grids` (Step 2); if `continent_obstacle_features`
        is given, distances are routed *around* continents instead of a great-circle straight
        line -- see that function's docstring.
    pybacktrack : bool, default: False
        If true, additionally run Step 5 (requires the optional `pybacktrack` package, plus
        `output_directory`, `static_polygon_filename` and `present_day_age_grid_filename`).
        pyBacktrack generates its output at a single increment rather than at a list of
        times, so that increment is taken from the spacing of the times in
        `age_grid_filenames_and_times`, which must therefore be evenly spaced.
    static_polygon_filename : str, optional
        Static polygons for Step 5 (pyBacktrack uses these to assign plate IDs); required if
        `pybacktrack` is true.
    present_day_age_grid_filename : str, optional
        The seafloor-age grid at 0 Ma, for Step 5 (regardless of the times being computed --
        pyBacktrack backtracks from the present day); required if `pybacktrack` is true.
    pybacktrack_kwargs : dict, optional
        Extra keyword arguments passed to
        :func:`gplately.grids.pybacktrack_paleobathymetry.merge_pybacktrack_paleobathymetry`.

    Returns
    -------
    dict
        Maps each time (as given in `age_grid_filenames_and_times`) to a ``(lon, lat, grid)``
        tuple: 1-D longitude/latitude coordinate arrays and a 2-D ``(lat, lon)`` array of
        paleobathymetry in metres, negative downwards. Reflects Step 4's output even when
        `pybacktrack` is true (Step 5's merged grids are written to
        ``<output_directory>/PaleobathymetryPyBacktrack/`` but not re-loaded into the returned
        dict, to avoid a second, potentially large, disk read).
    """
    # Imported here, not at module level, to avoid a circular import: sediment_thickness
    # imports dutkiewicz_2017_sediment_thickness from this module.
    from .sediment_thickness import (
        generate_distance_grids,
        generate_sediment_thickness_grids,
    )

    distance_output_dir = (
        os.path.join(output_directory, "Distances") if output_directory else None
    )
    sediment_thickness_output_dir = (
        os.path.join(output_directory, "SedimentThickness")
        if output_directory
        else None
    )

    distance_grid_kwargs = {}
    if continent_obstacle_features is not None:
        distance_grid_kwargs["continent_obstacle_features"] = (
            continent_obstacle_features
        )
        distance_grid_kwargs["shortest_path_grid_subdivision_depth"] = (
            shortest_path_grid_subdivision_depth
        )
        if plate_boundary_obstacle_feature_types is not None:
            distance_grid_kwargs["plate_boundary_obstacle_feature_types"] = (
                plate_boundary_obstacle_feature_types
            )

    age_grid_filenames_and_times = list(age_grid_filenames_and_times)

    # Probed now, and the result thrown away, so that an unusable model name or an
    # unreadable lookup table is reported before Step 2 runs rather than from Step 4 at the
    # end of it. (The table is cached, so Step 4 does not pay for this twice.)
    age_to_basement_depth(
        0.0, model=age_depth_model, richards_table_filename=richards_table_filename
    )

    # Every Step 5 prerequisite is checked here rather than at the Step 5 call at the end
    # of this function, so that a request Step 5 cannot satisfy is rejected before Steps
    # 1-4 spend hours producing grids it will not be able to merge.
    pybacktrack_time_increment = None
    if pybacktrack:
        if not output_directory:
            raise ValueError(
                "pybacktrack=True requires output_directory (Step 5 merges the Step 4 "
                "grids by reading them back off disk)."
            )
        if not static_polygon_filename:
            raise ValueError("pybacktrack=True requires static_polygon_filename.")
        if not present_day_age_grid_filename:
            raise ValueError(
                "pybacktrack=True requires present_day_age_grid_filename (the age grid at 0 Ma)."
            )
        if isinstance(rotation_model, pygplates.RotationModel):
            raise TypeError(
                "pybacktrack=True requires rotation_model to be filename(s), not an "
                "already-constructed pygplates.RotationModel -- pyBacktrack builds its own "
                "rotation model internally and needs the raw file path(s)."
            )
        # Imported here for the same reason the Step 5 call below does: pybacktrack is an
        # optional dependency. This function does not need it, only the model name.
        # pyBacktrack generates output at one increment between two times, rather than at a
        # list of times, so the times it is asked for must be describable that way.
        pybacktrack_time_increment = uniform_time_step(
            [time for _, time in age_grid_filenames_and_times]
        )

    # Forward the caller's value unresolved to each step, so that when it is None each
    # output keeps the default that reproduces its own original workflow -- in particular
    # the distance grids carry one decimal place of time where these grids carry none.
    # Resolving it here instead would rename the distance grids that the
    # 'generate-sediment-grids' subcommand later goes looking for.
    requested_decimal_places_in_time = decimal_places_in_time
    decimal_places_in_time = resolve_decimal_places_in_time(
        decimal_places_in_time, DEFAULT_DECIMAL_PLACES_IN_TIME
    )
    if output_directory:
        check_times_are_distinct_in_filenames(
            [time for _, time in age_grid_filenames_and_times],
            decimal_places_in_time,
            "paleobathymetry_{}Ma.nc",
        )

    distance_grids = generate_distance_grids(
        rotation_model=rotation_model,
        proximity_features=proximity_features,
        topological_features=topological_features,
        age_grid_filenames_and_times=age_grid_filenames_and_times,
        grid_spacing=grid_spacing,
        time_increment=time_increment,
        max_reconstruction_time=max_reconstruction_time,
        anchor_plate_id=anchor_plate_id,
        clamp_distance_km=clamp_distance_km,
        output_directory=distance_output_dir,
        decimal_places_in_time=requested_decimal_places_in_time,
        **distance_grid_kwargs,
    )
    # clamp_distance_km is the default for max_distance_km, but an explicit
    # sediment_thickness_kwargs["max_distance_km"] (a documented override) must win rather
    # than colliding with it as a duplicate keyword argument.
    # generate_distance_grids() skips an age grid it found nothing usable in, so the later
    # steps work from what it actually produced rather than from what was asked for.
    skipped = [
        filename
        for filename, time in age_grid_filenames_and_times
        if time not in distance_grids
    ]
    if skipped:
        if pybacktrack:
            raise ValueError(
                "Step 5 cannot run because Step 2 skipped "
                f"{len(skipped)} age grid(s) with no usable data ({skipped}): it merges the "
                "Step 4 grids by filename, and those times will have no grid. Remove those "
                "age grids from age_grid_filenames_and_times, or set pybacktrack=False."
            )
        logger.warning(
            "Continuing without %d age grid(s) that Step 2 skipped: %s",
            len(skipped),
            skipped,
        )
        age_grid_filenames_and_times = [
            (filename, time)
            for filename, time in age_grid_filenames_and_times
            if time in distance_grids
        ]

    sediment_thickness_grid_kwargs = {"max_distance_km": clamp_distance_km}
    sediment_thickness_grid_kwargs.update(sediment_thickness_kwargs or {})
    sediment_thickness_grids = generate_sediment_thickness_grids(
        age_grid_filenames_and_times,
        distance_grids,
        output_directory=sediment_thickness_output_dir,
        decimal_places_in_time=requested_decimal_places_in_time,
        **sediment_thickness_grid_kwargs,
    )

    if output_directory:
        os.makedirs(output_directory, exist_ok=True)

    results = {}
    for age_grid_filename, time in age_grid_filenames_and_times:
        lon, lat, sediment_thickness_m = sediment_thickness_grids[time]

        age_grid, grid_lon, grid_lat = read_netcdf_grid(
            age_grid_filename, return_grids=True
        )
        lon_2d, lat_2d = np.meshgrid(lon, lat)
        ages_on_output_grid = sample_grid(
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

        basement_depth_m = age_to_basement_depth(
            ages_on_output_grid,
            model=age_depth_model,
            richards_table_filename=richards_table_filename,
        )
        paleobathymetry_m = paleobathymetry(basement_depth_m, sediment_thickness_m)
        results[time] = (lon, lat, paleobathymetry_m)

        if output_directory:
            output_path = os.path.join(
                output_directory,
                "paleobathymetry_{}Ma.nc".format(
                    format_time_in_filename(time, decimal_places_in_time)
                ),
            )
            write_netcdf_grid(output_path, paleobathymetry_m)

    if pybacktrack:
        # Imported here, not at module level: pybacktrack is an optional dependency, and this
        # module (gplately.grids.pybacktrack_paleobathymetry) only imports it lazily too.
        from .pybacktrack_paleobathymetry import merge_pybacktrack_paleobathymetry

        times = [time for _, time in age_grid_filenames_and_times]
        pybacktrack_output_dir = os.path.join(
            output_directory, "PaleobathymetryPyBacktrack"
        )
        os.makedirs(pybacktrack_output_dir, exist_ok=True)

        merge_pybacktrack_paleobathymetry(
            output_file_prefix=os.path.join(
                pybacktrack_output_dir, "paleobathymetry_${time}Ma.nc"
            ),
            merge_paleobathymetry_filename_format=os.path.join(
                output_directory, "paleobathymetry_${time}Ma.nc"
            ),
            rotation_filenames=rotation_model,
            static_polygon_filename=static_polygon_filename,
            present_day_age_grid_filename=present_day_age_grid_filename,
            grid_spacing_degrees=grid_spacing,
            oldest_time=max(times),
            youngest_time=min(times),
            # pyBacktrack's time_increment is the increment its *output* is generated at,
            # between youngest_time and oldest_time -- not the increment Step 2 steps the
            # reconstruction by. Passing the latter makes pyBacktrack generate a different
            # set of times than Steps 1-4 did, and then look for Step 4 grids at times that
            # were never written.
            time_increment=pybacktrack_time_increment,
            age_depth_model=age_depth_model,
            richards_table_filename=richards_table_filename,
            anchor_plate_id=anchor_plate_id,
            decimal_places_in_time=requested_decimal_places_in_time,
            **(pybacktrack_kwargs or {}),
        )

    return results
