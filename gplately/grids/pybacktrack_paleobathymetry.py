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

"""Merge pyBacktrack's paleobathymetry with :mod:`gplately.grids.paleobathymetry`'s output.

This is Step 5 of EarthByte's `simple_paleobathymetry
<https://github.com/EarthByte/simple_paleobathymetry>`__ workflow (`gplately#447
<https://github.com/GPlates/gplately/issues/447>`__). Steps 1-4
(:func:`gplately.grids.paleobathymetry.simple_paleobathymetry`) only reconstruct **oceanic**
crust that is defined in the seafloor-age grids; they leave submerged continental crust, and
ocean crust that has since subducted (and so no longer appears in a present-day age grid), as
NaN.

`pyBacktrack <https://github.com/EarthByte/pyBacktrack>`__ (Muller, Cannon, Williams &
Dutkiewicz, 2018) fills that gap: it backtracks/backstrips a uniform grid of synthetic drill
sites on *present-day* crust -- including submerged continental crust -- back through time. It
cannot generate paleobathymetry for crust that no longer exists today (already subducted), so
this module merges in the Steps 1-4 grids to fill those regions (pyBacktrack's own
``merge_paleo_bathymetry_filename_format`` support: its reconstructed values take precedence on
crust that still exists today, the merged-in grids fill in the rest).

This module requires the optional `pybacktrack <https://github.com/EarthByte/pyBacktrack>`__
package (not a gplately dependency -- ``pip install pybacktrack`` or
``conda install -c conda-forge pybacktrack``), imported lazily inside
:func:`merge_pybacktrack_paleobathymetry` so that importing gplately never requires it.
"""

from ._utils import (
    DEFAULT_DECIMAL_PLACES_IN_TIME,
    resolve_decimal_places_in_time,
)

# Map this module's/`gplately.grids.paleobathymetry`'s `age_depth_model` names onto pyBacktrack's
# equivalent ocean age -> depth model constant (see gplately.grids.paleobathymetry.AGE_DEPTH_MODELS).
# Only models where the two implementations genuinely agree belong here; the others, and why
# they are excluded, are in _PYBACKTRACK_UNSUPPORTED_AGE_DEPTH_MODELS below.
_PYBACKTRACK_AGE_DEPTH_MODEL_ATTRS = {
    "gdh1": "AGE_TO_DEPTH_MODEL_GDH1",
    "rhcw18": "AGE_TO_DEPTH_MODEL_RHCW18",
}

# Models gplately offers for Steps 1-4 that Step 5 will not accept, and why. Spelling the
# reason out matters for crosby09 in particular, where a superficially matching pyBacktrack
# constant exists and mapping to it looks obviously right.
_PYBACKTRACK_UNSUPPORTED_AGE_DEPTH_MODELS = {
    "parsons_sclater": "pyBacktrack has no equivalent of this model.",
    "crosby09": (
        "pyBacktrack's AGE_TO_DEPTH_MODEL_CROSBY_2007 is a different model despite the "
        "similar name: it is the plate-cooling model of Crosby's 2007 thesis, whereas "
        "gplately's 'crosby09' is the empirical piecewise fit of Crosby & McKenzie (2009). "
        "They differ by up to ~218 m over 0-200 Ma (~51 m at the ridge crest), so Steps 1-4 "
        "computed with one and Step 5 with the other leave a step change at the merge "
        "boundary. Use 'gdh1' or 'rhcw18', whose pyBacktrack counterparts agree with "
        "gplately's to within a metre."
    ),
}


def check_age_depth_model_supported(age_depth_model, richards_table_filename=None):
    """Return the canonical name of `age_depth_model`, or raise if Step 5 cannot use it.

    Separate from :func:`merge_pybacktrack_paleobathymetry` so that callers running Steps
    1-4 first can find out before, rather than after, doing that work.

    Parameters
    ----------
    age_depth_model : str
        The model name, or any alias :func:`gplately.age_to_basement_depth` accepts.
    richards_table_filename : str, optional
        The lookup table Steps 1-4 would use for ``"rhcw18"``. A substituted table is a
        different age-depth relationship, and pyBacktrack has no way to be given it, so it
        is refused for the same reason a different model is.

    Raises
    ------
    ValueError
        If `age_depth_model` has no usable pyBacktrack equivalent, or if
        `richards_table_filename` would make the two sides disagree.
    """
    # Imported here rather than at module level to keep this module free of an import back
    # into gplately.grids.paleobathymetry, which imports this one.
    from .paleobathymetry import _AGE_DEPTH_MODEL_ALIASES

    # "richards" and "r18" are rhcw18; refusing them would refuse a model that does agree.
    key = _AGE_DEPTH_MODEL_ALIASES.get(
        str(age_depth_model).strip().lower(), str(age_depth_model).strip().lower()
    )
    if key in _PYBACKTRACK_AGE_DEPTH_MODEL_ATTRS:
        if key == "rhcw18" and richards_table_filename is not None:
            raise ValueError(
                "richards_table_filename cannot be combined with pyBacktrack: it replaces "
                "the RHCW18 age-depth table for Steps 1-4 only, and pyBacktrack has no way "
                "to be given the same table, so Step 5 would keep using its own. That is "
                "the same step change at the merge boundary that other mismatched models "
                "are refused for. Drop richards_table_filename, or do not run Step 5."
            )
        return key

    reason = _PYBACKTRACK_UNSUPPORTED_AGE_DEPTH_MODELS.get(key)
    raise ValueError(
        f"age_depth_model {age_depth_model!r} cannot be used with pyBacktrack. "
        + (reason + " " if reason else "")
        + "Supported models are: "
        + ", ".join(_PYBACKTRACK_AGE_DEPTH_MODEL_ATTRS)
        + "."
    )


__all__ = ["merge_pybacktrack_paleobathymetry"]


def merge_pybacktrack_paleobathymetry(
    output_file_prefix,
    merge_paleobathymetry_filename_format,
    rotation_filenames,
    static_polygon_filename,
    present_day_age_grid_filename,
    grid_spacing_degrees,
    oldest_time,
    youngest_time=0.0,
    time_increment=1,
    age_depth_model="gdh1",
    anchor_plate_id=0,
    use_all_cpus=False,
    *,
    decimal_places_in_time=None,
    **pybacktrack_kwargs,
):
    """Compute pyBacktrack paleobathymetry and merge in ``gplately``'s Steps 1-4 grids.

    A thin wrapper around ``pybacktrack.reconstruct_paleo_bathymetry_grids()``'s merge support,
    using the same plate model / ocean age-depth model / grid settings as Steps 1-4 so the two
    paleobathymetry sources stay aligned (see the *simple_paleobathymetry* README's Step 5).

    Parameters
    ----------
    output_file_prefix : str
        Passed straight through to pyBacktrack: either a plain path prefix (output files are
        named ``<output_file_prefix>_<time>.nc``), or a `Python Template string
        <https://docs.python.org/3/library/string.html#template-strings>`__ containing
        ``${time}`` (e.g. ``"output/paleobathymetry_${time}Ma.nc"``).
    merge_paleobathymetry_filename_format : str
        A Template string containing ``${time}`` identifying the Steps-1-4 paleobathymetry grid
        for each time -- e.g. the ``paleobathymetry_${time}Ma.nc`` files written by
        :func:`gplately.grids.paleobathymetry.simple_paleobathymetry` when given an
        `output_directory`.
    rotation_filenames : str or list of str
        Rotation file path(s). Unlike elsewhere in gplately, this must be filename(s) --
        pyBacktrack builds its own rotation model internally and does not accept an already
        constructed :class:`pygplates.RotationModel`.
    static_polygon_filename : str
        Static polygons, used by pyBacktrack to assign plate IDs to its synthetic drill sites.
    present_day_age_grid_filename : str
        The seafloor-age grid at 0 Ma (regardless of the time range being computed -- pyBacktrack
        backtracks from the present day).
    grid_spacing_degrees : float
        Should match the `grid_spacing` used for Steps 1-4, so the merge lines up.
    oldest_time, youngest_time : float
        Time range (Ma) to compute, inclusive.
    time_increment : float, default: 1
        The increment (Myr) that pyBacktrack generates its output at, between `youngest_time`
        and `oldest_time`. This is the *output* time step -- the spacing of the times Steps
        1-4 produced grids for -- and not the increment Step 2 stepped its reconstruction by
        (:func:`gplately.generate_distance_grids`'s `time_increment`), which is a different
        quantity and is usually finer.
    age_depth_model : str, default: "gdh1"
        One of ``"gdh1"`` or ``"rhcw18"`` (see
        :data:`gplately.grids.paleobathymetry.AGE_DEPTH_MODELS`). It should match the
        `age_depth_model` used for Steps 1-4, which is why only these two are accepted:
        ``"parsons_sclater"`` has no pyBacktrack counterpart, and pyBacktrack's
        ``CROSBY_2007`` is a different model from gplately's ``"crosby09"`` rather than the
        same one under another name, and the two differ by up to ~218 m over 0-200 Ma.
    anchor_plate_id : int, default: 0
        Should match the `anchor_plate_id` used for Steps 1-4.
    use_all_cpus : bool or int, default: False
        Passed to pyBacktrack: ``True`` to use all CPUs, or a specific number of CPUs to use.
    decimal_places_in_time : int, optional
        Decimal places of the reconstruction time in both the output filenames and the ``merge_paleobathymetry_filename_format`` names used to find the Steps 1-4 grids. Defaults to 0,
        reproducing the filenames of the workflow this was ported from. Times that are not
        distinct at this resolution would overwrite each other, so a clash raises
        `ValueError` rather than silently discarding grids -- raise this value when using a
        fractional time step. It is the same rule as pyBacktrack's
        ``output_file_decimal_places_in_time``.
    **pybacktrack_kwargs
        Extra keyword arguments passed to ``pybacktrack.reconstruct_paleo_bathymetry_grids()``
        (e.g. to override its bundled lithology/topography/sediment-thickness/crustal-thickness
        data -- see its docstring). By default pyBacktrack's own bundled global data is used, as
        recommended by pyBacktrack's documentation when swapping in a different plate model.

    Raises
    ------
    ImportError
        If the optional `pybacktrack` package is not installed.
    ValueError
        If `age_depth_model` has no usable pyBacktrack equivalent (see
        :func:`check_age_depth_model_supported`).

    References
    ----------
    Muller, R.D., Cannon, J., Williams, S. & Dutkiewicz, A. (2018). PyBacktrack 1.0: A tool for
    reconstructing paleobathymetry on oceanic and continental crust. *Geochemistry, Geophysics,
    Geosystems*, 19, 1898-1909, doi: 10.1029/2017GC007313.
    """
    try:
        import pybacktrack
    except ImportError as exc:
        raise ImportError(
            "merge_pybacktrack_paleobathymetry() requires the optional 'pybacktrack' package "
            "(`pip install pybacktrack` or `conda install -c conda-forge pybacktrack`); "
            "see https://github.com/EarthByte/pyBacktrack."
        ) from exc

    decimal_places_in_time = resolve_decimal_places_in_time(
        decimal_places_in_time, DEFAULT_DECIMAL_PLACES_IN_TIME
    )

    key = check_age_depth_model_supported(age_depth_model)
    ocean_age_to_depth_model = getattr(
        pybacktrack, _PYBACKTRACK_AGE_DEPTH_MODEL_ATTRS[key]
    )

    pybacktrack.reconstruct_paleo_bathymetry_grids(
        output_file_prefix,
        grid_spacing_degrees=grid_spacing_degrees,
        oldest_time=oldest_time,
        youngest_time=youngest_time,
        time_increment=time_increment,
        age_grid_filename=present_day_age_grid_filename,
        rotation_filenames=rotation_filenames,
        static_polygon_filename=static_polygon_filename,
        ocean_age_to_depth_model=ocean_age_to_depth_model,
        anchor_plate_id=anchor_plate_id,
        merge_paleo_bathymetry_filename_format=merge_paleobathymetry_filename_format,
        merge_paleo_bathymetry_file_decimal_places_in_time=decimal_places_in_time,
        merge_paleo_bathymetry_is_positive_below_sea_level=False,
        output_positive_bathymetry_below_sea_level=False,
        output_file_decimal_places_in_time=decimal_places_in_time,
        use_all_cpus=use_all_cpus,
        **pybacktrack_kwargs,
    )
