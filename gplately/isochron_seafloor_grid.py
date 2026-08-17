"""
Copyright (C) 2015-2026 The University of Sydney, Australia

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License, version 2, as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

import numpy as np
import pygplates
from enum import Enum, auto
from pathlib import Path
from .oceans import SeafloorGrid
from .reconstruction import PlateReconstruction
from .lib import isopolate


class OutputScalarType(Enum):
    AGE = auto()
    SPREADING_RATE = auto()
    SPREADING_ASYMMETRY = auto()
    FULL_SPREADING_RATE = auto()
    SPREADING_DIRECTION = auto()
    SPREADING_OBLIQUITY = auto()

    @property
    def gpml_name(self) -> str:
        return _OUTPUT_SCALAR_TYPE_NAMES[self]

    @property
    def gpml_column(self) -> str:
        return f"gpml:{self.gpml_name}"


_DEFAULT_OUTPUT_SCALAR_TYPES = set([OutputScalarType.AGE])  # for only age
# _DEFAULT_OUTPUT_SCALAR_TYPES = set(OutputScalarType) # for all
_OUTPUT_SCALAR_TYPE_SETTINGS = {
    OutputScalarType.AGE: "output_scalar_age",
    OutputScalarType.SPREADING_RATE: "output_scalar_spreading_rate",
    OutputScalarType.SPREADING_ASYMMETRY: "output_scalar_spreading_asymmetry",
    OutputScalarType.FULL_SPREADING_RATE: "output_scalar_full_spreading_rate",
    OutputScalarType.SPREADING_DIRECTION: "output_scalar_spreading_direction",
    OutputScalarType.SPREADING_OBLIQUITY: "output_scalar_spreading_obliquity",
}
_OUTPUT_SCALAR_TYPE_NAMES = {
    OutputScalarType.AGE: "Age",
    OutputScalarType.SPREADING_RATE: "SpreadingRate",
    OutputScalarType.SPREADING_ASYMMETRY: "SpreadingAsymmetry",
    OutputScalarType.FULL_SPREADING_RATE: "FullSpreadingRate",
    OutputScalarType.SPREADING_DIRECTION: "SpreadingDirection",
    OutputScalarType.SPREADING_OBLIQUITY: "SpreadingObliquity",
}


# https://github.com/EarthByte/presentday-agegridding/blob/master/compute_agegrid.sh
class IsochronSeafloorGrid(SeafloorGrid):
    """
    A SeafloorGrid that is constructed from a set of isochrons.
    """

    def __init__(
        self,
        plate_reconstruction: PlateReconstruction,
        time_steps: np.ndarray,
        *,
        ridges: pygplates.FeatureCollection,
        isochrons: pygplates.FeatureCollection,
        iso_cob: pygplates.FeatureCollection,
        **kwargs,
    ):
        kwargs["plate_reconstruction"] = plate_reconstruction
        self._time_steps = np.asarray(time_steps)
        kwargs["max_time"] = np.max(self._time_steps)
        kwargs["min_time"] = np.min(self._time_steps)
        if len(self._time_steps) > 1:
            kwargs["ridge_time_step"] = np.min(np.diff(self._time_steps))
        else:
            kwargs["ridge_time_step"] = 1

        self._ridges = ridges
        self._isochrons = isochrons
        self._iso_cob = iso_cob
        super().__init__(**kwargs)

    @property
    def time_steps(self) -> np.ndarray:
        return self._time_steps

    def generate(
        self,
        output_scalar_types: set[OutputScalarType] | None = None,
        *,
        grid_region: str = "d",
        grid_spacing: str = "0.1d",
        grid_output_dir: str | Path = "InterpolatedIsochrons_grids",
    ):
        """
        Generate the seafloor grids from the isochrons.
        """
        if output_scalar_types is None:
            output_scalar_types = _DEFAULT_OUTPUT_SCALAR_TYPES

        invalid_output_scalar_types = [
            output_scalar_type
            for output_scalar_type in output_scalar_types
            if not isinstance(output_scalar_type, OutputScalarType)
        ]
        if invalid_output_scalar_types:
            raise TypeError(
                "output_scalar_types must contain only OutputScalarType values: "
                f"{invalid_output_scalar_types!r}"
            )

        output_scalar_types = set(output_scalar_types)
        grid_output_dir = Path(grid_output_dir)

        try:
            import pygmt
        except ImportError as exc:
            raise ImportError(
                "PyGMT is required to generate NetCDF grids from interpolated isochrons."
            ) from exc

        for _time in self.time_steps:
            rotation_model = self.plate_reconstruction.rotation_model

            print("Current reconstruction time is", _time, "Ma")

            _valid_isochrons = pygplates.FeatureCollection(
                [
                    feature
                    for feature in self._isochrons
                    if feature.get_valid_time()[1] <= _time
                ]
            )

            _valid_ridges = pygplates.FeatureCollection(
                [
                    feature
                    for feature in self._ridges
                    if feature.get_valid_time()[1] <= _time
                ]
            )

            _valid_iso_cob = pygplates.FeatureCollection(
                [
                    feature
                    for feature in self._iso_cob
                    if feature.get_valid_time()[1] <= _time
                ]
            )

            # This is where isopolate is actually called - note that these one line interpolates the isochrons,
            # but does not reconstruct them

            # The spacing between lines (in radians) at which to generate interpolated isochrons.
            # Note that 'tessellate_threshold_radians' is the spacing ALONG lines.
            # interval_spacing_radians = isopolate.DEFAULT_INTERPOLATE_RESOLUTION_RADIANS  # Default spacing (0.1 degrees).
            interval_spacing_radians = 0.002

            # Keywords arguments for the interpolate_isochrons() function.
            interpolate_isochrons_kwargs = {}
            interpolate_isochrons_kwargs["print_debug_output"] = 1
            interpolate_isochrons_kwargs["tessellate_threshold_radians"] = (
                interval_spacing_radians  # Make spacing the same along *and* between lines.
            )
            # This is where we enable/disable outputing of various scalar types.
            for (
                output_scalar_type,
                setting_name,
            ) in _OUTPUT_SCALAR_TYPE_SETTINGS.items():
                interpolate_isochrons_kwargs[setting_name] = (
                    output_scalar_type in output_scalar_types
                )

            output_features = isopolate.interpolate_isochrons(
                rotation_model,
                (_valid_ridges, _valid_isochrons, _valid_iso_cob),
                interval_spacing_radians,
                # Unpack the keyword arguments dict into keyword arguments...
                **interpolate_isochrons_kwargs,
            )

            scalar_types = [
                pygplates.ScalarType.create_gpml(
                    _OUTPUT_SCALAR_TYPE_NAMES[output_scalar_type]
                )
                for output_scalar_type in OutputScalarType
                if output_scalar_type in output_scalar_types
            ]

            gdf = isopolate.get_geodataframe_from_coverage_features(
                output_features,
                output_scalar_types=scalar_types,
                rotation_model=rotation_model,
                reconstruction_time=_time,
            )
            if (
                _time >= 0
                and OutputScalarType.AGE in output_scalar_types
                and "gpml:Age" in gdf.columns
            ):
                gdf["gpml:Age"] = gdf["gpml:Age"] - _time
            print(gdf)

            # Build one NetCDF grid for each selected scalar type.
            for output_scalar_type in output_scalar_types:
                scalar_column = output_scalar_type.gpml_column
                if scalar_column not in gdf.columns:
                    continue

                if hasattr(gdf, "geometry"):
                    lons = gdf.geometry.x
                    lats = gdf.geometry.y
                else:
                    raise ValueError(
                        "Expected geometry column in interpolated GeoDataFrame."
                    )

                xyz_table = gdf.copy()
                xyz_table["lon"] = lons
                xyz_table["lat"] = lats
                xyz_table = xyz_table[["lon", "lat", scalar_column]].dropna()

                if xyz_table.empty:
                    continue

                print(
                    "Running sphinterpolate for",
                    output_scalar_type.gpml_name,
                    "at",
                    _time,
                    "Ma",
                )

                median_table = pygmt.blockmedian(
                    data=xyz_table,
                    region=grid_region,
                    spacing=grid_spacing,
                )

                scalar_dir = (
                    grid_output_dir / output_scalar_type.name.lower() / "NoMask"
                )
                scalar_dir.mkdir(parents=True, exist_ok=True)
                grid_path = scalar_dir / (
                    f"{output_scalar_type.name.lower()}grid_final_nomask_{_time}Ma.nc"
                )

                grid = pygmt.sphinterpolate(
                    data=median_table,
                    region=grid_region,
                    spacing=grid_spacing,
                    Q=0,
                )

                # Prevent small negative artifacts introduced by interpolation.
                clipped_grid = grid.clip(min=0.01)
                finite_values = clipped_grid.values[np.isfinite(clipped_grid.values)]
                if finite_values.size > 0:
                    min_value = float(np.min(finite_values))
                    max_value = float(np.max(finite_values))
                    clipped_grid.attrs["actual_range"] = np.array(
                        [min_value, max_value], dtype=np.float32
                    )
                    clipped_grid.attrs["valid_min"] = min_value
                    clipped_grid.attrs["valid_max"] = max_value
                else:
                    clipped_grid.attrs.pop("actual_range", None)
                clipped_grid.to_netcdf(grid_path)
