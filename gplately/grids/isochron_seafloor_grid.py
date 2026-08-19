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

import logging
import numpy as np
import pygplates
import xarray
import pygmt
from enum import Enum, auto
from pathlib import Path
from ..reconstruction import PlateReconstruction
from ..lib import isopolate
from rasterio.features import rasterize as _rasterize
from rasterio.transform import from_origin as _from_origin
from rasterio.enums import MergeAlg as _MergeAlg
from ..utils.io_utils import FeatureCollectionInput, load_feature_collection
from ..geometry import pygplates_to_shapely

logger = logging.getLogger("gplately")


def _parse_gmt_region(
    region: str | tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    """Convert a GMT-style region string (e.g. "d", "g", "minx/maxx/miny/maxy") to bounds."""
    if isinstance(region, tuple):
        return region
    if region == "d":
        return -180.0, 180.0, -90.0, 90.0
    if region == "g":
        return 0.0, 360.0, -90.0, 90.0
    minx, maxx, miny, maxy = (float(value) for value in region.split("/"))
    return minx, maxx, miny, maxy


def _parse_gmt_spacing(spacing: str | float) -> float:
    """Convert a GMT-style spacing string (e.g. "0.1d", "5m", "30s") to degrees."""
    if not isinstance(spacing, str):
        return float(spacing)
    unit = spacing[-1]
    if unit.isalpha():
        value = float(spacing[:-1])
        if unit == "m":
            return value / 60.0
        if unit == "s":
            return value / 3600.0
        if unit != "d":
            raise ValueError(f"Unsupported grid spacing unit: {unit!r}")
        return value
    return float(spacing)


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


# This implementation was based on https://github.com/EarthByte/presentday-agegridding/blob/master/compute_agegrid.sh
class IsochronSeafloorGrid:
    """
    Using isochrons interpolation to genereate seafloor grids.
    """

    def __init__(
        self,
        plate_reconstruction: PlateReconstruction,
        time_steps: np.ndarray | list,
        *,
        ridges: FeatureCollectionInput,
        isochrons: FeatureCollectionInput,
        iso_cob: FeatureCollectionInput,
        continental_polygons: FeatureCollectionInput | None = None,
        grid_region: str | tuple[float, float, float, float] = "d",
        grid_spacing: str | float = "0.1d",  # degrees
        interval_spacing_degrees: float = 0.1,
        grid_output_dir: str | Path = "seafloor_grids_by_isochron_interpolation",
    ):
        self._plate_reconstruction = plate_reconstruction
        self._time_steps = np.asarray(time_steps)
        self._ridges: pygplates.FeatureCollection = (
            ridges
            if isinstance(ridges, pygplates.FeatureCollection)
            else load_feature_collection(ridges)
        )
        self._isochrons: pygplates.FeatureCollection = (
            isochrons
            if isinstance(isochrons, pygplates.FeatureCollection)
            else load_feature_collection(isochrons)
        )
        self._iso_cob: pygplates.FeatureCollection = (
            iso_cob
            if isinstance(iso_cob, pygplates.FeatureCollection)
            else load_feature_collection(iso_cob)
        )
        self._continental_polygons: pygplates.FeatureCollection | None
        if continental_polygons is None or isinstance(
            continental_polygons, pygplates.FeatureCollection
        ):
            self._continental_polygons = continental_polygons
        else:
            self._continental_polygons = load_feature_collection(continental_polygons)
        self._interval_spacing_radians = np.radians(interval_spacing_degrees)
        if isinstance(grid_region, str) and grid_region.lower() == "global":
            grid_region = "d"
        self._grid_region = grid_region
        self._grid_spacing = grid_spacing
        self._grid_output_dir = Path(grid_output_dir)

    @property
    def time_steps(self) -> np.ndarray:
        return self._time_steps

    def generate(
        self,
        output_scalar_types: set[OutputScalarType] | None = None,
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

        for _time in self.time_steps:
            rotation_model = self._plate_reconstruction.rotation_model

            logger.info("Current reconstruction time is %s Ma", _time)

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

            # The spacing between lines (in radians) at which to generate interpolated isochrons.
            # Note that 'tessellate_threshold_radians' is the spacing ALONG lines.
            interval_spacing_radians = self._interval_spacing_radians

            # Keywords arguments for the interpolate_isochrons() function.
            interpolate_isochrons_kwargs = {}
            interpolate_isochrons_kwargs["print_debug_output"] = 0
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
            logger.debug("Interpolated GeoDataFrame:\n%s", gdf)

            # Computed once per time step and reused for every scalar type below —
            # it doesn't depend on the scalar type, and rasterizing/writing it
            # repeatedly per scalar type was wasted work.
            if self._continental_polygons is not None:
                mask = self._get_continent_mask(_time)
            else:
                mask = None

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

                logger.info(
                    "Running sphinterpolate for %s at %s Ma",
                    output_scalar_type.gpml_name,
                    _time,
                )

                median_table = pygmt.blockmedian(
                    data=xyz_table,
                    region=self._grid_region,
                    spacing=(
                        self._grid_spacing
                        if isinstance(self._grid_spacing, str)
                        else f"{self._grid_spacing}d"
                    ),
                )

                if mask is None:
                    scalar_dir = (
                        self._grid_output_dir
                        / output_scalar_type.name.lower()
                        / "NoMask"
                    )
                    scalar_dir.mkdir(parents=True, exist_ok=True)
                    grid_path = scalar_dir / (
                        f"seafloor_{output_scalar_type.name.lower()}_grid_nomask_{_time}Ma.nc"
                    )
                else:
                    scalar_dir = (
                        self._grid_output_dir
                        / output_scalar_type.name.lower()
                        / "Masked"
                    )
                    scalar_dir.mkdir(parents=True, exist_ok=True)
                    grid_path = scalar_dir / (
                        f"seafloor_{output_scalar_type.name.lower()}_grid_{_time}Ma.nc"
                    )

                grid = pygmt.sphinterpolate(
                    data=median_table,
                    region=self._grid_region,
                    spacing=(
                        self._grid_spacing
                        if isinstance(self._grid_spacing, str)
                        else f"{self._grid_spacing}d"
                    ),
                    Q=0,
                )
                if grid is None:
                    raise RuntimeError(
                        f"pygmt.sphinterpolate returned no grid for "
                        f"{output_scalar_type.gpml_name} at {_time} Ma."
                    )
                # GMT may not honour the requested spacing/region exactly, so resample
                # onto the exact target grid to guarantee the expected output size.
                target_lon, target_lat = self._target_lon_lat()
                grid = grid.interp(lon=target_lon, lat=target_lat)
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

                if mask is None:
                    clipped_grid.to_netcdf(grid_path)
                else:
                    masked_grid = clipped_grid.where(mask == 0)
                    masked_grid.to_netcdf(grid_path)

    def _target_lon_lat(self) -> tuple[np.ndarray, np.ndarray]:
        """Compute the exact lon/lat node coordinates for grid_region/grid_spacing."""
        minx, maxx, miny, maxy = _parse_gmt_region(self._grid_region)
        spacing = _parse_gmt_spacing(self._grid_spacing)
        nx = round((maxx - minx) / spacing) + 1
        ny = round((maxy - miny) / spacing) + 1
        return np.linspace(minx, maxx, nx), np.linspace(miny, maxy, ny)

    def _get_continent_mask(self, time: float, save: bool = False) -> xarray.DataArray:
        if self._continental_polygons is None:
            raise ValueError(
                "Cannot generate continent mask grid because continental polygons are not provided."
            )

        logger.info("Generating continent mask grid...")

        reconstructed_polygons = pygplates_to_shapely(
            self._plate_reconstruction.reconstruct(self._continental_polygons, time)
        )
        if not isinstance(reconstructed_polygons, list):
            reconstructed_polygons = [reconstructed_polygons]

        lon, lat = self._target_lon_lat()
        nx, ny = lon.size, lat.size
        if nx < 2 or ny < 2:
            raise ValueError(
                f"Grid must have at least 2 nodes per axis to rasterize a mask "
                f"(got nx={nx}, ny={ny})."
            )
        minx, maxx = float(lon.min()), float(lon.max())
        miny, maxy = float(lat.min()), float(lat.max())

        if not reconstructed_polygons:
            # No continental polygons at this reconstruction time.
            logger.info(
                "No continental polygons found at %s Ma; mask is all-ocean.", time
            )
            mask_array = np.zeros((ny, nx), dtype=np.uint8)
        else:
            # lon/lat are node (pixel-center) coordinates from linspace, spanning
            # nx and ny *points*, i.e. (nx-1)/(ny-1) intervals. rasterio's
            # from_bounds() instead treats nx/ny as cell *counts* spanning the
            # bounds, which would shift pixel centers off the lon/lat nodes
            # (worse the larger the grid). Build the affine transform directly so
            # pixel centers land exactly on the lon/lat node coordinates instead.
            dx = (maxx - minx) / (nx - 1)
            dy = (maxy - miny) / (ny - 1)
            transform = _from_origin(minx - dx / 2, maxy + dy / 2, dx, dy)
            mask_array = _rasterize(
                shapes=zip(reconstructed_polygons, [1] * len(reconstructed_polygons)),
                out_shape=(ny, nx),
                fill=0,
                dtype=np.uint8,
                merge_alg=_MergeAlg.replace,
                transform=transform,
            )
            if mask_array is None:
                raise RuntimeError(
                    f"Failed to rasterize continental polygons at {time} Ma."
                )
        # rasterize's row 0 is the top (maxy); flip to match ascending lat coordinates.
        mask_array = mask_array[::-1, :]

        mask_grid = xarray.DataArray(
            mask_array,
            dims=["lat", "lon"],
            coords={"lon": lon, "lat": lat},
            name="continent_mask",
        )
        mask_grid.lon.attrs.update(
            {"standard_name": "longitude", "units": "degrees_east", "axis": "X"}
        )
        mask_grid.lat.attrs.update(
            {"standard_name": "latitude", "units": "degrees_north", "axis": "Y"}
        )
        if save:
            continent_mask_dir = self._grid_output_dir / "continent_mask"
            continent_mask_dir.mkdir(parents=True, exist_ok=True)
            continent_mask_path = continent_mask_dir / f"continent_mask_{time}_Ma.nc"
            mask_grid.to_netcdf(continent_mask_path)
            logger.info(f"Continent mask grid saved to {continent_mask_path}")
        return mask_grid
