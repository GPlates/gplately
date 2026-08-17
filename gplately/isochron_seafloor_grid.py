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
from .oceans import SeafloorGrid
from .reconstruction import PlateReconstruction
from .lib import isopolate


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

    def generate(self):
        """
        Generate the seafloor grids from the isochrons.
        """
        for _time in self.time_steps:
            rotation_model = self.plate_reconstruction.rotation_model

            print("     Current reconstruction time is", _time, "Ma")

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
            # This determines what gets written to the '.xy' file.
            # If they're all False then only the geometry (lat, lon) is written (to '.xy' file).
            interpolate_isochrons_kwargs["output_scalar_age"] = True
            interpolate_isochrons_kwargs["output_scalar_spreading_rate"] = True
            interpolate_isochrons_kwargs["output_scalar_spreading_asymmetry"] = True
            interpolate_isochrons_kwargs["output_scalar_full_spreading_rate"] = True
            interpolate_isochrons_kwargs["output_scalar_spreading_direction"] = True
            interpolate_isochrons_kwargs["output_scalar_spreading_obliquity"] = True

            output_features = isopolate.interpolate_isochrons(
                rotation_model,
                [_valid_ridges, _valid_isochrons, _valid_iso_cob],
                interval_spacing_radians,
                # Unpack the keyword arguments dict into keyword arguments...
                **interpolate_isochrons_kwargs,
            )

            # The scalar types exported to the '.xy' file - depends on 'interpolate_isochrons_kwargs' settings above.
            output_scalar_types = []
            if interpolate_isochrons_kwargs["output_scalar_age"]:
                output_scalar_types.append(pygplates.ScalarType.create_gpml("Age"))
            if interpolate_isochrons_kwargs["output_scalar_spreading_rate"]:
                output_scalar_types.append(
                    pygplates.ScalarType.create_gpml("SpreadingRate")
                )
            if interpolate_isochrons_kwargs["output_scalar_spreading_asymmetry"]:
                output_scalar_types.append(
                    pygplates.ScalarType.create_gpml("SpreadingAsymmetry")
                )
            if interpolate_isochrons_kwargs["output_scalar_full_spreading_rate"]:
                output_scalar_types.append(
                    pygplates.ScalarType.create_gpml("FullSpreadingRate")
                )
            if interpolate_isochrons_kwargs["output_scalar_spreading_direction"]:
                output_scalar_types.append(
                    pygplates.ScalarType.create_gpml("SpreadingDirection")
                )
            if interpolate_isochrons_kwargs["output_scalar_spreading_obliquity"]:
                output_scalar_types.append(
                    pygplates.ScalarType.create_gpml("SpreadingObliquity")
                )

            # We are outputting scalar coverages to '.xy' format.
            #
            # We handle this as a special case so we can write out the scalar values
            # after the xy (lat/lon) values.
            #
            # We write the scalar types in the same order as they appear in 'output_scalar_types'.
            #
            isopolate.write_coverage_features_to_xy_file(
                "InterpolatedIsochrons_" + str(_time) + "Ma.xy",
                output_features,
                output_scalar_types,
                rotation_model,
                reconstruction_time=_time,
                print_debug_output=1,
            )

            # the output features only exist in memory, these lines save to a file
            OutFile = "./InterpolatedIsochrons.gpmlz"
            file_registry = pygplates.FeatureCollectionFileFormatRegistry()
            file_registry.write(pygplates.FeatureCollection(output_features), OutFile)

            # Here we do the last step that the isopolate command-line script does, reconstructing the interpolated
            # isochrons to the selected time. The output file format is picked up from the file name
            # extension (eg .gmt, .shp). The interpolated isochrons are created in a 'tmp' file
            # print("     Creating interpolated isochrons at", ReconTime, "Ma")
            # FNAME = './tmp/InterpolatedIsochrons_'+str(ReconTime)+'Ma.shp'
            # FNAME = 'InterpolatedIsochrons_'+str(ReconTime)+'Ma.shp'
            # pgp.reconstruct(OutFile, RotFile, FNAME, ReconTime, 0)
            # print("     Creating interpolated isochrons at", ReconTime, "Ma  --- done!")
