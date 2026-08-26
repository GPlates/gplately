#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
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

from ._grids import *
from ._utils import num_grid_points


# The Raster class is a top-level class in gplately, and should be imported from gplately directly,
# not from gplately.grids.
# For historical reasons, Raster is also available in gplately.grids, but this is deprecated and has been removed.
# This error will be raised if someone tries to import Raster from gplately.grids.
def __getattr__(name):
    if name in ["Raster"]:
        raise ImportError(
            f"The Raster class has been taken out of gplately.grids. Use `from gplately import Raster` instead of `from gplately.grids import Raster`"
        )
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "fill_raster",
    "read_netcdf_grid",
    "write_netcdf_grid",
    "default_netcdf_fill_value",
    "sample_grid",
    "reconstruct_grid",
    "rasterise",
    "rasterize",
    "num_grid_points",
]
