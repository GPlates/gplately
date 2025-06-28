#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
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

"""The "ptt" stands for Plate Tectonics Tools.

This "ptt" module provides a collection of common plate tectonic functionality that researchers can use in their workflows.
It is primarily built on top of the pyGPlates Python library.

"""

from . import (
    cleanup_topologies,
    continent_contours,
    convert_xy_to_gplates,
    remove_plate_rotations,
    resolve_topologies,
    ridge_spreading_rate,
    rotation_tools,
    separate_ridge_transform_segments,
    subduction_convergence,
    utils,
    velocity_tools,
)
