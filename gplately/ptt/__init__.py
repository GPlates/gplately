#
#    Copyright (C) 2024 The University of Sydney, Australia
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

"""
Plate Tectonic Tools
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
from .documentation import install_documentation

__pdoc__ = {
    "cleanup_topologies.add_arguments": False,
    "cleanup_topologies.iteritems": False,
    "cleanup_topologies.main": False,
    "cleanup_topologies.itervalues": False,
    "cleanup_topologies.listitems": False,
    "cleanup_topologies.listvalues": False,
    "fix_crossovers.main": False,
    "separate_ridge_transform_segments.main": False,
}
