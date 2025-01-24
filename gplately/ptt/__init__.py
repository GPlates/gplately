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

__pdoc__ = {
    "cleanup_topologies.add_arguments": False,
    "convert_xy_to_gplates.add_arguments": False,
    "convert_xy_to_gplates.main": False,
    "diagnose_rotations.add_arguments": False,
    "diagnose_rotations.main": False,
    "fix_crossovers.add_arguments": False,
    "gpmdb.add_arguments": False,
    "gpmdb.main": False,
    "gpmdb.ArgParser": False,
    "remove_plate_rotations.add_arguments": False,
    "remove_plate_rotations.main": False,
    "remove_plate_rotations.ArgParseAccuracyAction": False,
    "resolve_topologies.add_arguments": False,
    "resolve_topologies.main": False,
    "cleanup_topologies.iteritems": False,
    "cleanup_topologies.main": False,
    "cleanup_topologies.itervalues": False,
    "cleanup_topologies.listitems": False,
    "cleanup_topologies.listvalues": False,
    "fix_crossovers.main": False,
    "separate_ridge_transform_segments.main": False,
    "rotation_tools.add_arguments": False,
    "rotation_tools.iteritems": False,
    "rotation_tools.itervalues": False,
    "rotation_tools.listvalues": False,
    "rotation_tools.listitems": False,
    "rotation_tools.main": False,
    "rotation_tools.ArgParseAccuracyAction": False,
    "separate_ridge_transform_segments.add_arguments": False,
    "subduction_convergence.add_arguments": False,
    "subduction_convergence.main": False,
    "subduction_convergence.warning_format": False,
    "subduction_convergence.write_output_file": False,
}
