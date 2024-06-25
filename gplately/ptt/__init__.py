"""Plate Tectonic Tools"""

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
}
