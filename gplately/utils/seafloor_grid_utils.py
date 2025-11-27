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

"""Auxiliary functions for SeafloorGrid"""

import numpy as np
import pygplates

from .. import ptt
from ..lib.icosahedron import get_mesh, xyz2lonlat


def create_icosahedral_mesh(refinement_levels):
    """Return a Icospheres mesh as pygplates.MultiPointOnSphere.

    This global mesh will later be masked with a set of continental or COB terrane
    polygons to define the ocean basin at a given reconstruction time.
    The `refinement_levels` integer is proportional to the resolution of the
    mesh and the ocean/continent boundary.

    Parameters
    ----------
    refinement_levels : int
        Refine the number of points in the triangulation. The larger the
        refinement level, the sharper the ocean basin resolution.

    Returns
    -------
    multi_point : instance of <pygplates.MultiPointOnSphere>
        The longitues and latitudes that make up the icosahedral ocean mesh
        collated into a MultiPointOnSphere object.
    """

    # Create the ocean basin mesh (icosahedral spherical mesh)
    vertices, _ = get_mesh(refinement_levels, use_stripy_icosahedron=True)
    # return the mesh as MultiPointOnSphere
    return pygplates.MultiPointOnSphere(  # type:ignore
        [tuple(reversed(xyz2lonlat(v[0], v[1], v[2]))) for v in vertices]
    )
