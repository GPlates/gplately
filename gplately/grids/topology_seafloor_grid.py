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

import warnings

from .oceans import SeafloorGrid


class TopologySeafloorGrid(SeafloorGrid):
    """
    A class derived from :class:`SeafloorGrid` for generating seafloor grids using topologies.
    """

    def generate(self, use_topological_model=None):
        """
        Call :meth:`SeafloorGrid.reconstruct_by_topologies` to generate the seafloor grids using topologies.

        Parameters
        ----------
        use_topological_model : bool, optional
            Deprecated, and ignored. It used to choose between two ways of reconstructing the seed points,
            and there is now only one.

            .. deprecated:: 2.1

        """
        if use_topological_model is not None:
            warnings.warn(
                "`use_topological_model` keyword argument has been deprecated, it is no longer used",
                DeprecationWarning,
                stacklevel=2,
            )
        return super().reconstruct_by_topologies()
