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

from .oceans import SeafloorGrid


class TopologySeafloorGrid(SeafloorGrid):
    """
    A class derived from :class:`SeafloorGrid` for generating seafloor grids using topologies.
    """

    def generate(self, use_topological_model=True):
        """
        Call :meth:`SeafloorGrid.reconstruct_by_topological_model` or :meth:`SeafloorGrid.reconstruct_by_topologies` to generate the seafloor grids using topologies.

        Parameters
        ----------
        use_topological_model : bool, optional
            If True, use the `pygplates.TopologicalModel <https://www.gplates.org/docs/pygplates/generated/pygplates.TopologicalModel.html>`_  class to reconstruct seed points.
            If False, use the GPlately Python code to reconstruct. Default is True.

        """
        if use_topological_model:
            return super().reconstruct_by_topological_model()
        else:
            return super().reconstruct_by_topologies()
