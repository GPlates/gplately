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
from abc import ABC, abstractmethod

from geopandas.geodataframe import GeoDataFrame


class PlotEngine(ABC):
    """Abstract base class for map plotting.
    Do not use this base class directly. Use subclasses instead, such as :class:`CartopyPlotEngine` or :class:`PygmtPlotEngine`.
    """

    @abstractmethod
    def plot_geo_data_frame(self, ax_or_fig, gdf: GeoDataFrame, **kwargs):
        """Plot GeoPandas GeoDataFrame object (abstract method)"""
        pass  # This is an abstract method, no implementation here.

    @abstractmethod
    def plot_pygplates_features(self, ax_or_fig, features, **kwargs):
        """Plot one or more pygplates feature(s) (abstract method)"""
        pass  # This is an abstract method, no implementation here.

    @abstractmethod
    def plot_subduction_zones(
        self,
        ax_or_fig,
        gdf_subduction_left: GeoDataFrame,
        gdf_subduction_right: GeoDataFrame,
        color="blue",
        **kwargs,
    ):
        """Plot subduction zones with "teeth"(abstract method)"""
        pass  # This is an abstract method, no implementation here.
