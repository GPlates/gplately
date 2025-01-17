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

from enum import Enum
from abc import ABC, abstractmethod

from geopandas.geodataframe import GeoDataFrame


class PlotEngineType(Enum):
    CARTOPY = 1
    PYGMT = 2


class PlotEngine(ABC):
    @abstractmethod
    def plot_geo_data_frame(self, gdf: GeoDataFrame, **kwargs):
        pass  # This is an abstract method, no implementation here.

    @abstractmethod
    def plot_pygplates_features(self, features, **kwargs):
        pass  # This is an abstract method, no implementation here.

    @abstractmethod
    def plot_subduction_zones(
        self,
        gdf_subduction_left: GeoDataFrame,
        gdf_subduction_right: GeoDataFrame,
        color="blue",
        **kwargs,
    ):
        pass  # This is an abstract method, no implementation here.
