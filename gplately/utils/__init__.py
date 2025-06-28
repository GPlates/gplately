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

from .feature_utils import shapelify_features
from .io_utils import get_geometries, get_valid_geometries
from .plot_utils import plot_subduction_teeth
from .seafloor_grid_utils import (
    create_icosahedral_mesh,
    ensure_polygon_geometry,
    point_in_polygon_routine,
)

__all__ = [
    "shapelify_features",
    "get_geometries",
    "get_valid_geometries",
    "plot_subduction_teeth",
    "create_icosahedral_mesh",
    "ensure_polygon_geometry",
    "point_in_polygon_routine",
]

__pdoc__ = {
    "check_pmm": False,
    "dev_warning": False,
    "feature_utils": False,
    "io_utils": False,
    "log_utils": False,
    "plot_utils": False,
    "seafloor_grid_utils": False,
    "version": False,
    "get_geometries": False,
    "get_valid_geometries": False,
    # "crustal_production": False,
}
