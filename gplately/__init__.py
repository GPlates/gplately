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

from .utils import dev_warning
from .utils.check_pmm import ensure_plate_model_manager_compatible
from .utils.log_utils import setup_logging
from .utils.version import get_distribution_version

REQUIRED_PMM_VERSION = "1.3.0"  # TODO: get this from package meta
USING_DEV_VERSION = False  ## change this to False before official release

__version__ = get_distribution_version()

setup_logging()
del setup_logging

if USING_DEV_VERSION:
    dev_warning.print_dev_warning(__version__)
    dev_warning.print_using_source_code_warning(__version__)
del dev_warning

ensure_plate_model_manager_compatible(REQUIRED_PMM_VERSION)
del ensure_plate_model_manager_compatible

from plate_model_manager import PlateModel, PlateModelManager, PresentDayRasterManager

from . import auxiliary, ptt
from .download import DataServer
from .grids import Raster, read_netcdf_grid, reconstruct_grid
from .lib.reconstruct import (
    reconstruct_points,
    reconstruct_points_impl,
    reverse_reconstruct_points,
    reverse_reconstruct_points_impl,
)
from .mapping.cartopy_plot import CartopyPlotEngine
from .mapping.plot_engine import PlotEngine
from .mapping.pygmt_plot import PygmtPlotEngine
from .oceans import SeafloorGrid
from .plot import PlotTopologies
from .points import Points
from .ptt.resolve_topologies import (
    resolve_topological_snapshot,
    resolve_topological_snapshot_into_features,
    resolve_topologies,
    resolve_topologies_into_features,
)
from .ptt.ridge_spreading_rate import spreading_rates as ridge_spreading_rate
from .ptt.subduction_convergence import subduction_convergence
from .reconstruction import PlateReconstruction
from .tools import EARTH_RADIUS

__all__ = [
    # modules
    "auxiliary",
    "ptt",
    # main classes
    "DataServer",
    "PlateReconstruction",
    "PlotTopologies",
    "Points",
    "Raster",
    "SeafloorGrid",
    # other classes
    "PlateModel",
    "PlateModelManager",
    "PresentDayRasterManager",
    "PlotEngine",
    "CartopyPlotEngine",
    "PygmtPlotEngine",
    # functions
    "read_netcdf_grid",
    "reconstruct_grid",
    "reconstruct_points",
    "ridge_spreading_rate",
    "subduction_convergence",
    # constants
    "EARTH_RADIUS",
]
