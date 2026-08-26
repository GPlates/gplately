#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
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
from .utils.check_pmm import (
    ensure_plate_model_manager_compatible,
    get_required_pmm_version,
)
from .utils.log_utils import setup_logging
from .utils.version import get_distribution_version

setup_logging()
del setup_logging

REQUIRED_PMM_VERSION = get_required_pmm_version()

__version__ = get_distribution_version()


if any(s in __version__ for s in ["post", "git", "dirty"]):
    USING_DEV_VERSION = True
else:
    USING_DEV_VERSION = False

if USING_DEV_VERSION:
    dev_warning.print_dev_warning(__version__)
    dev_warning.print_using_source_code_warning(__version__)
del dev_warning

ensure_plate_model_manager_compatible(REQUIRED_PMM_VERSION)
del ensure_plate_model_manager_compatible
del get_required_pmm_version

from plate_model_manager import PlateModel, PlateModelManager, PresentDayRasterManager

from . import auxiliary, ptt
from .data_server import DataServer
from .raster import Raster
from .grids import (
    read_netcdf_grid,
    write_netcdf_grid,
    default_netcdf_fill_value,
    reconstruct_grid,
)
from .lib.reconstruct import (
    reconstruct_points,
    reconstruct_points_impl,
    reverse_reconstruct_points,
    reverse_reconstruct_points_impl,
)
from .plot.cartopy_plot import CartopyPlotEngine
from .plot.plot_engine import PlotEngine
from .plot.pygmt_plot import PygmtPlotEngine
from .plot.hillshade import get_topo_cmap
from .grids.oceans import SeafloorGrid
from .grids.topology_seafloor_grid import TopologySeafloorGrid
from .grids.isochron_seafloor_grid import IsochronSeafloorGrid, OutputScalarType

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
from .geometry import pygplates_to_shapely
from .utils.io_utils import load_feature_collection

# To make the `gplately.mapping` module available for backward compatibility, we import the `plot` module
# and assign it to `sys.modules["gplately.mapping"]`. This allows users to access the plotting functionalities through
# the `gplately.mapping` namespace, even though the actual implementation resides in the `plot` module.
# And also do the same thing for deprecated modules for backward compatibility.
import sys
from . import plot as _plot
from .grids import oceans as _oceans

# Import the deprecated modules for backward compatibility
from .deprecated import (
    pygplates as _pygplates,
    download as _download,
    data as _data,
    parallel as _parallel,
)

sys.modules["gplately.mapping"] = _plot
sys.modules["gplately.pygplates"] = _pygplates
sys.modules["gplately.download"] = _download
sys.modules["gplately.data"] = _data
sys.modules["gplately.parallel"] = _parallel
sys.modules["gplately.oceans"] = _oceans

# Clean up namespace
del _download
del _data
del _pygplates
del _plot
del _parallel
del _oceans
del sys

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
    "IsochronSeafloorGrid",
    "TopologySeafloorGrid",
    # other classes
    "PlateModel",
    "PlateModelManager",
    "PresentDayRasterManager",
    "PlotEngine",
    "CartopyPlotEngine",
    "PygmtPlotEngine",
    # functions
    "read_netcdf_grid",
    "write_netcdf_grid",
    "default_netcdf_fill_value",
    "reconstruct_grid",
    "reconstruct_points",
    "ridge_spreading_rate",
    "subduction_convergence",
    "load_feature_collection",
    # constants
    "EARTH_RADIUS",
]
