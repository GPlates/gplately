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

from .utils import dev_warning
from .utils.check_pmm import ensure_plate_model_manager_compatible
from .utils.log_utils import setup_logging
from .utils.version import get_distribution_version

REQUIRED_PMM_VERSION = "1.2.2"  # TODO: get this from package meta
USING_DEV_VERSION = True  ## change this to False before official release

__version__ = get_distribution_version()

setup_logging()
del setup_logging

if USING_DEV_VERSION:
    dev_warning.print_dev_warning(__version__)
    dev_warning.print_using_source_code_warning(__version__)
del dev_warning

ensure_plate_model_manager_compatible(REQUIRED_PMM_VERSION)
del ensure_plate_model_manager_compatible

from plate_model_manager import PlateModelManager, PresentDayRasterManager

from . import (
    auxiliary,
    data,
    download,
    geometry,
    gpml,
    grids,
    oceans,
    plot,
    ptt,
    pygplates,
    reconstruction,
    spatial,
)
from .data import DataCollection
from .download import DataServer
from .grids import Raster
from .oceans import SeafloorGrid
from .plot import PlotTopologies
from .reconstruction import (
    PlateReconstruction,
    Points,
    _ContinentCollision,
    _DefaultCollision,
    _ReconstructByTopologies,
)
from .tools import EARTH_RADIUS
from .utils import io_utils
from .utils.io_utils import get_geometries, get_valid_geometries

__all__ = [
    # Modules
    "auxiliary",
    "data",
    "download",
    "geometry",
    "gpml",
    "grids",
    "oceans",
    "plot",
    "pygplates",
    "io_utils",
    "reconstruction",
    "ptt",
    "spatial",
    # Classes
    "DataCollection",
    "PlateModelManager",
    "PresentDayRasterManager",
    "DataServer",
    "PlateReconstruction",
    "PlotTopologies",
    "Points",
    "Raster",
    "SeafloorGrid",
    "_ContinentCollision",
    "_DefaultCollision",
    "_ReconstructByTopologies",
    # Functions
    "get_geometries",
    "get_valid_geometries",
    # Constants
    "EARTH_RADIUS",
]

__pdoc__ = {
    "data": False,
    "download": False,
    "_DefaultCollision": False,
    "_ContinentCollision": False,
    "_ReconstructByTopologies": False,
    "examples": False,
    "notebooks": False,
    "commands": False,
    "decorators": False,
    "exceptions": False,
    "lib": False,
    "pygplates": False,
    "DataCollection": False,
    "get_geometries": False,
    "get_valid_geometries": False,
    "PlateModelManager": False,
    "PresentDayRasterManager": False,
    "mapping": False,  # this folder contains no public interface
}
