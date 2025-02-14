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

"""

![Intro GIF](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/docs_muller19_seed_points.gif)

## Version - latest dev

## Quick start🚀

### [PlateModelManager](https://pypi.org/project/plate-model-manager/)
The **PlateModelManager** module was introduced as a more efficient alternative to the **DataServer** class, 
designed specifically for downloading and managing plate reconstruction model files. 
More information about the PlateModelManager module can be found in [its GitHub repository](https://github.com/GPlates/plate-model-manager).

```python
from gplately import (
    PlateModelManager,
    PlateReconstruction,
    PlotTopologies,
    PresentDayRasterManager,
    Raster,
)

model = PlateModelManager().get_model(
    "Muller2019",  # model name
    data_dir="plate-model-repo",  # the folder to save the model files
)

recon_model = PlateReconstruction(
    model.get_rotation_model(),
    topology_features=model.get_layer("Topologies"),
    static_polygons=model.get_layer("StaticPolygons"),
)
gplot = PlotTopologies(
    recon_model,
    coastlines=model.get_layer("Coastlines"),
    COBs=model.get_layer("COBs"),
    time=55,
)
# get present-day topography raster
raster = Raster(PresentDayRasterManager().get_raster("topography"))
# get paleo-agegrid raster at 100Ma from Muller2019 model
agegrid = Raster(model.get_raster("AgeGrids", time=100))
```

For more example code, a [comprehensive example](https://github.com/GPlates/gplately/blob/master/Notebooks/Examples/introducing-plate-model-manager.py) 
on GitHub demonstrates how to use the PlateModelManager module in details. [Another example](https://github.com/GPlates/gplately/blob/master/Notebooks/Examples/working-with-plate-model-manager.py) 
shows how to use the PlateModelManager module with GPlately.

You may use the auxiliary functions to create the `PlateReconstruction` and `PlotTopologies` instances.

```python
from gplately.auxiliary import get_gplot, get_plate_reconstruction

# use the auxiliary function to create a PlateReconstruction instance
plate_reconstruction_instance = get_plate_reconstruction("Muller2019")

# use the auxiliary function to create a PlotTopologies instance
plot_topologies_instance = get_gplot("Muller2019", time=140)

# there is a PlateReconstruction instance inside a PlotTopologies instance.
# so, in most cases a single get_gplot() call is enough.
# you can get the PlateReconstruction instance from the PlotTopologies instance later.
another_plate_reconstruction_instance = plot_topologies_instance.plate_reconstruction
```

### [DataServer](https://gplates.github.io/gplately/download.html#gplately.download.DataServer)
The `DataServer` class allows users to automatically download and cache the necessary files for plate reconstructions to a designated folder on your system. 
These files include rotation models, topology features, and static geometries such as coastlines, continents, and continent-ocean boundaries. 
Additionally, it supports the retrieval of other data types, including rasters, grids, and feature data. 
(Use the newer **PlateModelManager** whenever possible.)

```python
from gplately.download import DataServer

gdownload = DataServer("Muller2019")

# Download plate reconstruction files and geometries from the Müller et al. 2019 model
rotation_model, topology_features, static_polygons = (
    gdownload.get_plate_reconstruction_files()
)
coastlines, continents, COBs = gdownload.get_topology_geometries()

# Download the Müller et al. 2019 100 Ma age grid
age_grid = gdownload.get_age_grid(times=100)

# Download the ETOPO1 geotiff raster
etopo = gdownload.get_raster("ETOPO1_tif")
```

Both `PlateModelManager` and `DataServer` support the following plate reconstruction models:

------------------

| **Model name string Identifier** | **Zenodo** | **Topology features** | **Static polygons** | **Coast-lines**  | **Cont-inents** | **COB**    | **Age grids**   | **SR grids**  |
|:--------------------------------:|:---------:|:--------------------:|:--------------------:|:-----------------:|:---------------:|:----------:|:--------------:|:--------------:|
|  Alfonso2024                     |     ✅     |          ✅           |          ✅          |        ✅        |        ❌       |     ❌    |       ❌       |       ❌      |
|  Cao2024                         |     ✅     |          ✅           |          ✅          |        ✅        |        ✅       |     ✅    |       ❌       |       ❌      |
|  Muller2022                      |     ✅     |          ✅           |          ✅          |        ✅        |        ✅       |     ✅    |       ❌       |       ❌      |
|  Zahirovic2022                   |     ✅     |          ✅           |          ✅          |        ✅        |        ✅       |     ❌    |       ✅       |       ✅      |
|  Merdith2021                     |     ✅     |          ✅           |          ✅          |        ✅        |        ✅       |     ❌    |       ❌       |       ❌      |
|  Clennett2020                    |     ✅     |          ✅           |          ✅          |        ✅        |        ❌       |     ❌    |       ✅       |       ✅      |
|  Clennett2020_M2019              |     ✅     |          ✅           |          ✅          |        ✅        |        ❌       |     ❌    |       ✅       |       ✅      |
|  Clennett2020_S2013              |     ✅     |          ✅           |          ✅          |        ✅        |        ❌       |     ❌    |       ❌       |       ❌      |
|  Muller2019                      |     ✅     |          ✅           |          ✅          |        ✅        |        ✅       |     ✅    |       ✅       |       ❌      |
|  Young2018                       |     ✅     |          ✅           |          ✅          |        ✅        |        ✅       |     ❌    |       ❌       |       ❌      |
|  TorsvikCocks2017                |     ❌     |          ❌           |          ✅          |        ✅        |        ❌       |     ❌    |       ❌       |       ❌      |
|  Matthews2016                    |     ✅     |          ✅           |          ✅          |        ✅        |        ✅       |     ❌    |       ❌       |       ❌      |
|  Matthews2016_pmag_ref           |     ❌     |          ❌           |          ✅          |        ✅        |        ✅       |     ❌    |       ❌       |       ❌      |
|  Muller2016                      |     ✅     |          ✅           |          ✅          |        ✅        |        ❌       |     ✅    |       ✅       |       ❌      |
|  Scotese2016                     |     ✅     |          ❌           |          ✅          |        ✅        |        ❌       |     ❌    |       ❌       |       ❌      |
|  Zahirovic2016                   |     ✅     |          ✅           |          ✅          |        ✅        |        ✅       |     ❌    |       ❌       |       ❌      |
|  Gibbons2015                     |     ✅     |          ✅           |          ✅          |        ✅        |        ❌       |     ❌    |       ❌       |       ❌      |
|  Zahirovic2014                   |     ✅     |          ❌           |          ✅          |        ✅        |        ❌       |     ❌    |       ❌       |       ❌      |
|  Shephard2013                    |     ✅     |          ✅           |          ✅          |        ✅        |        ❌       |     ❌    |       ❌       |       ❌      |
|  Gurnis2012                      |     ✅     |          ✅           |          ✅          |        ✅        |        ❌       |     ❌    |       ❌       |       ❌      |
|  Seton2012                       |     ✅     |          ✅           |          ✅          |        ✅        |        ✅       |     ✅    |       ✅       |       ❌      |
|  Muller2008                      |     ❌     |          ❌           |          ✅          |        ❌        |        ❌       |     ❌    |       ❌       |       ❌      |

**Please note that all models have rotation files. The "Zenodo" column indicates whether the model files are available on [Zenodo](https://zenodo.org/).**

------------------

### [PlateReconstruction](https://gplates.github.io/gplately/reconstruction.html#gplately.reconstruction.PlateReconstruction)
The `PlateReconstruction` class contains tools to reconstruct geological features like tectonic plates and plate boundaries,
and to interrogate plate kinematic data like plate motion velocities, and rates of subduction and seafloor spreading.

```python
from gplately import PlateReconstruction, PlateModelManager

model = PlateModelManager().get_model("Muller2019")

# Build a plate reconstruction model using a rotation model, a set of topology features and static polygons
recon_model = PlateReconstruction(
    model.get_rotation_model(),
    topology_features=model.get_layer("Topologies"),
    static_polygons=model.get_layer("StaticPolygons"),
)
```

Alternatively, you may use the auxiliary functions to create a `PlateReconstruction` instance.

```python
from gplately.auxiliary import get_plate_reconstruction

# use the auxiliary function to create a PlateReconstruction instance
plate_reconstruction_instance = get_plate_reconstruction("Muller2019")
```

This [02-PlateReconstructions.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb) demonstrates in details 
how to use the `PlateReconstruction` class.

### [Points](https://gplates.github.io/gplately/reconstruction.html#gplately.reconstruction.Points)
The methods in the `Points` class track the motion of a point (or group of points) represented by a latitude and longitude 
through geologic time. This motion can be visualised using flowlines or motion paths and quantified with point 
motion velocities.

```python
import numpy as np

from gplately import PlateModelManager, Points, auxiliary

model = PlateModelManager().get_model("Muller2019")

# Create a plate reconstruction model using a rotation model, a set of topology features and static polygons
recon_model = auxiliary.get_plate_reconstruction(model)

# Define some points using their latitude and longitude coordinates so we can track them though time!
pt_lons = np.array([140.0, 150.0, 160.0])
pt_lats = np.array([-30.0, -40.0, -50.0])

# Create a Points instance from these points
gpts = Points(recon_model, pt_lons, pt_lats)
```

The [03-WorkingWithPoints.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/03-WorkingWithPoints.ipynb) demonstrates in details 
how to use the `Points` class.

![PointsDemo](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/Hawaii_Emperor_motion_path.png)

### [Raster](https://gplates.github.io/gplately/grids.html#gplately.grids.Raster)
The `Raster` class contains methods to work with netCDF4 or MaskedArray gridded data. Grids may be filled, 
resized, resampled, and reconstructed back and forwards through geologic time. Other array data can also be 
interpolated onto `Raster` grids.  

```python
from gplately import PlateModelManager, PresentDayRasterManager, Raster, auxiliary

model_name = "Muller2019"
# Create a plate reconstruction model using a rotation model, a set of topology features and static polygons
recon_model = auxiliary.get_plate_reconstruction(model_name)

# Any numpy array can be turned into a Raster object!
raster = Raster(
    plate_reconstruction=recon_model,
    data=PresentDayRasterManager().get_raster("topography"),
    extent="global",  # equivalent to (-180, 180, -90, 90)
    origin="lower",  # or set extent to (-180, 180, -90, 90)
)

# Reconstruct the raster data to 50 million years ago!
reconstructed_raster = raster.reconstruct(
    time=50,
    partitioning_features=PlateModelManager()
    .get_model(model_name)
    .get_layer("ContinentalPolygons"),
)
```

The [06-Rasters.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/06-Rasters.ipynb) demonstrates in details 
how to use the `Raster` class.

![RasterDemo](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/etopo_reconstruction.png)

### [PlotTopologies](https://gplates.github.io/gplately/plot.html#gplately.plot.PlotTopologies)
The `PlotTopologies` class works with the aforementioned `PlateReconstruction` class to plot
geologic features of different types listed 
[here](https://gplates.github.io/gplately/plot.html#gplately.plot.PlotTopologies), as well as 
coastline, continent and continent-ocean boundary geometries reconstructed through time using pyGPlates. 

```python
from gplately import PlateModelManager, PlotTopologies, auxiliary

model = PlateModelManager().get_model("Muller2019")
recon_model = auxiliary.get_plate_reconstruction(model)

gplot = PlotTopologies(
    recon_model,
    coastlines=model.get_layer("Coastlines"),
    COBs=model.get_layer("COBs"),
    continents=model.get_layer("ContinentalPolygons"),
    time=55,
)
```

You may use the auxiliary functions to create a `PlotTopologies` instance.

```python
from gplately.auxiliary import get_gplot

# use the auxiliary function to create a PlotTopologies instance
plot_topologies_instance = get_gplot("Muller2019", time=55)
```

The [02-PlateReconstructions.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb) demonstrates in details 
how to use the `PlotTopologies` class.

![PlotTopologiesDemo](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/plottopologies.png)

### [SeafloorGrid](https://gplates.github.io/gplately/oceans.html#gplately.oceans.SeafloorGrid)
The `SeafloorGrid` class wraps an automatic workflow to grid seafloor ages and seafloor spreading rates
as encoded by a plate reconstruction model. 

```python
import os

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

from gplately import SeafloorGrid, auxiliary

if __name__ == "__main__":
    gplot = auxiliary.get_gplot("Muller2019")

    # Set up automatic gridding from 5Ma to present day
    seafloorgrid = SeafloorGrid(
        PlateReconstruction_object=gplot.plate_reconstruction,  # The PlateReconstruction object
        PlotTopologies_object=gplot,  # The PlotTopologies object
        max_time=5,  # start time (Ma)
        min_time=0,  # end time (Ma)
        ridge_time_step=1,  # time increment (Myr)
    )

    # Begin automatic gridding!
    seafloorgrid.reconstruct_by_topologies()
```

The [10-SeafloorGrids.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/10-SeafloorGrids.ipynb) is a tutorial notebook that demonstrates
how to set up and use the `SeafloorGrid` object, and shows a sample set of output grids. 

![SeafloorGridDemo](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/seafloorgrid.gif)


## Sample workflows

- [__01 - Getting Started__](01-GettingStarted.html): A brief overview of how to initialise GPlately's main objects
- [__02 - Plate Reconstructions__](02-PlateReconstructions.html): Setting up a `PlateReconstruction` object, reconstructing geological data through time 
- [__03 - Working with Points__](03-WorkingWithPoints.html): Setting up a `Points` object, reconstructing seed point locations through time with. This notebook uses point data from the Paleobiology Database (PBDB).
- [__04 - Velocity Basics__](04-VelocityBasics.html): Calculating plate velocities, plotting velocity vector fields
- [__05 - Working with Feature Geometries__](05-WorkingWithFeatureGeometries.html): Processing and plotting assorted polyline, polygon and point data from [GPlates 2.3's sample data sets](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/)
- [__06 - Rasters__](06-Rasters.html): Reading, resizing, resampling raster data, and linearly interpolating point data onto raster data
- [__07 - Plate Tectonic Stats__](07-WorkingWithPlateTectonicStats.html): Using [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools) to calculate and plot subduction zone and ridge data (convergence/spreading velocities, subduction angles, subduction zone and ridge lengths, crustal surface areas produced and subducted etc.) 
- [__08 - Predicting Slab Flux__](08-PredictingSlabFlux.html): Predicting the average slab dip angle of subducting oceanic lithosphere.
- [__09 - Motion Paths and Flowlines__](09-CreatingMotionPathsAndFlowlines.html): Using pyGPlates to create motion paths and flowines of points on a tectonic plate to illustrate the plate's trajectory through geological time.
- [__10 - SeafloorGrid__](10-SeafloorGrids.html): Defines the parameters needed to set up a `SeafloorGrid` object, and demonstrates how to produce age and spreading rate grids from a set of plate reconstruction model files.
- [__11 - AndesFluxes__](11-AndesFluxes.html): Demonstrates how the reconstructed subduction history along the Andean margin can be potentially used in the plate kinematics analysis and data mining.

## Examples

- [__PlateModelManager__](https://github.com/GPlates/gplately/tree/master/Notebooks/Examples/readme.md#plate-model-manager): Examples demonstrate how to use the PlateModelManager module

## Command-line interface

- list
- combine

"""
from .utils import dev_warning
from .utils.check_pmm import ensure_plate_model_manager_compatible
from .utils.log_utils import setup_logging
from .utils.version import get_distribution_version

REQUIRED_PMM_VERSION = "1.2.1"  # TODO: get this from package meta
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
}
