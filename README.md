<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/GPlately_White_logo.png">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/GPlately_Main_logo.png">
  <img alt="GPlately logo." src="https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/GPlately_Main_logo.png">
</picture>
</p>

![PyPI](https://img.shields.io/pypi/v/gplately?style=for-the-badge)
![PyPI - Downloads](https://img.shields.io/pypi/dm/gplately?style=for-the-badge)
![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/gplately?style=for-the-badge)
![downloads](https://img.shields.io/conda/d/conda-forge/gplately?style=for-the-badge)
![GitHub Unitttest Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/GPlates/gplately/build_and_test.yml?branch=master&style=for-the-badge&label=test)
![GitHub Build Doc Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/GPlates/gplately/deploy_documentation.yaml?branch=master&style=for-the-badge&label=doc)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/gplately?style=for-the-badge)


GPlately was created to accelerate spatio-temporal data analysis by leveraging [pyGPlates](https://www.gplates.org/docs/pygplates/index.html) and [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools) within a simplified Python interface. This object-oriented package enables the reconstruction of data through deep geologic time (such as points, lines, polygons, and rasters), the interrogation of plate kinematic information (plate velocities, rates of subduction and seafloor spreading), the rapid comparison of multiple plate motion models, and the plotting of reconstructed output data on maps. All tools are designed to be parallel-safe, accelerating spatio-temporal analysis over multiple CPU processors.

![SeedPointGIF](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/muller19_seedpoints.gif)

GPlately requires a working installation of pyGPlates, which is freely
available at https://www.gplates.org/download.
All major system architectures (e.g., Linux, macOS, Windows) are supported, and installation instructions
are [well documented](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation).
Sample data is also available from [EarthByte servers](https://www.earthbyte.org/category/resources/), which
include rasters, seafloor age grids, rotation files, and more to help you get started with plate reconstructions.

#### Citation

> Mather, B.R., Müller, R.D., Zahirovic, S., Cannon, J., Chin, M., Ilano, L., Wright, N.M., Alfonso, C., Williams, S., Tetley, M., Merdith, A. (2023) Deep time spatio-temporal data analysis using pyGPlates with PlateTectonicTools and GPlately. _Geoscience Data Journal_, 1–8. Available from: https://doi.org/10.1002/gdj3.185

```bib
@article{Mather2023,
author = {Mather, Ben R. and Müller, R. Dietmar and Zahirovic, Sabin and Cannon, John and Chin, Michael and Ilano, Lauren and Wright, Nicky M. and Alfonso, Christopher and Williams, Simon and Tetley, Michael and Merdith, Andrew},
title = {Deep time spatio-temporal data analysis using pyGPlates with PlateTectonicTools and GPlately},
year = {2023},
journal = {Geoscience Data Journal},
pages = {1-8},
keywords = {geospatial, plate reconstructions, pyGPlates, python, tectonics},
doi = {https://doi.org/10.1002/gdj3.185},
url = {https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/gdj3.185},
eprint = {https://rmets.onlinelibrary.wiley.com/doi/pdf/10.1002/gdj3.185},
}
```

## Dependencies

- [pyGPlates](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation)
- [plate-model-manager](https://pypi.org/project/plate-model-manager/) >= 1.2.2
- [Shapely](https://shapely.readthedocs.io/en/stable/project.html#installing-shapely)
- [NumPy](https://numpy.org/install/) > 1.16
- [SciPy](https://scipy.org/install/) > 1.0
- [Matplotlib](https://matplotlib.org/stable/users/installing/index.html)
- [Cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html#getting-started) (for mapping)
- [Shapely](https://shapely.readthedocs.io/en/stable/installation.html)
- [Pooch](https://github.com/fatiando/pooch)
- [GeoPandas](https://geopandas.org/en/stable/getting_started.html)
- [netCDF4](https://unidata.github.io/netcdf4-python/#quick-install)
- [pygmt](https://www.pygmt.org/latest/)
- [rioxarray](https://github.com/corteva/rioxarray)

## Installation

### 1. Using conda (recommended)

You can install the latest stable public release of `GPlately` and all its dependencies using conda.
This is the preferred method for installing `GPlately`, as it downloads binaries from the "conda-forge" channel.

```sh
conda install -c conda-forge gplately
```

#### Creating a new conda environment

We recommend creating a new conda environment in which to install `GPlately`. This avoids any potential conflicts within your base Python environment. In the example below, we create a new environment called `my-env`.

```sh
conda create -n my-env
conda activate my-env
conda install -c conda-forge gplately
```

The `my-env` environment needs to be activated before using `GPlately`, i.e., `conda activate my-env`.

### 2. Using pip

Alternatively, you can install the latest stable public release of `GPlately` using the pip package manager.

```sh
pip install gplately
```

or from this GitHub repository:

```sh
pip install git+https://github.com/GPlates/gplately.git
```

#### Pull from repository

**First-time installation:** To install the latest version of GPlately from a specific repository branch (e.g., `master`), run the following commands in your terminal:

```sh
cd /directory-to-keep-gplately-files # go into the directory in which you'd like to keep the GPlately repository
git clone https://github.com/GPlates/gplately.git gplately.git
cd gplately.git # go into the root directory of your cloned gplately repository
git checkout master # check out the "master" branch or the name of branch you want
git pull # fetch all recent code changes from the GitHub remote repository to your computer
pip install . # alternatively, you can use "pip install -e ." to install gplately in editable mode
```

**Update installation:** To update your installation of GPlately by fetching the latest code from a specific repository branch (e.g., `master`), run the following commands in your terminal:

```sh
cd /path-to-gplately-repository # go into the root directory of your cloned gplately repository, such as ./gplately.git
git checkout master # check out the "master" branch or the name of branch you want
git pull # fetch all recent code changes from the GitHub remote repository to your computer
pip install . # alternatively, you can use "pip install -e ." to install gplately in editable mode
```

### 3. Using Docker 🐳

👉 Run GPlately example notebooks with Docker

- `docker pull gplates/gplately`
- `docker run --rm -ti -p 8888:8888  gplates/gplately`
- http://localhost:8888

👉 Run GPlately command with Docker

- `docker run gplates/gplately gplately --version`
- `docker run gplates/gplately gplately --help`

👉 Run your Python script with Docker

- `docker run -it --rm -v "$PWD":/ws -w /ws gplates/gplately python my_script_to_run.py` (assume my_script_to_run.py is in the current working directory)

See details [docker/README.md](docker/README.md).

## Usage

- [Quick start](https://gplates.github.io/gplately/dev-doc/#quick-start) - a brief tutorial to help users get up to speed  
- [Sample workflows](#sample-workflows) - demonstrations of how GPlately can be used in various workflows
- [Examples](Notebooks/Examples/readme.md) - code snippets that demonstrate the usage of GPlately's classes and methods
- [command-line interface (CLI)](#command-line-interface-(CLI)) - use GPlately in a command-line environment

## Documentation

- [Documentation dev (latest unstable)](https://gplates.github.io/gplately/dev-doc)
- [Documentation v1.3.0 (**latest stable**)](https://gplates.github.io/gplately/v1.3.0/)
- [Documentation v1.2.0](https://gplates.github.io/gplately/v1.2.0/)
- [Documentation v1.1.0](https://gplates.github.io/gplately/v1.1.0/)
- [Documentation v1.0.0](https://gplates.github.io/gplately/v1.0.0/)

## Classes

- [`DataServer`](#the-dataserver-object) - download rotation files and topology features (use the newer **PlateModelManager** when possible)
- [`PlateModelManager`](#the-platemodelmanager-object) - download and manage the plate reconstruction model files
- [`PlateReconstruction`](#the-platereconstruction-object) - reconstruct features, tesselate mid ocean ridges, subduction zones
- [`Points`](#the-points-object) - partition points onto plates, rotate back through time
- [`Raster`](#the-raster-object) - read in NetCDF grids, interpolation, resampling
- [`PlotTopologies`](#the-plottopologies-object) - one stop shop for plotting ridges, trenches, subduction teeth

### The `PlateModelManager` class

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
    data_dir="plate-model-repo",  # the local folder where you would like to save the model files
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

You may use the auxiliary functions to create the PlateReconstruction and PlotTopologies instances.

```python
from gplately.auxiliary import get_gplot, get_plate_reconstruction

# use the auxiliary function to create a PlateReconstruction instance
plate_reconstruction_instance = get_plate_reconstruction("Muller2019")

# use the auxiliary function to create a PlotTopologies instance
plot_topologies_instance = get_gplot("Muller2019", age=140)

# there is a PlateReconstruction instance inside a PlotTopologies instance.
# so, in most cases a single get_gplot() call is enough.
# You can get the PlateReconstruction instance from a PlotTopologies instance using the one-line code below.
another_plate_reconstruction_instance = plot_topologies_instance.plate_reconstruction
```

### The `DataServer` class

The **DataServer** class can be used to download:

- rotation models
- topology features
- static polygons
- coastlines
- continents
- continent-ocean boundaries
- age grids and rasters
- geological feature data

from assorted plate reconstruction models. These files are needed to construct most of `GPlately`'s objects. For example,
we can download a `rotation model`, a set of `topology features` and some `static polygons` from the [Müller et al. 2019](https://www.earthbyte.org/muller-et-al-2019-deforming-plate-reconstruction-and-seafloor-age-grids-tectonics/)
global Mesozoic–Cenozoic deforming plate motion model. (Use the newer **PlateModelManager** class when it is possible.)

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

### The `PlateReconstruction` class

The **PlateReconstruction** class contains tools to reconstruct geological features like tectonic plates and plate boundaries, and to interrogate plate kinematic data like plate motion velocities, and rates of subduction and seafloor spreading.

A complete Jupyter notebook example is available at https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb .

```python
from plate_model_manager import PlateModelManager
from gplately import PlateReconstruction

model = PlateModelManager().get_model("Muller2019")

# Build a plate reconstruction model using a rotation model, a set of topology features and static polygons
recon_model = PlateReconstruction(
    model.get_rotation_model(),
    topology_features=model.get_layer("Topologies"),
    static_polygons=model.get_layer("StaticPolygons"),
)
```

### The `Points` class

The methods in the **Points** class track the motion of a point (or group of points) represented by a latitude and longitude through geologic time. This motion can be visualised using flowlines or motion paths and quantified with point motion velocities.

A complete Jupyter notebook example is available at https://github.com/GPlates/gplately/blob/master/Notebooks/03-WorkingWithPoints.ipynb .

```python
import numpy as np
import gplately
from plate_model_manager import PlateModelManager

model = PlateModelManager().get_model("Muller2019")

# Create a plate reconstruction model using a rotation model, a set of topology features and static polygons
recon_model = gplately.auxiliary.get_plate_reconstruction(model)

# Define some points using their latitude and longitude coordinates so we can track them though time!
pt_lons = np.array([140., 150., 160.])
pt_lats = np.array([-30., -40., -50.])

# Create a Points instance from these points
gpts = gplately.Points(recon_model, pt_lons, pt_lats)
```

![PointData](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/surface_hotspot_plumes.png)

### The `Raster` class

The **Raster** class contains methods to work with netCDF4 or MaskedArray gridded data. Grids may be filled, resized, resampled, and reconstructed back and forwards through geologic time. Other array data can also be interpolated onto Raster grids.

A complete Jupyter notebook example is available at https://github.com/GPlates/gplately/blob/master/Notebooks/06-Rasters.ipynb .

```python
import gplately
from plate_model_manager import PlateModelManager, PresentDayRasterManager

model = PlateModelManager().get_model("Muller2019")

# Create a plate reconstruction model using a rotation model, a set of topology features and static polygons
recon_model = gplately.auxiliary.get_plate_reconstruction(model)

# Any numpy array can be turned into a Raster object!
raster = gplately.Raster(
    plate_reconstruction=recon_model,
    data=PresentDayRasterManager().get_raster("topography"),
    extent="global",  # equivalent to (-180, 180, -90, 90)
    origin="lower",  # or set extent to (-180, 180, -90, 90)
)

# Reconstruct the raster data to 50 million years ago! 
reconstructed_raster = raster.reconstruct(
    time=50, 
    partitioning_features=model.get_layer("ContinentalPolygons")
    )
```

Below is a plot of the [ETOPO1 global relief raster](https://www.ncei.noaa.gov/products/etopo-global-relief-model) at present day, and reconstructed to 50Ma:

![RasterImg](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/etopo_reconstruction.png)

### The `PlotTopologies` class

The **PlotTopologies** class works with the PlateReconstruction class to plot geologic features of different types listed here, as well as coastline, continent and continent-ocean boundary geometries reconstructed through time using pyGPlates.

A complete Jupyter notebook example is available at https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb .

```python
import gplately
from plate_model_manager import PlateModelManager
from gplately import PlotTopologies

model = PlateModelManager().get_model("Muller2019")
recon_model = gplately.auxiliary.get_plate_reconstruction(model)

gplot = PlotTopologies(
        recon_model,
        coastlines=model.get_layer("Coastlines"),
        COBs=model.get_layer("COBs"),
        continents= model.get_layer("ContinentalPolygons"),
        time=55)
```

Below are some continents, coastlines, COBs, ridges and transforms, trenches, subduction teeth and
seafloor age grids plotted using `PlotTopologies`.

![ReconstructionImage](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/plot_topologies_img.png)

## Sample workflows

To see GPlately in action, launch a Jupyter Notebook environment and check out the [sample notebooks](./Notebooks):

- [**01 - Getting Started**](https://github.com/GPlates/gplately/blob/master/Notebooks/01-GettingStarted.ipynb): A brief overview of how to initialise GPlately's main objects
- [**02 - Plate Reconstructions**](https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb): Setting up a `PlateReconstruction` object, reconstructing geological data through time
- [**03 - Working with Points**](https://github.com/GPlates/gplately/blob/master/Notebooks/03-WorkingWithPoints.ipynb): Setting up a `Points` object, reconstructing seed point locations through time with. This notebook uses point data from the Paleobiology Database (PBDB).
- [**04 - Velocity Basics**](https://github.com/GPlates/gplately/blob/master/Notebooks/04-VelocityBasics.ipynb): Calculating plate velocities, plotting velocity vector fields
- [**05 - Working with Feature Geometries**](https://github.com/GPlates/gplately/blob/master/Notebooks/05-WorkingWithFeatureGeometries.ipynb): Processing and plotting assorted polyline, polygon and point data from [GPlates 2.3's sample data sets](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/)
- [**06 - Rasters**](https://github.com/GPlates/gplately/blob/master/Notebooks/06-Rasters.ipynb): Reading, resizing, resampling raster data, and linearly interpolating point data onto raster data
- [**07 - Plate Tectonic Stats**](https://github.com/GPlates/gplately/blob/master/Notebooks/07-WorkingWithPlateTectonicStats.ipynb): Using [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools) to calculate and plot subduction zone and ridge data (convergence/spreading velocities, subduction angles, subduction zone and ridge lengths, crustal surface areas produced and subducted etc.)
- [**08 - Predicting Slab Flux**](https://github.com/GPlates/gplately/blob/master/Notebooks/08-PredictingSlabFlux.ipynb): Predicting the average slab dip angle of subducting oceanic lithosphere.
- [**09 - Motion Paths and Flowlines**](https://github.com/GPlates/gplately/blob/master/Notebooks/09-CreatingMotionPathsAndFlowlines.ipynb): Using pyGPlates to create motion paths and flowlines of points on a tectonic plate to illustrate the plate's trajectory through geological time.
- [**10 - SeafloorGrid**](https://github.com/GPlates/gplately/blob/master/Notebooks/10-SeafloorGrids.ipynb): Defines the parameters needed to set up a `SeafloorGrid` object, and demonstrates how to produce age and spreading rate grids from a set of plate reconstruction model files.
- [**11 - AndesFluxes**](https://github.com/GPlates/gplately/blob/master/Notebooks/11-AndesFluxes.ipynb): Demonstrates how the reconstructed subduction history along the Andean margin can be potentially used in the plate kinematics analysis and data mining.

