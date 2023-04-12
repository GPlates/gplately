<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/GPlately_White_logo.png">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/GPlately_Main_logo.png">
  <img alt="GPlately logo." src="https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/GPlately_Main_logo.png">
</picture>
</p>

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/GPlates/gplately/build_and_test.yml?branch=master&style=for-the-badge)
![PyPI](https://img.shields.io/pypi/v/gplately?style=for-the-badge)
![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/gplately?style=for-the-badge)
![PyPI - Downloads](https://img.shields.io/pypi/dw/gplately?style=for-the-badge)


GPlately was created to accelerate spatio-temporal data analysis leveraging [pyGPlates](https://www.gplates.org/docs/pygplates/index.html) and [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools) within a simplified Python interface. This object-oriented package enables the reconstruction of data through deep geologic time (points, lines, polygons, and rasters), the interrogation of plate kinematic information (plate velocities, rates of subduction and seafloor spreading), the rapid comparison between multiple plate motion models, and the plotting of reconstructed output data on maps. All tools are designed to be parallel-safe to accelerate spatio-temporal analysis over multiple CPU processors.

![SeedPointGIF](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/muller19_seedpoints.gif)


GPlately requires a working installation of pyGPlates, which is freely
available at https://www.gplates.org/download.
All major system architectures (e.g. Linux, MacOS, Windows) are supported and installation instructions
are [well documented](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html\#installation).
Sample data is also available from [EarthByte servers](https://www.earthbyte.org/category/resources/), which
includes rasters, seafloor age grids, rotation files, and more to get started with plate reconstructions.

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
- [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools)
- [Shapely](https://shapely.readthedocs.io/en/stable/project.html#installing-shapely)
- [NumPy](https://numpy.org/install/) > 1.16
- [SciPy](https://scipy.org/install/) > 1.0
- [Matplotlib](https://matplotlib.org/stable/users/installing/index.html)
- [Cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html#getting-started) (for mapping)
- [Shapely](https://shapely.readthedocs.io/en/stable/installation.html)
- [Pooch](https://github.com/fatiando/pooch)
- [GeoPandas](https://geopandas.org/en/stable/getting_started.html)
- [netCDF4](https://unidata.github.io/netcdf4-python/#quick-install)


## Installation

### 1. Using conda (recommended)

You can install the latest stable public release of `GPlately` and all of its dependencies using conda.
This is the preferred method to install `GPlately` which downloads binaries from the conda-forge channel.

```sh
conda install -c conda-forge gplately
```

#### Creating a new conda environment

We recommend creating a new conda environment inside which to install `GPlately`. This avoids any potential conflicts in your base Python environment. In the example below we create a new environment called "`my-env`":

```sh
conda create -n my-env
conda activate my-env
conda install -c conda-forge gplately
```

`my-env` needs to be activated whenever you use `GPlately`: i.e. `conda activate my-env`.

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

**First-time installation:** To install the latest version of GPlately from a specific repository branch (e.g. `master`), copy the following commands into your terminal:

```sh
cd /path/to/desired/directory #Change your command directory to where you'd like to clone GPlately
git clone https://github.com/GPlates/gplately.git
cd gplately # navigate within the gplately folder
git checkout master # or the name of whichever branch you need
git pull # fetch all recent changes from this branch
pip install .
```

**Update installation from cloned repo:** To update your installation of GPlately by fetching the latest pushes from a specific repository branch (e.g. `master`), copy the following commands into your terminal:

```sh
cd /path/to/gplately/directory #Should be where gplately is cloned - must end in /.../gplately
git checkout master # or the name of whichever branch you need
git pull # fetch all recent changes from this branch
pip install .
```


## Usage

GPlately uses objects to accomplish a variety of common tasks. The common objects include:

- [`DataServer`](#the-dataserver-object) - download rotation files and topology features from plate models on EarthByte's webDAV server
- [`PlateReconstruction`](#the-platereconstruction-object) - reconstruct features, tesselate mid ocean ridges, subduction zones
- [`Points`](#the-points-object) - partition points onto plates, rotate back through time
- [`Raster`](#the-raster-object) - read in NetCDF grids, interpolation, resampling
- [`PlotTopologies`](#the-plottopologies-object) - one stop shop for plotting ridges, trenches, subduction teeth


### The `DataServer` object

`GPlately`'s `DataServer` object can be used to download:

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
global Mesozoic–Cenozoic deforming plate motion model. 

```python
gDownload = gplately.DataServer("Muller2019")
rotation_model, topology_features, static_polygons = gDownload.get_plate_reconstruction_files()
```

### The `PlateReconstruction` object

... contains methods to reconstruct the positions of present-day feature data back through geological time. You can also use
it to calculate plate model data like topological plate velocities, or total trench and ridge lengths per Ma! You can create
the object by passing a `rotation model`, a set of `topology features` and some `static polygons`: 

```python
model = gplately.PlateReconstruction(rotation_model, topology_features, static_polygons)
```
Launch the [Plate Reconstruction](https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb) notebook to see more.

### The `Points` object

... can be used to reconstruct the positions of geological point features and calculate their underlying plate velocities 
through geological time. 

```python
pt_lon = np.array([-107.662152, -58.082792, 17.483189, 133.674590, 80.412876])
pt_lat = np.array([48.797807, -12.654857, 11.884395, -26.415630, 31.368509])

# Call the Points object: pass the PlateReconstruction object, and the latitudes and longitudes of the seed points!
gpts = gplately.Points(model, pt_lon, pt_lat)
```
![PointData](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/surface_hotspot_plumes.png)

### The `Raster` object

...can be used to read, resample and resize assorted raster data like `netCDF4` seafloor age grids, continental grids and ETOPO
relief rasters. You can also reconstruct raster data back through geological time!

```python
etopo = gdownload.get_raster("ETOPO1_tif")

raster = gplately.Raster(
    model,
    data=etopo,
    time=0,
    origin="upper",
)
white_rgb = (255, 255, 255)  # RGB code for white, to fill gaps in output

reconstructed = raster.reconstruct(
    time=50,
    fill_value=white_rgb,
    threads=4,
)
```

Below is a plot of the [ETOPO1 global relief raster](https://www.ncei.noaa.gov/products/etopo-global-relief-model) at present day, and reconstructed to 50Ma:

![RasterImg](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/etopo_reconstruction.png)

### The `PlotTopologies` object

... can be used to visualise reconstructed feature geometries through time. To call the object, pass a set of `continents`, 
`coastlines` and `COBs` (either as file paths or as `<pyGPlates.FeatureCollection>` objects), as well as a `PlateReconstruction`
object, and a reconstruction `time`. 

```python
coastlines, continents, COBs = gDownload.get_topology_geometries()
time = 50 #Ma
gPlot = gplately.plot.PlotTopologies(model, time, coastlines, continents, COBs)
```
Below are some continents, coastlines, COBs, ridges and transforms, trenches, subduction teeth and
seafloor age grids plotted using `PlotTopologies`!

![ReconstructionImage](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/plot_topologies_img.png)


## Sample workflows

To see GPlately in action, launch a Jupyter Notebook environment and check out the [sample notebooks](./Notebooks):

- [__01 - Getting Started__](https://github.com/GPlates/gplately/blob/master/Notebooks/01-GettingStarted.ipynb): A brief overview of how to initialise GPlately's main objects
- [__02 - Plate Reconstructions__](https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb): Setting up a `PlateReconstruction` object, reconstructing geological data through time 
- [__03 - Working with Points__](https://github.com/GPlates/gplately/blob/master/Notebooks/03-WorkingWithPoints.ipynb): Setting up a `Points` object, reconstructing seed point locations through time with. This notebook uses point data from the Paleobiology Database (PBDB).
- [__04 - Velocity Basics__](https://github.com/GPlates/gplately/blob/master/Notebooks/04-VelocityBasics.ipynb): Calculating plate velocities, plotting velocity vector fields
- [__05 - Working with Feature Geometries__](https://github.com/GPlates/gplately/blob/master/Notebooks/05-WorkingWithFeatureGeometries.ipynb): Processing and plotting assorted polyline, polygon and point data from [GPlates 2.3's sample data sets](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/)
- [__06 - Rasters__](https://github.com/GPlates/gplately/blob/master/Notebooks/06-Rasters.ipynb): Reading, resizing, resampling raster data, and linearly interpolating point data onto raster data
- [__07 - Plate Tectonic Stats__](https://github.com/GPlates/gplately/blob/master/Notebooks/07-WorkingWithPlateTectonicStats.ipynb): Using [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools) to calculate and plot subduction zone and ridge data (convergence/spreading velocities, subduction angles, subduction zone and ridge lengths, crustal surface areas produced and subducted etc.) 
- [__08 - Predicting Slab Flux__](https://github.com/GPlates/gplately/blob/master/Notebooks/08-PredictingSlabFlux.ipynb): Predicting the average slab dip angle of subducting oceanic lithosphere.
- [__09 - Motion Paths and Flowlines__](https://github.com/GPlates/gplately/blob/master/Notebooks/09-CreatingMotionPathsAndFlowlines.ipynb): Using pyGPlates to create motion paths and flowines of points on a tectonic plate to illustrate the plate's trajectory through geological time.
- [__10 - SeafloorGrid__](https://github.com/GPlates/gplately/blob/master/Notebooks/10-SeafloorGrids.ipynb): Defines the parameters needed to set up a `SeafloorGrid` object, and demonstrates how to produce age and spreading rate grids from a set of plate reconstruction model files.

## API Documentation

Documentation of GPlately's objects and methods can be found [here](https://gplates.github.io/gplately/)!
