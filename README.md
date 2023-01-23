# GPlately

GPlately is an object-oriented interface to common pyGPlates and PlateTectonicTools routines.

![SeedPointGIF](./Notebooks/NotebookFiles/ReadMe_Files/muller19_seedpoints.gif)

## Dependencies

- [pyGPlates](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation)
- [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools)
- [Shapely <2.0](https://shapely.readthedocs.io/en/stable/project.html#installing-shapely)
- NumPy
- SciPy 1.10
- Matplotlib
- [Cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html#getting-started) (for mapping)
- Shapely
- [Pooch](https://github.com/fatiando/pooch)
- GeoPandas
- netCDF4


## Installation

You can install `GPlately` using the pip package manager,

```python
pip install gplately
```

... you can also install the most updated version of the `GPlately` repository with pip:

```python
pip install git+https://github.com/GPlates/gplately.git 
```


## Usage

GPlately uses objects to accomplish a variety of common tasks. The common objects include:

- `DataServer` - download rotation files and topology features from plate models on EarthByte's webDAV server
- `PlateReconstruction` - reconstruct features, tesselate mid ocean ridges, subduction zones
- `Points` - partition points onto plates, rotate back through time
- `Raster` - read in NetCDF grids, interpolation, resampling
- `PlotTopologies` - one stop shop for plotting ridges, trenches, subduction teeth


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
Launch the [Plate Reconstruction](./Notebooks/02-PlateReconstructions.ipynb) notebook to see more.


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

![ReconstructionImage](./Notebooks/NotebookFiles/ReadMe_Files/plot_topologies_img.png)

### The `Points` object

... can be used to reconstruct the positions of geological point features and calculate their underlying plate velocities 
through geological time. 

```python
pt_lon = np.array([-107.662152, -58.082792, 17.483189, 133.674590, 80.412876])
pt_lat = np.array([48.797807, -12.654857, 11.884395, -26.415630, 31.368509])

# Call the Points object: pass the PlateReconstruction object, and the latitudes and longitudes of the seed points!
gpts = gplately.Points(model, pt_lon, pt_lat)
```
![PointData](./Notebooks/NotebookFiles/ReadMe_Files/surface_hotspot_plumes.png)


### The `Raster` object

...can be used to read, resample and resize assorted raster data like `netCDF4` seafloor age grids, continental grids and ETOPO
relief rasters. You can also reconstruct raster data back through geological time!

```python
time = 0
agegrid = gdownload.get_age_grid(time)
graster = gplately.Raster(model, array=agegrid, extent=[-180,180,-90,90])
```

Below is a plot of the original Müller et al. (2019) age grid (L) and an age-grid resampled with the `Raster` object (R):

![RasterImg](./Notebooks/NotebookFiles/ReadMe_Files/muller19_raster_resample.png)

## Sample workflows

To see GPlately in action, launch a Jupyter Notebook environment and check out the [sample notebooks](./Notebooks):

- [__01 - Getting Started__](./Notebooks/01-GettingStarted.ipynb): A brief overview of how to initialise GPlately's main objects
- [__02 - Plate Reconstructions__](./Notebooks/02-PlateReconstructions.ipynb): Setting up a `PlateReconstruction` object, reconstructing geological data through time 
- [__03 - Working with Points__](./Notebooks/03-WorkingWithPoints.ipynb): Setting up a `Points` object, reconstructing seed point locations through time with. This notebook uses point data from the Paleobiology Database (PBDB).
- [__04 - Velocity Basics__](./Notebooks/04-VelocityBasics.ipynb): Calculating plate velocities, plotting velocity vector fields
- [__05 - Working with Feature Geometries__](./Notebooks/05-WorkingWithFeatureGeometries.ipynb): Processing and plotting assorted polyline, polygon and point data from [GPlates 2.3's sample data sets](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/)
- [__06 - Rasters__](./Notebooks/06-Rasters.ipynb): Reading, resizing, resampling raster data, and linearly interpolating point data onto raster data
- [__07 - Plate Tectonic Stats__](./Notebooks/07-WorkingWithPlateTectonicStats.ipynb): Using [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools) to calculate and plot subduction zone and ridge data (convergence/spreading velocities, subduction angles, subduction zone and ridge lengths, crustal surface areas produced and subducted etc.) 
- [__08 - Predicting Slab Dip__](./Notebooks/08-PredictingSlabDip.ipynb): Predicting the average slab dip angle of subducting oceanic lithosphere.
- [__09 - Motion Paths and Flowlines__](./Notebooks/09-CreatingMotionPathsAndFlowlines.ipynb): Using pyGPlates to create motion paths and flowines of points on a tectonic plate to illustrate the plate's trajectory through geological time.


## API Documentation

Documentation of GPlately's objects and methods can be found [here](https://gplates.github.io/gplately/)!

