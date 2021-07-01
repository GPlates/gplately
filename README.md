# GPlately

GPlately is an object-oriented interface to common pyGPlates and PlateTectonicTools routines.

___This repository is under active developement and the API is liable to change!___

## Dependencies

- [pyGPlates](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation)
- [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools)
- [Shapely](https://shapely.readthedocs.io/en/stable/project.html#installing-shapely)
- NumPy
- SciPy
- Matplotlib
- [Cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html#getting-started) (for mapping)

## Installation

You can install `gplately` using the pip package manager,

```python
pip3 install --user gplately
```

## Usage

GPlately uses objects to accomplish a variety of common tasks. The common objects include:

- `PlateReconstruction` - reconstruct features, tesselate mid ocean ridges, subduction zones
- `Points` - partition points onto plates, rotate back through time
- `Raster` - read in NetCDF grids, interpolation, resampling.
- `PlotTopologies` - one stop shop for plotting ridges, trenches, subduction teeth

Launch a Jupyter environment in the [Notebooks](./Notebooks) directory to get started.


### To-do list

__PlateReconstruction__

- [ ] Properly implement `from_time` in `reconstruct` method
- [ ] Options for saving to GPML files
- [ ] Present subduction data and MOR data (perhaps `Point` objects with multiple attributes

__TimeRaster__

- [ ] Convert raster into points and partition into plates
- [ ] Gridding at reconstruction timesteps
- [ ] Load in an entire timeseries as numpy arrays
- [ ] Save timeseries to NetCDF

__Points__

- [ ] Plate ID attributes
- [ ] Save to csv / xls file

__Parallel__

- [ ] Implement queue system with ordering
- [ ] Flag at initialisation (e.g. `nprocs=4`)
- [ ] Pass in custom functions (in particular for making figures)

__Download__

- [ ] On-the-fly downloader for plate reconstructions (similar to cartopy
- [ ] Download rasters e.g. age grids, spreading rate grids, etc.
