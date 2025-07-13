## GPlately 2.0.0 Release Notes (2025-07-14)


### How to Upgrade 
* `pip install gplately --upgrade`
* `conda update gplately`
* `docker pull gplates/gplately`

### Breaking Changes
* The `plot_ridges()` method in the `PlotTopologies` class now plots all features labeled as `gpml:MidOceanRidge` in the reconstruction model. Prior to GPlately 2.0.0, this method included an algorithm that attempted to distinguish between "ridges" and "transforms" within gpml:MidOceanRidge features. However, this algorithm had known limitations and could produce inaccurate results. To ensure scientific accuracy, we have removed the algorithm in version 2.0.0. Going forward, it is the responsibility of model creators to explicitly label mid-ocean ridges in their data.
* The `get_ridges()` method in the `PlotTopologies` class now returns all features labeled as `gpml:MidOceanRidge`.
* The `plot_transforms()` now plots all features labeled as `gpml:Transform`.
* The `get_transforms()` now returns all features labeled as `gpml:Transform`.
* The `get_misc_transforms()` method, `plot_misc_transforms()` method and `misc_transforms` property are deprecated. 
* The `plot_plate_id()` method has been renamed as `plot_plate_polygon_by_id()`.
* In the "Muller2016" model, all COB (Continent-Ocean Boundary) features now consist exclusively of polylines. Any polygons previously included have been converted to polylines.
* The "setup.py" has been removed. Instead of using `python setup.py install`, now you need to use `pip install`.
* The `gplately.pygplates` is deprecated. Users now need to `import pygplates` directly.
* Python 3.5, 3.6 and 3.7 support has been dropped.
* Removed Stripy dependency.

### New Features
* The [online documentation](https://gplates.github.io/gplately/) has been overhauled to use "Sphinx", replacing pdoc3. As a result, the online documentation website now features a completely redesigned look.
* A new [new plate-model-manager](https://pypi.org/project/plate-model-manager/) model has been introduced to facilitate the management of reconstruction models.
* A new logging framework has been introduced to assist with debugging and enhance usability.
* An initial implementation leveraging PyGMT for map plotting is now available.
* A new [command line interface (CLI)](https://gplates.github.io/gplately/latest/sphinx/html/command_line_interface.html) has been introduced to enhance usability.
* [New example workflows](https://gplates.github.io/gplately/latest/sphinx/html/examples.html#workflows).
* [New basic examples](https://gplates.github.io/gplately/latest/sphinx/html/examples.html#basics).
* New Raster query and clip functions.
* Add Python 3.12 and 3.13 support.
* Introduced GPlately Docker container images.
* New function to plot poles.
* Added support for boolean grids.
* Added "resize" parameter for reading netcdf grids.
* New `plate_surface_depth()` function.
* New `write_netcdf4()` function for unstructured data.

### Bug Fixes
* [#269](https://github.com/GPlates/gplately/pull/296) Incorrect trench normal angle calculation.
* [#270](https://github.com/GPlates/gplately/issues/270) Incorrect Column Name in tessellate_subduction_zones Function. 
* [#275](https://github.com/GPlates/gplately/issues/275) Agegrid cmdline - continents file inclusion instructions.
* [#286](https://github.com/GPlates/gplately/issues/286) Notify user if continent contouring is activated in the age gridding workflow.
* [#287](https://github.com/GPlates/gplately/issues/287) Cannot find variable z data in netcdf.
* [#289](https://github.com/GPlates/gplately/issues/289) The `raster.reconstruct()` does not work on global netcdfs in 0-360 longitude format.
* [#290](https://github.com/GPlates/gplately/issues/290) The `raster.reconstruct()` ignores anchor_plate_id of PlateReconstruction.  
* [#299](https://github.com/GPlates/gplately/issues/299) The 04-VelocityBasics results in a velocity of 0. 
* [#301](https://github.com/GPlates/gplately/issues/301) Rotating flowlines with gplately.
* [#326](https://github.com/GPlates/gplately/issues/326) Continent masks created by old gplately cannot work with new gplately code.
* [#341](https://github.com/GPlates/gplately/pull/341) Fix spurious lines of polylines/polygons intersecting the dateline. 
* [#346](https://github.com/GPlates/gplately/issues/346) Bug in Raster reconstruction with anchor plate ID.
* [#352](https://github.com/GPlates/gplately/issues/352) Fix slow multi-processing due to pickling pyGPlates objects in PlateReconstruction (and Points and PlotTopologies). 
* [#354](https://github.com/GPlates/gplately/issues/354) Fix subduction teeth occasionally in wrong direction. 

**Note: This "Bug Fixes" list is not exhaustive. See "Other Changes" below.**


### Performance Improvements 
* Improved performance of pickling.
* Improved performance of generating age grids.
* Improved performance of `rotate_reference_frames`.

### Other Changes
See [all changes](https://github.com/GPlates/gplately/compare/v1.3.0...v2.0.0) since GPlately 1.3.0.