Common Use Cases
================

.. contents::
   :local:
   :depth: 2
   
DataServer
----------

The `DataServer` class allows users to automatically download and cache the necessary files for plate reconstructions to a designated folder on your system.
These files include rotation models, topology features, and static geometries such as coastlines, continents, and continent-ocean boundaries.
Additionally, it supports the retrieval of other data types, including rasters, grids, and feature data.
(Use the newer **PlateModelManager** whenever possible.)

.. code-block:: python
   :linenos:

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


Both `PlateModelManager` and `DataServer` support the following plate reconstruction models:

.. list-table:: Plate Reconstruction Models
   :header-rows: 1
   :align: center
   :width: 100%
   :widths: 30 10 10 10 10 10 10 10 

   * - Model Name
     - Topology
     - Static Polygons
     - Coastlines
     - Continents
     - COB
     - Age Grids
     - SR Grids
   * - Alfonso2024_
     - ✅
     - ✅
     - ✅
     - ✅
     - ❌
     - ❌
     - ❌
   * - Cao2024_
     - ✅
     - ✅
     - ✅
     - ✅
     - ✅
     - ❌
     - ❌

.. _Cao2024: https://doi.org/10.5281/zenodo.11536686
.. _Alfonso2024: https://doi.org/10.5281/zenodo.11392268

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

**Please note that all models have rotation files.**

------------------

PlateReconstruction
-------------------

The `PlateReconstruction` class contains tools to reconstruct geological features like tectonic plates and plate boundaries,
and to interrogate plate kinematic data like plate motion velocities, and rates of subduction and seafloor spreading.

.. code-block:: python
   :linenos:

   from gplately import PlateReconstruction, PlateModelManager

   model = PlateModelManager().get_model("Muller2019")

   # Build a plate reconstruction model using a rotation model, a set of topology features and static polygons
   recon_model = PlateReconstruction(
      model.get_rotation_model(),
      topology_features=model.get_layer("Topologies"),
      static_polygons=model.get_layer("StaticPolygons"),
   )


Alternatively, you may use the auxiliary functions to create a `PlateReconstruction` instance.

.. code-block:: python
   :linenos:

   from gplately.auxiliary import get_plate_reconstruction

   # use the auxiliary function to create a PlateReconstruction instance
   plate_reconstruction_instance = get_plate_reconstruction("Muller2019")


This [02-PlateReconstructions.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb) demonstrates in details
how to use the `PlateReconstruction` class.

Points
------

The methods in the `Points` class track the motion of a point (or group of points) represented by a latitude and longitude
through geologic time. This motion can be visualised using flowlines or motion paths and quantified with point
motion velocities.

.. code-block:: python
   :linenos:

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


The [03-WorkingWithPoints.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/03-WorkingWithPoints.ipynb) demonstrates in details
how to use the `Points` class.

![PointsDemo](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/Reconstructed-Jurassic-Foraminifera-locations-min.png)

The [09-CreatingMotionPathsAndFlowlines.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/09-CreatingMotionPathsAndFlowlines.ipynb) demonstrates how to create motion paths and flowlines.

![motion paths and flowlines](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/Hawaii_Emperor_motion_path.png)

Raster
------

The `Raster` class contains methods to work with netCDF4 or MaskedArray gridded data. Grids may be filled,
resized, resampled, and reconstructed back and forwards through geologic time. Other array data can also be
interpolated onto `Raster` grids.

.. code-block:: python
   :linenos:

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


The [06-Rasters.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/06-Rasters.ipynb) demonstrates in details
how to use the `Raster` class.

![RasterDemo](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/etopo_reconstruction.png)

PlotTopologies
--------------

The `PlotTopologies` class works with the aforementioned `PlateReconstruction` class to plot
geologic features of different types listed
[here](https://gplates.github.io/gplately/plot.html#gplately.plot.PlotTopologies), as well as
coastline, continent and continent-ocean boundary geometries reconstructed through time using pyGPlates.

.. code-block:: python
   :linenos:

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


You may use the auxiliary functions to create a `PlotTopologies` object.

.. code-block:: python
   :linenos:

   from gplately.auxiliary import get_gplot

   # use the auxiliary function to create a PlotTopologies object
   plot_topologies_obj = get_gplot("Muller2019", time=55)


The [02-PlateReconstructions.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb) demonstrates in details
how to use the `PlotTopologies` class.

![PlotTopologiesDemo](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/plottopologies.png)

SeafloorGrid
------------

The `SeafloorGrid` class wraps an automatic workflow to grid seafloor ages and seafloor spreading rates
as encoded by a plate reconstruction model.

.. code-block:: python
   :linenos:

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


The [10-SeafloorGrids.ipynb](https://github.com/GPlates/gplately/blob/master/Notebooks/10-SeafloorGrids.ipynb) is a tutorial notebook that demonstrates
how to set up and use the `SeafloorGrid` object, and shows a sample set of output grids.

![SeafloorGridDemo](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/seafloorgrid.gif)