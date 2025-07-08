Common Use Cases
================

.. contents::
   :local:
   :depth: 2
   
DataServer
----------

The :py:class:`gplately.DataServer` class allows users to automatically download and cache the plate reconstruction model files to a designated folder 
on your system. These files include rotation models, topology features, and static geometries such as coastlines, continents, and 
continent-ocean boundaries. Additionally, it supports the retrieval of other data types, including rasters, grids, and feature data.

The plate-model-manager_ Python package is a newly developed alternative to the :py:class:`gplately.DataServer` class. 
Users are encouraged to use the newer plate-model-manager_ whenever possible.

.. _plate-model-manager: https://github.com/michaelchin/plate-model-manager

.. seealso::

   `plate-model-manager documentation`_

.. _`plate-model-manager documentation`: https://michaelchin.github.io/plate-model-manager/latest/

.. code-block:: python
   :linenos:
   :emphasize-lines: 3

   from gplately import DataServer

   data_server = DataServer("Muller2019")

   # Download plate reconstruction files and geometries from the Müller et al. 2019 model
   rotation_model, topology_features, static_polygons = (
      data_server.get_plate_reconstruction_files()
   )
   coastlines, continents, COBs = data_server.get_topology_geometries()

   # Download the age grid at 100Ma from the Müller et al. 2019 model
   age_grid = data_server.get_age_grid(times=100)

   # Download the ETOPO1 geotiff raster
   etopo = data_server.get_raster("ETOPO1_tif")


The table below provides a list of all available plate reconstruction models.

.. note::

      - Use horizontal scrolling to see all the columns in the table.
      - All models include rotation data.

.. list-table:: Plate Reconstruction Models
   :header-rows: 1
   :align: left
   :width: 100%
   :widths: 30 10 10 10 10 10 10 10 

   * - Model Name
     - Topology
     - Static Polygons
     - Coasts
     - Continents
     - COB
     - Age Grids
     - SR Grids
   * - Shirmard2025_
     - ✅
     - ✅
     - ✅
     - ✅
     - ✅
     - ❌
     - ❌
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
   * - Muller2022_ 
     - ✅
     - ✅ 
     - ✅ 
     - ✅
     - ✅
     - ❌ 
     - ❌
   * - Zahirovic2022_ 
     - ✅
     - ✅ 
     - ✅ 
     - ✅ 
     - ❌ 
     - ✅ 
     - ✅ 
   * - Merdith2021_ 
     - ✅ 
     - ✅ 
     - ✅  
     - ✅ 
     - ❌ 
     - ❌
     - ❌  
   * - Clennett2020_ 
     - ✅ 
     - ✅
     - ✅ 
     - ❌
     - ❌
     - ✅ 
     - ✅
   * - Clennett2020_M2019_
     - ✅ 
     - ✅ 
     - ✅ 
     - ❌ 
     - ❌ 
     - ✅ 
     - ✅ 
   * - Clennett2020_S2013_
     - ✅
     - ✅ 
     - ✅ 
     - ❌ 
     - ❌ 
     - ❌ 
     - ❌ 
   * - Muller2019_
     - ✅
     - ✅
     - ✅ 
     - ✅ 
     - ✅  
     - ✅ 
     - ❌ 
   * - Young2018_
     - ✅ 
     - ✅ 
     - ✅ 
     - ✅ 
     - ❌ 
     - ❌ 
     - ❌ 
   * - TorsvikCocks2017
     - ❌ 
     - ✅ 
     - ✅  
     - ❌ 
     - ❌ 
     - ❌
     - ❌
   * - Matthews2016_ 
     - ✅ 
     - ✅ 
     - ✅ 
     - ✅ 
     - ❌
     - ❌
     - ❌ 
   * - Matthews2016_pmag_ref
     - ❌ 
     - ✅ 
     - ✅ 
     - ✅ 
     - ❌
     - ❌
     - ❌ 
   * - Muller2016_
     - ✅ 
     - ✅
     - ✅ 
     - ❌ 
     - ✅ 
     - ✅ 
     - ❌ 
   * - Scotese2016_ 
     - ❌ 
     - ✅ 
     - ✅ 
     - ❌ 
     - ❌ 
     - ❌
     - ❌ 
   * - Zahirovic2016_
     - ✅ 
     - ✅ 
     - ✅ 
     - ✅ 
     - ❌ 
     - ❌ 
     - ❌ 
   * - Gibbons2015_ 
     - ✅
     - ✅ 
     - ✅ 
     - ❌ 
     - ❌
     - ❌ 
     - ❌ 
   * - Zahirovic2014_
     - ❌
     - ✅
     - ✅ 
     - ❌ 
     - ❌ 
     - ❌ 
     - ❌ 
   * - Shephard2013_
     - ✅
     - ✅ 
     - ✅ 
     - ❌ 
     - ❌ 
     - ❌
     - ❌
   * - Gurnis2012_
     - ✅
     - ✅ 
     - ✅ 
     - ❌
     - ❌ 
     - ❌ 
     - ❌
   * - Seton2012_
     - ✅ 
     - ✅ 
     - ✅ 
     - ✅ 
     - ✅ 
     - ✅ 
     - ❌ 
   * - Muller2008
     - ❌ 
     - ✅ 
     - ❌ 
     - ❌ 
     - ❌ 
     - ❌
     - ❌ 

.. _Shirmard2025: https://zenodo.org/records/15233548
.. _Cao2024: https://doi.org/10.5281/zenodo.11536686
.. _Alfonso2024: https://doi.org/10.5281/zenodo.11392268
.. _Muller2022: https://doi.org/10.5281/zenodo.10297173
.. _Zahirovic2022: https://zenodo.org/records/4729045
.. _Merdith2021: https://doi.org/10.5281/zenodo.10346399
.. _Clennett2020: https://doi.org/10.5281/zenodo.10348270
.. _Clennett2020_M2019: https://doi.org/10.5281/zenodo.10348270
.. _Clennett2020_S2013: https://doi.org/10.5281/zenodo.10348270
.. _Muller2019: https://doi.org/10.5281/zenodo.10525286
.. _Young2018: https://doi.org/10.5281/zenodo.10525369
.. _Matthews2016: https://doi.org/10.5281/zenodo.10526156
.. _Muller2016: https://doi.org/10.5281/zenodo.10565444
.. _Scotese2016: https://doi.org/10.5281/zenodo.10596609
.. _Zahirovic2016: https://doi.org/10.5281/zenodo.10531296
.. _Gibbons2015: https://doi.org/10.5281/zenodo.10595658
.. _Zahirovic2014: https://doi.org/10.5281/zenodo.10595658
.. _Shephard2013: https://doi.org/10.5281/zenodo.10595888
.. _Gurnis2012: https://doi.org/10.5281/zenodo.10596349
.. _Seton2012: https://doi.org/10.5281/zenodo.10596049

.. note::

   - ``Topology``: topological plate polygons and plate boundaries 
   - ``Coasts``: coastlines. The Greenland coastline uses the Danish Geological Survey dataset. The rest of the world uses the World Vector Shoreline from the Global Self-consistent Hierarchical High-resolution Geography (GSHHG) dataset.
   - ``Continents``: continental polygons, derived from the Static Polygons, containing continental crust and volcanically-modified oceanic crust (including island arcs). 
   - ``COB``: continent-ocean boundary. The COBs are represented as lines along passive margins and do not include data from active margins.
   - ``Age Grids``: numerical grid of seafloor age
   - ``SR Grids``: numerical grid of seafloor spreading rate


PlateReconstruction
-------------------

The :py:class:`gplately.PlateReconstruction` class contains tools to reconstruct geological features like tectonic plates and plate boundaries,
and to interrogate plate kinematic data like plate motion velocities, and rates of subduction and seafloor spreading.

.. code-block:: python
   :linenos:
   :emphasize-lines: 6

   from gplately import PlateReconstruction, PlateModelManager

   model = PlateModelManager().get_model("Muller2019")

   # Build a plate reconstruction model using a rotation model, a set of topology features and static polygons
   recon_model = PlateReconstruction(
      model.get_rotation_model(),
      topology_features=model.get_layer("Topologies"),
      static_polygons=model.get_layer("StaticPolygons"),
   )


Alternatively, you may use the auxiliary functions to create a :py:class:`gplately.PlateReconstruction` object.

.. code-block:: python
   :linenos:
   :emphasize-lines: 4

   from gplately.auxiliary import get_plate_reconstruction

   # use the auxiliary function to create a PlateReconstruction object
   plate_reconstruction_instance = get_plate_reconstruction("Muller2019")


The `PlateReconstructions example`_ demonstrates in detail how to use the :py:class:`gplately.PlateReconstruction` class.
The `02-PlateReconstructions.ipynb`_ Jupyter Notebook is available in the GitHub GPlately repository.

.. _`02-PlateReconstructions.ipynb`: https://github.com/GPlates/gplately/blob/master/Notebooks/02-PlateReconstructions.ipynb
.. _`PlateReconstructions example`: https://gplates.github.io/gplately/stable/02-PlateReconstructions.html

Points
------

The methods in the :py:class:`gplately.Points` class track the motion of a point (or group of points) represented by a latitude and longitude
through geologic time. This motion can be visualised using flowlines or motion paths and quantified with point motion velocities.

.. code-block:: python
   :linenos:
   :emphasize-lines: 15

   import numpy as np

   from gplately import PlateModelManager, Points, auxiliary

   model = PlateModelManager().get_model("Muller2019")

   # Create a plate reconstruction model using a rotation model, a set of topology features and static polygons
   recon_model = auxiliary.get_plate_reconstruction(model)

   # Define some points using their latitude and longitude coordinates so we can track them through time!
   pt_lons = np.array([140.0, 150.0, 160.0])
   pt_lats = np.array([-30.0, -40.0, -50.0])

   # Create a Points object from these points
   gpts = Points(recon_model, pt_lons, pt_lats)


The `WorkingWithPoints example`_ demonstrates in detail how to use the Points class. 
The `03-WorkingWithPoints.ipynb`_ Jupyter Notebook is available in the GitHub GPlately repository.

.. _`WorkingWithPoints example`: https://gplates.github.io/gplately/stable/03-WorkingWithPoints.html
.. _`03-WorkingWithPoints.ipynb`: https://github.com/GPlates/gplately/blob/master/Notebooks/03-WorkingWithPoints.ipynb

.. image:: images/Reconstructed-Jurassic-Foraminifera-locations-min.png
      :width: 600
      :alt: PointsDemo

The `CreatingMotionPathsAndFlowlines example`_ demonstrates how to create motion paths and flowlines.
The `09-CreatingMotionPathsAndFlowlines.ipynb`_ Jupyter Notebook is available in the GitHub GPlately repository.

.. _`CreatingMotionPathsAndFlowlines example`:
.. _`09-CreatingMotionPathsAndFlowlines.ipynb`: https://github.com/GPlates/gplately/blob/master/Notebooks/09-CreatingMotionPathsAndFlowlines.ipynb

.. image:: images/Hawaii_Emperor_motion_path.png
      :width: 600
      :alt: motion paths and flowlines

Raster
------

The :py:class:`gplately.Raster` class contains methods to work with netCDF4 or MaskedArray gridded data. Grids may be filled,
resized, resampled, and reconstructed back and forwards through geologic time. Other array data can also be
interpolated onto Raster grids.

.. code-block:: python
   :linenos:
   :emphasize-lines: 8, 16

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


The `Rasters example`_ demonstrates in detail how to use the :py:class:`gplately.Raster` class. 
The `06-Rasters.ipynb`_ Jupyter Notebook is available in the GitHub GPlately repository.

.. _`06-Rasters.ipynb`: https://github.com/GPlates/gplately/blob/master/Notebooks/06-Rasters.ipynb
.. _`Rasters example`: https://gplates.github.io/gplately/stable/06-Rasters.html

.. image:: images/etopo_reconstruction.png
      :width: 600
      :alt: RasterDemo

PlotTopologies
--------------

The :py:class:`gplately.PlotTopologies` class works with the aforementioned :py:class:`gplately.PlateReconstruction` class to plot
geologic features of different types, such as coastlines, continents and continent-ocean boundaries reconstructed through time using pyGPlates.

.. code-block:: python
   :linenos:
   :emphasize-lines: 6

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


You may use the auxiliary functions to create a :py:class:`gplately.PlotTopologies` object.

.. code-block:: python
   :linenos:
   :emphasize-lines: 4

   from gplately.auxiliary import get_gplot

   # use the auxiliary function to create a PlotTopologies object
   plot_topologies_obj = get_gplot("Muller2019", time=55)

The `PlateReconstructions example`_ demonstrates in detail how to use the :py:class:`gplately.PlotTopologies` class.
The `02-PlateReconstructions.ipynb`_ Jupyter Notebook is available in the GitHub GPlately repository.

.. image:: images/plottopologies.png
      :width: 600
      :alt: PlotTopologiesDemo

SeafloorGrid
------------

The :py:class:`gplately.SeafloorGrid` class wraps an automatic workflow to grid seafloor ages and seafloor spreading rates
as encoded by a plate reconstruction model.

.. code-block:: python
   :linenos:
   :emphasize-lines: 11, 20

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

The `SeafloorGrids example`_ is a tutorial notebook that demonstrates
how to set up and use the :py:class:`gplately.SeafloorGrid` object, and shows a sample set of output grids. 
The `10-SeafloorGrids.ipynb`_ Jupyter Notebook is available in the GitHub GPlately repository.

.. _`SeafloorGrids example`: https://gplates.github.io/gplately/dev-doc/10-SeafloorGrids.html
.. _`10-SeafloorGrids.ipynb`: https://github.com/GPlates/gplately/blob/master/Notebooks/10-SeafloorGrids.ipynb

.. image:: images/seafloorgrid.gif
      :width: 600
      :alt: SeafloorGridDemo