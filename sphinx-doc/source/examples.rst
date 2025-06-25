.. _gplately-examples:

Examples
========

.. contents::
   :local:
   :depth: 2

Workflows
---------

- `01 - Getting Started`_
   A brief overview of how to initialise GPlately's main objects
- `02 - Plate Reconstructions`_ 
   Setting up a :py:class:`gplately.PlateReconstruction` object, reconstructing geological data through time.
- `03 - Working with Points`_ 
   Setting up a :py:class:`gplately.Points` object, reconstructing seed point locations through time with. 
   This notebook uses point data from the Paleobiology Database (PBDB).
- `04 - Velocity Basics`_ 
   Calculating plate velocities, plotting velocity vector fields.
- `05 - Working with Feature Geometries`_ 
   Processing and plotting assorted polyline, polygon and point data from `GPlates 2.3's sample data sets`_.
- `06 - Rasters`_ 
   Reading, resizing, resampling raster data, and linearly interpolating point data onto raster data.
- `07 - Plate Tectonic Stats`_ 
   Calculating and plotting subduction zone and ridge data (convergence/spreading velocities, subduction angles, 
   subduction zone and ridge lengths, crustal surface areas produced and subducted etc.).
- `08 - Predicting Slab Flux`_ 
   Predicting the average slab dip angle of subducting oceanic lithosphere.
- `09 - Motion Paths and Flowlines`_ 
   Using pyGPlates to create motion paths and flowines of points on a tectonic plate to illustrate the plate's 
   trajectory through geological time.
- `10 - Seafloor Grid`_   
   Defines the parameters needed to set up a :py:class:`gplately.SeafloorGrid` object, and demonstrates 
   how to produce age and spreading rate grids from a set of plate reconstruction model files.
- `11 - Andes Fluxes`_ 
   Demonstrates how the reconstructed subduction history along the Andean margin can be potentially 
   used in the plate kinematics analysis and data mining.
- `12 - Mutschler World Porphyry Copper Deposits Regional Plots`_ 
   Generates regional plots for Mutschler world porphyry copper deposits.

.. _`01 - Getting Started`: ../../notebook-html/01-GettingStarted.html
.. _`02 - Plate Reconstructions`: ../../notebook-html/02-PlateReconstructions.html
.. _`03 - Working with Points`: ../../notebook-html/03-WorkingWithPoints.html
.. _`04 - Velocity Basics`: ../../notebook-html/04-VelocityBasics.html
.. _`05 - Working with Feature Geometries`: ../../notebook-html/05-WorkingWithFeatureGeometries.html
.. _`06 - Rasters`: ../../notebook-html/06-Rasters.html
.. _`07 - Plate Tectonic Stats`: ../../notebook-html/07-WorkingWithPlateTectonicStats.html
.. _`08 - Predicting Slab Flux`: ../../notebook-html/08-PredictingSlabFlux.html
.. _`09 - Motion Paths and Flowlines`: ../../notebook-html/09-CreatingMotionPathsAndFlowlines.html
.. _`10 - Seafloor Grid`: ../../notebook-html/10-SeafloorGrids.html
.. _`11 - Andes Fluxes`: ../../notebook-html/11-AndesFluxes.html
.. _`12 - Mutschler World Porphyry Copper Deposits Regional Plots`: ../notebook-html/12-MutschlerWorldPorphyryCopperDepositsRegionalPlots.html
.. _`GPlates 2.3's sample data sets`: https://www.earthbyte.org/gplates-2-3-software-and-data-sets/

.. note::

   All the `Jupyter Notebook <https://docs.jupyter.org/en/latest/#what-is-a-notebook>`__ files of these sample workflows 
   are available `here <https://github.com/GPlates/gplately/tree/master/Notebooks>`__ in the GPlately GitHub repository.


Basics
------

- `Hello World <../../notebook-html/hello_world.html>`__ 
   A minimal working example of GPlately.
- `Use Plate Model Manager <../../notebook-html/introducing_plate_model_manager.html>`__
   How to use Plate Model Manager to download plate reconstruction models.
- `Plot with Cartopy <../../notebook-html/plot_map_with_cartopy.html>`__
   Plot paleo-map using Cartopy.
- `Plot with PyGMT <../../notebook-html/plot_map_with_pygmt.html>`__
   Plot paleo-map using PyGMT.
- `Save Reconstructed Geometries to Files <../../notebook-html/save_reconstructed_data.html>`__
   How to save the reconstructed data to shapefiles.
- `Shortcut to create PlateReconstruction and PlotTopologies objects <../../notebook-html/use_auxiliary_functions.html>`__
   Easier way to get PlateReconstruction and PlotTopologies objects from the name of a plate reconstruction model.