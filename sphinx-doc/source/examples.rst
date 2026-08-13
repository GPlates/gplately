.. _gplately-examples:

Examples
========

It is recommended to use a Conda environment to run these examples.
Download `this YAML file <https://github.com/GPlates/gplately/blob/master/docker/env.yaml>`__ and
run the commands below to start Jupyter Notebook.

.. code:: console

   $ conda env create --name my-gplately-env --file=env.yaml
   $ conda activate my-gplately-env
   $ jupyter notebook

Alternatively, you can also use Docker to run these examples. Use the ``-v`` option to access the local directory on the host machine.
Visit `this page <https://docs.docker.com/engine/storage/bind-mounts/#options-for---volume>`__ for details.

.. code:: console

   $ docker pull gplates/gplately
   $ docker run --rm -ti -v .:/ws -w /ws -p 8888:8888 gplates/gplately

.. contents::
   :local:
   :depth: 2

Workflows
---------

- `01 - Getting Started`_
   A brief overview of how to initialize GPlately's main objects.
- `02 - Plate Reconstructions`_
   Setting up a :py:class:`gplately.PlateReconstruction` object and reconstructing geological data through time.
- `03 - Working with Points`_
   Setting up a :py:class:`gplately.Points` object and reconstructing seed point locations through time.
   This notebook uses point data from the Paleobiology Database (PBDB).
- `04 - Velocity Basics`_
   Calculating plate velocities and plotting velocity vector fields.
- `05 - Working with Feature Geometries`_
   Processing and plotting assorted polyline, polygon, and point data from `GPlates 2.3's sample data sets`_.
- `06 - Rasters`_
   Reading, resizing, and resampling raster data, and linearly interpolating point data onto raster data.
- `07 - Plate Tectonic Stats`_
   Calculating and plotting subduction zone and ridge data (convergence/spreading velocities, subduction angles,
   subduction zone and ridge lengths, crustal surface areas produced and subducted, etc.).
- `08 - Predicting Slab Flux`_
   Predicting the average slab dip angle of subducting oceanic lithosphere.
- `09 - Motion Paths and Flowlines`_
   Using pyGPlates to create motion paths and flowlines for points on a tectonic plate to illustrate the plate's
   trajectory through geological time.
- `10 - Seafloor Grid`_
   Defines the parameters needed to set up a :py:class:`gplately.SeafloorGrid` object and demonstrates
   how to produce age and spreading rate grids from a set of plate reconstruction model files.
- `11 - Andes Fluxes`_
   Demonstrates how the reconstructed subduction history along the Andean margin can potentially be used in
   plate kinematics analysis and data mining.
- `12 - Mutschler World Porphyry Copper Deposits Regional Plots`_
   Generates regional plots for Mutschler world porphyry copper deposits.
- `13 - Reconstructing Zircon Data`_
   Demonstrates how to reconstruct and plot zircon data on a global map through geological time.
- `14 - Rule Based GPML Processing Pipeline`_
   Demonstrates how to use the rule-based GPML processing pipeline to filter and transform geological data.
- `15 - Convert Grid Reference Frame`_
   Demonstrates how to convert the reference frame of grids.

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
.. _`12 - Mutschler World Porphyry Copper Deposits Regional Plots`: ../../notebook-html/12-MutschlerWorldPorphyryCopperDepositsRegionalPlots.html
.. _`13 - Reconstructing Zircon Data`: ../../notebook-html/13-ReconstructingZirconData.html
.. _`14 - Rule Based GPML Processing Pipeline`: ../../notebook-html/14-RuleBasedGPMLProcessingPipeline.html
.. _`15 - Convert Grid Reference Frame`: ../../notebook-html/15-ConvertGridReferenceFrame.html
.. _`GPlates 2.3's sample data sets`: https://www.earthbyte.org/gplates-2-3-software-and-data-sets/

.. note::

   All the `Jupyter Notebook <https://docs.jupyter.org/en/latest/#what-is-a-notebook>`__ files of these sample workflows
   are available `here <https://github.com/GPlates/gplately/tree/master/Notebooks>`__ in the GPlately GitHub repository.


Basics
------

- `Hello World <../../notebook-html/Examples/01-HelloWorld.html>`__
   A minimal working example of GPlately.
- `Use Plate Model Manager <../../notebook-html/Examples/02-PlateModelManager.html>`__
   Use Plate Model Manager to download plate reconstruction models.
- `Plot with Cartopy <../../notebook-html/Examples/03-PlotWithCartopy.html>`__
   Plot a paleo-map using Cartopy.
- `Plot with PyGMT <../../notebook-html/Examples/04-PlotWithPyGMT.html>`__
   Plot a paleo-map using PyGMT.
- `Reconstruct Files <../../notebook-html/Examples/05-ReconstructFiles.html>`__
   Reconstruct and plot shapefiles and `other supported files <https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection>`__.
- `Use Your Own Plate Model <../../notebook-html/Examples/06-LoadPlateModelFromFiles.html>`__
   Use your own plate model to reconstruct points.
- `Save Reconstructed Geometries to Files <../../notebook-html/Examples/07-SaveReconstructedData.html>`__
   Save the reconstructed data to shapefiles.
- `Shortcut to Create PlateReconstruction and PlotTopologies Objects <../../notebook-html/Examples/08-UseAuxiliaryFunctions.html>`__
   An easier way to get PlateReconstruction and PlotTopologies objects from the name of a plate reconstruction model.
- `Generate Icosahedron Mesh <../../notebook-html/Examples/09-IcosahedronMesh.html>`__
   Generate and visualize an icosahedron mesh.

.. note::

   The `Jupyter Notebook <https://docs.jupyter.org/en/latest/#what-is-a-notebook>`__ files of these basic examples
   are available `here <https://github.com/GPlates/gplately/tree/master/Notebooks/Examples>`__.