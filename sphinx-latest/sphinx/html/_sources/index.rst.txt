GPlately v2.0.0 documentation
=============================

GPlately_ was created to accelerate spatio-temporal data analysis by leveraging pyGPlates_ within a simplified Python interface. 
This object-oriented package enables the reconstruction of data through deep geologic time (such as points, lines, polygons, and rasters), 
the interrogation of plate kinematic information (plate velocities, rates of subduction and seafloor spreading), 
the rapid comparison of multiple plate motion models, and the plotting of reconstructed output data on maps. 
All tools are designed to be parallel-safe, accelerating spatio-temporal analysis over multiple CPU processors.

GPlately can be installed using either pip_ or conda_ (via the `conda-forge channel`_). 
See the `installation section`_ for details. Additionally, |Docker images| are available for your convenience.

Sample data is available from |EarthByte servers|, which include rasters, seafloor age grids, rotation files, 
and more to help you get started with plate reconstructions.

.. toctree::
   :caption: User Guide
   :maxdepth: 2

   installation
   basic_usages
   use_cases
   faq
   examples
   command_line_interface

.. toctree::
   :caption: Reference
   :maxdepth: 2

   api
   functions
   secondaries

Index
=====

* :ref:`genindex`

.. |EarthByte servers| raw:: html
   
   <a href="https://www.earthbyte.org/category/resources/" target="_blank">EarthByte servers</a>

.. _installation section: installation.html
.. _GPlately: https://github.com/GPlates/gplately
.. _pyGPlates: https://www.gplates.org/docs/pygplates/
.. _pip: https://pip.pypa.io/en/stable/
.. _conda: https://docs.conda.io/projects/conda/en/latest/index.html
.. _`conda-forge channel`: https://conda-forge.org/
.. |Docker images| raw:: html

   <a href="https://github.com/GPlates/gplately/blob/master/docker/README.md" target="_blank">Docker images</a>
