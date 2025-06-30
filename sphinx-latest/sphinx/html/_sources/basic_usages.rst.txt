Minimal working example
=======================

.. contents::
   :local:
   :depth: 2

Follow the steps below to create your own first paleogeographic map.

.. note::

   You need to know `how to use a terminal`_ to follow the instructions below.

.. _`how to use a terminal`: https://www.freecodecamp.org/news/command-line-for-beginners/

Option 1: Python script
-----------------------

- **Step 1**: use micromamba_ to create a GPlately environment
   
   .. code:: console

    $ micromamba create -n my-gplately-env
    $ micromamba activate my-gplately-env
    $ micromamba install -c conda-forge gplately

.. seealso::
    
   |How to install micromamba?|

- **Step 2**: create a gplately_hello_world.py file and copy & paste the Python code below

   .. code-block:: python
      :linenos:
      :emphasize-lines: 4,10,13

      import cartopy.crs as ccrs
      import matplotlib.pyplot as plt

      import gplately 

      # create a basemap using Mollweide projection
      ax = plt.figure(figsize=(8, 4)).add_subplot(111, projection=ccrs.Mollweide(180))

      # get a PlotTopologies object
      gplot = gplately.auxiliary.get_gplot("Muller2019", time=100) 

      # use the PlotTopologies object to plot a paleo-coastlines map
      gplot.plot_coastlines(ax, color="lightblue",facecolor="0.8")

      # add title for the map
      plt.title(f"{int(gplot.time)} Ma")

      # save the map to a .png file
      plt.gcf().savefig("gplately-hello-world.png")

- **Step 3**: run the Python file and check the map

   .. code:: console

      $ python3 gplately_hello_world.py

   Open the file ``gplately-hello-world.png`` in your current working directory.
   The paleogeographic map created by the code above shows the coastlines at 100 Million years ago.

   .. image:: images/gplately-helloworld.png
      :width: 400
      :alt: GPlately Hello World Map


Option 2: Jupyter Notebook
--------------------------

- **Step 1**: use micromamba_ to create a GPlately environment
   
   .. code:: console

    $ micromamba create -n my-gplately-env
    $ micromamba activate my-gplately-env
    $ micromamba install -c conda-forge gplately jupyter

.. seealso::
    
   |How to install micromamba?|
    
- **Step 2**: start a Jupyter Notebook server

   .. code:: console

      $ jupyter notebook

.. seealso::

   `Jupyter Notebook Documentation`_

.. _`Jupyter Notebook Documentation`: https://jupyter-notebook.readthedocs.io/en/latest/ 

- **Step 3**: create an empty notebook and copy & paste the Python code below

   .. code-block:: python
      :linenos:
      :emphasize-lines: 4,10,13

      import cartopy.crs as ccrs
      import matplotlib.pyplot as plt

      import gplately 

      # create a basemap using Mollweide projection
      ax = plt.figure(figsize=(8, 4)).add_subplot(111, projection=ccrs.Mollweide(180))

      # get a PlotTopologies object
      gplot = gplately.auxiliary.get_gplot("Muller2019", time=100) 

      # use the PlotTopologies object to plot a paleo-coastlines map
      gplot.plot_coastlines(ax, color="lightblue",facecolor="0.8")

      # add title for the map
      plt.title(f"{int(gplot.time)} Ma")

      # display the map
      plt.show()

- **Step 4**: run the notebook and check the map

   The map will be displayed inline within the Jupyter Notebook's web interface.
   The paleogeographic map created by the code above shows the coastlines at 100 Million years ago.

   .. image:: images/gplately-helloworld-notebook.png
      :width: 400
      :alt: GPlately Hello World Map

.. seealso::

   - `GPlately Hello World Python Script`_
   - `GPlately Hello World Jupyter Notebook`_

.. _`GPlately Hello World Python Script`: https://github.com/GPlates/gplately/blob/master/Notebooks/Examples/hello_world.py
.. _`GPlately Hello World Jupyter Notebook`: https://github.com/GPlates/gplately/blob/master/Notebooks/Examples/hello_world.ipynb
.. _micromamba: https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html
.. |How to install micromamba?| raw:: html
   
   <a href="https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html" target="_blank">How to install micromamba?</a>