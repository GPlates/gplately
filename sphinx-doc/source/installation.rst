.. _gplately-installation:

Installation
============

.. contents::
   :local:
   :depth: 2
   
Use Conda (recommended)
-----------------------

The latest stable public release of **GPlately** can be installed using conda_ from the `conda-forge channel`_. 
The following commands will create a new conda environment called **my-gplately-conda-env** and install GPlately within that environment.

.. code:: console

    $ conda create -n my-gplately-conda-env
    $ conda activate my-gplately-conda-env
    $ conda install -c conda-forge gplately

.. note::
    
    If conda gets **stuck while solving the environment** during the installation of GPlately, you can try to use micromamba_ instead.

Use Pip
-------

GPlately can also be installed using pip_.

👉 Install the latest stable public release from PyPI_.

.. code:: console

    $ pip install gplately


👉 Install from the `GitHub GPlately repository`_ (if you need the latest code changes on GitHub).

.. code:: console

    $ pip install git+https://github.com/GPlates/gplately.git


👉 Install from a local folder (if you need local code changes).

.. code:: console

    $ git clone https://github.com/GPlates/gplately.git gplately.git
    $ cd gplately.git 
    $ git checkout master 
    $ git pull 
    $ MAKE YOUR LOCAL CODE CHANGES HERE ...
    $ pip install -e . 
    
.. note::

    ✏️ The ``pip install -e .`` command installs GPlately in `editable mode`_.

.. _`editable mode`: https://pip.pypa.io/en/stable/topics/local-project-installs/#editable-installs


Use Docker
----------

👉 Run GPlately notebooks within a Docker container.

.. code:: console

    $ docker pull gplates/gplately
    $ docker run --rm -ti -p 8888:8888 gplates/gplately

The commands above will start a `Jupyter Notebook`_ server on port 8888. Open this link http://localhost:8888 in a web browser.

👉 Run `GPlately commands`_ within a Docker container.

.. code:: console

    $ docker run gplates/gplately gplately --version
    $ docker run gplates/gplately gplately --help

👉 Run your Python scripts within a Docker container.

.. code:: console

    $ docker run -it --rm -v THE_FULL_PATH_TO_YOUR_SCRIPT_FOLDER:/ws -w /ws gplates/gplately python my_script_to_run.py

.. note::

    Replace ``THE_FULL_PATH_TO_YOUR_SCRIPT_FOLDER`` with the full path to the folder containing your script file. 
    In **PowerShell**, you can use "$PWD" if your script is in the current working directory. On **Linux** or **macOS**, you can use \`pwd\` instead.

In certain shell environments, using a relative path may also work. You can try the following command to mount your current working directory to ``/ws``.

.. code:: console

    $ docker run -it --rm -p 8888:8888 -v .:/ws -w /ws gplates/gplately

GPlately Docker images are available at both `Docker Hub <https://hub.docker.com/>`__ and `GitHub Container Registry <https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry>`__.

- https://hub.docker.com/r/gplates/gplately/tags
- https://github.com/GPlates/gplately/pkgs/container/gplately 

Visit this `Docker README page`_ for more details about using Docker with GPlately.

.. _`conda-forge channel`: https://conda-forge.org/
.. _conda: https://docs.conda.io/projects/conda/en/latest/index.html
.. _micromamba: https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html
.. _pip: https://pip.pypa.io/en/stable/
.. _PyPI: https://pypi.org/project/gplately/
.. _`GitHub GPlately repository`: https://github.com/GPlates/gplately.git
.. _`Docker README page`: https://github.com/GPlates/gplately/tree/master/docker/README.md 
.. _`GPlately commands`: command_line_interface.html
.. _`Jupyter Notebook`: https://jupyter-notebook.readthedocs.io/en/latest/ 