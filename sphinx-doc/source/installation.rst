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

If you are planning to *develop* GPlately rather than just use a local checkout,
see `Install from source code`_ below for the prerequisites and the commands to
run the tests and build the documentation.


.. _install-from-source-code:

Install from source code
------------------------

This section covers setting up GPlately for development. GPlately itself is pure
Python, so there is no compiler or native build step — the only native piece is
its ``pygplates`` dependency (see `Prerequisites`_).

Prerequisites
~~~~~~~~~~~~~

- Python 3.10 or newer.
- ``pygplates`` available in the environment **before** installing GPlately.
  ``pyproject.toml`` declares it as ``pygplates>=1.0.0`` and PyPI ships wheels
  for CPython 3.8–3.13, so on those versions a plain ``pip install`` resolves it.
  There is no CPython 3.14 wheel yet (expected with the pyGPlates 1.1 release),
  so on 3.14 install it from the `conda-forge channel`_, from the project's
  `Docker image`_ (see `Use Docker`_), or from a local pyGPlates build.

Clone and install in editable mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: console

    $ git clone https://github.com/GPlates/gplately.git
    $ cd gplately
    $ pip install -e .
    $ pip install -e ".[dev]"   # adds black, isort, bumpver, pip-tools, pytest

.. note::

    ✏️ The ``.[dev]`` extra only brings in the development tools. If your
    ``pygplates`` comes from conda-forge, activate that environment (or
    ``micromamba activate gplately``) first so the editable install lands in it.

Run the tests
~~~~~~~~~~~~~

The suite that CI runs is ``tests-dir/pytestcases``; the folder is called
``tests-dir`` rather than ``tests`` to avoid clashes with Python packaging.

.. code:: console

    $ python -m pytest -vv tests-dir/pytestcases

Some raster cases are gated behind ``GPLATELY_TEST_LEVEL`` and download large
files, so they are skipped by default. To include them:

.. code:: console

    $ GPLATELY_TEST_LEVEL=100 python -m pytest -vv tests-dir/pytestcases

There are further entry points that need a working environment and network
access (see ``tests-dir/readme.md``): ``./tests-dir/test-cli.sh`` for CLI smoke
tests and ``./scripts/run_all_notebooks.sh`` to execute the example notebooks.
The scripts under ``tests-dir/unittest/`` need human/visual verification and are
not part of the regular suite.

The tests pull the plate models and rasters they need at run time via
``plate_model_manager`` instead of using committed fixtures, so expect a
``data-cache/`` and a ``plate-model-repo/`` directory to appear next to your
clone on the first run.

Build the documentation
~~~~~~~~~~~~~~~~~~~~~~~

The Sphinx sources live in ``sphinx-doc/source``:

.. code:: console

    $ pip install -U sphinx sphinx_rtd_theme
    $ sphinx-autogen -o sphinx-doc/source/generated sphinx-doc/source/*.rst
    $ cd sphinx-doc && make html

``scripts/build-sphinx-doc.sh`` wraps the same steps, but it assumes micromamba
and an environment named ``gplately`` — the three commands above are the portable
version.

.. note::

    ✏️ An editable install reports an inexact version number (for example
    ``2.1.0.post13+...``) because the version is derived from the Git state. That
    is expected — GPlately warns about it on import — and not a sign of a broken
    install.

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
.. _`Docker image`: https://hub.docker.com/r/gplates/gplately/tags
.. _`Docker README page`: https://github.com/GPlates/gplately/tree/master/docker/README.md 
.. _`GPlately commands`: command_line_interface.html
.. _`Jupyter Notebook`: https://jupyter-notebook.readthedocs.io/en/latest/ 