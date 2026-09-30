.. _gplately-installation:

Installation
============

.. contents::
   :local:
   :depth: 2

Use Micromamba or Conda (recommended)
--------------------------------------

The latest stable public release of **GPlately** can be installed using micromamba_ or conda_ from the `conda-forge channel`_.
**micromamba** is preferred — it resolves environments faster and is much less likely to get stuck than conda — but the same
commands work with conda too; just swap ``conda`` for ``micromamba`` (or vice versa).

The following commands will create a new environment called **my-gplately-env** and install GPlately within it.

.. code:: console

    $ micromamba create -n my-gplately-env
    $ micromamba activate my-gplately-env
    $ micromamba install -c conda-forge gplately

.. note::

    ✏️ Using conda instead? Run the same three commands with ``conda`` in place of ``micromamba``:

    .. code:: console

        $ conda create -n my-gplately-env
        $ conda activate my-gplately-env
        $ conda install -c conda-forge gplately

    If conda gets **stuck while solving the environment**, switching to micromamba_ usually fixes it.

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

----

If you are planning to *develop* GPlately or contribute to the project, you may want to install it from source code. See the next section for instructions.

Install from source code
------------------------

This section covers setting up GPlately for development.

👉 Install GPlately via micromamba first.

.. code:: console

    $ micromamba create -n my-gplately-env
    $ micromamba activate my-gplately-env
    $ micromamba install -c conda-forge gplately

👉 Clone and install GPlately in editable mode.

.. code:: console

    $ git clone https://github.com/GPlates/gplately.git
    $ cd gplately
    $ pip install -e .
    $ pip install -e ".[dev]"   # adds black, isort, bumpver, pip-tools, pytest

.. note::

    ✏️ The ``.[dev]`` extra only brings in the development tools specified in the ``pyproject.toml`` file.

👉 Run the tests.

.. code:: console

    $ python -m pytest -vv tests-dir/pytestcases

Some raster cases are gated behind ``GPLATELY_TEST_LEVEL`` and download large
files, so they are skipped by default. To include them:

.. code:: console

    $ GPLATELY_TEST_LEVEL=100 python -m pytest -vv tests-dir/pytestcases

See ``tests-dir/readme.md`` for more details.

👉 Build the documentation.

The Sphinx sources live in ``sphinx-doc/source``:

.. code:: console

    $ pip install -U sphinx sphinx_rtd_theme
    $ sphinx-autogen -o sphinx-doc/source/generated sphinx-doc/source/*.rst
    $ cd sphinx-doc && make html

.. _`conda-forge channel`: https://conda-forge.org/
.. _conda: https://docs.conda.io/projects/conda/en/latest/index.html
.. _micromamba: https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html
.. _pip: https://pip.pypa.io/en/stable/
.. _PyPI: https://pypi.org/project/gplately/
.. _`GitHub GPlately repository`: https://github.com/GPlates/gplately.git
.. _`Docker README page`: https://github.com/GPlates/gplately/tree/master/docker/README.md
.. _`GPlately commands`: command_line_interface.html
.. _`Jupyter Notebook`: https://jupyter-notebook.readthedocs.io/en/latest/
