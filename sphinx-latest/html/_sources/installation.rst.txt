Installation
============

.. contents::
   :local:
   :depth: 2
   
conda (recommended)
-------------------

The latest stable public release of **GPlately** can be installed using conda_ from the `conda-forge channel`_. 
The following commands will create a new conda environment called **my-gplately-conda-env** and install GPlately within that environment.

.. code:: console

    $ conda create -n my-gplately-conda-env
    $ conda activate my-gplately-conda-env
    $ conda install -c conda-forge gplately

.. note::
    
    ✏️ If conda gets **stuck while solving the environment** during the installation of GPlately, you can try to use micromamba_ instead.

pip
---

GPlately can also be installed using pip_.

👉 Install the latest stable public release from PyPI_

.. code:: console

    $ pip install gplately


👉 Install from the `GPlately GitHub repository`_ (if you need the latest code changes on GitHub)

.. code:: console

    $ pip install git+https://github.com/GPlates/gplately.git


👉 Install from a local folder (if you need local code changes)

.. code:: console

    $ git clone https://github.com/GPlates/gplately.git gplately.git
    $ cd gplately.git 
    $ git checkout master 
    $ git pull 
    $ MAKE YOUR LOCAL CODE CHANGES HERE ...
    $ pip install -e . 
    
.. note::

    ✏️ The "pip install -e ." command installs GPlately in `editable mode`_.

.. _`editable mode`: https://pip.pypa.io/en/stable/topics/local-project-installs/#editable-installs


docker 🐳
---------

👉 Run GPlately notebooks with Docker

.. code:: console

    $ docker pull gplates/gplately
    $ docker run --rm -ti -p 8888:8888 gplates/gplately

The commands above will start a `Jupyter Notebook`_ server on port 8888. Open this link http://localhost:8888 in a web browser.

👉 Run `GPlately commands`_ with Docker

.. code:: console

    $ docker run gplates/gplately gplately --version
    $ docker run gplates/gplately gplately --help

👉 Run your Python scripts with Docker

.. code:: console

    $ docker run -it --rm -v THE_FULL_PATH_TO_YOUR_SCRIPT_FOLDER:/ws -w /ws gplates/gplately python my_script_to_run.py

.. note::

    ✏️ Replace ``THE_FULL_PATH_TO_YOUR_SCRIPT_FOLDER`` with the full path to the folder containing your script file. 
    In **PowerShell**, you can use "$PWD" if your script is in the current working directory. On **Linux** or **macOS**, you can use \`pwd\` instead.

Visit this `Docker README page`_ for more details about using Docker with GPlately.

.. _`conda-forge channel`: https://conda-forge.org/
.. _conda: https://docs.conda.io/projects/conda/en/latest/index.html
.. _micromamba: https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html
.. _pip: https://pip.pypa.io/en/stable/
.. _PyPI: https://pypi.org/project/gplately/
.. _`GPlately GitHub repository`: https://github.com/GPlates/gplately.git
.. _`Docker README page`: https://github.com/GPlates/gplately/tree/master/docker/README.md 
.. _`GPlately commands`: command_line_interface.html
.. _`Jupyter Notebook`: https://jupyter-notebook.readthedocs.io/en/latest/ 