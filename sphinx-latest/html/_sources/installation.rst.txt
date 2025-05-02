Installation
============

.. contents::
   :local:
   :depth: 2
   
Using conda (recommended)
-------------------------

The latest stable public release of `GPlately` can be installed using conda from the "conda-forge" channel. The following commands will create a new conda environment called "my-gplately-conda-env" and install GPlately within that environment.

.. code:: console

    $ conda create -n my-gplately-conda-env
    $ conda activate my-gplately-conda-env
    $ conda install -c conda-forge gplately


✏️ If `conda` gets __stuck while solving the environment__ during the installation of `GPlately`, you can try using [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) instead.

Using pip
---------

`GPlately` can also be installed using `pip`.

🟢 Install the latest stable public release from [PyPI](https://pypi.org/project/gplately/)

.. code:: console

    $ pip install gplately


🟢 Install from [GitHub repository](https://github.com/GPlates/gplately.git) (if you need the latest code changes on GitHub)

.. code:: console

    $ pip install git+https://github.com/GPlates/gplately.git


🟢 Install from a local folder (if you need local code changes)

.. code:: console

    $ git clone https://github.com/GPlates/gplately.git gplately.git
    $ cd gplately.git # go into the folder created by "git clone" command
    $ git checkout master # check out the "master" branch or the name of branch you want
    $ git pull # fetch all recent code changes from the GitHub remote repository
    # make your local code changes
    $ pip install . # alternatively, you can use "pip install -e ." to install gplately in editable mode


Using Docker 🐳
---------------

👉 Run GPlately notebooks with Docker

.. code:: console

    $ docker pull gplates/gplately
    $ docker run --rm -ti -p 8888:8888  gplates/gplately

http://localhost:8888

👉 Run GPlately command with Docker

.. code:: console

    $ docker run gplates/gplately gplately --version
    $ docker run gplates/gplately gplately --help

👉 Run your Python script with Docker

.. code:: console

    $ docker run -it --rm -v THE_FULL_PATH_TO_YOUR_SCRIPT_FOLDER:/ws -w /ws gplates/gplately python my_script_to_run.py

✏️ Replace ``THE_FULL_PATH_TO_YOUR_SCRIPT_FOLDER`` with the full path to the folder containing your script file. In PowerShell, you can use "$PWD"  if your script is in the current working directory. On Linux or macOS, you can use \`pwd\` instead.

Visit [this page](https://github.com/GPlates/gplately/tree/master/docker/README.md) for more details about using Docker with GPlately.