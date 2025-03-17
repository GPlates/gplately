<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/GPlately_White_logo.png">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/GPlately_Main_logo.png">
  <img alt="GPlately logo." src="https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/GPlately_Main_logo.png">
</picture>
</p>

![PyPI](https://img.shields.io/pypi/v/gplately?style=for-the-badge)
![PyPI - Downloads](https://img.shields.io/pypi/dm/gplately?style=for-the-badge)
![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/gplately?style=for-the-badge)
![downloads](https://img.shields.io/conda/d/conda-forge/gplately?style=for-the-badge)
![GitHub Unitttest Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/GPlates/gplately/build_and_test.yml?branch=master&style=for-the-badge&label=test)
![GitHub Build Doc Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/GPlates/gplately/deploy_documentation.yaml?branch=master&style=for-the-badge&label=doc)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/gplately?style=for-the-badge)


GPlately was created to accelerate spatio-temporal data analysis by leveraging [pyGPlates](https://www.gplates.org/docs/pygplates/index.html) and [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools) within a simplified Python interface. This object-oriented package enables the reconstruction of data through deep geologic time (such as points, lines, polygons, and rasters), the interrogation of plate kinematic information (plate velocities, rates of subduction and seafloor spreading), the rapid comparison of multiple plate motion models, and the plotting of reconstructed output data on maps. All tools are designed to be parallel-safe, accelerating spatio-temporal analysis over multiple CPU processors.

![SeedPointGIF](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/muller19_seedpoints.gif)
 
GPlately can be installed using either `pip` or `conda` (via the conda-forge channel). For detailed installation instructions, please refer to the [Installation](https://github.com/GPlates/gplately/tree/update-doc-examples-tests?tab=readme-ov-file#installation) section. Additionally, [Docker images](https://github.com/GPlates/gplately/tree/update-doc-examples-tests?tab=readme-ov-file#3-using-docker-) are available for your convenience.

Sample data is available from [EarthByte servers](https://www.earthbyte.org/category/resources/), which
include rasters, seafloor age grids, rotation files, and more to help you get started with plate reconstructions.

#### Citation

> Mather, B.R., Müller, R.D., Zahirovic, S., Cannon, J., Chin, M., Ilano, L., Wright, N.M., Alfonso, C., Williams, S., Tetley, M., Merdith, A. (2023) Deep time spatio-temporal data analysis using pyGPlates with PlateTectonicTools and GPlately. _Geoscience Data Journal_, 1–8. Available from: https://doi.org/10.1002/gdj3.185

```bib
@article{Mather2023,
author = {Mather, Ben R. and Müller, R. Dietmar and Zahirovic, Sabin and Cannon, John and Chin, Michael and Ilano, Lauren and Wright, Nicky M. and Alfonso, Christopher and Williams, Simon and Tetley, Michael and Merdith, Andrew},
title = {Deep time spatio-temporal data analysis using pyGPlates with PlateTectonicTools and GPlately},
year = {2023},
journal = {Geoscience Data Journal},
pages = {1-8},
keywords = {geospatial, plate reconstructions, pyGPlates, python, tectonics},
doi = {https://doi.org/10.1002/gdj3.185},
url = {https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/gdj3.185},
eprint = {https://rmets.onlinelibrary.wiley.com/doi/pdf/10.1002/gdj3.185},
}
```

## Usage

- [Quick start](https://gplates.github.io/gplately/dev-doc/#quick-start) - a brief tutorial to help users get up to speed  
- [Sample workflows](https://github.com/GPlates/gplately?tab=readme-ov-file#sample-workflows) - demonstrations of how GPlately can be used in various workflows
- [Examples](https://github.com/GPlates/gplately/blob/master/Notebooks/Examples/readme.md) - reusable code snippets demonstrating how to use GPlately for specific tasks
- [command-line interface (CLI)](https://github.com/GPlates/gplately/blob/master/gplately/commands/readme.md) - use GPlately in a command-line environment

## API Reference

- [v1.3.0 (**latest stable**)](https://gplates.github.io/gplately/v1.3.0/)
- [latest in dev](https://gplates.github.io/gplately/dev-doc)
- [v1.2.0](https://gplates.github.io/gplately/v1.2.0/)
- [v1.1.0](https://gplates.github.io/gplately/v1.1.0/)
- [v1.0.0](https://gplates.github.io/gplately/v1.0.0/)

## Installation

### 1. Using conda (recommended)

The latest stable public release of `GPlately` can be installed using conda from the "conda-forge" channel. The following commands will create a new conda environment called "my-gplately-conda-env" and install GPlately within that environment.

```sh
conda create -n my-gplately-conda-env
conda activate my-gplately-conda-env
conda install -c conda-forge gplately
```

✏️ If `conda` gets __stuck while solving the environment__ during the installation of `GPlately`, you can try using [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) instead.

### 2. Using pip

`GPlately` can also be installed using `pip`.

🟢 Install the latest stable public release from [PyPI](https://pypi.org/project/gplately/)

```sh
pip install gplately
```

🟢 Install from [GitHub repository](https://github.com/GPlates/gplately.git) (if you need the latest code changes on GitHub)

```sh
pip install git+https://github.com/GPlates/gplately.git
```

🟢 Install from a local folder (if you need local code changes)

```sh
git clone https://github.com/GPlates/gplately.git gplately.git
cd gplately.git # go into the folder created by "git clone" command
git checkout master # check out the "master" branch or the name of branch you want
git pull # fetch all recent code changes from the GitHub remote repository
# make your local code changes 
pip install . # alternatively, you can use "pip install -e ." to install gplately in editable mode
```

### 3. Using Docker 🐳

👉 Run GPlately notebooks with Docker

- `docker pull gplates/gplately`
- `docker run --rm -ti -p 8888:8888  gplates/gplately`
- http://localhost:8888

👉 Run GPlately command with Docker

- `docker run gplates/gplately gplately --version`
- `docker run gplates/gplately gplately --help`

👉 Run your Python script with Docker

- `docker run -it --rm -v THE_FULL_PATH_TO_YOUR_SCRIPT_FOLDER:/ws -w /ws gplates/gplately python my_script_to_run.py` 

✏️ Replace __THE_FULL_PATH_TO_YOUR_SCRIPT_FOLDER__ with the full path to the folder containing your script file. In PowerShell, you can use "$PWD"  if your script is in the current working directory. On Linux or macOS, you can use \`pwd\` instead.

Visit [this page](https://github.com/GPlates/gplately/tree/master/docker/README.md) for more details about using Docker with GPlately.

## Dependencies

- [pyGPlates](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation) 
- [plate-model-manager](https://pypi.org/project/plate-model-manager/) 
- [Shapely](https://shapely.readthedocs.io/en/stable/project.html#installing-shapely)
- [NumPy](https://numpy.org/install/) 
- [SciPy](https://scipy.org/install/) 
- [Matplotlib](https://matplotlib.org/stable/users/installing/index.html)
- [Cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html#getting-started)
- [Shapely](https://shapely.readthedocs.io/en/stable/installation.html)
- [Pooch](https://github.com/fatiando/pooch)
- [GeoPandas](https://geopandas.org/en/stable/getting_started.html)
- [netCDF4](https://unidata.github.io/netcdf4-python/#quick-install)
- [pygmt](https://www.pygmt.org/latest/)
- [rioxarray](https://github.com/corteva/rioxarray)

