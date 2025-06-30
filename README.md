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


GPlately was created to accelerate spatio-temporal data analysis by leveraging [pyGPlates](https://www.gplates.org/docs/pygplates/index.html) within a simplified Python interface. This object-oriented package enables the reconstruction of data through deep geologic time (such as points, lines, polygons, and rasters), the interrogation of plate kinematic information (plate velocities, rates of subduction and seafloor spreading), the rapid comparison of multiple plate motion models, and the plotting of reconstructed output data on maps. All tools are designed to be parallel-safe, accelerating spatio-temporal analysis over multiple CPU processors.

![SeedPointGIF](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/ReadMe_Files/muller19_seedpoints.gif)
 
GPlately can be installed using either `pip` or `conda` (via the conda-forge channel). For detailed installation instructions, please refer to the [Installation](https://gplates.github.io/gplately/latest/sphinx/html/installation.html) section. Additionally, [Docker images](https://gplates.github.io/gplately/latest/sphinx/html/installation.html#use-docker) are available for your convenience.

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

## Documentation

- [latest stable](https://gplates.github.io/gplately/stable/sphinx/html/index.html)
- [latest development](https://gplates.github.io/gplately/latest/sphinx/html/index.html)
- [v2.0.0](https://gplates.github.io/gplately/v2.0.0/sphinx/html/index.html)

### Older versions
- [v1.3.0](https://gplates.github.io/gplately/v1.3.0/)
- [v1.2.0](https://gplates.github.io/gplately/v1.2.0/)
- [v1.1.0](https://gplates.github.io/gplately/v1.1.0/)
- [v1.0.0](https://gplates.github.io/gplately/v1.0.0/)

## Dependencies

- [pyGPlates](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation) 
- [plate-model-manager](https://pypi.org/project/plate-model-manager/) 
- [Shapely](https://shapely.readthedocs.io/en/stable/manual.html)
- [NumPy](https://numpy.org/install/) 
- [SciPy](https://scipy.org/install/) 
- [Matplotlib](https://matplotlib.org/stable/users/installing/index.html)
- [Cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html#getting-started)
- [Pooch](https://github.com/fatiando/pooch)
- [GeoPandas](https://geopandas.org/en/stable/getting_started.html)
- [netCDF4](https://unidata.github.io/netcdf4-python/#quick-install)
- [pygmt](https://www.pygmt.org/latest/)
- [rioxarray](https://github.com/corteva/rioxarray)

