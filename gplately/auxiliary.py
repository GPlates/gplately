#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

"""A set of helper functions designed to streamline the use of GPlately’s functionalities,
minimizing the coding effort required from users."""

from typing import Union

import pygmt
from plate_model_manager import PlateModel, PlateModelManager

from .download import path_to_cache
from .mapping.cartopy_plot import CartopyPlotEngine
from .mapping.plot_engine import PlotEngine
from .plot import PlotTopologies
from .reconstruction import PlateReconstruction


def get_plate_reconstruction(model: Union[str, PlateModel], model_repo_dir: str = "./"):
    """Return a :py:class:`gplately.PlateReconstruction` object for a given model name or :class:`gplately.PlateModel` object.

    Parameters
    ----------
    model : str or PlateModel
        model name or a :class:`gplately.PlateModel` object
    model_repo_dir: str, default="./"
        the folder in which you would like to keep the model files

    Returns
    -------
    PlateReconstruction
        a :class:`gplately.PlateReconstruction` object


    .. seealso::

        `usage example <https://github.com/GPlates/gplately/blob/master/Notebooks/Examples/use_auxiliary_functions.py>`__
    """
    if isinstance(model, str):
        model_name: str = model
        try:
            plate_model = PlateModelManager().get_model(
                model_name, data_dir=model_repo_dir
            )
        except:
            plate_model = PlateModel(model_name, data_dir=model_repo_dir, readonly=True)

        if plate_model is None:
            raise Exception(f"Unable to get model ({model_name})")
    else:
        plate_model = model

    topology_features = None
    static_polygons = None

    all_layers = plate_model.get_avail_layers()

    if "Topologies" in all_layers:
        topology_features = plate_model.get_layer("Topologies")
    if "StaticPolygons" in all_layers:
        static_polygons = plate_model.get_layer("StaticPolygons")

    return PlateReconstruction(
        plate_model.get_rotation_model(),
        topology_features=topology_features,
        static_polygons=static_polygons,
        plate_model=plate_model,
    )


def get_gplot(
    model: Union[str, PlateModel],
    model_repo_dir: str = "./",
    time: Union[int, float] = 0,
    plot_engine: PlotEngine = CartopyPlotEngine(),
) -> PlotTopologies:
    """Return a :py:class:`gplately.PlotTopologies` object for a given model name or :class:`gplately.PlateModel` object.

    Parameters
    ----------
    model : str or PlateModel
        model name or a :class:`gplately.PlateModel` object
    model_repo_dir: str, default="./"
        the folder in which you would like to keep the model files
    time: int or float, default=0
        the reconstruction age/time
    plot_engine: PlotEngine, default=CartopyPlotEngine()
        two choices - CartopyPlotEngine() or PygmtPlotEngine()

    Returns
    -------
    PlotTopologies
        a :class:`gplately.PlotTopologies` object


    .. seealso::

        `usage example <https://github.com/GPlates/gplately/blob/master/Notebooks/Examples/use_auxiliary_functions.py>`__
    """
    if isinstance(model, str):
        model_name: str = model
        try:
            plate_model = PlateModelManager().get_model(
                model_name, data_dir=model_repo_dir
            )
        except:
            plate_model = PlateModel(model_name, data_dir=model_repo_dir, readonly=True)

        if plate_model is None:
            raise Exception(f"Unable to get model ({model_name})")
    else:
        plate_model = model

    m = get_plate_reconstruction(plate_model)

    coastlines = None
    COBs = None
    continents = None

    all_layers = plate_model.get_avail_layers()

    if "Coastlines" in all_layers:
        coastlines = plate_model.get_layer("Coastlines")
    if "COBs" in all_layers:
        COBs = plate_model.get_layer("COBs")
    if "ContinentalPolygons" in all_layers:
        continents = plate_model.get_layer("ContinentalPolygons")

    return PlotTopologies(
        m,
        coastlines=coastlines,
        COBs=COBs,
        continents=continents,
        time=time,
        plot_engine=plot_engine,
    )


def get_pygmt_basemap_figure(projection="N180/10c", region="d"):
    """A helper function to return a ``pygmt.Figure()`` object

    Parameters
    ----------
    projection: str, default="N180/10c"
        string to define the map projection in GMT style
    region: str, default="d"
        string to define the map extent in GMT style


    Returns
    -------
    pygmt.Figure()
       a ``pygmt.Figure()`` object for map plotting

    """
    fig = pygmt.Figure()
    fig.basemap(region=region, projection=projection, frame="lrtb")
    return fig


def get_data_server_cache_path():
    """Return the path to the :class:`gplately.DataServer` cache as a ``os.PathLike`` object.

    .. seealso::

        :py:attr:`gplately.DataServer.cache_path`
    """
    return path_to_cache()
