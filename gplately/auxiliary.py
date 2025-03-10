from typing import Union

from plate_model_manager import PlateModel, PlateModelManager

from .mapping.cartopy_plot import CartopyPlotEngine
from .mapping.plot_engine import PlotEngine
from .plot import PlotTopologies
from .reconstruction import PlateReconstruction


def get_plate_reconstruction(model: Union[str, PlateModel], model_repo_dir: str = "./"):
    """Convenient function to return a PlateReconstruction object. Check out the [usage example](https://gplates.github.io/gplately/dev-doc/#platemodelmanager).

    Parameters
    ----------
    model : str or PlateModel
        model name or a PlateModel object
    model_repo_dir: str, default="./"
        the folder in which you would like to keep the model files

    Returns
    -------
    PlateReconstruction
        a PlateReconstruction object

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
    )


def get_gplot(
    model: Union[str, PlateModel],
    model_repo_dir: str = "./",
    time: Union[int, float] = 0,
    plot_engine: PlotEngine = CartopyPlotEngine(),
) -> PlotTopologies:
    """Convenient function to return a PlotTopologies object. Check out the [usage example](https://gplates.github.io/gplately/dev-doc/#platemodelmanager).

    Parameters
    ----------
    model : str or PlateModel
        model name or a PlateModel object
    model_repo_dir: str, default="./"
        the folder in which you would like to keep the model files
    time: int or float, default=0
        the reconstruction age/time
    plot_engine: PlotEngine, default=CartopyPlotEngine()
        two choices - CartopyPlotEngine() or PygmtPlotEngine()

    Returns
    -------
    PlotTopologies
        a PlotTopologies object
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
