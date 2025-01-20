from typing import Union

from plate_model_manager import PlateModel, PlateModelManager

from .mapping.plot_engine import PlotEngine
from .mapping.cartopy_plot import CartopyPlotEngine
from .plot import PlotTopologies
from .reconstruction import PlateReconstruction


def get_gplot(
    model_name: str,
    model_repo_dir: str,
    age: Union[int, float],
    plot_engine: PlotEngine = CartopyPlotEngine(),
) -> PlotTopologies:
    """auxiliary function to get gplot object"""
    try:
        model = PlateModelManager().get_model(model_name, data_dir=model_repo_dir)
    except:
        model = PlateModel(model_name, data_dir=model_repo_dir, readonly=True)

    if model is None:
        raise Exception(f"Unable to get model ({model_name})")

    topology_features = None
    static_polygons = None
    coastlines = None
    COBs = None
    continents = None

    all_layers = model.get_avail_layers()

    if "Topologies" in all_layers:
        topology_features = model.get_layer("Topologies")
    if "StaticPolygons" in all_layers:
        static_polygons = model.get_layer("StaticPolygons")
    if "Coastlines" in all_layers:
        coastlines = model.get_layer("Coastlines")
    if "COBs" in all_layers:
        COBs = model.get_layer("COBs")
    if "ContinentalPolygons" in all_layers:
        continents = model.get_layer("ContinentalPolygons")

    m = PlateReconstruction(
        model.get_rotation_model(),
        topology_features=topology_features,
        static_polygons=static_polygons,
    )
    return PlotTopologies(
        m,
        coastlines=coastlines,
        COBs=COBs,
        continents=continents,
        time=age,
        plot_engine=plot_engine,
    )
