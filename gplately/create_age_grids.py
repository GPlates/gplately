import argparse
import warnings
from typing import Optional, Sequence, Union

import pygplates
from plate_model_manager import PlateModelManager

from gplately import PlateReconstruction, PlotTopologies, SeafloorGrid


def add_parser(parser: argparse.ArgumentParser):
    """add command line argument parser"""

    agegrid_cmd = parser.add_parser(
        "agegrid",
        aliases=("ag",),
        help=__description__,
        add_help=True,
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # agegrid command arguments
    agegrid_cmd.set_defaults(func=_run_create_agegrids)
    agegrid_cmd.add_argument(
        metavar="INPUT_FILE",
        nargs="*",
        help="input reconstruction files, including rotation files and topology files",
        dest="input_filenames",
    )
    agegrid_cmd.add_argument(
        metavar="OUTPUT_DIR",
        help="output directory",
        dest="output_dir",
    )
    agegrid_cmd.add_argument(
        "-m",
        "--model",
        metavar="MODEL_NAME",
        help="reconstruction model name",
        dest="model_name",
        default=None,
    )
    agegrid_cmd.add_argument(
        "-c",
        "--continents",
        metavar="CONTINENTS_FILE",
        nargs="*",
        help="input continental polygons files",
        dest="continents_filenames",
        default=None,
    )
    agegrid_cmd.add_argument(
        "-r",
        "--resolution",
        metavar="RESOLUTION",
        type=float,
        help="grid resolution (degrees); default: 0.1",
        default=0.1,
        dest="grid_spacing",
    )
    agegrid_cmd.add_argument(
        "--refinement-levels",
        metavar="LEVELS",
        type=int,
        help="mesh refinement levels; default: 5",
        default=5,
        dest="refinement_levels",
    )
    agegrid_cmd.add_argument(
        "--ridge-sampling",
        metavar="RESOLUTION",
        type=float,
        help="MOR sampling resolution (degrees); default: 0.5",
        default=0.5,
        dest="ridge_sampling",
    )
    agegrid_cmd.add_argument(
        "--initial-spreadrate",
        metavar="SPREADRATE",
        type=float,
        help="initial ocean spreading rate (km/Myr); default: 75",
        default=75,
        dest="initial_spreadrate",
    )
    agegrid_cmd.add_argument(
        "-e",
        "--min-time",
        metavar="MIN_TIME",
        type=float,
        help="minimum time (Ma); default: 0",
        default=0,
        dest="min_time",
    )
    agegrid_cmd.add_argument(
        "-s",
        "--max-time",
        metavar="MAX_TIME",
        type=float,
        help="maximum time (Ma); default: 0",
        default=0,
        dest="max_time",
    )
    agegrid_cmd.add_argument(
        "-j",
        "--n_jobs",
        help="number of processes to use; default: 1",
        metavar="N_JOBS",
        default=1,
        dest="n_jobs",
    )
    agegrid_cmd.add_argument(
        "-f",
        "--file-collection",
        help="file collection name (optional)",
        metavar="NAME",
        default=None,
        dest="file_collection",
    )
    agegrid_cmd.add_argument(
        "-u",
        "--include-unmasked",
        help="create unmasked grids in addition to masked ones",
        action="store_true",
        dest="unmasked",
    )


__description__ = """Create age grids for a plate model.

example usage: gplately ag plate-model-repo/muller2019/Rotations/Muller_etal_2019_CombinedRotations.rot plate-model-repo/muller2019/Topologies/Muller_etal_2019_PlateBoundaries_DeformingNetworks.gpmlz output -c plate-model-repo/muller2019/ContinentalPolygons/Global_EarthByte_GPlates_PresentDay_ContinentalPolygons_2019_v1.shp -e 0 -s 10
               gplately ag output -m muller2019 -e 0 -s 10

"""


def create_agegrids(
    model_name: str,
    input_filenames: Union[str, Sequence[str]],
    continents_filenames: Union[str, Sequence[str]],
    output_dir: str,
    min_time: float,
    max_time: float,
    ridge_time_step: float = 1,
    n_jobs: int = 1,
    refinement_levels: int = 5,
    grid_spacing: float = 0.1,
    ridge_sampling: float = 0.5,
    initial_spreadrate: float = 75,
    file_collection: Optional[str] = None,
    unmasked: bool = False,
) -> None:
    """Create age grids for a plate model."""

    if model_name:
        pm_manager = PlateModelManager()
        plate_model = pm_manager.get_model(model_name, data_dir="plate-model-repo")

        rotations = pygplates.FeatureCollection(
            pygplates.FeaturesFunctionArgument(
                plate_model.get_rotation_model()
            ).get_features()
        )
        topologies = pygplates.FeatureCollection(
            pygplates.FeaturesFunctionArgument(
                plate_model.get_topologies()
            ).get_features()
        )
        continents = pygplates.FeatureCollection(
            pygplates.FeaturesFunctionArgument(
                plate_model.get_layer("ContinentalPolygons")
            ).get_features()
        )
    else:
        features = pygplates.FeaturesFunctionArgument(input_filenames).get_features()
        rotations = []
        topologies = []
        for i in features:
            if (
                i.get_feature_type().to_qualified_string()
                == "gpml:TotalReconstructionSequence"
            ):
                rotations.append(i)
            else:
                topologies.append(i)
        topologies = pygplates.FeatureCollection(topologies)
        rotations = pygplates.RotationModel(rotations)

        if continents_filenames is None:
            continents = pygplates.FeatureCollection()
        else:
            continents = pygplates.FeatureCollection(
                pygplates.FeaturesFunctionArgument(continents_filenames).get_features()
            )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ImportWarning)
        reconstruction = PlateReconstruction(
            rotation_model=rotations,
            topology_features=topologies,
        )
        gplot = PlotTopologies(
            reconstruction,
            continents=continents,
        )

        grid = SeafloorGrid(
            reconstruction,
            gplot,
            min_time=min_time,
            max_time=max_time,
            save_directory=output_dir,
            ridge_time_step=ridge_time_step,
            refinement_levels=refinement_levels,
            grid_spacing=grid_spacing,
            ridge_sampling=ridge_sampling,
            initial_ocean_mean_spreading_rate=initial_spreadrate,
            file_collection=file_collection,
        )
        grid.reconstruct_by_topologies()
        for val in ("SEAFLOOR_AGE", "SPREADING_RATE"):
            grid.lat_lon_z_to_netCDF(val, unmasked=unmasked, nprocs=n_jobs)


def _run_create_agegrids(args):
    create_agegrids(
        model_name=args.model_name,
        input_filenames=args.input_filenames,
        continents_filenames=args.continents_filenames,
        output_dir=args.output_dir,
        min_time=args.min_time,
        max_time=args.max_time,
        n_jobs=args.n_jobs,
        refinement_levels=args.refinement_levels,
        grid_spacing=args.grid_spacing,
        ridge_sampling=args.ridge_sampling,
        initial_spreadrate=args.initial_spreadrate,
        file_collection=args.file_collection,
        unmasked=args.unmasked,
    )
