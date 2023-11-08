import argparse
import os
import sys
import warnings
from typing import List, Optional, Sequence, Union

import pygplates

from gplately import __version__, feature_filter

from .ptt import fix_crossovers, remove_plate_rotations


def combine_feature_collections(input_files: List[str], output_file: str):
    """Combine multiple feature collections into one."""
    feature_collection = pygplates.FeatureCollection()
    for file in input_files:
        if not os.path.isfile(file):
            raise Exception(f"{file} is not a file.")
        feature_collection.add(pygplates.FeatureCollection(file))

    feature_collection.write(output_file)

    print(f"Done! The combined feature collection has been saved to {output_file}.")


def _run_combine_feature_collections(args):
    combine_feature_collections(
        [args.combine_first_input_file] + args.combine_other_input_files,
        args.combine_output_file,
    )


def create_agegrids(
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
    from gplately import PlateReconstruction, PlotTopologies, SeafloorGrid

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


class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"error: {message}\n")
        self.print_help()
        sys.exit(1)


def main():
    parser = ArgParser()

    parser.add_argument("-v", "--version", action="store_true")

    # sub-commands
    subparser = parser.add_subparsers(
        dest="command",
        title="subcommands",
        description="valid subcommands",
    )

    # add combine feature sub-command
    combine_cmd = subparser.add_parser(
        "combine",
        help=combine_feature_collections.__doc__,
        description=combine_feature_collections.__doc__,
    )

    # add feature filter sub-command
    feature_filter.add_parser(subparser)

    agegrid_cmd = subparser.add_parser(
        "agegrid",
        aliases=("ag",),
        help=create_agegrids.__doc__,
        add_help=True,
        description=create_agegrids.__doc__,
    )

    # add fix crossovers sub-command
    fix_crossovers_cmd = subparser.add_parser(
        "fix_crossovers",
        help="fix crossovers",
        add_help=True,
    )
    fix_crossovers.add_arguments(fix_crossovers_cmd)

    # add remove plate rotations sub-command
    remove_plate_rotations_cmd = subparser.add_parser(
        "remove_rotations",
        help="remove plate rotations",
        add_help=True,
    )
    remove_plate_rotations.add_arguments(remove_plate_rotations_cmd)

    # combine command arguments
    combine_cmd.set_defaults(func=_run_combine_feature_collections)
    combine_cmd.add_argument("combine_first_input_file", type=str)
    combine_cmd.add_argument("combine_other_input_files", nargs="+", type=str)
    combine_cmd.add_argument("combine_output_file", type=str)

    # agegrid command arguments
    agegrid_cmd.set_defaults(func=_run_create_agegrids)
    agegrid_cmd.add_argument(
        metavar="INPUT_FILE",
        nargs="+",
        help="input reconstruction files",
        dest="input_filenames",
    )
    agegrid_cmd.add_argument(
        metavar="OUTPUT_DIR",
        help="output directory",
        dest="output_dir",
    )
    agegrid_cmd.add_argument(
        "-c",
        "--continents",
        metavar="CONTINENTS_FILE",
        nargs="+",
        help="input continent files",
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

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.version:
        print(__version__)
        sys.exit(0)

    args.func(args)


if __name__ == "__main__":
    main()
