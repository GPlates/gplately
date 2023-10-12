import argparse
import os
import sys
from typing import (
    Iterable,
    List,
    Optional,
    Union,
)

import pygplates

from gplately import __version__, feature_filter


def combine_feature_collections(input_files: List[str], output_file: str):
    """combine multiply feature collections into one"""
    feature_collection = pygplates.FeatureCollection()
    for file in input_files:
        if not os.path.isfile(file):
            raise Exception(f"{file} is not a file.")
        feature_collection.add(pygplates.FeatureCollection(file))

    feature_collection.write(output_file)

    print(f"Done! The combined feature collection has been saved to {output_file}.")


def filter_feature_collection(args):
    """filter the input feature collection according to command line arguments"""
    input_feature_collection = pygplates.FeatureCollection(args.filter_input_file)

    filters = []
    if args.names:
        filters.append(
            feature_filter.FeatureNameFilter(
                args.names,
                exact_match=args.exact_match,
                case_sensitive=args.case_sensitive,
            )
        )
    elif args.exclude_names:
        filters.append(
            feature_filter.FeatureNameFilter(
                args.exclude_names,
                exclude=True,
                exact_match=args.exact_match,
                case_sensitive=args.case_sensitive,
            )
        )

    if args.pids:
        filters.append(feature_filter.PlateIDFilter(args.pids))
    elif args.exclude_pids:
        filters.append(feature_filter.PlateIDFilter(args.exclude_pids, exclude=True))

    # print(args.max_birth_age)
    if args.max_birth_age is not None:
        filters.append(
            feature_filter.BirthAgeFilter(args.max_birth_age, keep_older=False)
        )
    elif args.min_birth_age is not None:
        filters.append(feature_filter.BirthAgeFilter(args.min_birth_age))

    new_fc = feature_filter.filter_feature_collection(
        input_feature_collection,
        filters,
    )

    new_fc.write(args.filter_output_file)
    print(
        f"Done! The filtered feature collection has been saved to {args.filter_output_file}."
    )


def create_agegrids(
    input_filenames: Union[str, Iterable[str]],
    continents_filenames: Union[str, Iterable[str]],
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
    from gplately import (
        PlateReconstruction,
        PlotTopologies,
        SeafloorGrid,
    )

    features = pygplates.FeaturesFunctionArgument(
        input_filenames
    ).get_features()
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
            pygplates.FeaturesFunctionArgument(
                continents_filenames
            ).get_features()
        )

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


class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"error: {message}\n")
        self.print_help()
        sys.exit(1)


def main():
    parser = ArgParser()

    parser.add_argument("-v", "--version", action="store_true")

    # sub-commands
    subparser = parser.add_subparsers(dest="command")
    combine_cmd = subparser.add_parser("combine")
    filter_cmd = subparser.add_parser("filter")
    agegrid_cmd = subparser.add_parser(
        "agegrid",
        aliases=("ag",),
        add_help=True,
        description="Create age grids for a plate model.",
    )

    # combine command arguments
    combine_cmd.add_argument("combine_first_input_file", type=str)
    combine_cmd.add_argument("combine_other_input_files", nargs="+", type=str)
    combine_cmd.add_argument("combine_output_file", type=str)

    # feature filter command arguments
    filter_cmd.add_argument("filter_input_file", type=str)
    filter_cmd.add_argument("filter_output_file", type=str)

    name_group = filter_cmd.add_mutually_exclusive_group()
    name_group.add_argument("-n", "--names", type=str, dest="names", nargs="+")
    name_group.add_argument(
        "--exclude-names", type=str, dest="exclude_names", nargs="+"
    )

    pid_group = filter_cmd.add_mutually_exclusive_group()
    pid_group.add_argument("-p", "--pids", type=int, dest="pids", nargs="+")
    pid_group.add_argument("--exclude-pids", type=int, dest="exclude_pids", nargs="+")

    birth_age_group = filter_cmd.add_mutually_exclusive_group()
    birth_age_group.add_argument(
        "-a", "--min-birth-age", type=float, dest="min_birth_age"
    )
    birth_age_group.add_argument("--max-birth-age", type=float, dest="max_birth_age")

    filter_cmd.add_argument(
        "--case-sensitive", dest="case_sensitive", action="store_true"
    )
    filter_cmd.add_argument("--exact-match", dest="exact_match", action="store_true")

    # agegrid command arguments
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
        help="minimum time (Ma)",
        default=0,
        dest="min_time",
    )
    agegrid_cmd.add_argument(
        "-s",
        "--max-time",
        metavar="MAX_TIME",
        type=float,
        help="maximum time (Ma)",
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
        help="file collection name; default: None",
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

    if args.command == "combine":
        combine_feature_collections(
            [args.combine_first_input_file] + args.combine_other_input_files,
            args.combine_output_file,
        )
    elif args.command == "filter":
        filter_feature_collection(args)
    elif args.command == "agegrid":
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
    else:
        print(f"Unknown command {args.command}!")
        parser.print_help(sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
