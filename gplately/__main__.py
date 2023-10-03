import argparse
import os
import sys

import pygplates

from gplately import __version__, feature_filter


def combine_feature_collections(input_files: list[str], output_file: str):
    """combine multiply feature collections into one"""
    feature_collection = pygplates.FeatureCollection()
    for file in input_files:
        if not os.path.isfile(file):
            raise Exception(f"{file} is not a file.")
        feature_collection.add(pygplates.FeatureCollection(file))

    feature_collection.write(output_file)

    print(f"Done! The result has been saved to {output_file}.")


def filter_feature_collection(args):
    input_feature_collection = pygplates.FeatureCollection(
        args.filter_input_file
        # "Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"
    )

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
        "-a", "--min-birth-age", type=int, dest="min_birth_age"
    )
    birth_age_group.add_argument("--max-birth-age", type=int, dest="max_birth_age")

    filter_cmd.add_argument(
        "--case-sensitive", dest="case_sensitive", action="store_true"
    )
    filter_cmd.add_argument("--exact-match", dest="exact_match", action="store_true")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

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
    else:
        print(f"Unknow command {args.command}!")
        parser.print_help(sys.stderr)
        sys.exit(1)

    # print(args)


if __name__ == "__main__":
    main()
