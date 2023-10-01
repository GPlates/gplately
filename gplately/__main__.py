import argparse
import os
import sys

import pygplates

from gplately import __version__


def combine_feature_collections(input_files: list[str], output_file: str):
    """combine multiply feature collections into one"""
    feature_collection = pygplates.FeatureCollection()
    for file in input_files:
        if not os.path.isfile(file):
            raise Exception(f"{file} is not a file.")
        feature_collection.add(pygplates.FeatureCollection(file))

    feature_collection.write(output_file)

    print(f"Done! The result has been saved to {output_file}.")


class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"error: {message}\n")
        self.print_help()
        sys.exit(1)


def main():
    parser = ArgParser()

    parser.add_argument("-v", "--version", action="store_true")

    subparser = parser.add_subparsers(dest="command")
    combine_cmd = subparser.add_parser("combine")
    filter_cmd = subparser.add_parser("filter")

    combine_cmd.add_argument("combine_first_input_file", type=str)
    combine_cmd.add_argument("combine_other_input_files", nargs="+", type=str)
    combine_cmd.add_argument("combine_output_file", type=str)

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
        print("Coming soon...")
    else:
        print(f"Unknow command {args.command}!")
        parser.print_help(sys.stderr)
        sys.exit(1)

    # print(args)


if __name__ == "__main__":
    main()
