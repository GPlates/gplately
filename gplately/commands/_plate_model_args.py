#
#    Copyright (C) 2026 The University of Sydney, Australia
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

"""Let a CLI command get its rotation and topology files from a plate model name (instead of '-r' and '-m')."""

DEFAULT_MODEL_REPO = "plate-model-repo"


def add_model_arguments(parser, run_func):
    """Add '--model' and '--model-repo' to a command that has '-r' (rotation files) and '-m' (topology files).

    The '-r' and '-m' arguments should not be required (by argparse), since either they or '--model' must be specified.
    When the command runs, the rotation and topology files are resolved (from '--model' if specified) into
    ``args.rotation_filenames`` and ``args.topology_filenames`` before calling ``run_func(args)``.
    """
    parser.add_argument(
        "--model",
        metavar="MODEL_NAME",
        dest="model_name",
        default=None,
        help="Name of a plate model (see 'gplately list') to use the rotation and topology files of, "
        "instead of specifying them with '-r' and '-m'. The model is downloaded if needed.",
    )
    parser.add_argument(
        "--model-repo",
        metavar="DIR",
        dest="model_repo",
        default=None,
        help="Folder to download plate models to (only used with '--model'). "
        f"Defaults to '{DEFAULT_MODEL_REPO}' in the current folder.",
    )

    def _run(args):
        if args.model_name:
            if args.rotation_filenames or args.topology_filenames:
                parser.error(
                    "argument --model: not allowed with -r/--rotation-filenames or -m/--topology-filenames"
                )
            try:
                args.rotation_filenames, args.topology_filenames = (
                    _get_rotation_and_topology_files(
                        args.model_name, args.model_repo or DEFAULT_MODEL_REPO
                    )
                )
            except ValueError as error:
                parser.error(f"argument --model: {error}")
        else:
            if args.model_repo:
                parser.error("argument --model-repo: only allowed with --model")
            if not (args.rotation_filenames and args.topology_filenames):
                message = (
                    "the following arguments are required: "
                    "-r/--rotation-filenames and -m/--topology-filenames (or --model)"
                )
                if args.topology_filenames and not args.rotation_filenames:
                    # Eg, 'gplately sc -m Muller2019' (in 'seafloor-grids', '-m' is a model name).
                    message += ". Note that -m is for topology files; use --model for a plate model name"
                parser.error(message)
        return run_func(args)

    parser.set_defaults(func=_run)


def _get_rotation_and_topology_files(model_name, model_repo):
    """Return the rotation and topology files of a plate model (raises ValueError if unavailable)."""
    # Import here so that 'gplately -h' (and other commands) don't pay for it.
    from ..auxiliary import get_plate_model

    try:
        plate_model = get_plate_model(model_name, model_repo)
    except Exception as error:
        # An unknown model name raises RuntimeError. If the model can't be downloaded (eg, offline)
        # then reading it from 'model_repo' instead can also fail (with a plain Exception).
        raise ValueError(
            f"unable to get plate model '{model_name}'. Check the model name "
            "(run 'gplately list' to see the available models) and your network connection"
        ) from error

    rotation_files = plate_model.get_rotation_model()
    topology_files = plate_model.get_layer("Topologies", return_none_if_not_exist=True)
    if not rotation_files or not topology_files:
        missing = "rotations" if not rotation_files else "topologies"
        raise ValueError(
            f"plate model '{model_name}' has no {missing} "
            f"(or, if there's no network connection, they haven't been downloaded to '{model_repo}' yet). "
            f"Run 'gplately list -m {model_name}' to see the model's layers"
        )

    return rotation_files, topology_files
