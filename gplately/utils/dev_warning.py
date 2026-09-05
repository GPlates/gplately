#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
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

import logging, os, warnings
from . import settings

logger = logging.getLogger("gplately")

disable_dev_warning = settings.disable_dev_warning


def print_dev_warning(version: str):
    if not disable_dev_warning:
        _warn_if_dev_version(version)


class GPlatelyDevWarning(UserWarning):
    """Raised when importing/using a DEV build of GPlately."""


def _warn_if_dev_version(version: str) -> None:
    if os.environ.get("GPLATELY_DISABLE_DEV_WARNING", "").lower() == "true":
        return

    message = f"""\
    

        You are using a development version of GPlately ({version}).

        Some functionality in DEV builds is not yet fully tested and may break
        your code or produce incorrect results. Proceed with caution.

        You likely also need the DEV version of plate_model_manager to work
        with this GPlately build:
            https://github.com/gplates/plate-model-manager

        To silence this warning, set the GPLATELY_DISABLE_DEV_WARNING
        environment variable to "true" before importing GPlately:

            Python:      os.environ["GPLATELY_DISABLE_DEV_WARNING"] = "true"
            bash/zsh:    export GPLATELY_DISABLE_DEV_WARNING=true
            PowerShell:  $env:GPLATELY_DISABLE_DEV_WARNING = "true"

        To make this permanent, add the shell command above to your shell
        startup file (e.g. .bashrc, .zshrc, or your PowerShell profile).

        Alternatively, you can filter this warning in your code with:
            import warnings
            warnings.filterwarnings("ignore", category=UserWarning, module="gplately.utils.dev_warning")
        """
    warnings.warn(message, category=GPlatelyDevWarning, stacklevel=2)


def print_using_source_code_warning(version: str):
    if not disable_dev_warning:
        if os.path.isdir(
            f"{os.path.dirname(os.path.realpath(__file__))}/../.git"
        ) or not (
            os.path.isfile(
                f"{os.path.dirname(os.path.realpath(__file__))}/../../../../bin/gplately"
            )
            or os.path.isfile(
                f"{os.path.dirname(os.path.realpath(__file__))}/../../../../Scripts/gplately.exe"
            )
        ):
            logger.warning(
                f"It seems that you are using GPlately source code directly or installed editable package with `pip install -e .`, "
                + f"the version number({version}) may not be accurate in these cases."
            )

        logger.info(
            f"The location of GPlately currently in use is {os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}. "
        )
