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

import logging

from packaging.version import InvalidVersion, parse

from importlib.metadata import requires
from packaging.requirements import Requirement

logger = logging.getLogger("gplately")


def get_required_pmm_version(default="1.3.2"):
    reqs = requires("gplately") or []
    for req in reqs:
        r = Requirement(req)
        if r.name == "plate-model-manager":
            # pick the >= lower bound if present
            for spec in r.specifier:
                if spec.operator in (">=", "=="):
                    # logger.info("Found required plate-model-manager version: %s", spec.version)
                    return spec.version
    # logger.info("Using default plate-model-manager version: %s", default)
    return default


def install_and_update_pmm():
    import subprocess
    import sys

    subprocess.call(
        [sys.executable, "-m", "pip", "install", "plate-model-manager", "--upgrade"]
    )


def is_pmm_version_good_enough(installed_version, required_version):
    """return True if the version of installed pmm is good enough, otherwise False."""
    try:
        return parse(installed_version) >= parse(required_version)
    except InvalidVersion:
        logger.error(
            "Invalid version string while checking plate-model-manager compatibility: "
            "installed=%s required=%s",
            installed_version,
            required_version,
        )
        return False


def ensure_plate_model_manager_compatible(REQUIRED_PMM_VERSION: str):
    try:
        import plate_model_manager

        if not is_pmm_version_good_enough(
            plate_model_manager.__version__, REQUIRED_PMM_VERSION
        ):
            logger.fatal(
                "The installed version of plate-model-manager is not compatible with gplately. Please update it to version %s or later. https://pypi.org/project/plate-model-manager/",
                REQUIRED_PMM_VERSION,
            )
    except (ImportError, ModuleNotFoundError):
        logger.fatal(
            "The plate-model-manager package is required but not installed. Please install it https://pypi.org/project/plate-model-manager/."
        )
        # install_and_update_pmm()
        # import plate_model_manager
        # import importlib
        # importlib.reload(plate_model_manager)
