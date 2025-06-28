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

logger = logging.getLogger("gplately")


def install_and_update_pmm():
    import subprocess
    import sys

    subprocess.call(
        [sys.executable, "-m", "pip", "install", "plate-model-manager", "--upgrade"]
    )


def is_pmm_version_good_enough(installed_version, required_version):
    """return True if the version of installed pmm is good enough, otherwise False.
    assume the version string something like 1.2.0
    """
    installed_version_numbers = installed_version.split(".")
    required_version_numbers = required_version.split(".")
    if int(installed_version_numbers[0]) > int(required_version_numbers[0]):
        return True
    elif int(installed_version_numbers[0]) == int(required_version_numbers[0]):
        if int(installed_version_numbers[1]) > int(required_version_numbers[1]):
            return True
        elif int(installed_version_numbers[1]) == int(required_version_numbers[1]):
            if int(installed_version_numbers[2]) >= int(required_version_numbers[2]):
                return True
    return False


def ensure_plate_model_manager_compatible(REQUIRED_PMM_VERSION: str):
    try:
        import plate_model_manager
    except (ImportError, ModuleNotFoundError):
        logger.info("The plate_model_manager is not installed, installing it now!")
        install_and_update_pmm()
        import plate_model_manager

    if not is_pmm_version_good_enough(
        plate_model_manager.__version__, REQUIRED_PMM_VERSION
    ):
        logger.info("The plate_model_manager is outdated, updating it now!")
        install_and_update_pmm()
        import importlib

        importlib.reload(plate_model_manager)
