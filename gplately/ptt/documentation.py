#
#    Copyright (C) 2024 The University of Sydney, Australia
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

import os
from distutils import dir_util as _dir_util

import pkg_resources as _pkg_resources


def install_documentation(path="./PlateTectonicTools-Examples"):
    """
    Install the examples for PlateTectonicTools in the given location.

    WARNING: If the path exists, the files will be written into the path
    and will overwrite any existing files with which they collide. The default
    path ("./PlateTectonicTools-Examples") is chosen to make collision less likely/problematic

    The documentation for PlateTectonicTools is in the form of jupyter notebooks.
    """

    Notebooks_Path = _pkg_resources.resource_filename("ptt", os.path.join("Examples"))

    ct = _dir_util.copy_tree(
        Notebooks_Path,
        path,
        preserve_mode=1,
        preserve_times=1,
        preserve_symlinks=1,
        update=0,
        verbose=1,
        dry_run=0,
    )

    return
