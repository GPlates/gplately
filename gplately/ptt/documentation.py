

import pkg_resources as _pkg_resources
from distutils import dir_util as _dir_util
import os


def install_documentation(path="./PlateTectonicTools-Examples"):
    """
    Install the examples for PlateTectonicTools in the given location.

    WARNING: If the path exists, the files will be written into the path
    and will overwrite any existing files with which they collide. The default
    path ("./PlateTectonicTools-Examples") is chosen to make collision less likely/problematic

    The documentation for PlateTectonicTools is in the form of jupyter notebooks.
    """

    Notebooks_Path = _pkg_resources.resource_filename("ptt", os.path.join("Examples"))

    ct = _dir_util.copy_tree(Notebooks_Path, path, preserve_mode=1, preserve_times=1, preserve_symlinks=1, update=0, verbose=1, dry_run=0)

    return
