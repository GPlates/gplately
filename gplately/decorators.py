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
from functools import wraps

from .exceptions import ReconstructionTimeNotSet

logger = logging.getLogger("gplately")


def append_docstring(docstring_to_add):
    """append text to the end of the function's __doc__

    Parameters
    ----------
    docstring_to_add : str
        the text to append to the function's __doc__
    """

    def inner(func_pointer):
        if func_pointer.__doc__:
            func_pointer.__doc__ += docstring_to_add
        else:
            func_pointer.__doc__ = docstring_to_add

        @wraps(func_pointer)
        def wrapper(*args, **kwargs):
            return func_pointer(*args, **kwargs)

        return wrapper

    return inner


def validate_topology_availability(feature_name):
    """check if topology features exist before plotting feature. if not, do nothing and return

    Parameters
    ----------
    feature_name : str
        this parameter is used in logging warning message, indicating in which function the error happened.

    """

    def inner(func_pointer):
        @wraps(func_pointer)
        def wrapper(self, ax, **kwargs):
            if not self.plate_reconstruction.topology_features:
                logger.warning(
                    f"Plate model does not have topology features. Unable to plot {feature_name}."
                )
                return ax
            return func_pointer(self, ax, **kwargs)

        return wrapper

    return inner


def validate_reconstruction_time(func_pointer):
    """check if reconstruction time is None. If so, raise ReconstructionTimeNotSet exception"""

    @wraps(func_pointer)
    def wrapper(self, *args, **kwargs):
        validate_reconstruction_time = kwargs.pop("validate_reconstruction_time", True)
        if validate_reconstruction_time and self._time is None:
            raise ReconstructionTimeNotSet()
        return func_pointer(self, *args, **kwargs)

    return wrapper
