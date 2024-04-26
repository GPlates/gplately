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
        func_pointer.__doc__ += docstring_to_add

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
                logger.warn(
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
