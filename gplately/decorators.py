def append_docstring(doctring_to_add):
    def inner(func_pointer):
        def wrapper(*args, **kwargs):
            return func_pointer(*args, **kwargs)

        wrapper.__doc__ = func_pointer.__doc__ + doctring_to_add
        return wrapper

    return inner
