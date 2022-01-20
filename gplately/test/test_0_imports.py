import pytest

## ==========================

def test_numpy_import():
    import numpy
    return

def test_scipy_import():
    import scipy
    print("\t\t You have scipy version {}".format(scipy.__version__))


def test_cartopy_import():
    import cartopy


def test_pygplates_import():
    import pygplates


def test_pooch_import():
    import pooch


def test_gplately_modules():
    import gplately
    from gplately import plot
    from gplately import download
    from gplately import tools
    from gplately import grids


def test_jupyter_available():
    from subprocess import check_output
    try:
        result = str(check_output(['which', 'jupyter']))[2:-3]
    except:
        print("Jupyter not installed")
        print("Jupyter is needed to run the example documentation")
