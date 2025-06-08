import logging

logger = logging.getLogger("gplately-pytest-logger")

## ==========================


def test_numpy_import():
    import numpy

    return


def test_scipy_import():
    import scipy

    logger.info("\t\t You have scipy version {}".format(scipy.__version__))


def test_cartopy_import():
    import cartopy


def test_pygplates_import():
    import pygplates


def test_pooch_import():
    import pooch


def test_gplately_modules():
    import gplately
    from gplately import download, grids, plot, ptt, pygplates, tools

    logger.info("\t\t You have gplately version {}".format(gplately.__version__))
    logger.info(
        "\t\t GPlately is using pygplates version {}".format(pygplates.__version__)
    )


def test_jupyter_available():
    from subprocess import check_output

    try:
        result = str(check_output(["which", "jupyter"]))[2:-3]
    except:
        logger.info("Jupyter not installed")
        logger.info("Jupyter is needed to run the example documentation")


def test_plate_model_manager_import():
    import plate_model_manager

    logger.info(plate_model_manager.__version__)
