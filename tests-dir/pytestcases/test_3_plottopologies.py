import pickle

import matplotlib.pyplot as plt
import pytest
from conftest import gplately_plot_topologies_object as gplot
from conftest import logger, reconstruction_times

import gplately

# ========================================= <gplately.PlotTopologies> =========================================

""" 
A series of automated tests that ensure GPlately's <DataServer> object collects the necessary
continent, coastline and continent-ocean boundary shapefiles/GPML files to initialise the 
<PlotTopologies> object with the Müller et al. (2019) plate reconstruction model. The following 
methods in the object are tested:

    - __init__ 
        This includes testing update_time; this ensures all geometries in the provided coastline,
        COB and/or continent files can be reconstructed to 0 and 100 Ma. It also ensures that the 
        outputs of Plate Tectonic Tools' resolve_topologies method:

            - resolved topology features (topological plates and networks)
            - ridge and transform boundary sections (resolved features)
            - ridge boundary sections (resolved features)
            - transform boundary sections (resolved features)
            - subduction boundary sections (resolved features)
            - left subduction boundary sections (resolved features)
            - right subduction boundary sections (resolved features)
            - other boundary sections (resolved features) that are not subduction zones or mid-ocean ridges (ridge/transform)

        are initialised in the <PlotTopologies> object for a given plate reconstruction model.

    - Plotting functions: (not implemented yet)
        This ensures that <cartopy.mpl.GeoAxis> objects are created with the attributes passed to
        <PlotTopologies>. Namely, the following plotting methods are tested:

            - plot_coastlines
            - plot_continents
            - plot_continent_ocean_boundaries
            - plot_ridges
            - plot_ridges_and_transforms
            - plot_transforms
            - plot_trenches
            - plot_misc_boundaries
            - plot_plate_id
            - plot_subduction_teeth
            - plot_grid
            - plot_grid_from_netCDF
            - plot_plate_motion_vectors
"""


def test_gplately_plotTopologies_object(
    gplately_plate_reconstruction_object,
    gplately_muller_static_geometries,
):
    gplot = gplately.PlotTopologies(
        gplately_plate_reconstruction_object,
        *gplately_muller_static_geometries,
    )
    assert gplot.time >= 0

    gplot = gplately.PlotTopologies(
        gplately_plate_reconstruction_object,
        *gplately_muller_static_geometries,
        time=100,
    )
    assert gplot.time == 100

    try:
        gplot = gplately.PlotTopologies(
            gplately_plate_reconstruction_object,
            *gplately_muller_static_geometries,
            time=None,  # type: ignore
        )
        assert False
    except ValueError:
        assert True

    try:
        gplot = gplately.PlotTopologies(
            gplately_plate_reconstruction_object,
            *gplately_muller_static_geometries,
            time=-1,
        )
        assert False
    except ValueError:
        assert True

    try:
        gplot = gplately.PlotTopologies(
            gplately_plate_reconstruction_object,
            *gplately_muller_static_geometries,
            time="abc",  # type: ignore
        )
        assert False
    except ValueError:
        assert True


@pytest.mark.parametrize("time", reconstruction_times)
def test_PlotTopologies_features(time, gplot):
    logger.info(f"test_PlotTopologies_features with time ({time}).")

    gplot.time = time

    assert gplot.ridges, f"No ridge features at {time} Ma in gplately.PlotTopologies."

    for ridge in gplot.ridges:
        assert (
            ridge.get_feature_type().to_qualified_string() == "gpml:MidOceanRidge"
        ), "ridge is not gpml:MidOceanRidge."

    assert (
        gplot.topologies
    ), f"No topological features at {time} Ma in gplately.PlotTopologies"

    assert (
        gplot.transforms
    ), f"No transforms features at {time} Ma in gplately.PlotTopologies."

    for transform in gplot.transforms:
        assert (
            transform.get_feature_type().to_qualified_string() == "gpml:Transform"
        ), "transform is not gpml:Transform."

    assert (
        gplot.trenches
    ), f"No trenches features at {time} Ma in gplately.PlotTopologies."

    for trench in gplot.trenches:
        assert (
            trench.get_feature_type().to_qualified_string() == "gpml:SubductionZone"
        ), "trench is not gpml:SubductionZone."

    assert (
        gplot.trench_left
    ), f"No trench (L) features at {time} Ma in gplately.PlotTopologies."

    for trench in gplot.trench_left:
        assert (
            trench.get_feature_type().to_qualified_string() == "gpml:SubductionZone"
        ), "trench (L) is not gpml:SubductionZone."

    assert (
        gplot.trench_right
    ), f"No trenches features at {time} Ma in gplately.PlotTopologies."

    for trench in gplot.trench_right:
        assert (
            trench.get_feature_type().to_qualified_string() == "gpml:SubductionZone"
        ), "trenchn (R) is not gpml:SubductionZone."


# Subduction teeth
def test_PlotTopologies_subduction_teeth(gplot):
    ax = plt.gca()
    gplot.plot_subduction_teeth(ax=ax)
    fig = plt.gcf()
    plt.close(fig)


def test_pickle_plotTopologies_object(gplot):
    gplot_dump = pickle.dumps(gplot)
    gplot_load = pickle.loads(gplot_dump)
    assert gplot_load.coastlines and len(gplot_load.coastlines) == len(gplot.coastlines)
    assert gplot_load.continents and len(gplot_load.continents) == len(gplot.continents)
    assert gplot_load.COBs and len(gplot_load.COBs) == len(gplot.COBs)


def test_set_invalid_time(gplot):
    try:
        gplot.time = None
    except ValueError as ex:
        logger.info(ex)
        assert True
    try:
        gplot.time = -1
    except ValueError as ex:
        logger.info(ex)
        assert True
    try:
        gplot.time = "abc"
    except ValueError as ex:
        logger.info(ex)
        assert True
