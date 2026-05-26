import matplotlib.pyplot as plt  # pyright: ignore[reportMissingModuleSource]
import pytest  # pyright: ignore[reportMissingImports]
from conftest import gplately_plot_topologies_object as gplot
from conftest import logger, muller_2019_model, reconstruction_times
from pygplates import FeatureCollection  # pyright: ignore[reportMissingModuleSource]

import gplately
from gplately.gpml import load_feature_collection_from_files  # type: ignore

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


def test_PlotTopologies_subduction_teeth(gplot):
    """Test that the plot_subduction_teeth method adds artists to the target axes when plotting subduction teeth."""
    ccrs = pytest.importorskip("cartopy.crs")

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    assert gplot.trench_left, "No left-trench features available for plotting."
    assert gplot.trench_right, "No right-trench features available for plotting."

    collections_before = len(ax.collections)
    gplot.plot_subduction_teeth(ax=ax, color="black")
    collections_after = len(ax.collections)

    assert (
        collections_after > collections_before
    ), "plot_subduction_teeth did not add any artists to the target axes."

    plt.close(fig)


def test_pickle_plotTopologies_object(gplot, muller_2019_model):
    import pickle

    # Also create a plot object from actual pygplates objects (instead of filenames).
    pygplates_coastline_features = load_feature_collection_from_files(
        muller_2019_model.get_layer("Coastlines")
    )
    pygplates_continent_features = load_feature_collection_from_files(
        muller_2019_model.get_layer("ContinentalPolygons")
    )
    pygplates_COB_features = load_feature_collection_from_files(
        muller_2019_model.get_layer("COBs")
    )

    pygplates_gplot = gplately.PlotTopologies(
        gplot.plate_reconstruction,
        coastlines=pygplates_coastline_features,
        continents=pygplates_continent_features,
        COBs=pygplates_COB_features,
    )
    pygplates_gplot.time = 0

    # Test both the plot object created using PlateModelManager (ie, using filenames) and
    # the plot object created using actual pygplates objects (ie, not filenames).
    #
    # Pickling of the former will be faster than the latter.
    for g in (gplot, pygplates_gplot):
        assert isinstance(g, gplately.PlotTopologies)
        g.time = 0
        pickled_plot = pickle.loads(pickle.dumps(g))

        assert pickled_plot is not g
        assert pickled_plot.time == pytest.approx(g.time)
        assert pickled_plot.anchor_plate_id == g.anchor_plate_id
        assert isinstance(pickled_plot._plot_engine, type(g._plot_engine))

        # Check the associated PlateReconstruction model got pickled properly.
        assert pickled_plot.plate_reconstruction.rotation_model.get_rotation(
            100.0, 701
        ) == g.plate_reconstruction.rotation_model.get_rotation(100.0, 701)
        assert pickled_plot.plate_reconstruction.topology_features and len(
            g.plate_reconstruction.topology_features
        ) == len(pickled_plot.plate_reconstruction.topology_features)
        assert pickled_plot.plate_reconstruction.static_polygons and len(
            g.plate_reconstruction.static_polygons
        ) == len(pickled_plot.plate_reconstruction.static_polygons)

        # Test pickled coastline/continent/COB.
        # Since these are handled specially when pickling.
        assert g.coastlines
        assert pickled_plot.coastlines
        assert len(pickled_plot.coastlines) == len(g.coastlines)

        assert g.continents
        assert pickled_plot.continents
        assert len(pickled_plot.continents) == len(g.continents)

        assert g.COBs
        assert pickled_plot.COBs
        assert len(pickled_plot.COBs) == len(g.COBs)

        # Reconstructed trench sections should be available after unpickling.
        assert len(pickled_plot.trench_left) == len(g.trench_left)
        assert len(pickled_plot.trench_right) == len(g.trench_right)


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
