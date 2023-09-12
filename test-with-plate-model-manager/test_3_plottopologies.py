import pytest
import matplotlib.pyplot as plt
import numpy as np
from conftest import reconstruction_times
from conftest import gplately_plot_topologies_object as gplot

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

    - Plotting functions:
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

# CALL THE PLOT TOPOLOGIES OBJECT
@pytest.mark.parametrize("time", reconstruction_times)
def test_gplately_plotTopologies_object(time, gplot):
    #assert gplot, "No <gplately.PlotTopologies> object made with {}.".format(model)
    assert gplot, "Unable to create a <gplately.PlotTopologies> object with Müller et al. (2019) at {} Ma.".format(time)

# ================================================================================================================================================================================================================

# ENSURE PLOT TOPOLOGIES' ATTRIBUTES EXIST AND ARE CORRECT
# Topologies
@pytest.mark.parametrize("time", reconstruction_times)
def test_PlotTopologies_topologies(time, gplot):
    assert gplot.topologies, "No topological features from Müller et al. (2019) at {} Ma are attributed to <gplately.PlotTopologies>.".format(time)  

# Ridge_transforms
@pytest.mark.parametrize("time", reconstruction_times)
def test_PlotTopologies_ridge_transforms(time, gplot):
    assert gplot.ridge_transforms, "No ridge/transform features from Müller et al. (2019) at {} Ma are attributed to <gplately.PlotTopologies>.".format(time)
    assert [ridge_transform.get_feature_type() == "gpml:MidOceanRidge" for ridge_transform in gplot.ridge_transforms], "<gplately.PlotTopologies> ridge/transforms are not all of type gpml:MidOceanRidge in Müller et al. (2019) at {} Ma.".format(time)

# Ridges
@pytest.mark.parametrize("time", reconstruction_times)
def test_PlotTopologies_ridges(time, gplot):
    assert gplot.ridges, "No ridge features from Müller et al. (2019) at {} Ma are attributed to <gplately.PlotTopologies>.".format(time)
    assert [ridge.get_feature_type() == "gpml:MidOceanRidge" for ridge in gplot.ridges], "<gplately.PlotTopologies> ridges are not all of type gpml:MidOceanRidge in Müller et al. (2019) at {} Ma.".format(time)

# Transforms
@pytest.mark.parametrize("time", reconstruction_times)
def test_PlotTopologies_transforms(time, gplot):
    assert gplot.transforms, "No transform features from Müller et al. (2019) at {} Ma are attributed to <gplately.PlotTopologies>.".format(time)
    assert [transform.get_feature_type() == "gpml:MidOceanRidge" for transform in gplot.transforms], "<gplately.PlotTopologies> transforms are not all of type gpml:MidOceanRidge in Müller et al. (2019) at {} Ma.".format(time)

# Trenches
@pytest.mark.parametrize("time", reconstruction_times)
def test_PlotTopologies_trenches(time, gplot):
    assert gplot.trenches, "No trench features from Müller et al. (2019) at {} Ma are attributed to <gplately.PlotTopologies>.".format(time)
    assert [trench.get_feature_type() == "gpml:SubductionZone" for trench in gplot.trenches], "<gplately.PlotTopologies> trenches are not all of type gpml:SubductionZone in Müller et al. (2019) at {} Ma.".format(time)

# Trench_left
@pytest.mark.parametrize("time", reconstruction_times)
def test_PlotTopologies_trench_L(time, gplot):
    assert gplot.trench_left, "No trench (L) features from Müller et al. (2019) at {} Ma are attributed to <gplately.PlotTopologies>.".format(time)
    assert [trench_l.get_feature_type() == "gpml:SubductionZone" for trench_l in gplot.trench_left], "<gplately.PlotTopologies> trenches (L) are not all of type gpml:SubductionZone in Müller et al. (2019) at {} Ma.".format(time)

# Trench_right
@pytest.mark.parametrize("time", reconstruction_times)
def test_PlotTopologies_trench_R(time, gplot):
    assert gplot.trench_right, "No trench (R) features from Müller et al. (2019) at {} Ma are attributed to <gplately.PlotTopologies>.".format(time)
    assert [trench_r.get_feature_type() == "gpml:SubductionZone" for trench_r in gplot.trench_right], "<gplately.PlotTopologies> trenches (R) are not all of type gpml:SubductionZone in Müller et al. (2019) at {} Ma.".format(time)


# Subduction teeth
def test_PlotTopologies_subduction_teeth(gplot):
    ax = plt.gca()
    gplot.plot_subduction_teeth(ax=ax)
    fig = plt.gcf()
    plt.close(fig)


def test_pickle_Points():
    import pickle
    gplot_dump = pickle.dumps(gplot)
    gplot_load = pickle.loads(gplot_dump)
