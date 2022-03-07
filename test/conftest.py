import pytest
import gplately
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from gplately import DataCollection

## ==========================

# We will test GPlately functionalities on the Müller et al. (2019) plate reconstruction
# model at 0 and 100 Ma.
reconstruction_times = [0, 100]
pt_lon = np.array([-155.4696, 164.3])
pt_lat = np.array([19.8202, 53.5])


# SET UP ALL MAIN GPLATELY OBJECTS FOR TESTING 
@pytest.fixture(scope="module")
def gplately_data_server_object():
    gdownload = gplately.download.DataServer("Muller2019")
    return gdownload


@pytest.fixture(scope="module")
def gplately_plate_reconstruction_object(gplately_data_server_object):
    gdownload = gplately_data_server_object
    rotation_model, topology_features, static_polygons = gdownload.get_plate_reconstruction_files()
    model = gplately.PlateReconstruction(
        rotation_model, 
        topology_features, 
        static_polygons
    )
    return model


@pytest.fixture(scope="module")
def gplately_plot_topologies_object(
    gplately_plate_reconstruction_object, 
    gplately_data_server_object
    ):
    model = gplately_plate_reconstruction_object
    time = 0 #Ma, will change to 100 when called in test_3.
    gdownload = gplately_data_server_object
    coastlines, continents, COBs = gdownload.get_topology_geometries()
    gplot = gplately.plot.PlotTopologies(model, time, coastlines, continents, COBs)
    return gplot


@pytest.fixture(scope="module")
def gplately_geo_axis_param():
    fig = plt.figure(figsize=(16,12), dpi=100)
    ax = fig.add_subplot(111, projection=ccrs.Mollweide(central_longitude = 20))
    return ax


@pytest.fixture(scope="module")
def gplately_points_object(request, gplately_plate_reconstruction_object):
    model = gplately_plate_reconstruction_object
    time = 0 #Ma, will change to 100 Ma in test 2.

    # For example: Longitude and latitude of the Hawaiian-Emperor Seamount chain seed points 
    pt_lon = np.array([-155.4696, 164.3])
    pt_lat = np.array([19.8202, 53.5])

    # Call the Points object: pass the PlateReconstruction object, and the latitudes and longitudes of the seed points.
    gpts = gplately.Points(model, pt_lon, pt_lat)
    return gpts


@pytest.fixture(scope="module")
def download_multifeature_gplates_sample_data(gplately_data_server_object):
    gdownload = gplately_data_server_object
    johansson_2018 = gdownload.get_feature_data("Johannson2018")
    assert isinstance(johansson_2018, pygplates.FeatureCollection), "Cached GPlates 2.3 sample data not an instance of <pygplates.FeatureCollection>."


@pytest.fixture(scope="module")
def download_multifeature_gplates_sample_data(gplately_data_server_object):
    gdownload = gplately_data_server_object
    seafloor_fabric = gdownload.get_feature_data("SeafloorFabric")
    assert [isinstance(sf, pygplates.FeatureCollection) for sf in seafloor_fabric], "Cached GPlates 2.3 sample data not an instance of <pygplates.FeatureCollection>."
