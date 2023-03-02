import pytest
import gplately
import numpy as np
import os
import tempfile
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from gplately import DataCollection

## ==========================

# We will test GPlately functionalities on the Müller et al. (2019) plate reconstruction
# model at 0 and 100 Ma.
reconstruction_times = [0, 100]
gridding_times = [249., 250.]
pt_lon = np.array([-155.4696, 164.3])
pt_lat = np.array([19.8202, 53.5])
test_geometry_n_points = 1000
test_geometry_origins = ((20, -10), (175, -40))  # (lon, lat)
test_geometry_radii = (100, 500, 1000, 2000)  # km
test_geometry_azimuths = (45, -100)  # degrees


@pytest.fixture(scope="module")
def gplately_muller_server():
    return gplately.DataServer("Muller2019")


@pytest.fixture(scope="module")
def gplately_merdith_server():
    return gplately.DataServer("Merdith2021")


@pytest.fixture(scope="module")
def gplately_muller_static_geometries(gplately_muller_server):
    return gplately_muller_server.get_topology_geometries()


@pytest.fixture(scope="module")
def gplately_merdith_static_geometries(gplately_merdith_server):
    coastlines, continents, _ = gplately_merdith_server.get_topology_geometries()
    return coastlines, continents


@pytest.fixture(scope="module")
def gplately_muller_reconstruction_files(gplately_muller_server):
    return gplately_muller_server.get_plate_reconstruction_files()


@pytest.fixture(scope="module")
def gplately_merdith_reconstruction_files(gplately_merdith_server):
    return gplately_merdith_server.get_plate_reconstruction_files()


@pytest.fixture(scope="module")
def gplately_plate_reconstruction_object(gplately_muller_reconstruction_files):
    return gplately.PlateReconstruction(
        *gplately_muller_reconstruction_files
    )


@pytest.fixture(scope="module")
def gplately_merdith_reconstruction(gplately_merdith_reconstruction_files):
    return gplately.PlateReconstruction(
        *gplately_merdith_reconstruction_files
    )


@pytest.fixture(scope="module")
def gplately_plot_topologies_object(
    gplately_plate_reconstruction_object,
    gplately_muller_static_geometries,
):
    time = 0
    return gplately.PlotTopologies(
        gplately_plate_reconstruction_object,
        time,
        *gplately_muller_static_geometries,
    )


@pytest.fixture(scope="module")
def gplately_points_object(gplately_plate_reconstruction_object):
    model = gplately_plate_reconstruction_object
    time = 0 #Ma, will change to 100 Ma in test 2.

    # For example: Longitude and latitude of the Hawaiian-Emperor Seamount chain seed points 
    pt_lon = np.array([-155.4696, 164.3])
    pt_lat = np.array([19.8202, 53.5])

    # Call the Points object: pass the PlateReconstruction object, and the latitudes and longitudes of the seed points.
    gpts = gplately.Points(model, pt_lon, pt_lat)
    return gpts


@pytest.fixture(scope="module")
def gplately_raster_object(
    gplately_muller_server,
    gplately_plate_reconstruction_object,
):
    model = gplately_plate_reconstruction_object
    time = 0
    masked_age_grid = gplately_muller_server.get_age_grid(time)

    graster = gplately.Raster(model, data=masked_age_grid, extent=[-180,180,-90,90])
    return graster


@pytest.fixture(scope="module")
def gplately_merdith_raster(
    gplately_merdith_server,
    gplately_merdith_reconstruction,
):
    etopo = gplately_merdith_server.get_raster("ETOPO1_grd").astype("float")
    downsampled = etopo[::15, ::15]
    raster = gplately.Raster(
        plate_reconstruction=gplately_merdith_reconstruction,
        data=downsampled,
        origin="lower",
    )
    return raster


# Create a temporary directory for testing seafloorgrid
#@pytest.fixture(scope="module")
#def temp_save_directory():
#    tmpdir = tempfile.mkdtemp()
#    os.mkdir(tmpdir, exist_ok=True)
#    return tmpdir


@pytest.fixture(scope="session")
def test_save_directory(tmp_path_factory):
    parent_dir = tmp_path_factory.mktemp("seafloorgrid_test")
    return parent_dir


@pytest.fixture(scope="module")
def gplately_seafloorgrid_object(
    gplately_plate_reconstruction_object, 
    gplately_plot_topologies_object, 
    test_save_directory
):
    
    model = gplately_plate_reconstruction_object
    gplot = gplately_plot_topologies_object

    seafloorgrid = gplately.SeafloorGrid(
        model, 
        gplot, 
        max_time = 250,
        min_time = 249,
        ridge_time_step=1.,
        save_directory=str(test_save_directory),
        file_collection = "Muller2019",
        spacing_degrees = 0.25
    )
    return seafloorgrid

