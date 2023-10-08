import gzip
import os
import shutil

import numpy as np
import pygplates
import pytest
import gplately
from plate_model_manager import (
    PlateModelManager,
    PresentDayRasterManager,
    network_requests,
)


## ==========================

# We will test GPlately functionalities on the Müller et al. (2019) plate reconstruction
# model at 0 and 100 Ma.
reconstruction_times = [0, 100]
gridding_times = [249.0, 250.0]
pt_lon = np.array([-155.4696, 164.3])
pt_lat = np.array([19.8202, 53.5])
test_geometry_n_points = 1000
test_geometry_origins = ((20, -10), (175, -40))  # (lon, lat)
test_geometry_radii = (100, 500, 1000, 2000)  # km
test_geometry_azimuths = (45, -100)  # degrees


@pytest.fixture(scope="module")
def muller_2019_model():
    pm_manger = PlateModelManager()
    return pm_manger.get_model("Muller2019")


@pytest.fixture(scope="module")
def merdith_2021_model():
    pm_manger = PlateModelManager()
    return pm_manger.get_model("Merdith2021")


@pytest.fixture(scope="module")
def gplately_muller_static_geometries(muller_2019_model):
    coastlines = pygplates.FeatureCollection()
    for file in muller_2019_model.get_layer("Coastlines"):
        coastlines.add(pygplates.FeatureCollection(file))
    continental_polygons = pygplates.FeatureCollection()
    for file in muller_2019_model.get_layer("ContinentalPolygons"):
        continental_polygons.add(pygplates.FeatureCollection(file))
    COBs = pygplates.FeatureCollection()
    for file in muller_2019_model.get_layer("COBs"):
        COBs.add(pygplates.FeatureCollection(file))
    return coastlines, continental_polygons, COBs


@pytest.fixture(scope="module")
def gplately_merdith_static_geometries(merdith_2021_model):
    coastlines = pygplates.FeatureCollection()
    for file in merdith_2021_model.get_layer("Coastlines"):
        coastlines.add(pygplates.FeatureCollection(file))
    continental_polygons = pygplates.FeatureCollection()
    for file in merdith_2021_model.get_layer("ContinentalPolygons"):
        continental_polygons.add(pygplates.FeatureCollection(file))

    return coastlines, continental_polygons


@pytest.fixture(scope="module")
def gplately_muller_reconstruction_files(muller_2019_model):
    rotation_model = pygplates.RotationModel(muller_2019_model.get_rotation_model())
    topology_features = pygplates.FeatureCollection()
    for file in muller_2019_model.get_layer("Topologies"):
        topology_features.add(pygplates.FeatureCollection(file))
    static_polygons = pygplates.FeatureCollection()
    for file in muller_2019_model.get_layer("StaticPolygons"):
        static_polygons.add(pygplates.FeatureCollection(file))
    return rotation_model, topology_features, static_polygons


@pytest.fixture(scope="module")
def gplately_merdith_reconstruction_files(merdith_2021_model):
    rotation_model = pygplates.RotationModel(merdith_2021_model.get_rotation_model())
    topology_features = pygplates.FeatureCollection()
    for file in merdith_2021_model.get_layer("Topologies"):
        topology_features.add(pygplates.FeatureCollection(file))
    static_polygons = pygplates.FeatureCollection()
    for file in merdith_2021_model.get_layer("StaticPolygons"):
        static_polygons.add(pygplates.FeatureCollection(file))
    return rotation_model, topology_features, static_polygons


@pytest.fixture(scope="module")
def gplately_plate_reconstruction_object(muller_2019_model):
    return gplately.PlateReconstruction(
        rotation_model=muller_2019_model.get_rotation_model(),
        topology_features=muller_2019_model.get_layer("Topologies"),
        static_polygons=muller_2019_model.get_layer("StaticPolygons"),
    )


@pytest.fixture(scope="module")
def gplately_merdith_reconstruction(gplately_merdith_reconstruction_files):
    return gplately.PlateReconstruction(*gplately_merdith_reconstruction_files)


@pytest.fixture(scope="module")
def gplately_plot_topologies_object(
    gplately_plate_reconstruction_object,
    gplately_muller_static_geometries,
):
    gplot = gplately.PlotTopologies(
        gplately_plate_reconstruction_object,
        *gplately_muller_static_geometries,
    )
    gplot.time = 0
    return gplot


@pytest.fixture(scope="module")
def gplately_points_object(gplately_plate_reconstruction_object):
    model = gplately_plate_reconstruction_object
    time = 0  # Ma, will change to 100 Ma in test 2.

    # For example: Longitude and latitude of the Hawaiian-Emperor Seamount chain seed points
    pt_lon = np.array([-155.4696, 164.3])
    pt_lat = np.array([19.8202, 53.5])

    # Call the Points object: pass the PlateReconstruction object, and the latitudes and longitudes of the seed points.
    gpts = gplately.Points(model, pt_lon, pt_lat)
    return gpts


@pytest.fixture(scope="module")
def gplately_raster_object(
    muller_2019_model,
    gplately_plate_reconstruction_object,
):
    model = gplately_plate_reconstruction_object
    time = 0
    agegrid_path = muller_2019_model.get_raster("AgeGrids", time)

    graster = gplately.Raster(
        data=agegrid_path,
        plate_reconstruction=model,
        extent=[-180, 180, -90, 90],
    )

    return graster


@pytest.fixture(scope="module")
def gplately_merdith_raster(
    gplately_merdith_reconstruction,
):
    raster_manager = PresentDayRasterManager()
    etopo = gplately.Raster(data=raster_manager.get_raster("ETOPO1_grd"))
    etopo = etopo.data.astype("float")
    downsampled = etopo[::15, ::15]
    raster = gplately.Raster(
        plate_reconstruction=gplately_merdith_reconstruction,
        data=downsampled,
        origin="lower",
    )
    return raster


@pytest.fixture(scope="module")
def gplately_seafloorgrid_object(
    gplately_plate_reconstruction_object, gplately_plot_topologies_object
):
    model = gplately_plate_reconstruction_object
    gplot = gplately_plot_topologies_object

    seafloorgrid = gplately.SeafloorGrid(
        model,
        gplot,
        max_time=250,
        min_time=249,
        ridge_time_step=1.0,
        save_directory="test-seafloor-grid",
        file_collection="Muller2019",
        grid_spacing=0.25,
    )
    return seafloorgrid
