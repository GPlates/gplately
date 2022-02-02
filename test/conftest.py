import pytest
import gplately

@pytest.fixture(scope="module")
def generic_reconstuction_object():

    # Call GPlately's DataServer from the download.py module
    gdownload = gplately.download.DataServer("Muller2019")

    # Obtain all rotation files, topology features and static polygons from Muller et al. 2019
    rotation_model, topology_features, static_polygons = gdownload.get_plate_reconstruction_files()

    model = gplately.PlateReconstruction(rotation_model, topology_features, static_polygons)

    return model

