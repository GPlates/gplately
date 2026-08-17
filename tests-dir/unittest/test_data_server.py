#!/usr/bin/env python3
import os

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
from common import get_logger

# pyright: reportMissingImports=false
test_new_parameters = True
try:
    from plate_model_manager import (
        ReferenceFrame,
        GenerationMethod,
    )
except ImportError:
    test_new_parameters = False
    from gplately.data_server import (
        ReferenceFrame,  # pyright: ignore[reportAttributeAccessIssue]
        GenerationMethod,  # pyright: ignore[reportAttributeAccessIssue]
    )
logger = get_logger()

import gplately

print(gplately.__file__)


def test_data_server_1():
    logger.info("Start test_data_server .......")

    logger.info(gplately.auxiliary.get_data_server_cache_path())

    data_server = gplately.DataServer("Muller2019")

    rotation_model, topology_features, static_polygons = (
        data_server.get_plate_reconstruction_files()
    )
    logger.info(rotation_model)
    logger.info(topology_features)
    logger.info(static_polygons)

    coastlines, continents, COBs = data_server.get_topology_geometries()
    logger.info(coastlines)
    logger.info(continents)
    logger.info(COBs)

    r1 = data_server.get_raster("ETOPO1_tif")
    assert isinstance(r1, gplately.Raster)
    logger.info(r1)

    r2 = data_server.get_age_grid(times=100)
    assert isinstance(r2, gplately.Raster)
    logger.info(r2)

    r3 = data_server.get_age_grid(times=[100, 50, 0])
    assert len(r3) == 3
    logger.info(r3)


def test_data_server_2():
    logger.info("Start test_data_server_2 .......")

    data_server = gplately.DataServer("Zahirovic2022")

    r4 = data_server.get_age_grid(
        times=100,
        reference_frame=ReferenceFrame.PmagReferenceFrame,
        generated_from=GenerationMethod.Isochrons,
    )
    assert isinstance(r4, gplately.Raster)
    logger.info(r4)

    r7 = data_server.get_spreading_rate_grid(
        times=100,
        reference_frame=ReferenceFrame.PmagReferenceFrame,
        generated_from=GenerationMethod.Topologies,
    )
    assert isinstance(r7, gplately.Raster)
    logger.info(r7)

    r5 = data_server.get_age_grid(times=[10, 50, 100])
    assert len(r5) == 3
    logger.info(r5)

    data_server = gplately.DataServer("Muller2025")

    r6 = data_server.get_age_grid(times=100)
    assert isinstance(r6, gplately.Raster)
    logger.info(r6)

    r8 = data_server.get_spreading_rate_grid(times=100)
    assert isinstance(r8, gplately.Raster)
    logger.info(r8)


if __name__ == "__main__":

    test_data_server_1()
    if test_new_parameters:
        test_data_server_2()
    print(gplately.DataServer.get_feature_data("SeafloorFabric"))
    print(gplately.DataServer.get_feature_data("Johansson2018"))
    print(gplately.DataServer.get_feature_data("Whittaker2015"))
    print(gplately.DataServer.get_feature_data("Hotspots"))
    logger.info("test_data_server has finished successfully!")
