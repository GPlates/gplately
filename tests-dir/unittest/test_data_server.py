#!/usr/bin/env python3
import os

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
from common import get_logger

import gplately

print(gplately.__file__)

if __name__ == "__main__":

    logger = get_logger()

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

    logger.info("test_data_server has finished successfully!")
