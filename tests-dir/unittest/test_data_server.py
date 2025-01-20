#!/usr/bin/env python3
import gplately

if __name__ == "__main__":
    gDownload = gplately.DataServer("Muller2019")
    rotation_model, topology_features, static_polygons = (
        gDownload.get_plate_reconstruction_files()
    )
    print(rotation_model)
    print(topology_features)
    print(static_polygons)
    coastlines, continents, COBs = gDownload.get_topology_geometries()
    print(coastlines)
    print(continents)
    print(COBs)
    r1 = gDownload.get_raster("ETOPO1_tif")
    print(r1)
    r2 = gDownload.get_age_grid(times=100)
    print(r2)
