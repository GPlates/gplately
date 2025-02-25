#!/usr/bin/env python3

from gplately.download import DataServer

# the DataServer class has been superseded by a newer module https://pypi.org/project/plate-model-manager/.
# although we will keep the DataServer interface indefinitely,
# you are encouraged to use the newer plate-model-manager instead.

# create a DataServer object/instance for the reconstruction model Muller2019
data_server = DataServer("Muller2019")

# now download the plate reconstruction files and geometries from the Müller et al. 2019 model
rotation_model, topology_features, static_polygons = (
    data_server.get_plate_reconstruction_files()
)
coastlines, continents, COBs = data_server.get_topology_geometries()
print(rotation_model)
print(topology_features)
print(static_polygons)
print(coastlines)
print(continents)
print(COBs)

# download the age grid at 100Ma
age_grid = data_server.get_age_grid(times=100)
print(age_grid)

# download the ETOPO1 geotiff raster
etopo = data_server.get_raster("ETOPO1_tif")
print(etopo)

# as we can see from the printout, the methods in the DataServer class return Python objects
# if you would like to know where the files are saved, or
# you would like to save the files in a folder that you choose,
# you might want to check out the new plate-model-manager module.
