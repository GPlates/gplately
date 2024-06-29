#!/usr/bin/env python3

import os

import pygplates as _pygplates
from common import *

import gplately

print(gplately.__file__)
from plate_model_manager import PlateModel, PlateModelManager

from gplately import PlateReconstruction, PlotTopologies, SeafloorGrid
from gplately.pygplates import FeatureCollection as gFeatureCollection
from gplately.pygplates import RotationModel as gRotationModel


def main():
    model_name = "merdith2021"
    try:
        plate_model = PlateModelManager().get_model(
            model_name, data_dir="plate-model-repo"
        )
    except:
        plate_model = PlateModel(model_name, data_dir="plate-model-repo", readonly=True)

    rotation_files = plate_model.get_rotation_model()
    topology_files = plate_model.get_topologies()
    continent_files = plate_model.get_layer("ContinentalPolygons")
    if "Cratons" in plate_model.get_avail_layers():
        continent_files += plate_model.get_layer("Cratons")

    rotations = gRotationModel(rotation_files)
    topologies = gFeatureCollection.from_file_list(topology_files)
    continents = gFeatureCollection.from_file_list(continent_files)

    reconstruction = PlateReconstruction(
        rotation_model=rotations,
        topology_features=topologies,
    )
    gplot = PlotTopologies(
        reconstruction,
        continents=continents,
    )

    grid = SeafloorGrid(
        reconstruction,
        gplot,
        min_time=400,
        max_time=410,
        save_directory="test-age-grid-output-0627",
        ridge_time_step=1,
        refinement_levels=5,
        grid_spacing=0.1,
        ridge_sampling=0.5,
        initial_ocean_mean_spreading_rate=75,
        file_collection=model_name,
        resume_from_checkpoints=True,
    )
    test_new = 1
    if test_new:
        grid.reconstruct_by_topological_model()
        for val in ("SEAFLOOR_AGE", "SPREADING_RATE"):
            grid.save_netcdf_files(val)
    else:
        grid.reconstruct_by_topologies()
        for val in ("SEAFLOOR_AGE", "SPREADING_RATE"):
            grid.lat_lon_z_to_netCDF(val, unmasked=False, nprocs=5)


if __name__ == "__main__":
    main()
