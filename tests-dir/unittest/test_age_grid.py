#!/usr/bin/env python3

import os

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

from common import *
from plate_model_manager import PlateModel, PlateModelManager

from gplately import PlateReconstruction, PlotTopologies, SeafloorGrid

# test using SeafloorGrid class to generate age grids.


def main():
    model_name = "merdith2021"
    try:
        plate_model = PlateModelManager().get_model(
            model_name, data_dir="plate-model-repo"
        )
    except:
        plate_model = PlateModel(model_name, data_dir="plate-model-repo", readonly=True)

    assert plate_model

    rotation_files = plate_model.get_rotation_model()
    topology_files = plate_model.get_topologies()
    continent_files = plate_model.get_layer("ContinentalPolygons")
    if "Cratons" in plate_model.get_avail_layers():
        cratons_files = plate_model.get_layer("Cratons")
        assert cratons_files
        assert continent_files
        continent_files += cratons_files

    reconstruction = PlateReconstruction(
        rotation_model=rotation_files,
        topology_features=topology_files,
    )
    gplot = PlotTopologies(
        reconstruction,
        continents=continent_files,
    )

    if True:
        # test the new reconstruct_by_topological_model method
        grid = SeafloorGrid(
            reconstruction,
            gplot,
            min_time=400,
            max_time=410,
            save_directory="./output/test-age-grid-reconstruct-by-topological-model",
            ridge_time_step=1,
            refinement_levels=5,
            grid_spacing=0.1,
            ridge_sampling=0.5,
            initial_ocean_mean_spreading_rate=75,
            file_collection=model_name,
            resume_from_checkpoints=True,
            use_continent_contouring=False,
        )
        grid.reconstruct_by_topological_model()
        for val in (
            SeafloorGrid.SEAFLOOR_AGE_KEY,
            SeafloorGrid.SPREADING_RATE_KEY,
        ):
            grid.lat_lon_z_to_netCDF(val)

    if True:
        # test the old reconstruct_by_topologies method (*without* continent contouring)
        grid = SeafloorGrid(
            reconstruction,
            gplot,
            min_time=400,
            max_time=410,
            save_directory="./output/test-age-grid-reconstruct-by-topologies",
            ridge_time_step=1,
            refinement_levels=5,
            grid_spacing=0.1,
            ridge_sampling=0.5,
            initial_ocean_mean_spreading_rate=75,
            file_collection=model_name,
            resume_from_checkpoints=True,
            use_continent_contouring=False,
        )

        grid.reconstruct_by_topologies()
        for val in (
            SeafloorGrid.SEAFLOOR_AGE_KEY,
            SeafloorGrid.SPREADING_RATE_KEY,
        ):
            grid.lat_lon_z_to_netCDF(val, unmasked=False, nprocs=5)

    if True:
        # test the old reconstruct_by_topologies method (*with* continent contouring)
        grid = SeafloorGrid(
            reconstruction,
            gplot,
            min_time=400,
            max_time=410,
            save_directory="./output/test-age-grid-use_continent_contouring",
            ridge_time_step=1,
            refinement_levels=5,
            grid_spacing=0.1,
            ridge_sampling=0.5,
            initial_ocean_mean_spreading_rate=75,
            file_collection=model_name,
            resume_from_checkpoints=True,
            use_continent_contouring=True,
        )

        grid.reconstruct_by_topological_model()
        for val in (
            SeafloorGrid.SEAFLOOR_AGE_KEY,
            SeafloorGrid.SPREADING_RATE_KEY,
        ):
            grid.lat_lon_z_to_netCDF(val)


if __name__ == "__main__":
    main()
