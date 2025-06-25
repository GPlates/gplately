#!/usr/bin/env python3
import os

from gplately import PlateModelManager, PresentDayRasterManager
from gplately.commands.list_models import get_model_names


def main():
    pm_manager = PlateModelManager()

    #
    # get all available models(names) in the PlateModelManager
    #
    print("Available models: ")
    print("*****************")
    gplately_model_names = get_model_names()
    for name in pm_manager.get_available_model_names():
        # the pm_manager.get_available_model_names() returns a superset of GPlately models
        # we need to check the the model name agaist a list of GPlately models
        if name in gplately_model_names:
            print(name)
    print()

    #
    # download model "Muller2019" and put the files in folder "plate-model-repo"
    #
    print("Folders/files in 'Muller2019' folder: ")
    print("**************************************")
    model = pm_manager.get_model("Muller2019")
    assert model
    model.set_data_dir("plate-model-repo")
    for layer in model.get_avail_layers():
        model.get_layer(layer)

    # now let's see what are inside the "plate-model-repo/muller2019" folder
    print(os.listdir("plate-model-repo/muller2019"))
    print()

    #
    # List all vailable layers in model Muller2019
    #
    print("Available layers in model Muller2019:")
    print("*************************************")
    for layer in model.get_avail_layers():
        print(layer)
    print()

    #
    # download rotation files
    #
    print("Rotation file(s):")
    print("*****************")
    print(model.get_rotation_model())
    print()

    #
    # download static polygons
    #
    print("StaticPolygons file(s):")
    print("***********************")
    print(model.get_layer("StaticPolygons"))
    print()

    #
    # download Coastlines
    #
    print("Coastlines file(s):")
    print("******************")
    print(model.get_layer("Coastlines"))
    print()

    #
    # download all layers
    #
    print("All layer files:")
    print("***************")
    for layer in model.get_avail_layers():
        print(model.get_layer(layer))
    print()

    #
    # get a list of time dependent rasters
    #
    print("Available time dependent rasters:")
    print("*********************************")
    for n in model.get_avail_time_dependent_raster_names():
        print(n)
    print()

    #
    # download AgeGrids rasters
    #
    print(model.get_rasters("AgeGrids", times=[10, 20, 30]))
    print(model.get_raster("AgeGrids", time=100))
    print()

    #
    # download AgeGrids rasters for all times
    # this function will take some time and download a large volume data
    #
    # model.download_time_dependent_rasters("AgeGrids")

    #
    # list the names of all present-day rasters
    #
    print("Available present-day rasters:")
    print("******************************")
    print(PresentDayRasterManager().list_present_day_rasters())
    print()

    #
    # get "topography" present-day raster
    #
    print("Present-day topography raster:")
    print("******************************")
    print(PresentDayRasterManager().get_raster("topography"))
    print()


if __name__ == "__main__":
    main()
