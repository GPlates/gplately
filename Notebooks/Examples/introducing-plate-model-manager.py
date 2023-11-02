#!/usr/bin/env python3

from plate_model_manager import PlateModelManager


def main():
    pm_manager = PlateModelManager()

    print("available models: ")
    print("*****************")
    for name in pm_manager.get_available_model_names():
        print(name)
    print()

    model = pm_manager.get_model("Muller2019")
    model.set_data_dir("plate-model-repo")

    print("available layers in model Muller2019:")
    print("*************************************")
    for layer in model.get_avail_layers():
        print(layer)
    print()

    print("rotation file(s):")
    print("*****************")
    print(model.get_rotation_model())
    print()

    print("StaticPolygons file(s):")
    print("***********************")
    print(model.get_layer("StaticPolygons"))
    print()

    # for layer in model.get_avail_layers():
    #    print(model.get_layer(layer))


if __name__ == "__main__":
    main()
