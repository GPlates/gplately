#!/usr/bin/env python

from common import MODEL_REPO_DIR
from plate_model_manager import PlateModelManager

import gplately

print(gplately.__version__)
print(gplately.__file__)


def main():
    pm_manger = PlateModelManager()
    model = pm_manger.get_model("Muller2019", data_dir=MODEL_REPO_DIR)

    print(model.get_avail_layers())

    print(model.get_rotation_model())

    print(model.get_layer("Coastlines"))

    print(model.get_COBs())

    print(model.get_topologies())

    model.download_all_layers()

    model.download_time_dependent_rasters("AgeGrids", times=[1, 2])

    print(model.get_data_dir())


if __name__ == "__main__":
    main()
