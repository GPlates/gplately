import sys

sys.path.insert(0, "../")
from gplately import plate_model

model = plate_model.PlateModel("Muller2019", data_dir="test-plate-model-folder")

print(model.get_avail_layers())

print(model.get_rotation_model())

print(model.get_layer("Coastlines"))

print(model.get_COBs())

print(model.get_topologies())

model.download_all_layers()

model.download_time_dependent_rasters("AgeGrids", times=[1, 2])

print(model.get_data_dir())
