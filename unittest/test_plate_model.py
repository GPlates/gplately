import sys

sys.path.insert(0, "../")
from gplately import plate_model

model = plate_model.PlateModel("Muller2019")

# print(model.get_avail_layers())

# print(model.get_rotation_model())

# print(model.get_layer("Coastlines"))

model.download(dst_path="test-download-folder")

# model.download_time_dependent_rasters("AgeGrids", "test-age-agrids-download")
