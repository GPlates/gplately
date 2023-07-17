import sys

sys.path.insert(0, "../")
from gplately import plate_model

model = plate_model.PlateModel("Muller2019")

print(model.get_avail_layers())

print(model.get_rotation_model())
