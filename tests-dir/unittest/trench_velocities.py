import gplately
import numpy as np
from plate_model_manager import PlateModelManager

zahirovic2022_model = PlateModelManager().get_model(
    "Zahirovic2022", data_dir="plate-model-repo"
)
assert zahirovic2022_model

z22_reconstruction_model = gplately.PlateReconstruction(
    zahirovic2022_model.get_rotation_model(),
    zahirovic2022_model.get_topologies(),
    zahirovic2022_model.get_static_polygons(),
)

# change this variable for other reconstruction time
reconstruction_time = 0

subduction_data = z22_reconstruction_model.tessellate_subduction_zones(
    reconstruction_time
)

# Latitudes and longitudes of points along trench segments
subduction_lon = subduction_data[:, 0]
subduction_lat = subduction_data[:, 1]
print(f"The longitudes of trench sample locations({len(subduction_lon)}):")
print(subduction_lon)
print()
print(f"The latitudes of trench sample locations({len(subduction_lat)}):")
print(subduction_lat)
print()

# Absolute velocity(cm/year) of the trench’s motion in the orthogonal direction towards the overriding plate.
# Negative if moving towards the overriding plate (Trench advance)
# Positive if moving away from the overriding plate (Trench retreat)
# The purpose of "-np.fabs()" is to allow "np.cos()" to produce the correct sign of value(negative or positive).
trench_absolute_vel = -np.fabs(subduction_data[:, 4]) * np.cos(
    np.radians(subduction_data[:, 5])
)
print(
    f"The trench absolute velocities({len(trench_absolute_vel)}) at {reconstruction_time} Ma:"
)
print(trench_absolute_vel)
print()

# Global mean and standard deviation of trench velocity(cm/year) of all trench segments
trench_abs_vel_mean = np.mean(trench_absolute_vel)
trench_abs_vel_std = np.std(trench_absolute_vel)
print(f"The mean of trench absolute velocities at {reconstruction_time} Ma:")
print(trench_abs_vel_mean)
print()
print(
    f"The standard deviation of trench absolute velocities at {reconstruction_time} Ma:"
)
print(trench_abs_vel_std)
