#!/usr/bin/env python3

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pygplates
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModel, PlateModelManager

MODEL_NAME = "Muller2019"

try:
    model = PlateModelManager().get_model(MODEL_NAME, data_dir=MODEL_REPO_DIR)
except:
    model = PlateModel(MODEL_NAME, data_dir=MODEL_REPO_DIR, readonly=True)

age = 52

resolved_topologies = []
shared_boundary_sections = []
pygplates.resolve_topologies(
    model.get_layer("Topologies"),
    model.get_rotation_model(),
    resolved_topologies,
    52,
    shared_boundary_sections,
)

fig = plt.figure(figsize=(10, 5), dpi=96)
ax = fig.add_subplot(
    111, projection=ccrs.Orthographic(central_longitude=-165, central_latitude=15)
)
for t in resolved_topologies:
    lat, lon = zip(*(t.get_resolved_boundary().to_lat_lon_list()))
    plt.plot(
        lon,
        lat,
        color="blue",
        linewidth=0.5,  # the topological plates boundaries in blue
        transform=ccrs.Geodetic(),
    )


plt.show()
