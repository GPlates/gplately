# %% [markdown]
# ### This example demonstrates how to use your own plate model and reconstruct points.

# %% [markdown]
# **⚠️🚫 Don't commit changes to this notebook directly into GitHub repository. 🚫⚠️**

# This notebook is generated from use_your_own_plate_model.py using the command
# `jupytext --to notebook Notebooks/Examples/use_your_own_plate_model.py -o Notebooks/Examples/use_your_own_plate_model.ipynb`.
# If you need to commit changes to this notebook to the GPlately repository, make your edits in use_your_own_plate_model.py and a GitHub workflow will regenerate this Jupyter Notebook file automatically.
# The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.
# %% [markdown]
# #### Step 0: Download the test files for this example.**
#
# Don't worry about the `FileDownloader`. If you'd like, you may download the zip file manually.

# %%
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pygplates
from plate_model_manager.utils import download

import gplately
from gplately.auxiliary import get_plate_reconstruction

time = 50
data_dir = "gplately-example-data"

downloader = download.FileDownloader(
    "https://repo.gplates.org/webdav/gplately-test-data/test_model.zip",
    f"{data_dir}/.metadata.json",
    f"{data_dir}",
)
# only re-download when necessary
if downloader.check_if_file_need_update():
    downloader.download_file_and_update_metadata()
else:
    print(f"The local files are still good. No need to download again!")

# %% [markdown]
# #### Step 1: Create `PlateReconstruction` and `PlotTopologies` objects with your own plate model


# %%
model = gplately.PlateReconstruction(
    rotation_model=f"{data_dir}/test_model/test_rotations.rot",
    topology_features=f"{data_dir}/test_model/test_topology.gpmlz",
    static_polygons=f"{data_dir}/test_model/test_static_polygons.gpmlz",
)
gplot = gplately.PlotTopologies(
    model,
    coastlines=f"{data_dir}/test_model/test_coastlines.gpmlz",
    continents=f"{data_dir}/test_model/test_continental_polygons.shp",
    COBs=f"{data_dir}/test_model/test_cobs.gpmlz",
    time=time,
)

# %% [markdown]
# #### Step 2: Reconstruct points in shapefile


# %%
points = [
    p.get_geometry().to_lat_lon()
    for p in pygplates.FeatureCollection(f"{data_dir}/test_model/Australia_Points.shp")
]

g_points = gplately.Points(
    plate_reconstruction=model, lons=[p[1] for p in points], lats=[p[0] for p in points]
)
reconstructed_points = g_points.reconstruct(time=time, return_array=True)

# %% [markdown]
# #### Step 3: Plot the map


# %%
ax = plt.figure(figsize=(8, 6)).add_subplot(
    111, projection=ccrs.Robinson(central_longitude=180)
)

ax.scatter(
    reconstructed_points[0],  # type: ignore
    reconstructed_points[1],  # type: ignore
    transform=ccrs.PlateCarree(),
    marker="o",
    color="blue",
)
gplot.plot_coastlines(ax, color="grey")
ax.set_global()  # type: ignore
plt.title(f"{time} Ma")
plt.show()
