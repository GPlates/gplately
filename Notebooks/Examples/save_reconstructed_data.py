# %% [markdown]
# ### This example demonstrates how to save the reconstructed data to shapefiles.

# The "get_" methods of PlotTopologies class return GeoDataFrame objects whose `to_file()` method can be used to save the data to other file formats, check out the [GeoDataFrame doc](https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_file.html).

# %% [markdown]
# **⚠️🚫 Don't commit changes to this notebook directly into GitHub repository. 🚫⚠️**

# This notebook is generated from save_reconstructed_data.py using the command
# `jupytext --to notebook Notebooks/Examples/save_reconstructed_data.py -o Notebooks/Examples/save_reconstructed_data.ipynb`.
# If you need to commit changes to this notebook to the GPlately repository, make your edits in save_reconstructed_data.py and a GitHub workflow will regenerate this Jupyter Notebook file automatically.
# The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.
# %%
import os

from gplately import auxiliary

# use the auxiliary function to create a PlotTopologies instance
gplot = auxiliary.get_gplot("Muller2019", time=100)  # 100Ma

# get GeoDataFrame objects and save them to shapefiles
data_dir = "gplately-example-data"
gplot.get_continents().to_file(f"{data_dir}/continents_100Ma.shp")
gplot.get_coastlines().to_file(f"{data_dir}/coastlines_100Ma.shp")
gplot.get_topological_plate_boundaries().to_file(
    f"{data_dir}/topological_plate_boundaries_100Ma.shp"
)
print(f"The files have been saved in folder {data_dir}.")
# %%
# list the created files in the `gplately-example-data` folder and print their names.
# sort the files by name to make sure the order is consistent across different operating systems.
for filename in sorted(os.listdir(data_dir)):
    if (
        filename.startswith("continents_100Ma")
        or filename.startswith("coastlines_100Ma")
        or filename.startswith("topological_plate_boundaries_100Ma")
    ):
        print(f" - {filename}")
