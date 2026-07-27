# %% [markdown]
# ### This notebook demonstrate how to use the PlateModelManager to access plate model files.
# %% [markdown]
# ⚠️ This notebook is generated from introducing_plate_model_manager.py using the command
# `jupytext --to notebook Notebooks/Examples/introducing_plate_model_manager.py -o Notebooks/Examples/introducing_plate_model_manager.ipynb`.
# If you need to commit changes to this notebook to the GPlately repository, make your edits in introducing_plate_model_manager.py and a GitHub workflow will regenerate this Jupyter Notebook file automatically.
# The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.
# %%
import os

from gplately import PlateModelManager, PresentDayRasterManager
from gplately.commands.list_models import get_model_names

pm_manager = PlateModelManager()

# %% [markdown]
# #### Get the names of all available models in the PlateModelManager

# %%
gplately_model_names = get_model_names()
for name in pm_manager.get_available_model_names():
    # the pm_manager.get_available_model_names() returns a superset of GPlately models.
    # we need to check the the model name agaist a list of GPlately officially supported models.
    # https://gplates.github.io/gplately/latest/sphinx/html/plate_models.html
    if name in gplately_model_names:
        print(name)

# %% [markdown]
# #### Download model "Muller2019" and put the files in folder "plate-model-repo"

# %%
model = pm_manager.get_model("Muller2019")
assert model
model.set_data_dir("plate-model-repo")
for layer in model.get_avail_layers():
    model.get_layer(layer)

# now let's see what are inside the "plate-model-repo/muller2019" folder
print(os.listdir("plate-model-repo/muller2019"))

# %% [markdown]
# #### List all vailable layers in model Muller2019

# %%
for layer in model.get_avail_layers():
    print(layer)

# %% [markdown]
# #### Download rotation files

# %%
rotation_files = model.get_rotation_model()
print(rotation_files)

# %% [markdown]
# #### Download static polygons

# %%
static_polygon_files = model.get_layer("StaticPolygons")
print(static_polygon_files)

# %% [markdown]
# #### Download Coastlines

# %%
coasts_files = model.get_layer("Coastlines")
print(coasts_files)

# %% [markdown]
# #### Download all layers

# %%
for layer in model.get_avail_layers():
    print(model.get_layer(layer))

# %% [markdown]
# #### Get a list of time dependent rasters

# %%
for raster in model.get_avail_time_dependent_raster_names():
    print(raster)

# %% [markdown]
# #### Download AgeGrids rasters

# %%
import warnings

warnings.filterwarnings("ignore")
print(model.get_rasters("AgeGrids", times=[10, 20, 30]))
print(model.get_raster("AgeGrids", time=100))


# %% [markdown]
# #### Download AgeGrids rasters for all available times
#
# This function will take a while to finish and download a large volume data. Uncomment the code in the code cell below to try it.

# %%
# model.download_time_dependent_rasters("AgeGrids")

# %% [markdown]
# #### List the names of all present-day rasters

# %%
print(PresentDayRasterManager().list_present_day_rasters())

# %% [markdown]
# #### Get "topography" present-day raster

# %%
print(PresentDayRasterManager().get_raster("topography"))
