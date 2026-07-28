# %% [markdown]
# ### This is a minimal working example of GPlately to demonstrate how easy it is to create a paleo-map using GPlately.
# %% [markdown]
# **⚠️🚫 Don't commit changes to this notebook directly into GitHub repository. 🚫⚠️**

# This notebook is generated from hello_world.py using the command
# `jupytext --to notebook Notebooks/Examples/hello_world.py -o Notebooks/Examples/hello_world.ipynb`.
# If you need to commit changes to this notebook to the GPlately repository, make your edits in hello_world.py and a GitHub workflow will regenerate this Jupyter Notebook file automatically.
# The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.

# %%
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import gplately

# create a basemap using Mollweide projection
ax = plt.figure(figsize=(8, 4)).add_subplot(111, projection=ccrs.Mollweide(180))

# get a PlotTopologies object
gplot = gplately.auxiliary.get_gplot("Muller2019", time=100)

# use the PlotTopologies object to plot a paleo-coastlines map
gplot.plot_coastlines(ax, color="lightblue", facecolor="0.8")

# add title for the map
plt.title(f"Paleo-coastlines at {int(gplot.time)} Ma")  # pyright: ignore

# save the map to a .png file
plt.gcf().savefig("gplately-hello-world.png")

# display the map
plt.show()
