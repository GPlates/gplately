# %% [markdown]
# ### This notebook demonstrate how to create a more interesting paleo-map using GPlately.
# %% [markdown]
# **⚠️🚫 Don't commit changes to this notebook directly into GitHub repository. 🚫⚠️**

# This notebook is generated from plot_map_with_cartopy.py using the command
# `jupytext --to notebook Notebooks/Examples/plot_map_with_cartopy.py -o Notebooks/Examples/plot_map_with_cartopy.ipynb`.
# If you need to commit changes to this notebook to the GPlately repository, make your edits in plot_map_with_cartopy.py and a GitHub workflow will regenerate this Jupyter Notebook file automatically.
# The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.
# %%
import warnings, os

warnings.filterwarnings("ignore", category=UserWarning, module="gplately")
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from gplately import PlateModelManager, Raster, auxiliary
from gplately.plot.gmt_cpt import get_cmap_from_gmt_cpt 

data_dir = "./gplately-example-data"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
# download age grid CPT file from https://raw.githubusercontent.com/GPlates/gplately/refs/heads/master/tests-dir/unittest/create-age-grids-video/agegrid.cpt
cpt_file = f"{data_dir}/agegrid.cpt"
if not os.path.isfile(cpt_file):
    import urllib.request

    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/GPlates/gplately/refs/heads/master/tests-dir/unittest/create-age-grids-video/agegrid.cpt",
        cpt_file,
    )
# %%
# use the auxiliary function to create a PlotTopologies instance
gplot = auxiliary.get_gplot("Muller2019", time=100)  # 100Ma
assert gplot.time is not None

# download the age grid at the reconstruction time
agegrid = Raster(
    data=PlateModelManager()
    .get_model("Muller2019")
    .get_raster("AgeGrids", int(gplot.time))  # type: ignore
)


fig = plt.figure(figsize=(8, 4))
ax1 = fig.add_subplot(111, projection=ccrs.Mollweide(190))

# plot something for fun
gplot.plot_continents(ax1, facecolor="0.8")
gplot.plot_coastlines(ax1, color="0.5")
gplot.plot_all_topological_sections(
    ax1,
    plot_subduction_teeth=True,
    other_kwargs={
        "color": "grey",
        "linewidth": 0.5,
    },
    ridge_kwargs={
        "color": "black",
        "linewidth": 0.7,
    },
    transform_kwargs={
        "color": "green",
        "linewidth": 0.7,
    },
    trench_kwargs={
        "color": "blue",
        "linewidth": 0.7,
    },
)
im = gplot.plot_grid(
    ax1, agegrid.data, cmap=get_cmap_from_gmt_cpt(cpt_file), vmin=0, vmax=200
)
gplot.plot_plate_motion_vectors(
    ax1, spacingX=10, spacingY=10, normalise=True, zorder=10, alpha=0.5
)
assert im is not None

fig.colorbar(im, orientation="horizontal", shrink=0.4, pad=0.05, label="Age (Ma)")
plt.title(f"{int(gplot.time)} Ma")


# save the map as a .png file
output_file = f"{data_dir}/plot_map_with_cartopy.png"
fig.savefig(output_file, dpi=120, bbox_inches="tight")  # transparent=True)
print(f"Done! The {output_file} has been saved successfully.")

plt.show()
plt.close(fig)
