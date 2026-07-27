# %% [markdown]
# ### This example demonstrates how to use the `PlateModelManager` with GPlately.

# %% [markdown]
# ⚠️ This notebook is generated from working_with_plate_model_manager.py using the command
# `jupytext --to notebook Notebooks/Examples/working_with_plate_model_manager.py -o Notebooks/Examples/working_with_plate_model_manager.ipynb`.
# If you need to commit changes to this notebook to the GPlately repository, make your edits in working_with_plate_model_manager.py and a GitHub workflow will regenerate this Jupyter Notebook file automatically.
# The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.

# %%
import warnings

warnings.filterwarnings("ignore", category=UserWarning, module="gplately")
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from plate_model_manager import PlateModelManager

from gplately import PlateReconstruction, PlotTopologies

# use `PlateModelManager` to create `PlateReconstruction` and `PlotTopologies` objects
model = PlateModelManager().get_model("Muller2019")
model.set_data_dir("plate-model-repo")  # type: ignore

age = 55
test_model = PlateReconstruction(
    model.get_rotation_model(),  # type: ignore
    topology_features=model.get_layer("Topologies"),  # type: ignore
    static_polygons=model.get_layer("StaticPolygons"),  # type: ignore
)
gplot = PlotTopologies(
    test_model,
    coastlines=model.get_layer("Coastlines"),  # type: ignore
    COBs=model.get_layer("COBs"),  # type: ignore
    time=age,
)

fig = plt.figure(figsize=(12, 6), dpi=72)
ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=180))
ax.set_global()  # type: ignore

# now use PlotTopologies object to plot some model data
gplot.plot_continent_ocean_boundaries(ax, color="cornflowerblue")
gplot.plot_coastlines(ax, color="black")
gplot.plot_all_topological_sections(
    ax,
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

plt.title(f"{age} Ma")

# save the map as a .png file
output_file = f"./gplately-example-data/working_with_pmm.png"
fig.savefig(output_file, dpi=120, bbox_inches="tight")  # transparent=True)
print(f"Done! The output file {output_file} has been saved.")

plt.show()
plt.close(fig)
