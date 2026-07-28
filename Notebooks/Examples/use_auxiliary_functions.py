# %% [markdown]
# ### This example demonstrates how to use the auxiliary functions to create the `PlateReconstruction` and `PlotTopologies` objects quickly.

# %% [markdown]
# ⚠️ This notebook is generated from use_auxiliary_functions.py using the command
# `jupytext --to notebook Notebooks/Examples/use_auxiliary_functions.py -o Notebooks/Examples/use_auxiliary_functions.ipynb`.
# If you need to commit changes to this notebook to the GPlately repository, make your edits in use_auxiliary_functions.py and a GitHub workflow will regenerate this Jupyter Notebook file automatically.
# The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.
# %%
from gplately.auxiliary import get_gplot, get_plate_reconstruction

# use the auxiliary function to create a PlateReconstruction object
plate_reconstruction_obj = get_plate_reconstruction("Muller2019")
print(plate_reconstruction_obj)

# use the auxiliary function to create a PlotTopologies object
plot_topologies_obj = get_gplot("Muller2019", time=140)
print(plot_topologies_obj)

# %% [markdown]
# There is a `PlateReconstruction` object inside the `PlotTopologies` object.
# So, in most cases, a single get_gplot() call is enough.
# You can retrieve the `PlateReconstruction` object from the `PlotTopologies` object later.
# %%
plate_reconstruction_obj_within_plot_topologies_obj = (
    plot_topologies_obj.plate_reconstruction
)
print(plate_reconstruction_obj_within_plot_topologies_obj)
