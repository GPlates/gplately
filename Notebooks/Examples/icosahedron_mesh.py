# %% [markdown]
# ### This example demonstrates how to generate an Icosahedron mesh and plot it with [plot_trisurf](https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.plot_trisurf.html).
# %% [markdown]
# **⚠️🚫 Don't commit changes to this notebook directly into GitHub repository. 🚫⚠️**

# This notebook is generated from icosahedron_mesh.py using the command
# `jupytext --to notebook Notebooks/Examples/icosahedron_mesh.py -o Notebooks/Examples/icosahedron_mesh.ipynb`.
# If you need to commit changes to this notebook to the GPlately repository, make your edits in icosahedron_mesh.py and a GitHub workflow will regenerate this Jupyter Notebook file automatically.
# The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.
# %%
import os

import matplotlib.pyplot as plt
import numpy as np

from gplately.lib.icosahedron import get_mesh, xyz2lonlat

mesh_resolution = 3
vertices_0, faces_0 = get_mesh(mesh_resolution)
# print(vertices_0.shape, faces_0.shape)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 8))

ax.plot_trisurf(  # type: ignore
    vertices_0[:, 0],
    vertices_0[:, 1],
    vertices_0[:, 2],
    triangles=faces_0,
    cmap="jet",
    linewidths=1,
)

ax.view_init(elev=-160.0, azim=45)  # type: ignore

ax.set_box_aspect((1, 1, 0.9))  # type: ignore
plt.show()

# %% [markdown]
# #### Plot the  Icosahedron mesh with [Poly3DCollection](https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.art3d.Poly3DCollection.html).
# %%
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

data = np.array(
    [
        [vertices_0[face[0]], vertices_0[face[1]], vertices_0[face[2]]]
        for face in faces_0
    ]
)
# print(data.shape)

ax.add_collection3d(
    Poly3DCollection(
        data,
        facecolors=[
            np.random.rand(
                3,
            )
            for _ in data
        ],
        linewidths=1,
    )
)

ax.view_init(elev=-160.0, azim=45)
ax.set_box_aspect((1, 1, 0.9))
plt.show()

# %% [markdown]
# #### Save the Icosahedron mesh vertices into a .gmt file.

# You may open the file icosahedron_mesh.gmt in [GPlates](https://www.gplates.org/).
# %%
data_dir = "gplately-example-data"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
seen = set()
with open(f"{data_dir}/icosahedron_mesh.gmt", "w+") as f:
    for v in vertices_0:
        lon, lat = xyz2lonlat(v[0], v[1], v[2])
        line = f"{lon:0.2f} {lat:0.2f}\n"
        if line in seen:
            continue
        f.write(line)
        seen.add(line)
print(
    f"Now you can open the file {data_dir}/icosahedron_mesh.gmt in GPlates to see the Icosahedron mesh vertices."
)
print(
    "If you don't have GPlates installed yet, you may download it at https://www.gplates.org/."
)
