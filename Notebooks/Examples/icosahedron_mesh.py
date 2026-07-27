#!/usr/bin/env python
# coding: utf-8

# **Generate Icosahedron mesh and plot with [plot_trisurf](https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.plot_trisurf.html).**

# In[ ]:


import matplotlib.pyplot as plt
import numpy as np

import gplately
from gplately.lib.icosahedron import get_mesh, xyz2lonlat

mesh_resolution=3
vertices_0, faces_0 = get_mesh(mesh_resolution)
#print(vertices_0.shape, faces_0.shape)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(8,8))

ax.plot_trisurf(vertices_0[:,0], vertices_0[:,1], vertices_0[:,2], triangles=faces_0, cmap='jet', linewidths=1)

ax.view_init(elev=-160., azim=45)

ax.set_box_aspect((1, 1, 0.9)) 
plt.show()

# **Plot the  Icosahedron mesh with [Poly3DCollection](https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.art3d.Poly3DCollection.html).**

# In[ ]:


from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')

data = np.array([[vertices_0[face[0]], vertices_0[face[1]], vertices_0[face[2]]] for face in faces_0 ])
#print(data.shape)

ax.add_collection3d(Poly3DCollection(data, facecolors=[np.random.rand(3,) for _ in data], linewidths=1))

ax.view_init(elev=-160., azim=45)
ax.set_box_aspect((1, 1, 0.9)) 
plt.show()

# **Save the Icosahedron mesh vertices into a .gmt file, and then you may open the file in [GPlates](https://www.gplates.org/).**

# In[ ]:


seen = set()
with open("icosahedron_mesh.gmt", "w+") as f:
    for v in vertices_0:
        lon, lat = xyz2lonlat(v[0], v[1], v[2])
        line = f"{lon:0.2f} {lat:0.2f}\n"
        if line in seen:
            continue
        f.write(line)
        seen.add(line)
print("Now you can open the file icosahedron_mesh.gmt with GPlates to see the Icosahedron mesh vertices.")
print("If you don't have GPlates installed yet, you may download it at https://www.earthbyte.org/download-gplates-2-5/.")
