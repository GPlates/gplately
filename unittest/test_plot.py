import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys

sys.path.insert(0, "../")
from gplately import plate_model, PlateReconstruction, PlotTopologies

# test the plot function with the new PlateModel class

age = 55
model = plate_model.PlateModel("Muller2019")

rotation_model = model.get_rotation_model()

test_model = PlateReconstruction(
    rotation_model,
    topology_features=model.get_layer("Topologies"),
    static_polygons=model.get_layer("StaticPolygons"),
)
gplot = PlotTopologies(
    test_model,
    coastlines=model.get_layer("Coastlines"),
    COBs=model.get_layer("COBs"),
    time=age,
)

fig = plt.figure(figsize=(10, 10), dpi=100)
ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=0))

gplot.plot_continent_ocean_boundaries(ax, color="cornflowerblue")
gplot.plot_coastlines(ax, color="black")
gplot.plot_ridges_and_transforms(ax, color="red")
gplot.plot_trenches(ax, color="orange")
gplot.plot_subduction_teeth(ax, color="orange")
ax.set_global()

ids = set([f.get_reconstruction_plate_id() for f in gplot.topologies])
for id in ids:
    gplot.plot_plate_id(ax, id, facecolor="None", edgecolor="lightgreen")
plt.show()
