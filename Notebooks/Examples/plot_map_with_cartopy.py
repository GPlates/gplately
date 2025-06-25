#!/usr/bin/env python3
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from gplately import PlateModelManager, Raster, auxiliary

# use the auxiliary function to create a PlotTopologies instance
gplot = auxiliary.get_gplot("Muller2019", time=100)  # 100Ma

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
gplot.plot_ridges(ax1, color="red")
gplot.plot_trenches(ax1, color="k")
gplot.plot_subduction_teeth(ax1, color="k")
im = gplot.plot_grid(ax1, agegrid.data, cmap="YlGnBu", vmin=0, vmax=200)
gplot.plot_plate_motion_vectors(
    ax1, spacingX=10, spacingY=10, normalise=True, zorder=10, alpha=0.5
)

fig.colorbar(im, orientation="horizontal", shrink=0.4, pad=0.05, label="Age (Ma)")
plt.title(f"{int(gplot.time)} Ma")

# save the map as a .png file
output_file = f"plot_map_with_cartopy.png"
plt.gcf().savefig(output_file, dpi=120, bbox_inches="tight")  # transparent=True)
print(f"Done! The {output_file} has been saved.")
plt.close(plt.gcf())
