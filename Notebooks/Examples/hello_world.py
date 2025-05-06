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
plt.title(f"{int(gplot.time)} Ma")

# save the map to a .png file
plt.gcf().savefig("gplately-hello-world.png")
