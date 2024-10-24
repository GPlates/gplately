import cartopy.crs as ccrs
import cartopy.mpl.gridliner as grd
import geopandas as gpd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from plate_model_manager import PlateModelManager
from shapely.geometry import LineString, Point

import gplately

pm_manager = PlateModelManager()
muller2019_model = pm_manager.get_model("Muller2019", data_dir="plate-model-repo")
rotation_model = muller2019_model.get_rotation_model()
topology_features = muller2019_model.get_topologies()
static_polygons = muller2019_model.get_static_polygons()
coastlines = muller2019_model.get_coastlines()
continents = muller2019_model.get_continental_polygons()
COBs = muller2019_model.get_COBs()

## loading reconstruction
model = gplately.PlateReconstruction(rotation_model, topology_features, static_polygons)
time = 100  # Ma
subduction_data = model.tessellate_subduction_zones(time, return_geodataframe=True)


profiles = []
trench_lats = []
trench_lons = []
n_steps = 20
for i in range(0, len(subduction_data) - 1):
    dlon1 = n_steps * np.sin(
        np.radians(subduction_data.iloc[i]["trench normal angle (degrees)"])
    )
    dlat1 = n_steps * np.cos(
        np.radians(subduction_data.iloc[i]["trench normal angle (degrees)"])
    )
    x1 = subduction_data.iloc[i]["geometry"].x
    y1 = subduction_data.iloc[i]["geometry"].y
    ilon1 = x1 + dlon1
    ilat1 = y1 + dlat1
    start_point = Point(x1 + 0.1 * dlon1, y1 + 0.1 * dlat1)
    end_point = Point(ilon1, ilat1)
    profile = LineString([start_point, end_point])
    profiles.append(profile)

central_longitude = 180
Profile = gpd.GeoDataFrame(geometry=profiles)
mollweide_proj = f"+proj=moll +lon_0={central_longitude} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

Profile = Profile.set_crs(mollweide_proj)


cmap = "tab20c"
fig, axes = plt.subplots(
    2,
    1,
    figsize=(16, 12),
    dpi=100,
    subplot_kw={"projection": ccrs.Mollweide(central_longitude=central_longitude)},
)

# First subplot (left side)
ax1 = axes[0]
ax1.gridlines(
    color="0.7",
    linestyle="--",
    xlocs=np.arange(-180, 180, 15),
    ylocs=np.arange(-90, 90, 15),
)

gplot1 = gplately.PlotTopologies(
    model, coastlines=coastlines, continents=continents, COBs=COBs, time=time
)
# ax1.set_title('Subduction zones and mid-ocean ridges reconstructed to %i Ma' % (time))

# Plot shapefile features, subduction zones, and MOR boundaries at the specified time
gplot1.plot_continent_ocean_boundaries(ax1, color="b", alpha=0.05)
gplot1.plot_continents(ax1, facecolor="palegoldenrod", alpha=0.2)
gplot1.plot_coastlines(ax1, color="DarkKhaki")
# gplot1.plot_ridges_and_transforms(ax1, color="red")
gplot1.plot_trenches(ax1, color="k", alpha=0.3)
gplot1.plot_subduction_teeth(ax1, color="k", alpha=0.7)

# Add subduction data plot (with a colorbar)
cbar1 = subduction_data.plot(
    ax=ax1,
    column="trench normal angle (degrees)",
    transform=ccrs.PlateCarree(),
    alpha=0.2,
    cmap=cmap,
    legend=False,
)

sm = cm.ScalarMappable(cmap=cmap)
sm.set_array(subduction_data["trench normal angle (degrees)"])
sm.set_clim(0, 360)

# Add colorbar using the same Axes object used for plotting
colorbar = plt.colorbar(
    sm,
    ax=ax1,
    orientation="vertical",
    shrink=0.6,
    label="trench normal angle (degrees)",
)

# Second subplot (right side)
ax2 = axes[1]
ax2.gridlines(
    color="0.7",
    linestyle="--",
    xlocs=np.arange(-180, 180, 15),
    ylocs=np.arange(-90, 90, 15),
)

gplot2 = gplately.PlotTopologies(
    model, coastlines=coastlines, continents=continents, COBs=COBs, time=time
)
ax2.set_title("Profile View")

# Plot shapefile features
gplot2.plot_continent_ocean_boundaries(ax2, color="b", alpha=0.05)
gplot2.plot_continents(ax2, facecolor="palegoldenrod", alpha=0.2)
gplot2.plot_coastlines(ax2, color="DarkKhaki")
# gplot2.plot_ridges_and_transforms(ax2, color="red")
gplot2.plot_trenches(ax2, color="k", alpha=0.2)
gplot2.plot_subduction_teeth(ax2, color="k")

# Add profile plot (with colorbar)
cbar2 = Profile.plot(ax=ax2, transform=ccrs.PlateCarree(), alpha=0.2)
# fig.colorbar(cbar2, ax=ax2, fraction=0.04, pad=0.04)

# Set the global extents for both subplots
ax1.set_global()
ax2.set_global()

# Adjust layout
plt.tight_layout()
plt.show()
