#!/usr/bin/env python
# coding: utf-8

# **Step 0: Download the shapefiles for this example.**
# 
# Don't worry about the `FileDownloader`. If you'd like, you may download the zip file manually.

# In[ ]:


import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pygplates
from plate_model_manager.utils import download

import gplately
from gplately.auxiliary import get_plate_reconstruction

time = 100
data_dir = "gplately-example-data"

downloader = download.FileDownloader(
    "https://repo.gplates.org/webdav/gplately-test-data/Global_Paleogeography_Cao_etal.zip",
    f"{data_dir}/.metadata.json",
    f"{data_dir}",
)
# only re-download when necessary
if downloader.check_if_file_need_update():
    downloader.download_file_and_update_metadata()
else:
    print(f"The local files are still good. No need to download again!")

# **Step 1: Load files**

# In[ ]:


landmass_fc = pygplates.FeatureCollection(f"{data_dir}/Global_Cao_etal/lm_402_2.shp")
mountain_fc = pygplates.FeatureCollection(f"{data_dir}/Global_Cao_etal/m_402_2.shp")
icesheet_fc = pygplates.FeatureCollection(f"{data_dir}/Global_Cao_etal/i_402_2.shp")
sallow_marine_fc = pygplates.FeatureCollection(f"{data_dir}/Global_Cao_etal/sm_402_2.shp")

# **Step 2: Reconstruct files**

# In[ ]:


# use the auxiliary function to create a PlateReconstruction object
plate_reconstruction_obj = get_plate_reconstruction("cao2024")

reconstructed_landmass = plate_reconstruction_obj.reconstruct(landmass_fc, time)
reconstructed_mountain = plate_reconstruction_obj.reconstruct(mountain_fc, time)
reconstructed_icesheet = plate_reconstruction_obj.reconstruct(icesheet_fc, time)
reconstructed_sallow_marine = plate_reconstruction_obj.reconstruct(
    sallow_marine_fc, time
)

# **Step 3: Plot the map**

# In[ ]:


ax = plt.figure(figsize=(8, 6)).add_subplot(
    111, projection=ccrs.Robinson(central_longitude=0)
)

ax.add_geometries(
    gplately.plot.shapelify_features(reconstructed_sallow_marine),
    crs=ccrs.PlateCarree(),
    facecolor="#add8e6",
    edgecolor="none",
)
ax.add_geometries(
    gplately.plot.shapelify_features(reconstructed_landmass),
    crs=ccrs.PlateCarree(),
    facecolor="#fedf00",
    edgecolor="none",
)
ax.add_geometries(
    gplately.plot.shapelify_features(reconstructed_mountain),
    crs=ccrs.PlateCarree(),
    facecolor="#ffa500",
    edgecolor="none",
)
ax.add_geometries(
    gplately.plot.shapelify_features(reconstructed_icesheet),
    crs=ccrs.PlateCarree(),
    facecolor="#e9e9ff",
    edgecolor="none",
)

ax.set_global()
plt.title(f"Paleogeography at {time} Ma")
plt.show()

