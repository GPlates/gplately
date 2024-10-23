#!/usr/bin/env python3

# this debugging script was created from NW's Jupyter notebook. Thank NW for providing us the code.

import gplately
from pathlib import Path
import os
import requests
from plate_model_manager import PlateModelManager, PlateModel
import xarray as xr

workdir = "debug-folder-github-issue-287"
Path(workdir).mkdir(parents=True, exist_ok=True)

miocene_file = os.path.join(workdir, "miocene_topo_pollard_antscape_dolan_0.5x0.5.nc")

if not os.path.isfile(miocene_file):
    r = requests.get(
        "https://repo.gplates.org/webdav/gplately-test-data/issue-287/miocene_topo_pollard_antscape_dolan_0.5x0.5.nc",
        allow_redirects=True,
    )
    open(miocene_file, "wb").write(r.content)

model_name = "merdith2021"
try:
    plate_model = PlateModelManager().get_model(model_name, data_dir=workdir)
except:
    plate_model = PlateModel(model_name, data_dir=workdir, readonly=True)

if not plate_model:
    raise Exception(f"Unable to get model({model_name})")

rotation_files = plate_model.get_rotation_model()
topology_files = plate_model.get_topologies()
continent_files = plate_model.get_layer("ContinentalPolygons")

r_model = gplately.PlateReconstruction(
    rotation_model=rotation_files,
    topology_features=topology_files,
    static_polygons=plate_model.get_layer("StaticPolygons"),
)

raster = gplately.Raster(
    data=miocene_file, time=15, plate_reconstruction=r_model, realign=True
)
reconstructed_raster = raster.reconstruct(
    0, threads=4, partitioning_features=plate_model.get_layer("StaticPolygons")
)
outfile = os.path.join(workdir, "reconstructed_grid.nc")
gplately.grids.write_netcdf_grid(
    filename=outfile,
    grid=reconstructed_raster._data,
    extent=[
        reconstructed_raster._lons.min(),
        reconstructed_raster._lons.max(),
        reconstructed_raster._lats.min(),
        reconstructed_raster._lats.max(),
    ],
)
# mio_nc2 = xr.open_dataset(outfile)
# mio_nc2.z.plot()
