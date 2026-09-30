import os, io
import zipfile

import requests
from pathlib import Path

import pygplates
from gplately import IsochronSeafloorGrid, OutputScalarType, Raster
import gplately
from gplately.auxiliary import get_gplot

data_dir = Path("./test-isochron-seafloor-grids")
age_grid_input_dir = data_dir / "AgeGridInput"
age_grid_output_dir = data_dir / "output"
age_grid_output_dir.mkdir(parents=True, exist_ok=True)
if not os.path.isdir(age_grid_input_dir):
    # download the isochron file if not exists
    r = requests.get(
        "https://repo.gplates.org/webdav/gplately/AgeGridInput.zip",
        allow_redirects=True,
    )

    if r.status_code in [200]:
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall(data_dir)

gplot = get_gplot("zahirovic2022", model_repo_dir="plate-model-repo")

times = [0, 50, 100, 200]
IsochronSeafloorGrid(
    plate_reconstruction=gplot.plate_reconstruction,
    time_steps=times,
    ridges=pygplates.FeatureCollection(
        f"{age_grid_input_dir}/Global_EarthByte_GeeK07_Ridges.gpml"
    ),
    isochrons=pygplates.FeatureCollection(
        f"{age_grid_input_dir}/Global_EarthByte_GeeK07_Isochrons.gpml"
    ),
    iso_cob=pygplates.FeatureCollection(
        f"{age_grid_input_dir}/Global_EarthByte_GeeK07_IsoCOB.gpml"
    ),
    continental_polygons=gplately.load_feature_collection(
        gplot.plate_reconstruction.plate_model.get_layer("ContinentalPolygons")
    ),
    # continental_polygons=pygplates.FeatureCollection(
    #    f"{age_grid_input_dir}/Global_EarthByte_GeeK07_COB_Terranes_ContinentsOnly.gpml"
    # ),
    interval_spacing_degrees=0.5,
    grid_output_dir=age_grid_output_dir,
).generate(set([OutputScalarType.AGE, OutputScalarType.SPREADING_RATE]))

image_output_dir = age_grid_output_dir / "images"
image_output_dir.mkdir(parents=True, exist_ok=True)

for time in times:
    Raster(
        data=age_grid_output_dir
        / "isochron_seafloor_age"
        / "Masked"
        / f"seafloor_age_grid_masked_{time}Ma.nc"
    ).save_as_image(
        image_output_dir / f"seafloor_age_{time}Ma.png",
        cmap=gplately.plot.get_age_grid_cmap(),
        vmin=0,
        vmax=200,
        colorbar_label="Seafloor Age (Ma)",
    )
