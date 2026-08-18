import os, io
import zipfile

import requests

import pygplates
from gplately.isochron_seafloor_grid import IsochronSeafloorGrid, OutputScalarType
from gplately.auxiliary import get_gplot

age_grid_input_dir = "AgeGridInput"
if not os.path.isdir(age_grid_input_dir):
    # download the isochron file if not exists
    r = requests.get(
        "https://repo.gplates.org/webdav/mchin/data/AgeGridInput.zip",
        allow_redirects=True,
    )

    if r.status_code in [200]:
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall(".")

gplot = get_gplot("zahirovic2022", model_repo_dir="plate-model-repo")

o = IsochronSeafloorGrid(
    plate_reconstruction=gplot.plate_reconstruction,
    time_steps=[0, 100],
    ridges=pygplates.FeatureCollection(
        f"{age_grid_input_dir}/Global_EarthByte_GeeK07_Ridges.gpml"
    ),
    isochrons=pygplates.FeatureCollection(
        f"{age_grid_input_dir}/Global_EarthByte_GeeK07_Isochrons.gpml"
    ),
    iso_cob=pygplates.FeatureCollection(
        f"{age_grid_input_dir}/Global_EarthByte_GeeK07_IsoCOB.gpml"
    ),
    # continental_polygons=gplately.utils.io_utils.load_feature_collection(
    #    gplot.plate_reconstruction.plate_model.get_layer("ContinentalPolygons")
    # ),
    continental_polygons=pygplates.FeatureCollection(
        f"{age_grid_input_dir}/Global_EarthByte_GeeK07_COB_Terranes_ContinentsOnly.gpml"
    ),
    interval_spacing_degrees=0.5,
)
# o._get_continent_mask(0)
o.generate(set([OutputScalarType.AGE, OutputScalarType.SPREADING_RATE]))
