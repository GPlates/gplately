#!/usr/bin/env python3
import math

import pygplates
from common import MODEL_REPO_DIR
from plate_model_manager import PlateModel, PlateModelManager

import gplately
from gplately import PlateReconstruction, PlotTopologies
from gplately.ptt import continent_contours

print(gplately.__file__)

max_distance_of_subduction_zone_from_active_margin_kms = 500
max_distance_of_subduction_zone_from_active_margin_radians = (
    max_distance_of_subduction_zone_from_active_margin_kms
    / pygplates.Earth.mean_radius_in_kms
)
continent_contouring_point_spacing_degrees = 0.25
continent_contouring_area_threshold_square_kms = 0
continent_contouring_area_threshold_steradians = (
    continent_contouring_area_threshold_square_kms
    / (pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms)
)
continent_exclusion_area_threshold_square_kms = 800000
continent_exclusion_area_threshold_steradians = (
    continent_exclusion_area_threshold_square_kms
    / (pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms)
)
continent_separation_distance_threshold_radians = (
    continent_contours.DEFAULT_CONTINENT_SEPARATION_DISTANCE_THRESHOLD_RADIANS
)


def continent_contouring_buffer_and_gap_distance_radians(time, contoured_continent):
    # One distance for time interval [1000, 300] and another for time interval [250, 0].
    # And linearly interpolate between them over the time interval [300, 250].
    pre_pangea_distance_radians = math.radians(2.5)  # convert degrees to radians
    post_pangea_distance_radians = math.radians(0.0)  # convert degrees to radians
    if time > 300:
        buffer_and_gap_distance_radians = pre_pangea_distance_radians
    elif time < 250:
        buffer_and_gap_distance_radians = post_pangea_distance_radians
    else:
        # Linearly interpolate between 250 and 300 Ma.
        interp = float(time - 250) / (300 - 250)
        buffer_and_gap_distance_radians = (
            interp * pre_pangea_distance_radians
            + (1 - interp) * post_pangea_distance_radians
        )

    # Area of the contoured continent.
    area_steradians = contoured_continent.get_area()

    # Linearly reduce the buffer/gap distance for contoured continents with area smaller than 1 million km^2.
    area_threshold_square_kms = 500000
    area_threshold_steradians = area_threshold_square_kms / (
        pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms
    )
    if area_steradians < area_threshold_steradians:
        buffer_and_gap_distance_radians *= area_steradians / area_threshold_steradians

    return buffer_and_gap_distance_radians


if __name__ == "__main__":
    MODEL_NAME = "merdith2021"
    try:
        model = PlateModelManager().get_model(MODEL_NAME, data_dir=MODEL_REPO_DIR)
    except:
        model = PlateModel(MODEL_NAME, data_dir=MODEL_REPO_DIR, readonly=True)

    rotation_model = pygplates.RotationModel(model.get_rotation_model())
    continent_files = model.get_continental_polygons() + model.get_layer("Cratons")
    continent_features = [pygplates.FeatureCollection(f) for f in continent_files]
    print(continent_features)
    continent_contouring = continent_contours.ContinentContouring(
        rotation_model,
        continent_features,
        continent_contouring_point_spacing_degrees,
        continent_contouring_area_threshold_steradians,
        continent_contouring_buffer_and_gap_distance_radians,
        continent_exclusion_area_threshold_steradians,
        continent_separation_distance_threshold_radians,
    )
    time = 0
    continent_mask, contoured_continents = (
        continent_contouring.get_continent_mask_and_contoured_continents(time)
    )
    continent_mask_filename = "continent_mask_{}.nc".format(time)
    # Note that we need to convert the boolean mask grid to a non-boolean number type for NetCDF (and it seems floating-point for gplately).
    continent_mask_grid = continent_mask.astype("float")
    gplately.grids.write_netcdf_grid(continent_mask_filename, continent_mask_grid)
