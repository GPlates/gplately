#!/usr/bin/env python3
import os
import sys

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pygplates
from common import MODEL_REPO_DIR, save_fig
from plate_model_manager import PlateModelManager

import gplately

print(gplately.__file__)


def main(show=True):
    # Call GPlately's PlateModelManager object and request data from the Müller et al. 2019 study
    pm_manager = PlateModelManager()
    muller2019_model = pm_manager.get_model("Muller2019", data_dir=MODEL_REPO_DIR)
    assert muller2019_model
    rotation_model = muller2019_model.get_rotation_model()
    topology_features = muller2019_model.get_topologies()
    static_polygons = muller2019_model.get_static_polygons()

    model = gplately.PlateReconstruction(
        rotation_model, topology_features, static_polygons
    )

    time = 50  # Ma

    use_ptt = False
    tessellation_threshold_radians = np.radians(1)
    spreading_feature_types = pygplates.FeatureType.gpml_mid_ocean_ridge
    transform_segment_deviation_in_radians = (
        gplately.ptt.separate_ridge_transform_segments.DEFAULT_TRANSFORM_SEGMENT_DEVIATION_RADIANS
    )
    include_network_boundaries = False
    convergence_threshold_in_cm_per_yr = None
    divergence_threshold_in_cm_per_yr = None
    output_obliquity_and_normal_and_left_right_plates = True
    subduction_output_distance_to_nearest_edge_of_trench = True
    subduction_output_distance_to_start_edge_of_trench = True
    subduction_output_convergence_velocity_components = True
    subduction_output_trench_absolute_velocity_components = True
    subduction_output_subducting_absolute_velocity = True
    subduction_output_subducting_absolute_velocity_components = True
    subduction_output_trench_normal = True

    # Subduction zones.
    subduction_data = model.tessellate_subduction_zones(
        time,
        tessellation_threshold_radians=tessellation_threshold_radians,
        ignore_warnings=True,
        use_ptt=use_ptt,
        include_network_boundaries=include_network_boundaries,
        convergence_threshold_in_cm_per_yr=convergence_threshold_in_cm_per_yr,
        return_geodataframe=True,
        output_distance_to_nearest_edge_of_trench=subduction_output_distance_to_nearest_edge_of_trench,
        output_distance_to_start_edge_of_trench=subduction_output_distance_to_start_edge_of_trench,
        output_convergence_velocity_components=subduction_output_convergence_velocity_components,
        output_trench_absolute_velocity_components=subduction_output_trench_absolute_velocity_components,
        output_subducting_absolute_velocity=subduction_output_subducting_absolute_velocity,
        output_subducting_absolute_velocity_components=subduction_output_subducting_absolute_velocity_components,
        output_trench_normal=subduction_output_trench_normal,
    )
    total_subduction_zone_length_in_kms = model.total_subduction_zone_length(
        time,
        use_ptt=use_ptt,
        ignore_warnings=True,
        include_network_boundaries=include_network_boundaries,
        convergence_threshold_in_cm_per_yr=convergence_threshold_in_cm_per_yr,
    )
    print("total subduction zone length (kms):", total_subduction_zone_length_in_kms)

    # Ridges.
    ridge_data = model.tessellate_mid_ocean_ridges(
        time,
        tessellation_threshold_radians=tessellation_threshold_radians,
        ignore_warnings=True,
        use_ptt=use_ptt,
        spreading_feature_types=spreading_feature_types,
        transform_segment_deviation_in_radians=transform_segment_deviation_in_radians,
        include_network_boundaries=include_network_boundaries,
        divergence_threshold_in_cm_per_yr=divergence_threshold_in_cm_per_yr,
        output_obliquity_and_normal_and_left_right_plates=output_obliquity_and_normal_and_left_right_plates,
        return_geodataframe=True,
    )
    total_ridge_length_in_kms = model.total_ridge_length(
        time,
        ignore_warnings=True,
        use_ptt=use_ptt,
        spreading_feature_types=spreading_feature_types,
        transform_segment_deviation_in_radians=transform_segment_deviation_in_radians,
        include_network_boundaries=include_network_boundaries,
        divergence_threshold_in_cm_per_yr=divergence_threshold_in_cm_per_yr,
    )
    print("total ridge length (kms):", total_ridge_length_in_kms)

    # Plate velocity stats.
    diverging_data, converging_data = model.divergent_convergent_plate_boundaries(
        time,
        uniform_point_spacing_radians=0.001,
        divergence_velocity_threshold=0.2,  # cm/yr
        convergence_velocity_threshold=0.2,  # cm/yr
        velocity_units=pygplates.VelocityUnits.cms_per_yr,
    )
    total_crustal_production_rate, total_crustal_destruction_rate = (
        model.crustal_production_destruction_rate(
            time,
            uniform_point_spacing_radians=0.001,
            divergence_velocity_threshold_in_cms_per_yr=0.2,  # cm/yr
            convergence_velocity_threshold_in_cms_per_yr=0.2,  # cm/yr
        )
    )
    print(
        "total crustal production/destruction rate (km^2/yr):",
        total_crustal_production_rate,
        total_crustal_destruction_rate,
    )

    # Non-geometry columns of subduction data.
    subduction_data_columns = [
        column
        for column in subduction_data.columns  # type: ignore
        if column != subduction_data.geometry.name  # type: ignore
    ]

    # Non-geometry columns of ridge data.
    ridge_data_columns = [
        column for column in ridge_data.columns if column != ridge_data.geometry.name  # type: ignore
    ]

    # Subduction and ridge subplots followed by a single plate boundary statistics subplot.
    num_subplots = len(subduction_data_columns) + len(ridge_data_columns) + 1

    nun_subplot_columns = 4
    num_subplot_rows = num_subplots // nun_subplot_columns + 1

    # First plots are subduction data followed by ridge data followed by boundary statistics.
    #
    # Note: You'll need to zoom to 100% if viewing figure saved to a file.
    fig = plt.figure(
        figsize=(6 * nun_subplot_columns, 2 * num_subplot_rows), layout="constrained"
    )
    fig.suptitle(f"Time {time} Ma", fontsize="x-large")
    axes = fig.subplots(
        num_subplot_rows,
        nun_subplot_columns,
        subplot_kw=dict(projection=ccrs.Mollweide()),
    )

    subplot_index = 0

    # Subduction data subplots.
    for subduction_data_column in subduction_data_columns:
        ax = axes.flatten()[subplot_index]
        ax.set_global()
        im_subduction = ax.scatter(
            subduction_data.geometry.x,  # longitude # type: ignore
            subduction_data.geometry.y,  # latitude # type: ignore
            c=subduction_data.loc[:, subduction_data_column],  # type: ignore
            s=2,
            transform=ccrs.PlateCarree(),
        )
        ax.set_title(f"Subduction: {subduction_data_column}")
        fig.colorbar(im_subduction, ax=ax)

        subplot_index += 1

    # Ridge data subplots.
    for ridge_data_column in ridge_data_columns:
        ax = axes.flatten()[subplot_index]
        ax.set_global()
        im_ridge = ax.scatter(
            ridge_data.geometry.x,  # longitude # type: ignore
            ridge_data.geometry.y,  # latitude # type: ignore
            c=ridge_data.loc[:, ridge_data_column],  # type: ignore
            s=2,
            transform=ccrs.PlateCarree(),
        )
        ax.set_title(f"Ridge: {ridge_data_column}")
        fig.colorbar(im_ridge, ax=ax)

        subplot_index += 1

    # Plate boundary statistics subplot.
    converging_lat_lon_locations = np.fromiter(
        (data.boundary_point.to_lat_lon() for data in converging_data),
        dtype=np.dtype((float, 2)),
    )
    diverging_lat_lon_locations = np.fromiter(
        (data.boundary_point.to_lat_lon() for data in diverging_data),
        dtype=np.dtype((float, 2)),
    )
    ax = axes.flatten()[subplot_index]
    ax.set_global()
    ax.scatter(
        converging_lat_lon_locations[:, 1],  # longitude
        converging_lat_lon_locations[:, 0],  # latitude
        color="blue",
        s=2,
        transform=ccrs.PlateCarree(),
    )
    ax.scatter(
        diverging_lat_lon_locations[:, 1],  # longitude
        diverging_lat_lon_locations[:, 0],  # latitude
        color="red",
        s=2,
        transform=ccrs.PlateCarree(),
    )
    ax.set_title("Plate Boundary Stats: converging(blue)/diverging(red)")

    subplot_index += 1

    # Remove unused subplots.
    while subplot_index < len(axes.flatten()):
        axes.flatten()[subplot_index].remove()
        subplot_index += 1

    if show:
        # LOOK HERE! 👀👀 👇👇
        # If the figure did not show up, you need to set your matplotlib plotting backend properly.
        # On Windows, you may install PyQt and do
        # import matplotlib
        # matplotlib.use('QtAgg')

        # if you are interested in finding what backends available on your computer and what is your current backend, do the following
        # import matplotlib.rcsetup as rcsetup
        # print(rcsetup.all_backends) # get all available backends
        # import matplotlib
        # matplotlib.get_backend() # your current backend
        #
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
