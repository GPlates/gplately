#!/bin/bash
set -e # Exit immediately if any command exits with a non-zero status

BASEDIR=$(dirname "$0")

$BASEDIR/test_age_grid.py

$BASEDIR/test_anchor_plate_id.py save

$BASEDIR/test_continent_contouring.py

$BASEDIR/test_crustal_production.py

$BASEDIR/test_data_server.py

$BASEDIR/test_discretize_polyline.py 

$BASEDIR/test_feature_filter.sh

$BASEDIR/test_issue_250.py save

$BASEDIR/test_motion_path.py save

$BASEDIR/test_plate_model.py

$BASEDIR/test_plot_models.py save

$BASEDIR/test_plot_with_raster.py save

$BASEDIR/test_plot.py save 

$BASEDIR//test_pygmt_plot.py

$BASEDIR/test_raster_reconstruction.py 701 save

$BASEDIR/test_raster_reconstruction.py save

$BASEDIR/test_raster.py save

$BASEDIR/test_reconstruct_points.py save

$BASEDIR/test_reset_feature_type.sh

$BASEDIR/test_ridges_and_transforms.py save

$BASEDIR/test_subduction_teeth.py save

$BASEDIR/test_tessellate.py save

