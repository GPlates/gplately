#!/bin/bash

set -euo pipefail

export GPLATELY_DISABLE_DEV_WARNING=true

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

$BASEDIR//test_pygmt_plot.py save

$BASEDIR/test_reconstruct_points.py save

$BASEDIR/test_reset_feature_type.sh

$BASEDIR/test_ridges_and_transforms.py save

$BASEDIR/test_subduction_teeth.py save

$BASEDIR/test_tessellate.py save

jupyter nbconvert \
    --to notebook \
    --execute \
    --ExecutePreprocessor.force_raise_errors=True \
    --output test_raster_tmp.ipynb \
    --output-dir $BASEDIR \
    $BASEDIR/test_raster.ipynb

jupyter nbconvert \
    --to notebook \
    --execute \
    --ExecutePreprocessor.force_raise_errors=True \
    --output test_pygmt_tmp.ipynb \
    --output-dir $BASEDIR \
    $BASEDIR/test_pygmt.ipynb

  $BASEDIR/test_feature_filter.sh

echo "All tests passed!"
