#!/bin/bash

BASEDIR=$(dirname "$0")

$BASEDIR/test_plate_model.py

$BASEDIR/test_subduction_teeth.py save

$BASEDIR/test_plot.py save 

$BASEDIR/test_anchor_plate_id.py save

$BASEDIR/test_plot_with_raster.py save

$BASEDIR/test_reconstruct_points.py save

$BASEDIR/test_tessellate.py save

$BASEDIR/test_motion_path.py save