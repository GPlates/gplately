#!/usr/bin/env bash

notebooks=(
  01-GettingStarted
  02-PlateReconstructions
  03-WorkingWithPoints
  04-VelocityBasics
  05-WorkingWithFeatureGeometries
  06-Rasters
  07-WorkingWithPlateTectonicStats
  08-PredictingSlabFlux
  09-CreatingMotionPathsAndFlowlines
  10-SeafloorGrids
  11-AndesFluxes
  12-MutschlerWorldPorphyryCopperDepositsRegionalPlots
  13-ReconstructingZirconData
  14-RuleBasedGPMLProcessingPipeline
  15-ConvertGridReferenceFrame
  Examples/plot_map_with_pygmt
  Examples/hello_world
  Examples/reconstruct_files
  Examples/icosahedron_mesh
  Examples/save_reconstructed_data
  Examples/introducing_plate_model_manager
  Examples/use_auxiliary_functions
  Examples/plot_map_with_cartopy
  Examples/use_your_own_plate_model
)

./scripts/run_notebooks.sh "${notebooks[@]}"