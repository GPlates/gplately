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
  Examples/01-HelloWorld
  Examples/02-PlateModelManager
  Examples/03-PlotWithCartopy
  Examples/04-PlotWithPyGMT
  Examples/05-ReconstructFiles
  Examples/06-LoadPlateModelFromFiles
  Examples/07-SaveReconstructedData
  Examples/08-UseAuxiliaryFunctions
  Examples/09-IcosahedronMesh
  Examples/10-ColorMapAndColorPaletteTable
)
export GPLATELY_NOTEBOOK_QUICK_RUN=true
./scripts/run_notebooks.sh "${notebooks[@]}"