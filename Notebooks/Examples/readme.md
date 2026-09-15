# GPlately Example Notebooks

This directory contains focused examples for common GPlately workflows.
Each notebook demonstrates one task so you can quickly copy, adapt, and run it.

## Available Examples

### [01-HelloWorld.ipynb](01-HelloWorld.ipynb)
The minimal working example — create a paleo-map with GPlately in a few lines.

### [02-PlateModelManager.ipynb](02-PlateModelManager.ipynb)
Use `PlateModelManager` to access plate model files.

### [03-PlotWithCartopy.ipynb](03-PlotWithCartopy.ipynb)
Build a more detailed paleo-map with Cartopy.

### [04-PlotWithPyGMT.ipynb](04-PlotWithPyGMT.ipynb)
Plot maps using GPlately's PyGMT integration.

### [05-ReconstructFiles.ipynb](05-ReconstructFiles.ipynb)
Reconstruct shapefiles with GPlately.

### [06-LoadPlateModelFromFiles.ipynb](06-LoadPlateModelFromFiles.ipynb)
Use your own local plate model and reconstruct points.

### [07-SaveReconstructedData.ipynb](07-SaveReconstructedData.ipynb)
Save reconstructed data to shapefiles and other formats.

### [08-UseAuxiliaryFunctions.ipynb](08-UseAuxiliaryFunctions.ipynb)
Quickly create `PlateReconstruction` and `PlotTopologies` objects with `gplately.auxiliary`.

### [09-IcosahedronMesh.ipynb](09-IcosahedronMesh.ipynb)
Generate and plot an icosahedron mesh.

### [10-ColorMapAndColorPaletteTable.ipynb](10-ColorMapAndColorPaletteTable.ipynb)
Work with Matplotlib colormaps and GMT Color Palette Tables (CPT).

### [PNG_reconstruction_copper_deposits.ipynb](PNG_reconstruction_copper_deposits.ipynb)
A larger, standalone worked example: reconstructing Papua New Guinea copper deposit locations
through time.

## Notes

- Some examples require external datasets and may download inputs on first run.
- PyGMT examples require a working GMT/PyGMT environment (see [`docker/env.yaml`](../../docker/env.yaml)).
