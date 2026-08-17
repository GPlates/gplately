# GPlately Example Notebooks

This directory contains focused examples for common GPlately workflows.
Each notebook demonstrates one task so you can quickly copy, adapt, and run it.

## Available Examples

### [hello_world.ipynb](hello_world.ipynb)
Minimal first example to verify your setup and run a simple reconstruction workflow.

### [introducing_plate_model_manager.ipynb](introducing_plate_model_manager.ipynb)
Overview of PlateModelManager and how to access model components.

### [plot_map_with_cartopy.ipynb](plot_map_with_cartopy.ipynb)
Map plotting with Cartopy using reconstructed features and grids.

### [plot_map_with_pygmt.ipynb](plot_map_with_pygmt.ipynb)
Map plotting with PyGMT using GPlately plotting helpers.

### [save_reconstructed_data.ipynb](save_reconstructed_data.ipynb)
How to save reconstructed outputs to files for later analysis.

### [icosahedron_mesh.ipynb](icosahedron_mesh.ipynb)
Create and work with an icosahedral mesh workflow.

### [reconstruct_files.ipynb](reconstruct_files.ipynb)
Reconstruct external files and export the results.

### [use_auxiliary_functions.ipynb](use_auxiliary_functions.ipynb)
Practical examples for helper utilities in gplately.auxiliary.

### [use_your_own_plate_model.ipynb](use_your_own_plate_model.ipynb)
Run reconstructions with your own local plate model data.

### [working_with_plate_model_manager.ipynb](working_with_plate_model_manager.ipynb)
Use PlateModelManager with GPlately.

## Editing Workflow

Notebooks in this folder are generated from paired `.py` source files when available.

1. Edit the `.py` file (preferred for version control and code review).
2. Regenerate notebooks:

```bash
./convert_py_to_ipynb.sh
```

The conversion script uses git status checks and skips unchanged `.py` files.

## Notes

- Some examples require external datasets and may download inputs on first run.
- PyGMT examples require a working GMT/PyGMT environment.