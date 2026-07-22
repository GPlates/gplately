# %% [markdown]

# This is a test file for the `raster.py` module. It contains unit tests for various functions
# and methods defined in the `raster.py` module.


# Put the test functions in the `test_raster_impl.py` file. Use this file to generate a Jupyter notebook for interactive testing and visualization.

# `jupytext --to notebook test_raster.py -o test_raster.ipynb`

# %%
import time

from test_raster_impl import (
    test_raster_reconstruct_lon_0_360,
    test_raster_query,
    test_raster_clip_by_extent,
    test_upper_origin_raster_reconstruction,
    test_raster_reconstruction,
)

all_flag = False
test_raster_reconstruct_lon_0_360_flag = False
test_raster_query_flag = False
test_raster_clip_by_extent_flag = False
test_upper_origin_raster_reconstruction_flag = False
test_raster_reconstruction_flag = True

# %%
if all_flag or test_raster_reconstruct_lon_0_360_flag:
    test_raster_reconstruct_lon_0_360()
# %%
if all_flag or test_raster_query_flag:
    test_raster_query()
# %%
if all_flag or test_raster_clip_by_extent_flag:
    test_raster_clip_by_extent()
# %%
if all_flag or test_upper_origin_raster_reconstruction_flag:
    test_upper_origin_raster_reconstruction()
# %%
if all_flag or test_raster_reconstruction_flag:
    start = time.perf_counter()
    test_raster_reconstruction()
    print(f"Elapsed: {time.perf_counter() - start:.4f} seconds")
# %%
if True:
    start = time.perf_counter()
    test_raster_reconstruction(use_spatial_tree=True)
    print(f"Elapsed: {time.perf_counter() - start:.4f} seconds")
# %%
if True:
    start = time.perf_counter()
    test_raster_reconstruction(use_old_implementation=True)
    print(f"Elapsed: {time.perf_counter() - start:.4f} seconds")
