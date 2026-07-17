# %% [markdown]

# `jupytext --to notebook test_raster.py -o test_raster.ipynb`

# %%
from test_raster_impl import (
    test_raster_longitude_conversion,
    test_raster_query,
    test_raster_clip_by_extent,
    test_upper_origin_raster_reconstruction,
)

all_flag = False
test_raster_longitude_conversion_flag = False
test_raster_query_flag = False
test_raster_clip_by_extent_flag = False
test_upper_origin_raster_reconstruction_flag = True

# %%
if all_flag or test_raster_longitude_conversion_flag:
    test_raster_longitude_conversion()
# %%
if all_flag or test_raster_query_flag:
    test_raster_query()
# %%
if all_flag or test_raster_clip_by_extent_flag:
    test_raster_clip_by_extent()
# %%
if all_flag or test_upper_origin_raster_reconstruction_flag:
    test_upper_origin_raster_reconstruction()
