# %% [markdown]

## Search and Filter Feature Collection


# This notebook demonstrates how to use GPlately to query and filter feature collections based on criteria such as feature name, plate ID, birth age, disappearance age, and region of interest. 

# The search and filter functionality is implemented in the `gplately.utils.feature_filter` module. The main function is `filter_feature_collection`, which takes a feature collection and a list of filters, and returns a new feature collection that contains only the features that satisfy all the filters. Each filter is a subclass of the `FeatureFilter` class, which defines a `should_keep` method that takes a feature and returns True if the feature should be kept, and False otherwise.


# %%
from os import makedirs
from os.path import exists
from urllib.request import urlretrieve

import pygplates # type: ignore

from gplately.utils.feature_filter import (
    FeatureNameFilter,
    PlateIDFilter,
    BirthAgeFilter,
    FeatureTypeFilter,
    filter_feature_collection,
)

DATA_DIR="./feature-search-notebook-data"
makedirs(DATA_DIR, exist_ok=True)

# Download the input feature collection file if it does not exist.
# The input file we use in this notbook is a global coastline feature collection at present day. 
# You can also use your own feature collection as input and apply the same search and filter methods as shown in this notebook.
IN_FILE=f"{DATA_DIR}/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"
FILE_URL="https://repo.gplates.org/webdav/mchin/data/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"

if not exists(IN_FILE):
    print(f"Downloading test file to {IN_FILE}...")
    urlretrieve(FILE_URL, IN_FILE)

input_feature_collection = pygplates.FeatureCollection(IN_FILE) # type: ignore

# %% [markdown]
#### Search and filter feature collection by feature names

# Search feature collection by feature names (case sensitive, partial match)

# Features whose name contains "Australia", "New Zealand" or "Tasmania" will be kept.

# %%
filters = []
filters.append(
    FeatureNameFilter(
        ["Australia", "New Zealand","Tasmania"],
        exact_match=False,
        case_sensitive=True,
    )
)
features = filter_feature_collection(
    input_feature_collection,
    filters,
)
if len(features) == 0:
    print(
        "Warning: No features matched the search criteria. "
        "The output feature collection will be empty."
    )

features.write(f"{DATA_DIR}/Australia_and_New_Zealand_coastlines.gpmlz")

# %% [markdown]
# Filter feature collection by feature names (case sensitive, partial match)

# Features whose name contains "Australia", "New Zealand" or "Tasmania" will be discarded.

# %%
filters = []
filters.append(
    FeatureNameFilter(
        ["Australia", "new zealand","Tasmania"],
        exclude=True,
        exact_match=False,
        case_sensitive=False,
    )
)
features = filter_feature_collection(
    input_feature_collection,
    filters,
)
if len(features) == 0:
     print(
        "Warning: No features matched the search criteria. "
        "The output feature collection will be empty."
    )

features.write(f"{DATA_DIR}/coastlines_without_Australia_and_New_Zealand.gpmlz")

'''
filters.append(PlateIDFilter(args.pids))

filters.append(PlateIDFilter(args.exclude_pids, exclude=True))


filters.append(BirthAgeFilter(args.max_birth_age, keep_older=False))

filters.append(BirthAgeFilter(args.min_birth_age))


filters.append(FeatureTypeFilter(args.feature_type_re))

new_fc = filter_feature_collection(
    input_feature_collection,
    filters,
)

if len(new_fc) == 0:
    print(
        "No feature is left after the filtering. The output feature collection will be empty."
    )

new_fc.write(args.filter_output_file)
'''

