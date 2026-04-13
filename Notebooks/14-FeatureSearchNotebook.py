# %% [markdown]

## Search and Filter Feature Collection

# This notebook demonstrates how to use GPlately to query and filter feature collections based on criteria such as feature name, plate ID, birth age, disappearance age, and region of interest.

# The search and filter functionality is implemented in the `gplately.utils.feature_filter` module. The main function is `filter_feature_collection`, which takes a feature collection and a list of filters, and returns a new feature collection that contains only the features that satisfy all the filters. Each filter is a subclass of the `FeatureFilter` class, which defines a `should_keep` method that takes a feature and returns True if the feature should be kept, and False otherwise.


# %% [markdown]

# 📒 This notebook is generated from 14-FeatureSearchNotebook.py using the command `jupytext --to notebook Notebooks/14-FeatureSearchNotebook.py -o Notebooks/14-FeatureSearch.ipynb`. If you need to commit changes to this notebook to the GPlately repository, make your edits in 14-FeatureSearchNotebook.py and then regenerate this Jupyter Notebook file. The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.

# %%
from os import makedirs
from os.path import exists
from urllib.request import urlretrieve

import pygplates  # type: ignore

from gplately.utils.feature_filter import (
    EndTimeFilter,
    FeatureNameFilter,
    PlateIDFilter,
    BirthAgeFilter,
    FeatureTypeFilter,
    filter_feature_collection,
)

DATA_DIR = "./feature-search-notebook-data"
makedirs(DATA_DIR, exist_ok=True)

# Download the input feature collection file if it does not exist.
# The input file we use in this notbook is a global coastline feature collection at present day.
# You can also use your own feature collection as input and apply the same search and filter methods as shown in this notebook.
IN_FILE = f"{DATA_DIR}/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"
FILE_URL = "https://repo.gplates.org/webdav/mchin/data/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"

if not exists(IN_FILE):
    print(f"Downloading test file to {IN_FILE}...")
    urlretrieve(FILE_URL, IN_FILE)

input_feature_collection = pygplates.FeatureCollection(IN_FILE)  # type: ignore

# %% [markdown]
#### Search and filter feature collection by feature names

# Search feature collection by feature names (case sensitive, partial match)

# Features whose name contains "Australia", "New Zealand" or "Tasmania" will be kept.

# %%
filters = []
filters.append(
    FeatureNameFilter(
        ["Australia", "New Zealand", "Tasmania"],
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
print(
    f"Coastlines of Australia and New Zealand have been written to {DATA_DIR}/Australia_and_New_Zealand_coastlines.gpmlz"
)

# %% [markdown]
#### Filter feature collection by feature names (case sensitive, partial match)

# Features whose name contains "Australia", "New Zealand" or "Tasmania" will be discarded.

# %%
filters = []
filters.append(
    FeatureNameFilter(
        ["Australia", "new zealand", "Tasmania"],
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

features.write(f"{DATA_DIR}/coastlines_exclude_Australia_and_New_Zealand.gpmlz")
print(
    f"Coastlines excluding Australia and New Zealand have been written to {DATA_DIR}/coastlines_exclude_Australia_and_New_Zealand.gpmlz"
)

# %% [markdown]
#### Search and filter feature collection by plate IDs

# Search feature collection by plate IDs

# Features whose plate ID matches the specified IDs will be kept. In this example, we keep features with plate IDs from 701 to 715, which correspond to the Africa and its neighboring plates.

# %%
filters = []
filters.append(
    PlateIDFilter(
        list(
            range(701, 716, 1)
        ),  # This is the list of plate IDs for the Africa and its neighboring plates
        exclude=False,
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

features.write(f"{DATA_DIR}/coastlines_with_plate_id_from_701_to_715.gpmlz")
print(
    f"Features with plate IDs from 701 to 715 have been written to {DATA_DIR}/coastlines_with_plate_id_from_701_to_715.gpmlz"
)

# %% [markdown]

##### Filter feature collection by plate IDs

# Features whose plate ID matches the specified IDs will be discarded. In this example, we keep features with plate IDs from 701 to 715, which correspond to the Africa and its neighboring plates.

# %%
filters = []
filters.append(
    PlateIDFilter(
        list(
            range(701, 716, 1)
        ),  # This is the list of plate IDs for the Africa and its neighboring plates
        exclude=True,  # exclude features with the specified plate IDs
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

features.write(f"{DATA_DIR}/coastlines_exclude_plate_id_from_701_to_715.gpmlz")
print(
    f"Features with plate IDs from 701 to 715 have been written to {DATA_DIR}/coastlines_exclude_plate_id_from_701_to_715.gpmlz"
)

# %% [markdown]

#### Search features in a feature collection by their birth ages

# Search all the features whose birth ages are older than 500 million years


# %%
filters = []
filters.append(
    BirthAgeFilter(
        500, keep_older=True
    )  # This filter will keep features that were born more than 500 million years ago.
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

features.write(f"{DATA_DIR}/coastlines_older_than_500_million_years.gpmlz")
print(
    f"Features with birth ages older than 500 million years have been written to {DATA_DIR}/coastlines_older_than_500_million_years.gpmlz"
)

# %% [markdown]

##### Keep all the features whose birth ages are younger than 500 million years


# %%
filters = []
filters.append(
    BirthAgeFilter(
        500, keep_older=False
    )  # This filter will keep features that were born less than 500 million years ago.
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

features.write(f"{DATA_DIR}/coastlines_younger_than_500_million_years.gpmlz")
print(
    f"Features with birth ages younger than 500 million years have been written to {DATA_DIR}/coastlines_younger_than_500_million_years.gpmlz"
)

# %% [markdown]

##### Keep all the features which had disappeared before 100 million years ago

# %%
# we need to use a different input file for this example.
topology_features = f"{DATA_DIR}/Feature_Geometries.gpmlz"
topology_features_url = (
    "https://repo.gplates.org/webdav/mchin/data/Feature_Geometries.gpmlz"
)

if not exists(topology_features):
    print(f"Downloading test file to {topology_features}...")
    urlretrieve(topology_features_url, topology_features)

topology_feature_collection = pygplates.FeatureCollection(topology_features)  # type: ignore

filters = []
filters.append(
    EndTimeFilter(
        100, disappear_before=True
    )  # This filter will keep features that had disappeared before 100 million years ago.
)

features = filter_feature_collection(
    topology_feature_collection,
    filters,
)
if len(features) == 0:
    print(
        "Warning: No features matched the search criteria. "
        "The output feature collection will be empty."
    )

features.write(
    f"{DATA_DIR}/topology_features_disappeared_before_100_million_years.gpmlz"
)
print(
    f"Topology features that had disappeared before 100 million years ago have been written to {DATA_DIR}/topology_features_disappeared_before_100_million_years.gpmlz"
)

# %% [markdown]

##### Keep all the features that had not yet disappeared as of 100 million years ago

# %%
filters = []
filters.append(
    EndTimeFilter(
        100, disappear_before=False
    )  # This filter will keep features that had not yet disappeared as of 100 million years ago
)

features = filter_feature_collection(
    topology_feature_collection,
    filters,
)
if len(features) == 0:
    print(
        "Warning: No features matched the search criteria. "
        "The output feature collection will be empty."
    )

features.write(
    f"{DATA_DIR}/topology_features_disappeared_after_100_million_years.gpmlz"
)
print(
    f"Topology features that had not yet disappeared as of 100 million years ago have been written to {DATA_DIR}/topology_features_disappeared_after_100_million_years.gpmlz"
)

# %% [markdown]

#### Search features in a feature collection by feature type

# keep all the features whose feature type is gpml:Basin or gpml:IslandArc


# %%
filters = []
filters.append(
    FeatureTypeFilter(
        "gpml:Basin|gpml:IslandArc",
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

features.write(f"{DATA_DIR}/coastlines_basins_and_island_arcs.gpmlz")
print(
    f"Features with feature type Basin or IslandArc have been written to {DATA_DIR}/coastlines_basins_and_island_arcs.gpmlz"
)


# %%
