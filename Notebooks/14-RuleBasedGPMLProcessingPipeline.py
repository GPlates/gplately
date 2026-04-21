# %% [markdown]

## Rule Based GPML Processing Pipeline

# This notebook demonstrates how to use GPlately to query and manipulate GPML files based on user-defined rules,
# which can be based on various criteria such as feature name, plate ID, birth age, disappearance age, region of interest, and more.
# This allows users to create customized GPML processing pipelines that can be applied to a wide range of use cases in geoscience research.

# The search and filter functionality is implemented in the `gplately.utils.feature_filter` module.
# The main function is `filter_feature_collection`, which takes a feature collection and a list of filters,
# and returns a new feature collection that contains only the features that satisfy all the filters.
# Each filter is a subclass of the `FeatureFilter` class, which defines a `should_keep` method that takes a feature and
# returns True if the feature should be kept, and False otherwise.


# %% [markdown]

# 📒 This notebook is generated from 14-RuleBasedGPMLProcessingPipeline.py using the command `jupytext --to notebook Notebooks/14-RuleBasedGPMLProcessingPipeline.py -o Notebooks/14-RuleBasedGPMLProcessingPipeline.ipynb`.
# If you need to commit changes to this notebook to the GPlately repository, make your edits in 14-RuleBasedGPMLProcessingPipeline.py and then regenerate this Jupyter Notebook file.
# The reason that a .py file is used is to allow for easier version control and collaboration. And it is also more Copilot and code auto-formatting friendly.

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
    PropertyExistsFilter,
    PropertyValueFilter,
    RegionOfInterestFilter,
    filter_feature_collection,
)

DATA_DIR = "./feature-search-notebook-data"
makedirs(DATA_DIR, exist_ok=True)

# Download the input feature collection file if it does not exist.
# The input file we use in this notbook is a global coastline feature collection at present day.
# You can also use your own feature collection as input and apply the same search and filter methods as shown in this notebook.
coastlines_file = f"{DATA_DIR}/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"
coastlines_url = "https://repo.gplates.org/webdav/mchin/data/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"

if not exists(coastlines_file):
    print(f"Downloading test file to {coastlines_file}...")
    urlretrieve(coastlines_url, coastlines_file)

coastlines_feature_collection = pygplates.FeatureCollection(coastlines_file)  # type: ignore

# %% [markdown]
#### Search coastlines by feature names (case sensitive, partial match)

# Features whose name contains "Australia", "New Zealand" or "Tasmania" will be saved in the output file.

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
    coastlines_feature_collection,
    filters,
)

features.write(f"{DATA_DIR}/Australia_and_New_Zealand_coastlines.gpmlz")
print(
    f"The coastlines of Australia and New Zealand have been written to {DATA_DIR}/Australia_and_New_Zealand_coastlines.gpmlz"
)

# plot the output feature collection to check if the search worked as expected.
# You should see the coastlines of Australia, New Zealand and Tasmania in the plot below.
from gplately.mapping.cartopy_plot import _plot_feature_collection

_plot_feature_collection(
    pygplates.FeatureCollection(  # type: ignore
        f"{DATA_DIR}/Australia_and_New_Zealand_coastlines.gpmlz",
    ),
    title="Coastlines of Australia and New Zealand",
)


# %% [markdown]
#### Filter out "Australia", "New Zealand" and "Tasmania" from the global coastlines (case sensitive, partial match)

# Features whose name contains "Australia", "New Zealand" or "Tasmania" will be taken out of the global coastlines.

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
    coastlines_feature_collection,
    filters,
)

features.write(f"{DATA_DIR}/coastlines_exclude_Australia_and_New_Zealand.gpmlz")
print(
    f"The coastlines excluding Australia and New Zealand have been written to {DATA_DIR}/coastlines_exclude_Australia_and_New_Zealand.gpmlz"
)

# plot the output feature collection to check if the filter worked as expected.
# You should see the coastlines of Australia, New Zealand and Tasmania were taken out of the global coastlines in the plot below.
from gplately.mapping.cartopy_plot import _plot_feature_collection

_plot_feature_collection(
    pygplates.FeatureCollection(  # type: ignore
        f"{DATA_DIR}/coastlines_exclude_Australia_and_New_Zealand.gpmlz",
    ),
    title="Coastlines excluding Australia and New Zealand",
    figsize=(10, 5),
)

# %% [markdown]
#### Search Africa coastlines by plate IDs

# Features whose plate ID matches the specified IDs will be kept. In this example, we keep features with plate IDs from 701 to 715,
# which correspond to the Africa and its neighboring plates.

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
    coastlines_feature_collection,
    filters,
)

features.write(f"{DATA_DIR}/coastlines_with_plate_id_from_701_to_715.gpmlz")
print(
    f"Features with plate IDs from 701 to 715 have been written to {DATA_DIR}/coastlines_with_plate_id_from_701_to_715.gpmlz"
)

# plot the output feature collection to check if the search worked as expected.
# You should see the coastlines of Africa in the plot below.
import cartopy.crs as ccrs  # type: ignore
from gplately.mapping.cartopy_plot import _plot_feature_collection

_plot_feature_collection(
    pygplates.FeatureCollection(  # type: ignore
        f"{DATA_DIR}/coastlines_with_plate_id_from_701_to_715.gpmlz",
    ),
    title="Coastlines of Africa",
    projection=ccrs.Robinson(central_longitude=0),
)

# %% [markdown]

##### Filter coastlines by plate IDs

# Features whose plate ID matches the specified IDs will be discarded. In this example, we filter out features with plate IDs from 701 to 715,
# which correspond to the Africa and its neighboring plates.

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
    coastlines_feature_collection,
    filters,
)
features.write(f"{DATA_DIR}/coastlines_exclude_plate_id_from_701_to_715.gpmlz")
print(
    f"Features with plate IDs from 701 to 715 have been written to {DATA_DIR}/coastlines_exclude_plate_id_from_701_to_715.gpmlz"
)
# plot the output feature collection to check if the filter worked as expected.
# You should see that the coastlines of Africa are excluded in the plot below.
import cartopy.crs as ccrs  # type: ignore
from gplately.mapping.cartopy_plot import _plot_feature_collection

_plot_feature_collection(
    pygplates.FeatureCollection(  # type: ignore
        f"{DATA_DIR}/coastlines_exclude_plate_id_from_701_to_715.gpmlz",
    ),
    title="Coastlines excluding Africa",
    projection=ccrs.Robinson(central_longitude=0),
)
# %% [markdown]

#### Search features by birth ages

# Find all the features whose birth ages are older than 500 million years

# %%
filters = []
filters.append(
    BirthAgeFilter(
        500, keep_older=True
    )  # This filter will keep features that were born more than 500 million years ago.
)

features = filter_feature_collection(
    coastlines_feature_collection,
    filters,
)

features.write(f"{DATA_DIR}/coastlines_older_than_500_million_years.gpmlz")
print(
    f"Features with birth ages older than 500 million years have been written to {DATA_DIR}/coastlines_older_than_500_million_years.gpmlz"
)

# plot the output feature collection to check if the search worked as expected.
# You should see that the coastlines older than 500 million years in the plot below.
from gplately.mapping.cartopy_plot import _plot_feature_collection

_plot_feature_collection(
    pygplates.FeatureCollection(  # type: ignore
        f"{DATA_DIR}/coastlines_older_than_500_million_years.gpmlz",
    ),
    title="Coastlines older than 500 million years",
)

# %% [markdown]

##### Find all the coastlines younger than 500 million years

# %%
filters = []
filters.append(
    BirthAgeFilter(
        500, keep_older=False
    )  # This filter will keep features that were born less than 500 million years ago.
)

features = filter_feature_collection(
    coastlines_feature_collection,
    filters,
)

features.write(f"{DATA_DIR}/coastlines_younger_than_500_million_years.gpmlz")
print(
    f"Features with birth ages younger than 500 million years have been written to {DATA_DIR}/coastlines_younger_than_500_million_years.gpmlz"
)

# plot the output feature collection to check if the search worked as expected.
# You should see that the coastlines younger than 500 million years in the plot below.
from gplately.mapping.cartopy_plot import _plot_feature_collection

_plot_feature_collection(
    pygplates.FeatureCollection(  # type: ignore
        f"{DATA_DIR}/coastlines_younger_than_500_million_years.gpmlz",
    ),
    title="Coastlines younger than 500 million years",
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

features.write(
    f"{DATA_DIR}/topology_features_disappeared_after_100_million_years.gpmlz"
)
print(
    f"Topology features that had not yet disappeared as of 100 million years ago have been written to {DATA_DIR}/topology_features_disappeared_after_100_million_years.gpmlz"
)

# %% [markdown]

#### Search by feature type

# find all the features whose feature type is gpml:Basin or gpml:IslandArc

# %%
filters = []
filters.append(
    FeatureTypeFilter(
        "gpml:Basin|gpml:IslandArc",
    )
)

features = filter_feature_collection(
    coastlines_feature_collection,
    filters,
)

features.write(f"{DATA_DIR}/coastlines_basins_and_island_arcs.gpmlz")
print(
    f"Features with feature type Basin or IslandArc have been written to {DATA_DIR}/coastlines_basins_and_island_arcs.gpmlz"
)

# plot the output feature collection to check if the search worked as expected.
# You should see basin and island arcs in the plot below.
from gplately.mapping.cartopy_plot import _plot_feature_collection

_plot_feature_collection(
    pygplates.FeatureCollection(  # type: ignore
        f"{DATA_DIR}/coastlines_basins_and_island_arcs.gpmlz",
    ),
    title="Coastlines with feature type Basin or IslandArc",
)

# %% [markdown]

#### Find features with gpml:subductionPolarity property

# %%
filters = []
filters.append(PropertyExistsFilter("gpml:subductionPolarity", not_exists=False))

features = filter_feature_collection(
    topology_feature_collection,
    filters,
)

features.write(f"{DATA_DIR}/topology_features_with_subduction_polarity.gpmlz")
print(
    f"Topology features with subduction polarity have been written to {DATA_DIR}/topology_features_with_subduction_polarity.gpmlz"
)

# plot the output feature collection to check if the search worked as expected.
# You should see features with gpml:subductionPolarity property in the plot below.
from gplately.mapping.cartopy_plot import _plot_feature_collection

_plot_feature_collection(
    pygplates.FeatureCollection(  # type: ignore
        f"{DATA_DIR}/topology_features_with_subduction_polarity.gpmlz",
    ),
    title="Topology features with subduction polarity",
)

# %% [markdown]

#### Find features whose gpml:subductionPolarity is "Left"

# %%
# open the .gpml file to find the property value name, such as this SubductionPolarityEnumeration
subduction_polarity_left = pygplates.Enumeration(  # type: ignore
    pygplates.EnumerationType.create_gpml("SubductionPolarityEnumeration"), "Left"  # type: ignore
)

filters = []
filters.append(
    PropertyValueFilter(
        "gpml:subductionPolarity", subduction_polarity_left, not_match=False
    )
)

features = filter_feature_collection(
    topology_feature_collection,
    filters,
)

features.write(f"{DATA_DIR}/topology_features_with_subduction_polarity_left.gpmlz")
print(
    f"Topology features with subduction polarity of Left have been written to {DATA_DIR}/topology_features_with_subduction_polarity_left.gpmlz"
)

# plot the output feature collection to check if the search worked as expected.
# You should see features whose gpml:subductionPolarity is Left in the plot below.
from gplately.mapping.cartopy_plot import _plot_feature_collection

_plot_feature_collection(
    pygplates.FeatureCollection(  # type: ignore
        f"{DATA_DIR}/topology_features_with_subduction_polarity_left.gpmlz",
    ),
    title="Topology features with subduction polarity of Left",
)

# %% [markdown]
#### Find points inside a region of interest

# Firstly, we create a feature collection for the vertices of an icosahedron mesh.
# Then we search for the vertices that are located within a region of interest defined by a bounding box (left, right, bottom, top).
# Finally, we create a feature for the region of interest and add it to the output feature collection for visualization.


# %%
from gplately.lib.icosahedron import get_mesh, xyz2lonlat

# Create a feature collection for the vertices of an icosahedron mesh.
# We will use this feature collection as input for searching features within a region of interest in the next step.
mesh_resolution = 5
vertices_0, faces_0 = get_mesh(mesh_resolution)
seen = set()
mesh_feature_collection = pygplates.FeatureCollection()  # type: ignore
for v in vertices_0:
    lon, lat = xyz2lonlat(v[0], v[1], v[2])
    point = f"{lon:0.2f} {lat:0.2f}"
    if point in seen:
        continue
    else:
        seen.add(point)
    feature = pygplates.Feature()  # type: ignore
    feature.set_geometry(pygplates.PointOnSphere(lat, lon))  # type: ignore
    mesh_feature_collection.add(feature)  # type: ignore

mesh_feature_collection.write(f"{DATA_DIR}/icosahedron_mesh_{mesh_resolution}.gpmlz")
print(
    f"Icosahedron mesh have been written to {DATA_DIR}/icosahedron_mesh_{mesh_resolution}.gpmlz"
)

# %% [markdown]
# if you open the icosahedron_mesh_5.gpmlz file in GPlates, you will see the vertices of the icosahedron mesh,
# which are represented as point features.

# <figure>
#  <img src="https://raw.githubusercontent.com/GPlates/gplately/refs/heads/377-feature-request-gpml-file-management-workflow/Notebooks/NotebookFiles/Notebook14/icosahedron_mesh.png" width="300">
#  <figcaption><b>Icosahedron Mesh</b></figcaption>
# </figure>

# %% [markdown]
# Search for the vertices that are located within a region of interest defined by a bounding box (left, right, bottom, top) in longitude and latitude.


# %%
# define the bounding box for the region of interest
(left, right, bottom, top) = (120, 140, -10, 10)
filters = []
filters.append(RegionOfInterestFilter((left, right, bottom, top), exclude=False))

features = filter_feature_collection(
    mesh_feature_collection,
    filters,
)

# Create a feature for the region of interest and add it to the output feature collection for visualization.
region_of_interest_feature = pygplates.Feature()  # type: ignore
region_of_interest = pygplates.PolygonOnSphere((lat, lon) for lon, lat in [(left, bottom), (left, top), (right, top), (right, bottom)])  # type: ignore
region_of_interest_feature.set_geometry(region_of_interest)  # type: ignore
features.add(region_of_interest_feature)  # type: ignore

features.write(f"{DATA_DIR}/icosahedron_vertices_in_region.gpmlz")
print(
    f"Icosahedron vertices in the region of interest have been written to {DATA_DIR}/icosahedron_vertices_in_region.gpmlz"
)

# %% [markdown]
# If you open icosahedron_vertices_in_region.gpmlz in GPlates, you will see the vertices of the icosahedron mesh that
# are located within the bounding box defined by (120, 140, -10, 10),
# as well as a polygon feature that represents the region of interest defined by the bounding box.

# <figure>
#  <img src="https://raw.githubusercontent.com/GPlates/gplately/refs/heads/377-feature-request-gpml-file-management-workflow/Notebooks/NotebookFiles/Notebook14/icosahedron_vertices_in_region.png" width="300">
#  <figcaption><b>Icosahedron Vertices within (120, 140, -10, 10)</b></figcaption>
# </figure>

# %% [markdown]
##### Search for the vertices that are located within the mainland of Australia.


# %%
# Firstly, we search for the feature that corresponds to the mainland of Australia in the global coastline feature collection.
# Then we use the geometry of this feature to define the region of interest for filtering the vertices of the icosahedron mesh.
filters = []
filters.append(
    FeatureNameFilter(
        ["Australia"],
        exact_match=True,
        case_sensitive=True,
    )
)
features = filter_feature_collection(
    coastlines_feature_collection,
    filters,
)

australia_mainland_geometry = None
australia_mainland_feature = None
largest_area = 0
for feature in features:
    geometry = feature.get_geometry()  # type: ignore
    if isinstance(geometry, pygplates.PolygonOnSphere):  # type: ignore
        if geometry.get_area() > largest_area:  # type: ignore
            largest_area = geometry.get_area()  # type: ignore
            australia_mainland_geometry = geometry
            australia_mainland_feature = feature
    else:
        continue

##### search for the vertices that are located within the mainland of Australia
filters = []
filters.append(RegionOfInterestFilter(australia_mainland_geometry, exclude=False))

features = filter_feature_collection(
    mesh_feature_collection,
    filters,
)

features.add(australia_mainland_feature)  # type: ignore
features.write(f"{DATA_DIR}/icosahedron_vertices_within_australia.gpmlz")
print(
    f"Icosahedron vertices in the region of interest have been written to {DATA_DIR}/icosahedron_vertices_within_australia.gpmlz"
)
# %% [markdown]
# If you open icosahedron_vertices_within_australia.gpmlz in GPlates, you will see the vertices of the icosahedron mesh that
# are located within the mainland of Australia, as well as a polygon feature that represents the mainland of Australia

# <figure>
#  <img src="https://raw.githubusercontent.com/GPlates/gplately/refs/heads/377-feature-request-gpml-file-management-workflow/Notebooks/NotebookFiles/Notebook14/icosahedron_vertices_within_australia.png" width="300">
#  <figcaption><b>Icosahedron Vertices within Australia</b></figcaption>
# </figure>

# %%

from gplately.gpml import (
    get_unique_feature_names,
    GPML_to_GeoDataFrame,
    feature_type_getter,
    feature_name_getter,
    plate_id_getter,
)

print(
    GPML_to_GeoDataFrame(
        coastlines_feature_collection,
        property_getters={
            "name": feature_name_getter,
            "type": feature_type_getter,
            "plate_id": plate_id_getter,
        },
    )
)
