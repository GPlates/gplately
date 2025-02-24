#!/usr/bin/env python3
from gplately import auxiliary

# use the auxiliary function to create a PlotTopologies instance
gplot = auxiliary.get_gplot("Muller2019", time=100)  # 100Ma

# save the reconstructed data to shapefiles
# the "get_" methods return GeoDataFrame objects
# you may save the data to other file formats, check out the GeoDataFrame doc
# https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_file.html
gplot.get_continents().to_file("continents_100Ma.shp")
gplot.get_coastlines().to_file("coastlines_100Ma.shp")
gplot.get_topological_plate_boundaries().to_file(
    "topological_plate_boundaries_100Ma.shp"
)
print("The files have been saved.")
