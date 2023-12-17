"""Tools for reconstructing and plotting geological features and feature data through time.

Methods in this module reconstruct geological features using
[pyGPlates' `reconstruct` function](https://www.gplates.org/docs/pygplates/generated/pygplates.reconstruct.html),
turn them into plottable Shapely geometries, and plot them onto
Cartopy GeoAxes using Shapely and GeoPandas.

Classes
-------
 * `PlotTopologies`
 * `SubductionTeeth`

Functions
---------
 * `plot_subduction_teeth`
 * `shapelify_features` and its aliases:
   -  `shapelify_feature_lines`
   -  `shapelify_feature_polygons`
"""
from ..geometry import (
    shapelify_features,
    shapelify_feature_lines,
    shapelify_feature_polygons,
)
from .plot_topologies import PlotTopologies
from .subduction_teeth import (
    SubductionTeeth,
    plot_subduction_teeth,
)

__all__ = [
    "PlotTopologies",
    "SubductionTeeth",
    "plot_subduction_teeth",
    "shapelify_features",
    "shapelify_feature_lines",
    "shapelify_feature_polygons",
]
