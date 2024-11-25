#
#    Copyright (C) 2024 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

"""Tools for reconstructing and plotting geological features and feature data through time.

Methods in `plot.py` reconstruct geological features using 
[pyGPlates' `reconstruct` function](https://www.gplates.org/docs/pygplates/generated/pygplates.reconstruct.html),
turns them into plottable Shapely geometries, and plots them onto Cartopy GeoAxes using Shapely and GeoPandas.

Classes
-------
PlotTopologies
"""

import logging
import math
import warnings
from functools import partial

import cartopy.crs as ccrs
import geopandas as gpd
import numpy as np
import pygplates
from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
from shapely.ops import linemerge

from . import ptt
from .decorators import (
    append_docstring,
    validate_reconstruction_time,
    validate_topology_availability,
)
from .gpml import _load_FeatureCollection
from .pygplates import FeatureCollection as _FeatureCollection
from .reconstruction import PlateReconstruction as _PlateReconstruction
from .tools import EARTH_RADIUS
from .utils.feature_utils import shapelify_features as _shapelify_features
from .utils.plot_utils import _clean_polygons, _meridian_from_ax
from .utils.plot_utils import plot_subduction_teeth as _plot_subduction_teeth

logger = logging.getLogger("gplately")

PLOT_DOCSTRING = """

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the `{0}` lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc. for plotting `{0}` geometries.
            See `Matplotlib` keyword arguments [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with `{0}` features plotted onto the chosen map projection.
"""

GET_DATE_DOCSTRING = """

        Parameters
        ----------
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Returns
        -------
        geopandas.GeoDataFrame
            A pandas.DataFrame that has a column with `{0}` geometry.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to
            `PlotTopologies`. This is needed to construct `{0}`
            to the requested `time` and thus populate the GeoDataFrame.

        Notes
        -----
        The `{0}` needed to produce the GeoDataFrame are automatically
        constructed if the optional `time` parameter is passed to the
        `PlotTopologies` object before calling this function. `time` can be
        passed either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 # Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `{0}`
        are reconstructed, they are converted into Shapely lines whose
        coordinates are passed to a geopandas GeoDataFrame.

"""


def shapelify_features(*args, **kwargs):
    return _shapelify_features(*args, **kwargs)


def shapelify_feature_lines(*args, **kwargs):
    return _shapelify_features(*args, **kwargs)


def shapelify_feature_polygons(*args, **kwargs):
    return _shapelify_features(*args, **kwargs)


def plot_subduction_teeth(*args, **kwargs):
    return _plot_subduction_teeth(*args, **kwargs)


plot_subduction_teeth.__doc__ = _plot_subduction_teeth.__doc__
shapelify_features.__doc__ = _shapelify_features.__doc__
shapelify_feature_lines.__doc__ = _shapelify_features.__doc__
shapelify_feature_polygons.__doc__ = _shapelify_features.__doc__


class PlotTopologies(object):
    """A class with tools to read, reconstruct and plot topology features at specific
    reconstruction times.

    `PlotTopologies` is a shorthand for PyGPlates and Shapely functionalities that:

    * Read features held in GPlates GPML (GPlates Markup Language) files and
    ESRI shapefiles;
    * Reconstruct the locations of these features as they migrate through
    geological time;
    * Turn these reconstructed features into Shapely geometries for plotting
    on `cartopy.mpl.geoaxes.GeoAxes` or `cartopy.mpl.geoaxes.GeoAxesSubplot` map
    Projections.

    To call the `PlotTopologies` object, supply:

    * an instance of the GPlately `plate_reconstruction` object

    and optionally,

    * a `coastline_filename`
    * a `continent_filename`
    * a `COB_filename`
    * a reconstruction `time`
    * an `anchor_plate_id`

    For example:

        # Calling the PlotTopologies object
        gplot = gplately.plot.PlotTopologies(plate_reconstruction,
                                            coastline_filename=None,
                                            continent_filename=None,
                                            COB_filename=None,
                                            time=None,
                                            anchor_plate_id=0,
                )

        # Setting a new reconstruction time
        gplot.time = 20 # Ma

    The `coastline_filename`, `continent_filename` and `COB_filename` can be single
    strings to GPML and/or shapefiles, as well as instances of `pygplates.FeatureCollection`.
    If using GPlately's `DataServer` object to source these files, they will be passed as
    `pygplates.FeatureCollection` items.

    Some features for plotting (like plate boundaries) are taken from the `PlateReconstruction`
    object's`topology_features` attribute. They have already been reconstructed to the given
    `time` using [Plate Tectonic Tools](https://github.com/EarthByte/PlateTectonicTools).
    Simply provide a new reconstruction time by changing the `time` attribute, e.g.

        gplot.time = 20 # Ma

    which will automatically reconstruct all topologies to the specified time.
    You __MUST__ set `gplot.time` before plotting anything.

    A variety of geological features can be plotted on GeoAxes/GeoAxesSubplot maps
    as Shapely `MultiLineString` or `MultiPolygon` geometries, including:

    * subduction boundaries & subduction polarity teeth
    * mid-ocean ridge boundaries
    * transform boundaries
    * miscellaneous boundaries
    * coastline polylines
    * continental polygons and
    * continent-ocean boundary polylines
    * topological plate velocity vector fields
    * netCDF4 MaskedArray or ndarray raster data:
        - seafloor age grids
        - paleo-age grids
        - global relief (topography and bathymetry)
    * assorted reconstructable feature data, for example:
        - seafloor fabric
        - large igneous provinces
        - volcanic provinces

    Attributes
    ----------
    plate_reconstruction : instance of <gplately.reconstruction.PlateReconstruction>
        The GPlately `PlateReconstruction` object will be used to access a plate
        `rotation_model` and a set of `topology_features` which contains plate boundary
        features like trenches, ridges and transforms.

    anchor_plate_id : int, default 0
        The anchor plate ID used for reconstruction.

    base_projection : instance of <cartopy.crs.{transform}>, default <cartopy.crs.PlateCarree> object
        where {transform} is the map Projection to use on the Cartopy GeoAxes.
        By default, the base projection is set to cartopy.crs.PlateCarree. See the
        [Cartopy projection list](https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html)
        for all supported Projection types.

    coastlines : str, or instance of <pygplates.FeatureCollection>
        The full string path to a coastline feature file. Coastline features can also
        be passed as instances of the `pygplates.FeatureCollection` object (this is
        the case if these features are sourced from the `DataServer` object).

    continents : str, or instance of <pygplates.FeatureCollection>
        The full string path to a continent feature file. Continent features can also
        be passed as instances of the `pygplates.FeatureCollection` object (this is
        the case if these features are sourced from the `DataServer` object).

    COBs : str, or instance of <pygplates.FeatureCollection>
        The full string path to a COB feature file. COB features can also be passed
        as instances of the `pygplates.FeatureCollection` object (this is the case
        if these features are sourced from the `DataServer` object).

    coastlines : iterable/list of <pygplates.ReconstructedFeatureGeometry>
        A list containing coastline features reconstructed to the specified `time` attribute.

    continents : iterable/list of <pygplates.ReconstructedFeatureGeometry>
        A list containing continent features reconstructed to the specified `time` attribute.

    COBs : iterable/list of <pygplates.ReconstructedFeatureGeometry>
        A list containing COB features reconstructed to the specified `time` attribute.

    time : float
        The time (Ma) to reconstruct and plot geological features to.

    topologies : iterable/list of <pygplates.Feature>
        A list containing assorted topologies like:

        - pygplates.FeatureType.gpml_topological_network
        - pygplates.FeatureType.gpml_oceanic_crust
        - pygplates.FeatureType.gpml_topological_slab_boundary
        - pygplates.FeatureType.gpml_topological_closed_plate_boundary

    ridges : iterable/list of <pygplates.Feature>
        A list containing ridge and transform boundary sections of type
        pygplates.FeatureType.gpml_mid_ocean_ridge

    transforms : iterable/list of <pygplates.Feature>
        A list containing transform boundary sections of type pygplates.FeatureType.gpml_transforms

    trenches : iterable/list of <pygplates.Feature>
        A list containing trench boundary sections of type pygplates.FeatureType.gpml_subduction_zone

    trench_left : iterable/list of <pygplates.Feature>
        A list containing left subduction boundary sections of type pygplates.FeatureType.gpml_subduction_zone

    trench_right : iterable/list of <pygplates.Feature>
        A list containing right subduction boundary sections of type pygplates.FeatureType.gpml_subduction_zone

    other : iterable/list of <pygplates.Feature>
        A list containing other geological features like unclassified features, extended continental crusts,
        continental rifts, faults, orogenic belts, fracture zones, inferred paleo boundaries, terrane
        boundaries and passive continental boundaries.

    """

    def __init__(
        self,
        plate_reconstruction,
        coastlines=None,
        continents=None,
        COBs=None,
        time=None,
        anchor_plate_id=0,
    ):
        self.plate_reconstruction = plate_reconstruction

        if self.plate_reconstruction.topology_features is None:
            self.plate_reconstruction.topology_features = []
            logger.warn("Plate model does not have topology features.")

        self.base_projection = ccrs.PlateCarree()

        # store these for when time is updated
        # make sure these are initialised as FeatureCollection objects

        self._coastlines = _load_FeatureCollection(coastlines)
        self._continents = _load_FeatureCollection(continents)
        self._COBs = _load_FeatureCollection(COBs)

        self.coastlines = None
        self.continents = None
        self.COBs = None
        self._topological_plate_boundaries = None
        self._topologies = None

        self._anchor_plate_id = self._check_anchor_plate_id(anchor_plate_id)

        # store topologies for easy access
        # setting time runs the update_time routine
        self._time = None
        if time is not None:
            self.time = time

    def __getstate__(self):
        filenames = self.plate_reconstruction.__getstate__()

        # add important variables from Points object
        if self._coastlines:
            filenames["coastlines"] = self._coastlines.filenames
        if self._continents:
            filenames["continents"] = self._continents.filenames
        if self._COBs:
            filenames["COBs"] = self._COBs.filenames
        filenames["time"] = self.time
        filenames["plate_id"] = self._anchor_plate_id

        return filenames

    def __setstate__(self, state):
        plate_reconstruction_args = [state["rotation_model"], None, None]
        if "topology_features" in state:
            plate_reconstruction_args[1] = state["topology_features"]
        if "static_polygons" in state:
            plate_reconstruction_args[2] = state["static_polygons"]

        self.plate_reconstruction = _PlateReconstruction(*plate_reconstruction_args)

        self._coastlines = None
        self._continents = None
        self._COBs = None
        self.coastlines = None
        self.continents = None
        self.COBs = None

        # reinstate unpicklable items
        if "coastlines" in state:
            self._coastlines = _FeatureCollection()
            for feature in state["coastlines"]:
                self._coastlines.add(_FeatureCollection(feature))

        if "continents" in state:
            self._continents = _FeatureCollection()
            for feature in state["continents"]:
                self._continents.add(_FeatureCollection(feature))

        if "COBs" in state:
            self._COBs = _FeatureCollection()
            for feature in state["COBs"]:
                self._COBs.add(_FeatureCollection(feature))

        self._anchor_plate_id = state["plate_id"]
        self.base_projection = ccrs.PlateCarree()
        self._time = None

    @property
    def topological_plate_boundaries(self):
        if self._topological_plate_boundaries is None:
            self.update_time(self.time)
        return self._topological_plate_boundaries

    @property
    def topologies(self):
        self._resolve_both_boundaries_and_networks()
        return self._topologies

    @property
    def time(self):
        """The reconstruction time."""
        return self._time

    @time.setter
    def time(self, var):
        """Allows the time attribute to be changed. Updates all instances of the time attribute in the object (e.g.
        reconstructions and resolving topologies will use this new time).

        Raises
        ------
        ValueError
            If the chosen reconstruction time is <0 Ma.
        """
        if var < 0:
            raise ValueError("The 'time' property must be greater than 0.")

        if self.time is None or (not math.isclose(var, self.time)):
            self.update_time(var)

    @property
    def anchor_plate_id(self):
        """Anchor plate ID for reconstruction. Must be an integer >= 0."""
        return self._anchor_plate_id

    @anchor_plate_id.setter
    def anchor_plate_id(self, anchor_plate):
        self._anchor_plate_id = self._check_anchor_plate_id(anchor_plate)
        if self.time is not None:
            self.update_time(self.time)

    @staticmethod
    def _check_anchor_plate_id(id):
        id = int(id)
        if id < 0:
            raise ValueError("Invalid anchor plate ID: {}".format(id))
        return id

    def update_time(self, time):
        """Re-reconstruct features and topologies to the time specified by the `PlotTopologies` `time` attribute
        whenever it or the anchor plate is updated.

        Notes
        -----
        The following `PlotTopologies` attributes are updated whenever a reconstruction `time` attribute is set:

        - resolved topology features (topological plates and networks)
        - ridge and transform boundary sections (resolved features)
        - ridge boundary sections (resolved features)
        - transform boundary sections (resolved features)
        - subduction boundary sections (resolved features)
        - left subduction boundary sections (resolved features)
        - right subduction boundary sections (resolved features)
        - other boundary sections (resolved features) that are not subduction zones or mid-ocean ridges
        (ridge/transform)

        Moreover, coastlines, continents and COBs are reconstructed to the new specified `time`.
        """
        assert time is not None, "time must be set to a valid reconstruction time"

        self._time = float(time)
        (
            self._topological_plate_boundaries,
            self.ridges,
            self._ridges_do_not_use_for_now,
            self._transforms_do_not_use_for_now,
            self.trenches,
            self.trench_left,
            self.trench_right,
            self.other,
        ) = ptt.resolve_topologies.resolve_topologies_into_features(
            self.plate_reconstruction.rotation_model,
            self.plate_reconstruction.topology_features,
            self.time,
            anchor_plate_id=self.anchor_plate_id,
            # use ResolveTopologyType.boundary parameter to resolve rigid plate boundary only
            # because the Mid-ocean ridges(and transforms) should not contain lines from topological networks
            resolve_topology_types=pygplates.ResolveTopologyType.boundary,
        )
        self._topologies = None

        # miscellaneous boundaries
        self.continental_rifts = []
        self.faults = []
        self.fracture_zones = []
        self.inferred_paleo_boundaries = []
        self.terrane_boundaries = []
        self.transitional_crusts = []
        self.orogenic_belts = []
        self.sutures = []
        self.continental_crusts = []
        self.extended_continental_crusts = []
        self.passive_continental_boundaries = []
        self.slab_edges = []
        self.transforms = []
        self.unclassified_features = []

        for topol in self.other:
            if topol.get_feature_type() == pygplates.FeatureType.gpml_continental_rift:
                self.continental_rifts.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_fault:
                self.faults.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_fracture_zone:
                self.fracture_zones.append(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_inferred_paleo_boundary
            ):
                self.inferred_paleo_boundaries.append(topol)

            elif (
                topol.get_feature_type() == pygplates.FeatureType.gpml_terrane_boundary
            ):
                self.terrane_boundaries.append(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_transitional_crust
            ):
                self.transitional_crusts.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_orogenic_belt:
                self.orogenic_belts.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_suture:
                self.sutures.append(topol)

            elif (
                topol.get_feature_type() == pygplates.FeatureType.gpml_continental_crust
            ):
                self.continental_crusts.append(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_extended_continental_crust
            ):
                self.extended_continental_crusts.append(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_passive_continental_boundary
            ):
                self.passive_continental_boundaries.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_slab_edge:
                self.slab_edges.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_transform:
                self.transforms.append(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_unclassified_feature
            ):
                self.unclassified_features.append(topol)

        # reconstruct other important polygons and lines
        if self._coastlines:
            self.coastlines = self.plate_reconstruction.reconstruct(
                self._coastlines,
                self.time,
                from_time=0,
                anchor_plate_id=self.anchor_plate_id,
            )

        if self._continents:
            self.continents = self.plate_reconstruction.reconstruct(
                self._continents,
                self.time,
                from_time=0,
                anchor_plate_id=self.anchor_plate_id,
            )

        if self._COBs:
            self.COBs = self.plate_reconstruction.reconstruct(
                self._COBs, self.time, from_time=0, anchor_plate_id=self.anchor_plate_id
            )

    # subduction teeth
    def _tessellate_triangles(
        self, features, tesselation_radians, triangle_base_length, triangle_aspect=1.0
    ):
        """Places subduction teeth along subduction boundary line segments within a MultiLineString shapefile.

        Parameters
        ----------
        shapefilename  : str
            Path to shapefile containing the subduction boundary features

        tesselation_radians : float
            Parametrises subduction teeth density. Triangles are generated only along line segments with distances
            that exceed the given threshold tessellation_radians.

        triangle_base_length : float
            Length of teeth triangle base

        triangle_aspect : float, default=1.0
            Aspect ratio of teeth triangles. Ratio is 1.0 by default.

        Returns
        -------
        X_points : (n,3) array
            X points that define the teeth triangles
        Y_points : (n,3) array
            Y points that define the teeth triangles
        """

        tesselation_degrees = np.degrees(tesselation_radians)
        triangle_pointsX = []
        triangle_pointsY = []

        date_line_wrapper = pygplates.DateLineWrapper()

        for feature in features:
            cum_distance = 0.0

            for geometry in feature.get_geometries():
                wrapped_lines = date_line_wrapper.wrap(geometry)
                for line in wrapped_lines:
                    pts = np.array(
                        [
                            (p.get_longitude(), p.get_latitude())
                            for p in line.get_points()
                        ]
                    )

                    for p in range(0, len(pts) - 1):
                        A = pts[p]
                        B = pts[p + 1]

                        AB_dist = B - A
                        AB_norm = AB_dist / np.hypot(*AB_dist)
                        cum_distance += np.hypot(*AB_dist)

                        # create a new triangle if cumulative distance is exceeded.
                        if cum_distance >= tesselation_degrees:
                            C = A + triangle_base_length * AB_norm

                            # find normal vector
                            AD_dist = np.array([AB_norm[1], -AB_norm[0]])
                            AD_norm = AD_dist / np.linalg.norm(AD_dist)

                            C0 = A + 0.5 * triangle_base_length * AB_norm

                            # project point along normal vector
                            D = C0 + triangle_base_length * triangle_aspect * AD_norm

                            triangle_pointsX.append([A[0], C[0], D[0]])
                            triangle_pointsY.append([A[1], C[1], D[1]])

                            cum_distance = 0.0

        return np.array(triangle_pointsX), np.array(triangle_pointsY)

    @validate_reconstruction_time
    def get_feature(
        self,
        feature,
        central_meridian=0.0,
        tessellate_degrees=None,
        validate_reconstruction_time=True,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed features.

        Notes
        -----
        The feature needed to produce the GeoDataFrame should already be constructed to a `time`.
        This function converts the feature into a set of Shapely geometries whose coordinates are
        passed to a geopandas GeoDataFrame.

        Parameters
        ----------
        feature : instance of <pygplates.Feature>
            A feature reconstructed to `time`.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `feature` geometries.

        """

        if feature is None:
            raise ValueError(
                "The 'feature' parameter is None. Make sure a valid `feature` object has been provided."
            )
        shp = shapelify_features(
            feature,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

        return gpd.GeoDataFrame({"geometry": shp}, geometry="geometry")

    @append_docstring(PLOT_DOCSTRING.format("feature"))
    def plot_feature(self, ax, feature, feature_name="", color="black", **kwargs):
        """Plot pygplates.FeatureCollection or pygplates.Feature onto a map."""
        if not feature:
            logger.warn(
                f"The given feature({feature_name}:{feature}) in model:{self.plate_reconstruction.plate_model_name} is empty and will not be plotted."
            )
            return ax
        else:
            if "edgecolor" not in kwargs.keys():
                kwargs["edgecolor"] = color
            if "facecolor" not in kwargs.keys():
                kwargs["facecolor"] = "none"
            return self._plot_feature(ax, partial(self.get_feature, feature), **kwargs)

    def _plot_feature(self, ax, get_feature_func, **kwargs):
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", 1)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        if not callable(get_feature_func):
            raise Exception("The 'get_feature_func' parameter must be callable.")
        gdf = get_feature_func(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

        if len(gdf) == 0:
            logger.warn("No feature found for plotting. Do nothing and return.")
            return ax

        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection

        return gdf.plot(ax=ax, **kwargs)

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("coastlines"))
    def get_coastlines(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed coastline polygons."""
        return self.get_feature(
            self.coastlines,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("coastlines"))
    def plot_coastlines(self, ax, color="black", **kwargs):
        """Plot reconstructed coastline polygons onto a standard map Projection.

        Notes
        -----
        The `coastlines` for plotting are accessed from the `PlotTopologies` object's `coastlines` attribute.
        These `coastlines` are reconstructed to the `time` passed to the `PlotTopologies` object and converted into Shapely polylines.
        The reconstructed `coastlines` are added onto the GeoAxes or GeoAxesSubplot map `ax` using GeoPandas.
        Map resentation details (e.g. facecolor, edgecolor, alpha…) are permitted as keyword arguments.
        """
        return self.plot_feature(
            ax,
            self.coastlines,
            feature_name="coastlines",
            color=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("continents"))
    def get_continents(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed continental polygons."""
        return self.get_feature(
            self.continents,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("continents"))
    def plot_continents(self, ax, color="black", **kwargs):
        """Plot reconstructed continental polygons onto a standard map Projection.

        Notes
        -----
        The `continents` for plotting are accessed from the `PlotTopologies` object's `continents` attribute.
        These `continents` are reconstructed to the `time` passed to the `PlotTopologies` object and converted into Shapely polygons.
        The reconstructed `coastlines` are plotted onto the GeoAxes or GeoAxesSubplot map `ax` using GeoPandas.
        Map presentation details (e.g. facecolor, edgecolor, alpha…) are permitted as keyword arguments.
        """
        return self.plot_feature(
            ax,
            self.continents,
            feature_name="continents",
            color=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("COBs"))
    def get_continent_ocean_boundaries(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed continent-ocean boundary lines."""
        return self.get_feature(
            self.COBs,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("continent ocean boundaries"))
    def plot_continent_ocean_boundaries(self, ax, color="black", **kwargs):
        """Plot reconstructed continent-ocean boundary (COB) polygons onto a standard map Projection.

        Notes
        -----
        The `COBs` for plotting are accessed from the `PlotTopologies` object's
        `COBs` attribute. These `COBs` are reconstructed to the `time`
        passed to the `PlotTopologies` object and converted into Shapely polylines.
        The reconstructed `COBs` are plotted onto the GeoAxes or GeoAxesSubplot map
        `ax` using GeoPandas. Map presentation details (e.g. `facecolor`, `edgecolor`, `alpha`…)
        are permitted as keyword arguments.

        These COBs are transformed into shapely geometries and added onto the chosen map for a specific geological time
        (supplied to the PlotTopologies object). Map presentation details (e.g. facecolor, edgecolor, alpha…) are permitted.
        """
        return self.plot_feature(
            ax,
            self.COBs,
            feature_name="continent_ocean_boundaries",
            color=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("ridges"))
    def get_ridges(
        self,
        central_meridian=0.0,
        tessellate_degrees=1,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed mid-ocean ridge lines(gpml:MidOceanRidge)."""
        logger.debug("Getting reconstructed mid-ocean ridge lines")
        return self.get_feature(
            self.ridges,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_topology_availability("ridges")
    @append_docstring(PLOT_DOCSTRING.format("ridges"))
    def plot_ridges(self, ax, color="black", **kwargs):
        """Plot reconstructed mid-ocean ridge lines(gpml:MidOceanRidge) onto a map.

        Notes
        -----
        The `ridges` sections for plotting are accessed from the
        `PlotTopologies` object's `ridges` attribute. These `ridges`
        are reconstructed to the `time` passed to the `PlotTopologies` object and converted
        into Shapely polylines. The reconstructed `ridges` are plotted onto the
        GeoAxes or GeoAxesSubplot map `ax` using GeoPandas. Map presentation details
        (e.g. `facecolor`, `edgecolor`, `alpha`…) are permitted as keyword arguments.

        Note: The `ridges` geometries are wrapped to the dateline using
        pyGPlates' [DateLineWrapper](https://www.gplates.org/docs/pygplates/generated/pygplates.datelinewrapper)
        by splitting a polyline into multiple polylines at the dateline. This is to avoid
        horizontal lines being formed between polylines at longitudes of -180 and 180 degrees.
        Point features near the poles (-89 & 89 degree latitude) are also clipped to ensure
        compatibility with Cartopy.
        """
        logger.debug("Plotting reconstructed mid-ocean ridge lines")
        return self.plot_feature(
            ax,
            self.ridges,
            feature_name="ridges",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("trenches"))
    def get_trenches(self, central_meridian=0.0, tessellate_degrees=1):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed trench lines."""
        return self.get_feature(
            self.trenches,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_topology_availability("trenches")
    @append_docstring(PLOT_DOCSTRING.format("trenches"))
    def plot_trenches(self, ax, color="black", **kwargs):
        """Plot reconstructed subduction trench polylines onto a standard map
        Projection.

        Notes
        -----
        The trench sections for plotting are accessed from the
        `PlotTopologies` object's `trenches` attribute. These `trenches`
        are reconstructed to the `time` passed to the `PlotTopologies` object and converted
        into Shapely polylines. The reconstructed `trenches` are plotted onto the
        GeoAxes or GeoAxesSubplot map `ax` using GeoPandas. Map presentation details
        (e.g. `facecolor`, `edgecolor`, `alpha`…) are permitted as keyword arguments.

        Trench geometries are wrapped to the dateline using
        pyGPlates' [DateLineWrapper](https://www.gplates.org/docs/pygplates/generated/pygplates.datelinewrapper)
        by splitting a polyline into multiple polylines at the dateline. This is to avoid
        horizontal lines being formed between polylines at longitudes of -180 and 180 degrees.
        Point features near the poles (-89 & 89 degree latitude) are also clipped to ensure
        compatibility with Cartopy.
        """
        return self.plot_feature(
            ax,
            self.trenches,
            feature_name="trenches",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("other"))
    def get_misc_boundaries(self, central_meridian=0.0, tessellate_degrees=1):
        """Create a geopandas.GeoDataFrame object containing geometries of other reconstructed lines."""
        return self.get_feature(
            self.other,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("other"))
    def plot_misc_boundaries(self, ax, color="black", **kwargs):
        """Plot reconstructed miscellaneous plate boundary polylines onto a standard
        map Projection.

        Notes
        -----
        The miscellaneous boundary sections for plotting are accessed from the
        `PlotTopologies` object's `other` attribute. These `other` boundaries
        are reconstructed to the `time` passed to the `PlotTopologies` object and converted
        into Shapely polylines. The reconstructed `other` boundaries are plotted onto the
        GeoAxes or GeoAxesSubplot map `ax` using GeoPandas. Map presentation details
        (e.g. `facecolor`, `edgecolor`, `alpha`…) are permitted as keyword arguments.

        Miscellaneous boundary geometries are wrapped to the dateline using
        pyGPlates' [DateLineWrapper](https://www.gplates.org/docs/pygplates/generated/pygplates.datelinewrapper)
        by splitting a polyline into multiple polylines at the dateline. This is to avoid
        horizontal lines being formed between polylines at longitudes of -180 and 180 degrees.
        Point features near the poles (-89 & 89 degree latitude) are also clipped to ensure
        compatibility with Cartopy.
        """
        return self.plot_feature(
            ax,
            self.other,
            feature_name="misc_boundaries",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    def plot_subduction_teeth_deprecated(
        self, ax, spacing=0.1, size=2.0, aspect=1, color="black", **kwargs
    ):
        """Plot subduction teeth onto a standard map Projection.

        Notes
        -----
        Subduction teeth are tessellated from `PlotTopologies` object attributes `trench_left` and
        `trench_right`, and transformed into Shapely polygons for plotting.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        spacing : float, default=0.1
            The tessellation threshold (in radians). Parametrises subduction tooth density.
            Triangles are generated only along line segments with distances that exceed
            the given threshold ‘spacing’.

        size : float, default=2.0
            Length of teeth triangle base.

        aspect : float, default=1
            Aspect ratio of teeth triangles. Ratio is 1.0 by default.

        color : str, default=’black’
            The colour of the teeth. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as ‘alpha’, etc. for
            plotting subduction tooth polygons.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with subduction teeth plotted onto the chosen map projection.
        """
        import shapely

        # add Subduction Teeth
        subd_xL, subd_yL = self._tessellate_triangles(
            self.trench_left,
            tesselation_radians=spacing,
            triangle_base_length=size,
            triangle_aspect=-aspect,
        )
        subd_xR, subd_yR = self._tessellate_triangles(
            self.trench_right,
            tesselation_radians=spacing,
            triangle_base_length=size,
            triangle_aspect=aspect,
        )

        teeth = []
        for tX, tY in zip(subd_xL, subd_yL):
            triangle_xy_points = np.c_[tX, tY]
            shp = shapely.geometry.Polygon(triangle_xy_points)
            teeth.append(shp)

        for tX, tY in zip(subd_xR, subd_yR):
            triangle_xy_points = np.c_[tX, tY]
            shp = shapely.geometry.Polygon(triangle_xy_points)
            teeth.append(shp)

        return ax.add_geometries(teeth, crs=self.base_projection, color=color, **kwargs)

    def get_subduction_direction(self):
        """Create a geopandas.GeoDataFrame object containing geometries of trench directions.

        Notes
        -----
        The `trench_left` and `trench_right` geometries needed to produce the GeoDataFrame are automatically
        constructed if the optional `time` parameter is passed to the `PlotTopologies` object before calling
        this function. `time` can be passed either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `other` geometries are reconstructed, they are
        converted into Shapely features whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf_left : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `trench_left` geometry.
        gdf_right : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `trench_right` geometry.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `trench_left` or `trench_right` geometries to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No miscellaneous topologies have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.trench_left is None or self.trench_right is None:
            raise ValueError(
                "No trench_left or trench_right topologies passed to PlotTopologies."
            )

        trench_left_features = shapelify_feature_lines(self.trench_left)
        trench_right_features = shapelify_feature_lines(self.trench_right)

        gdf_left = gpd.GeoDataFrame(
            {"geometry": trench_left_features}, geometry="geometry"
        )
        gdf_right = gpd.GeoDataFrame(
            {"geometry": trench_right_features}, geometry="geometry"
        )

        return gdf_left, gdf_right

    def plot_subduction_teeth(
        self, ax, spacing=0.07, size=None, aspect=None, color="black", **kwargs
    ):
        """Plot subduction teeth onto a standard map Projection.

        Notes
        -----
        Subduction teeth are tessellated from `PlotTopologies` object attributes `trench_left` and
        `trench_right`, and transformed into Shapely polygons for plotting.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        spacing : float, default=0.07
            The tessellation threshold (in radians). Parametrises subduction tooth density.
            Triangles are generated only along line segments with distances that exceed
            the given threshold `spacing`.

        size : float, default=None
            Length of teeth triangle base (in radians). If kept at `None`, then
            `size = 0.5*spacing`.

        aspect : float, default=None
            Aspect ratio of teeth triangles. If kept at `None`, then `aspect = 2/3*size`.

        color : str, default='black'
            The colour of the teeth. By default, it is set to black.

        **kwargs :
            Keyword arguments parameters such as `alpha`, etc.
            for plotting subduction tooth polygons.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).
        """
        if not self.plate_reconstruction.topology_features:
            logger.warn(
                "Plate model does not have topology features. Unable to plot_subduction_teeth."
            )
            return ax
        if self._time is None:
            raise ValueError(
                "No topologies have been resolved. Set `PlotTopologies.time` to construct them."
            )
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")

        central_meridian = _meridian_from_ax(ax)
        tessellate_degrees = np.rad2deg(spacing)

        try:
            projection = ax.projection
        except AttributeError:
            print(
                "The ax.projection does not exist. You must set projection to plot Cartopy maps, such as ax = plt.subplot(211, projection=cartopy.crs.PlateCarree())"
            )
            projection = None

        if isinstance(projection, ccrs.PlateCarree):
            spacing = math.degrees(spacing)
        else:
            spacing = spacing * EARTH_RADIUS * 1e3

        if aspect is None:
            aspect = 2.0 / 3.0
        if size is None:
            size = spacing * 0.5

        height = size * aspect

        trench_left_features = shapelify_feature_lines(
            self.trench_left,
            tessellate_degrees=tessellate_degrees,
            central_meridian=central_meridian,
        )
        trench_right_features = shapelify_feature_lines(
            self.trench_right,
            tessellate_degrees=tessellate_degrees,
            central_meridian=central_meridian,
        )

        plot_subduction_teeth(
            trench_left_features,
            size,
            "l",
            height,
            spacing,
            projection=projection,
            ax=ax,
            color=color,
            **kwargs,
        )
        plot_subduction_teeth(
            trench_right_features,
            size,
            "r",
            height,
            spacing,
            projection=projection,
            ax=ax,
            color=color,
            **kwargs,
        )

    def plot_plate_polygon_by_id(self, ax, plate_id, color="black", **kwargs):
        """Plot a plate polygon with an associated `plate_id` onto a standard map Projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        plate_id : int
            A plate ID that identifies the continental polygon to plot. See the
            [Global EarthByte plate IDs list](https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/SampleData/FeatureCollections/Rotations/Global_EarthByte_PlateIDs_20071218.pdf)
            for a full list of plate IDs to plot.

        **kwargs :
            Keyword arguments for map presentation parameters such as
            `alpha`, etc. for plotting the grid.
            See `Matplotlib`'s `imshow` keyword arguments
            [here](https://matplotlib.org/3.5.1/api/_as_gen/matplotlib.axes.Axes.imshow.html).

        """
        return self.plot_feature(
            ax,
            [
                feature
                for feature in self.topologies
                if feature.get_reconstruction_plate_id() == plate_id
            ],
            color=color,
            **kwargs,
        )

    # the old function name(plot_plate_id) is bad. we should change the name
    # for backward compatibility, we have to allow users to use the old name
    def plot_plate_id(self, *args, **kwargs):
        logger.warn(
            "The class method plot_plate_id will be deprecated in the future soon. Use plot_plate_polygon_by_id instead."
        )
        return self.plot_plate_polygon_by_id(*args, **kwargs)

    plot_plate_id.__doc__ = plot_plate_polygon_by_id.__doc__

    def plot_grid(self, ax, grid, extent=[-180, 180, -90, 90], **kwargs):
        """Plot a `MaskedArray` raster or grid onto a standard map Projection.

        Notes
        -----
        Uses Matplotlib's `imshow`
        [function](https://matplotlib.org/3.5.1/api/_as_gen/matplotlib.axes.Axes.imshow.html).

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        grid : MaskedArray or `gplately.grids.Raster`
            A `MaskedArray` with elements that define a grid. The number of rows in the raster
            corresponds to the number of latitudinal coordinates, while the number of raster
            columns corresponds to the number of longitudinal coordinates.

        extent : 1d array, default=[-180,180,-90,90]
            A four-element array to specify the [min lon, max lon, min lat, max lat] with
            which to constrain the grid image. If no extents are supplied, full global
            extent is assumed.

        **kwargs :
            Keyword arguments for map presentation parameters such as
            `alpha`, etc. for plotting the grid.
            See `Matplotlib`'s `imshow` keyword arguments
            [here](https://matplotlib.org/3.5.1/api/_as_gen/matplotlib.axes.Axes.imshow.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with the grid plotted onto the chosen map projection.
        """
        # Override matplotlib default origin ('upper')
        origin = kwargs.pop("origin", "lower")

        from .grids import Raster

        if isinstance(grid, Raster):
            # extract extent and origin
            extent = grid.extent
            origin = grid.origin
            data = grid.data
        else:
            data = grid

        return ax.imshow(
            data,
            extent=extent,
            transform=self.base_projection,
            origin=origin,
            **kwargs,
        )

    def plot_grid_from_netCDF(self, ax, filename, **kwargs):
        """Read a raster from a netCDF file, convert it to a `MaskedArray` and plot it
        onto a standard map Projection.

        Notes
        -----
        `plot_grid_from_netCDF` uses Matplotlib's `imshow`
        [function](https://matplotlib.org/3.5.1/api/_as_gen/matplotlib.axes.Axes.imshow.html).

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        filename : str
            Full path to a netCDF filename.

        **kwargs :
            Keyword arguments for map presentation parameters for
            plotting the grid. See `Matplotlib`'s `imshow` keyword arguments
            [here](https://matplotlib.org/3.5.1/api/_as_gen/matplotlib.axes.Axes.imshow.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with the netCDF grid plotted onto the chosen map projection.
        """
        # Override matplotlib default origin ('upper')
        origin = kwargs.pop("origin", "lower")

        from .grids import read_netcdf_grid

        raster, lon_coords, lat_coords = read_netcdf_grid(filename, return_grids=True)
        extent = [lon_coords[0], lon_coords[-1], lat_coords[0], lat_coords[-1]]

        if lon_coords[0] < lat_coords[-1]:
            origin = "lower"
        else:
            origin = "upper"

        return self.plot_grid(ax, raster, extent=extent, origin=origin, **kwargs)

    def plot_plate_motion_vectors(
        self, ax, spacingX=10, spacingY=10, normalise=False, **kwargs
    ):
        """Calculate plate motion velocity vector fields at a particular geological time
        and plot them onto a standard map Projection.

        Notes
        -----
        `plot_plate_motion_vectors` generates a MeshNode domain of point features using
        given spacings in the X and Y directions (`spacingX` and `spacingY`). Each point in
        the domain is assigned a plate ID, and these IDs are used to obtain equivalent stage
        rotations of identified tectonic plates over a 5 Ma time interval. Each point and
        its stage rotation are used to calculate plate velocities at a particular geological
        time. Velocities for each domain point are represented in the north-east-down
        coordinate system and plotted on a GeoAxes.

        Vector fields can be optionally normalised by setting `normalise` to `True`. This
        makes vector arrow lengths uniform.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        spacingX : int, default=10
            The spacing in the X direction used to make the velocity domain point feature mesh.

        spacingY : int, default=10
            The spacing in the Y direction used to make the velocity domain point feature mesh.

        normalise : bool, default=False
            Choose whether to normalise the velocity magnitudes so that vector lengths are
            all equal.

        **kwargs :
            Keyword arguments for quiver presentation parameters for plotting
            the velocity vector field. See `Matplotlib` quiver keyword arguments
            [here](https://matplotlib.org/3.5.1/api/_as_gen/matplotlib.axes.Axes.quiver.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with the velocity vector field plotted onto the chosen map projection.
        """

        lons = np.arange(-180, 180 + spacingX, spacingX)
        lats = np.arange(-90, 90 + spacingY, spacingY)
        lonq, latq = np.meshgrid(lons, lats)

        # create a feature from all the points
        velocity_domain_features = ptt.velocity_tools.make_GPML_velocity_feature(
            lonq.ravel(), latq.ravel()
        )

        rotation_model = self.plate_reconstruction.rotation_model
        topology_features = self.plate_reconstruction.topology_features

        delta_time = 5.0
        all_velocities = ptt.velocity_tools.get_plate_velocities(
            velocity_domain_features,
            topology_features,
            rotation_model,
            self.time,
            delta_time,
            "vector_comp",
        )

        X, Y, U, V = ptt.velocity_tools.get_x_y_u_v(lons, lats, all_velocities)

        if normalise:
            mag = np.hypot(U, V)
            mag[mag == 0] = 1
            U /= mag
            V /= mag

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            quiver = ax.quiver(X, Y, U, V, transform=self.base_projection, **kwargs)
        return quiver

    def plot_pole(self, ax, lon, lat, a95, **kwargs):
        """
        Plot pole onto a matplotlib axes.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        lon : float
            Longitudinal coordinate to place pole
        lat : float
            Latitudinal coordinate to place pole
        a95 : float
            The size of the pole (in degrees)

        Returns
        -------
        matplotlib.patches.Circle handle
        """
        from matplotlib import patches as mpatches

        # Define the projection used to display the circle:
        proj1 = ccrs.Orthographic(central_longitude=lon, central_latitude=lat)

        def compute_radius(ortho, radius_degrees):
            phi1 = lat + radius_degrees if lat <= 0 else lat - radius_degrees
            _, y1 = ortho.transform_point(lon, phi1, ccrs.PlateCarree())
            return abs(y1)

        r_ortho = compute_radius(proj1, a95)

        # adding a patch
        patch = ax.add_patch(
            mpatches.Circle(xy=[lon, lat], radius=r_ortho, transform=proj1, **kwargs)
        )
        return patch

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("continental rifts"))
    def get_continental_rifts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed contiental rift lines."""
        return self.get_feature(
            self.continental_rifts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("continental rifts"))
    def plot_continental_rifts(self, ax, color="black", **kwargs):
        """Plot continental rifts on a standard map projection."""
        return self.plot_feature(
            ax,
            self.continental_rifts,
            feature_name="continental_rifts",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("faults"))
    def get_faults(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed fault lines."""
        return self.get_feature(
            self.faults,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("faults"))
    def plot_faults(self, ax, color="black", **kwargs):
        """Plot faults on a standard map projection."""
        return self.plot_feature(
            ax,
            self.faults,
            feature_name="faults",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("fracture zones"))
    def get_fracture_zones(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed fracture zone lines."""
        return self.get_feature(
            self.fracture_zones,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("fracturezones"))
    def plot_fracture_zones(self, ax, color="black", **kwargs):
        """Plot fracture zones on a standard map projection."""
        return self.plot_feature(
            ax,
            self.fracture_zones,
            feature_name="fracture_zones",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("inferred paleo-boundaries"))
    def get_inferred_paleo_boundaries(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed inferred paleo boundary lines."""
        return self.get_feature(
            self.inferred_paleo_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("inferred paleo-boundaries"))
    def plot_inferred_paleo_boundaries(self, ax, color="black", **kwargs):
        """Plot inferred paleo boundaries on a standard map projection."""
        return self.plot_feature(
            ax,
            self.inferred_paleo_boundaries,
            feature_name="inferred_paleo_boundaries",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("terrane boundaries"))
    def get_terrane_boundaries(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed terrane boundary lines."""
        return self.get_feature(
            self.terrane_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("terrane boundaries"))
    def plot_terrane_boundaries(self, ax, color="black", **kwargs):
        """Plot terrane boundaries on a standard map projection."""
        return self.plot_feature(
            ax,
            self.terrane_boundaries,
            feature_name="terrane_boundaries",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("transitional crusts"))
    def get_transitional_crusts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed transitional crust lines."""
        return self.get_feature(
            self.transitional_crusts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("transitional crusts"))
    def plot_transitional_crusts(self, ax, color="black", **kwargs):
        """Plot transitional crust on a standard map projection."""
        return self.plot_feature(
            ax,
            self.transitional_crusts,
            feature_name="transitional_crusts",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("orogenic belts"))
    def get_orogenic_belts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed orogenic belt lines."""
        return self.get_feature(
            self.orogenic_belts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("orogenic belts"))
    def plot_orogenic_belts(self, ax, color="black", **kwargs):
        """Plot orogenic belts on a standard map projection."""
        return self.plot_feature(
            ax,
            self.orogenic_belts,
            feature_name="orogenic_belts",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("sutures"))
    def get_sutures(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed suture lines."""
        return self.get_feature(
            self.sutures,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("sutures"))
    def plot_sutures(self, ax, color="black", **kwargs):
        """Plot sutures on a standard map projection."""
        return self.plot_feature(
            ax,
            self.sutures,
            feature_name="sutures",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("continental crusts"))
    def get_continental_crusts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed continental crust lines."""
        return self.get_feature(
            self.continental_crusts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("continental crusts"))
    def plot_continental_crusts(self, ax, color="black", **kwargs):
        """Plot continental crust lines on a standard map projection."""
        return self.plot_feature(
            ax,
            self.continental_crusts,
            feature_name="continental_crusts",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("extended continental crusts"))
    def get_extended_continental_crusts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed extended continental crust lines."""
        return self.get_feature(
            self.extended_continental_crusts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("extended continental crusts"))
    def plot_extended_continental_crusts(self, ax, color="black", **kwargs):
        """Plot extended continental crust lines on a standard map projection."""
        return self.plot_feature(
            ax,
            self.extended_continental_crusts,
            feature_name="extended_continental_crusts",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("passive continental boundaries"))
    def get_passive_continental_boundaries(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed passive continental boundary lines."""
        return self.get_feature(
            self.passive_continental_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("passive continental boundaries"))
    def plot_passive_continental_boundaries(self, ax, color="black", **kwargs):
        """Plot passive continental boundaries on a standard map projection."""
        return self.plot_feature(
            ax,
            self.passive_continental_boundaries,
            feature_name="passive_continental_boundaries",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("slab edges"))
    def get_slab_edges(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed slab edge lines."""
        return self.get_feature(
            self.slab_edges,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("slab edges"))
    def plot_slab_edges(self, ax, color="black", **kwargs):
        """Plot slab edges on a standard map projection."""
        return self.plot_feature(
            ax,
            self.slab_edges,
            feature_name="slab_edges",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("transforms"))
    def get_transforms(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed transform lines(gpml:Transform)."""
        logger.debug("Getting reconstructed transform boundary lines")
        return self.get_feature(
            self.transforms,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("transforms"))
    def plot_transforms(self, ax, color="black", **kwargs):
        """Plot transform boundaries(gpml:Transform) onto a map."""
        logger.debug("Plotting transform")
        return self.plot_feature(
            ax,
            self.transforms,
            feature_name="transforms",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("unclassified features"))
    def get_unclassified_features(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed unclassified feature lines."""
        return self.get_feature(
            self.unclassified_features,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING.format("unclassified features"))
    def plot_unclassified_features(self, ax, color="black", **kwargs):
        """Plot GPML unclassified features on a standard map projection."""
        return self.plot_feature(
            ax,
            self.unclassified_features,
            feature_name="unclassified_features",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("topologies"))
    def get_all_topologies(
        self,
        central_meridian=0.0,
        tessellate_degrees=1,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed unclassified feature lines."""

        # get plate IDs and feature types to add to geodataframe
        plate_IDs = []
        feature_types = []
        feature_names = []
        all_topologies = []

        if self.topologies:
            all_topologies = shapelify_features(
                self.topologies,
                central_meridian=central_meridian,
                tessellate_degrees=tessellate_degrees,
            )
            for topo in self.topologies:
                ft_type = topo.get_feature_type()

                plate_IDs.append(topo.get_reconstruction_plate_id())
                feature_types.append(ft_type)
                feature_names.append(ft_type.get_name())

        gdf = gpd.GeoDataFrame(
            {
                "geometry": all_topologies,
                "reconstruction_plate_ID": plate_IDs,
                "feature_type": feature_types,
                "feature_name": feature_names,
            },
            geometry="geometry",
        )
        return gdf

    @validate_topology_availability("all topologies")
    @append_docstring(PLOT_DOCSTRING.format("topologies"))
    def plot_all_topologies(self, ax, color="black", **kwargs):
        """Plot topological polygons and networks on a standard map projection."""
        if "edgecolor" not in kwargs.keys():
            kwargs["edgecolor"] = color

        return self._plot_feature(
            ax,
            self.get_all_topologies,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("topologies"))
    def get_all_topological_sections(
        self,
        central_meridian=0.0,
        tessellate_degrees=1,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries ofresolved topological sections."""

        # get plate IDs and feature types to add to geodataframe
        geometries = []
        plate_IDs = []
        feature_types = []
        feature_names = []
        for topo in [
            *self.ridges,
            *self.trenches,
            *self.trench_left,
            *self.trench_right,
            *self.other,
        ]:
            converted = shapelify_features(
                topo,
                central_meridian=central_meridian,
                tessellate_degrees=tessellate_degrees,
            )
            if not isinstance(converted, BaseGeometry):
                if len(converted) > 1:
                    tmp = []
                    for i in converted:
                        if isinstance(i, BaseMultipartGeometry):
                            tmp.extend(list(i.geoms))
                        else:
                            tmp.append(i)
                    converted = tmp
                    del tmp
                    converted = linemerge(converted)
                elif len(converted) == 1:
                    converted = converted[0]
                else:
                    continue
            geometries.append(converted)
            plate_IDs.append(topo.get_reconstruction_plate_id())
            feature_types.append(topo.get_feature_type())
            feature_names.append(topo.get_name())

        gdf = gpd.GeoDataFrame(
            {
                "geometry": geometries,
                "reconstruction_plate_ID": plate_IDs,
                "feature_type": feature_types,
                "feature_name": feature_names,
            },
            geometry="geometry",
        )
        return gdf

    @validate_topology_availability("all topological sections")
    @append_docstring(PLOT_DOCSTRING.format("topologies"))
    def plot_all_topological_sections(self, ax, color="black", **kwargs):
        """Plot all topologies on a standard map projection."""

        return self._plot_feature(
            ax,
            self.get_all_topological_sections,
            color=color,
            **kwargs,
        )

    def _resolve_both_boundaries_and_networks(self):
        """need to resolve_topologies for both rigid boundaries and networks if not done yet"""
        if self._topologies is None:
            resolved_topologies = []
            shared_boundary_sections = []
            pygplates.resolve_topologies(
                self.plate_reconstruction.topology_features,
                self.plate_reconstruction.rotation_model,
                resolved_topologies,
                self.time,
                shared_boundary_sections,
                anchor_plate_id=self.anchor_plate_id,
            )
            self._topologies = [t.get_resolved_feature() for t in resolved_topologies]

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("topological plate boundaries"))
    def get_topological_plate_boundaries(
        self, central_meridian=0.0, tessellate_degrees=1
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed rigid topological plate boundaries."""
        return self.get_feature(
            self._topological_plate_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_topology_availability("topological plate boundaries")
    @append_docstring(PLOT_DOCSTRING.format("topological plate boundaries"))
    def plot_topological_plate_boundaries(self, ax, color="black", **kwargs):
        return self.plot_feature(
            ax,
            self._topological_plate_boundaries,
            feature_name="topological plate boundaries",
            color=color,
            **kwargs,
        )

    @property
    @append_docstring(GET_DATE_DOCSTRING.format("ridge transforms"))
    def ridge_transforms(self):
        """Combine ridges and transforms into one feature set."""
        if self._topological_plate_boundaries is None:
            logger.debug("Accessing ridge transforms - updating topological plate boundaries.")
            self.update_time(self.time)
        logger.debug("Returning ridge_transforms property.")
        return self.ridges + self.transforms

    def get_ridge_transforms(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a GeoDataFrame containing geometries of reconstructed ridge transforms."""
        logger.debug("Getting reconstructed ridge transforms.")
        return self.get_feature(
            self.ridge_transforms,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees
        )

    @validate_topology_availability("ridge transforms")
    @append_docstring(PLOT_DOCSTRING.format("ridge transforms"))
    def plot_ridge_transforms(self, ax, color="black", **kwargs):
        """Plot reconstructed ridge transforms onto a map."""
        logger.debug("Plotting reconstructed ridge transforms.")
        return self.plot_feature(
            ax,
            self.ridge_transforms,
            feature_name="ridge_transforms",
            edgecolor=color,
            **kwargs,
        )