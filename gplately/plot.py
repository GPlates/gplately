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

import cartopy.crs as ccrs
import geopandas as gpd
import numpy as np
import pygplates
from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
from shapely.ops import linemerge

from . import ptt
from .decorators import append_docstring
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
            The colour of the ridge lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc. for
            plotting ridge geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with ridge features plotted onto the chosen map projection.
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

    ridge_transforms : iterable/list of <pygplates.Feature>
        A list containing ridge and transform boundary sections of type
        pygplates.FeatureType.gpml_mid_ocean_ridge

    ridges : iterable/list of <pygplates.Feature>
        A list containing ridge boundary sections of type pygplates.FeatureType.gpml_mid_ocean_ridge

    transforms : iterable/list of <pygplates.Feature>
        A list containing transform boundary sections of type pygplates.FeatureType.gpml_mid_ocean_ridge

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
        if var == self.time:
            pass
        elif var >= 0:
            self.update_time(var)
        else:
            raise ValueError("Enter a valid time >= 0")

    @property
    def anchor_plate_id(self):
        """Anchor plate ID for reconstruction. Must be an integer >= 0."""
        return self._anchor_plate_id

    @anchor_plate_id.setter
    def anchor_plate_id(self, anchor_plate):
        self._anchor_plate_id = self._check_anchor_plate_id(anchor_plate)
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
        self._time = float(time)
        resolved_topologies = ptt.resolve_topologies.resolve_topologies_into_features(
            self.plate_reconstruction.rotation_model,
            self.plate_reconstruction.topology_features,
            self.time,
            anchor_plate_id=self.anchor_plate_id,
        )

        (
            self.topologies,
            self.ridge_transforms,
            self.ridges,
            self.transforms,
            self.trenches,
            self.trench_left,
            self.trench_right,
            self.other,
        ) = resolved_topologies

        self.ridges, self.transforms = (
            ptt.separate_ridge_transform_segments.separate_features_into_ridges_and_transforms(
                self.plate_reconstruction.rotation_model, self.ridge_transforms
            )
        )

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
        self.misc_transforms = []
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
                self.misc_transforms.append(topol)

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
        if validate_reconstruction_time and self._time is None:
            raise ValueError(
                "The reconstruction time has not been set yet. Set `PlotTopologies.time` before calling plotting functions."
            )
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

    def plot_feature(self, ax, feature, **kwargs):
        """Plot pygplates.FeatureCollection  or pygplates.Feature onto a map.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        **kwargs :
            Keyword arguments for parameters such as `facecolor`, `alpha`,
            etc. for plotting coastline geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with coastline features plotted onto the chosen map projection.
        """
        return self._plot_feature(ax, feature=feature, **kwargs)

    def _plot_feature(self, ax, feature=None, get_feature_func=None, **kwargs):
        if feature and get_feature_func:
            logger.warn(
                "Both 'feature' and 'get_feature_func' parameters are not None. Use 'feature' and ignore 'get_feature_func'."
            )
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
        if feature:
            gdf = self.get_feature(
                feature,
                central_meridian=central_meridian,
                tessellate_degrees=tessellate_degrees,
            )
        elif get_feature_func:
            gdf = get_feature_func(
                central_meridian=central_meridian,
                tessellate_degrees=tessellate_degrees,
            )
        else:
            raise Exception(
                "The caller must provide either a 'feature' or 'get_feature_func' parameter. Unable to plot the feature if both parameters are None."
            )

        if len(gdf) == 0:
            logger.warn("No feature found for plotting. Do nothing and return.")
            return ax

        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection

        return gdf.plot(ax=ax, **kwargs)

    def get_coastlines(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed coastline polygons.

        Notes
        -----
        The `coastlines` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `coastlines` are reconstructed, they are
        converted into Shapely polygons whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `coastlines` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `coastlines` to the requested `time` and thus populate the GeoDataFrame.

        """
        if self._time is None:
            raise ValueError(
                "No coastlines have been resolved. Set `PlotTopologies.time` to construct coastlines."
            )

        if self.coastlines is None:
            raise ValueError("Supply coastlines to PlotTopologies object")

        coastline_polygons = shapelify_feature_polygons(
            self.coastlines,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": coastline_polygons}, geometry="geometry")
        return gdf

    def plot_coastlines(self, ax, **kwargs):
        """Plot reconstructed coastline polygons onto a standard map Projection.

        Notes
        -----
        The `coastlines` for plotting are accessed from the `PlotTopologies` object's
        `coastlines` attribute. These `coastlines` are reconstructed to the `time`
        passed to the `PlotTopologies` object and converted into Shapely polylines. The
        reconstructed `coastlines` are added onto the GeoAxes or GeoAxesSubplot map `ax` using
        GeoPandas.
        Map resentation details (e.g. facecolor, edgecolor, alpha…) are permitted as keyword
        arguments.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        **kwargs :
            Keyword arguments for parameters such as `facecolor`, `alpha`,
            etc. for plotting coastline geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with coastline features plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_coastlines(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, **kwargs)

    def get_continents(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed continental polygons.

        Notes
        -----
        The `continents` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `continents` are reconstructed, they are
        converted into Shapely polygons whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `continents` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `continents` to the requested `time` and thus populate the GeoDataFrame.

        """
        if self._time is None:
            raise ValueError(
                "No continents have been resolved. Set `PlotTopologies.time` to construct continents."
            )

        if self.continents is None:
            raise ValueError("Supply continents to PlotTopologies object")

        continent_polygons = shapelify_feature_polygons(
            self.continents,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": continent_polygons}, geometry="geometry")
        return gdf

    def plot_continents(self, ax, **kwargs):
        """Plot reconstructed continental polygons onto a standard map Projection.

        Notes
        -----
        The `continents` for plotting are accessed from the `PlotTopologies` object's
        `continents` attribute. These `continents` are reconstructed to the `time`
        passed to the `PlotTopologies` object and converted into Shapely polygons.
        The reconstructed `coastlines` are plotted onto the GeoAxes or GeoAxesSubplot map `ax` using
        GeoPandas.
        Map presentation details (e.g. facecolor, edgecolor, alpha…) are permitted as
        keyword arguments.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        **kwargs :
            Keyword arguments for parameters such as `facecolor`, `alpha`,
            etc. for plotting continental geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with continent features plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_continents(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, **kwargs)

    def get_continent_ocean_boundaries(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed continent-ocean
        boundary lines.

        Notes
        -----
        The `COBs` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `COBs` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `COBs` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `COBs` to the requested `time` and thus populate the GeoDataFrame.

        """
        if self._time is None:
            raise ValueError(
                "No geometries have been resolved. Set `PlotTopologies.time` to construct topologies."
            )

        if self.COBs is None:
            raise ValueError("Supply COBs to PlotTopologies object")

        COB_lines = shapelify_feature_lines(
            self.COBs,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": COB_lines}, geometry="geometry")
        return gdf

    def plot_continent_ocean_boundaries(self, ax, **kwargs):
        """Plot reconstructed continent-ocean boundary (COB) polygons onto a standard
        map Projection.

        Notes
        -----
        The `COBs` for plotting are accessed from the `PlotTopologies` object's
        `COBs` attribute. These `COBs` are reconstructed to the `time`
        passed to the `PlotTopologies` object and converted into Shapely polylines.
        The reconstructed `COBs` are plotted onto the GeoAxes or GeoAxesSubplot map
        `ax` using GeoPandas. Map presentation details (e.g. `facecolor`, `edgecolor`, `alpha`…)
        are permitted as keyword arguments.

        These COBs are transformed into shapely
        geometries and added onto the chosen map for a specific geological time (supplied to the
        PlotTopologies object). Map presentation details (e.g. facecolor, edgecolor, alpha…)
        are permitted.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        **kwargs :
            Keyword arguments for parameters such as `facecolor`, `alpha`,
            etc. for plotting COB geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with COB features plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_continent_ocean_boundaries(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, **kwargs)

    def get_ridges_and_transforms(
        self,
        central_meridian=0.0,
        tessellate_degrees=1,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed ridge and transform lines.

        Notes
        -----
        The `ridge_transforms` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `ridge_transforms` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `ridges` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `ridge_transforms` to the requested `time` and thus populate the GeoDataFrame.

        """
        if self._time is None:
            raise ValueError(
                "No ridges and transforms have been resolved. Set `PlotTopologies.time` to construct ridges and transforms."
            )

        if self.ridge_transforms is None:
            raise ValueError(
                "No ridge and transform topologies passed to PlotTopologies."
            )

        ridge_transform_lines = shapelify_feature_lines(
            self.ridge_transforms,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": ridge_transform_lines}, geometry="geometry")
        return gdf

    def plot_ridges_and_transforms(self, ax, color="black", **kwargs):
        """Plot reconstructed ridge & transform boundary polylines onto a standard map
        Projection.

        Notes
        -----
        The ridge & transform sections for plotting are accessed from the
        `PlotTopologies` object's `ridge_transforms` attribute. These `ridge_transforms`
        are reconstructed to the `time` passed to the `PlotTopologies` object and converted
        into Shapely polylines. The reconstructed `ridge_transforms` are plotted onto the
        GeoAxes or GeoAxesSubplot map `ax` using GeoPandas. Map presentation details
        (e.g. `facecolor`, `edgecolor`, `alpha`…) are permitted as keyword arguments.

        Note: Ridge & transform geometries are wrapped to the dateline using
        pyGPlates' [DateLineWrapper](https://www.gplates.org/docs/pygplates/generated/pygplates.datelinewrapper)
        by splitting a polyline into multiple polylines at the dateline. This is to avoid
        horizontal lines being formed between polylines at longitudes of -180 and 180 degrees.
        Point features near the poles (-89 & 89 degree latitude) are also clipped to ensure
        compatibility with Cartopy.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the ridge & transform lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as ‘alpha’, etc. for
            plotting ridge & transform geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with ridge & transform features plotted onto the chosen map projection.
        """
        if not self.plate_reconstruction.topology_features:
            logger.warn(
                "Plate model does not have topology features. Unable to plot_ridges_and_transforms."
            )
            return ax

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

        gdf = self.get_ridges_and_transforms(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_transforms(self, central_meridian=0.0, tessellate_degrees=1):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed transform lines.

        Notes
        -----
        The `transforms` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `transforms` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `transforms` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `transforms` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No transforms have been resolved. Set `PlotTopologies.time` to construct transforms."
            )

        if self.transforms is None:
            raise ValueError("No transform topologies passed to PlotTopologies.")

        transform_lines = shapelify_feature_lines(
            self.transforms,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": transform_lines}, geometry="geometry")
        return gdf

    def plot_transforms(self, ax, color="black", **kwargs):
        """Plot reconstructed transform boundary polylines onto a standard map.

        Notes
        -----
        The transform sections for plotting are accessed from the
        `PlotTopologies` object's `transforms` attribute. These `transforms`
        are reconstructed to the `time` passed to the `PlotTopologies` object and converted
        into Shapely polylines. The reconstructed `transforms` are plotted onto the
        GeoAxes or GeoAxesSubplot map `ax` using GeoPandas. Map presentation details
        (e.g. `facecolor`, `edgecolor`, `alpha`…) are permitted as keyword arguments.

        Transform geometries are wrapped to the dateline using
        pyGPlates' [DateLineWrapper](https://www.gplates.org/docs/pygplates/generated/pygplates.datelinewrapper)
        by splitting a polyline into multiple polylines at the dateline. This is to avoid
        horizontal lines being formed between polylines at longitudes of -180 and 180 degrees.
        Point features near the poles (-89 & 89 degree latitude) are also clipped to ensure
        compatibility with Cartopy.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the transform lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting transform geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with transform features plotted onto the chosen map projection.
        """
        if not self.plate_reconstruction.topology_features:
            logger.warn(
                "Plate model does not have topology features. Unable to plot_transforms."
            )
            return ax

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

        gdf = self.get_transforms(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_trenches(self, central_meridian=0.0, tessellate_degrees=1):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed trench lines.

        Notes
        -----
        The `trenches` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `trenches` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `trenches` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `trenches` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No trenches have been resolved. Set `PlotTopologies.time` to construct trenches."
            )

        if self.trenches is None:
            raise ValueError("No trenches passed to PlotTopologies.")

        trench_lines = shapelify_feature_lines(
            self.trenches,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": trench_lines}, geometry="geometry")
        return gdf

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

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with transform features plotted onto the chosen map projection.
        """
        if not self.plate_reconstruction.topology_features:
            logger.warn(
                "Plate model does not have topology features. Unable to plot_trenches."
            )
            return ax
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

        gdf = self.get_trenches(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_misc_boundaries(self, central_meridian=0.0, tessellate_degrees=1):
        """Create a geopandas.GeoDataFrame object containing geometries of other reconstructed lines.

        Notes
        -----
        The `other` geometries needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `other` geometries are reconstructed, they are
        converted into Shapely features whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `other` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `other` geometries to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No miscellaneous topologies have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.other is None:
            raise ValueError("No miscellaneous topologies passed to PlotTopologies.")

        lines = shapelify_features(
            self.other,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": lines}, geometry="geometry")
        return gdf

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

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the boundary lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as ‘alpha’, etc. for
            plotting miscellaneous boundary geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with miscellaneous boundary features plotted onto the chosen map projection.
        """
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

        gdf = self.get_misc_boundaries(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

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

    def plot_plate_polygon_by_id(self, ax, plate_id, **kwargs):
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
        tessellate_degrees = kwargs.pop("tessellate_degrees", 1)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        for feature in self.topologies:
            if feature.get_reconstruction_plate_id() == plate_id:
                ft_plate = shapelify_feature_polygons(
                    [feature],
                    central_meridian=central_meridian,
                    tessellate_degrees=tessellate_degrees,
                )
                return ax.add_geometries(ft_plate, crs=self.base_projection, **kwargs)

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

    def get_continental_rifts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed contiental rift lines.

        Notes
        -----
        The `continental_rifts` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `continental_rifts` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `continental_rifts` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `continental_rifts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No continental rifts have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.continental_rifts is None:
            raise ValueError("No continental rifts passed to PlotTopologies.")

        continental_rift_lines = shapelify_feature_lines(
            self.continental_rifts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame(
            {"geometry": continental_rift_lines}, geometry="geometry"
        )
        return gdf

    def plot_continental_rifts(self, ax, color="black", **kwargs):
        """Plot continental rifts on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with continental rifts plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_continental_rifts(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_faults(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed fault lines.

        Notes
        -----
        The `faults` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `faults` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `faults` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `faults` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No faults have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.faults is None:
            raise ValueError("No faults passed to PlotTopologies.")

        fault_lines = shapelify_feature_lines(
            self.faults,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": fault_lines}, geometry="geometry")
        return gdf

    def plot_faults(self, ax, color="black", **kwargs):
        """Plot faults on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with faults plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_faults(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_fracture_zones(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed fracture zone lines.

        Notes
        -----
        The `fracture_zones` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `fracture_zones` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `fracture_zones` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `fracture_zones` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No fracture zones have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.fracture_zones is None:
            raise ValueError("No fracture zones passed to PlotTopologies.")

        fracture_zone_lines = shapelify_feature_lines(
            self.fracture_zones,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": fracture_zone_lines}, geometry="geometry")
        return gdf

    def plot_fracture_zones(self, ax, color="black", **kwargs):
        """Plot fracture zones on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with fracture zones plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_fracture_zones(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_inferred_paleo_boundaries(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed inferred paleo boundary lines.

        Notes
        -----
        The `inferred_paleo_boundaries` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `inferred_paleo_boundaries` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `inferred_paleo_boundaries` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `inferred_paleo_boundaries` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No inferred paleo boundaries have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.inferred_paleo_boundaries is None:
            raise ValueError("No inferred paleo boundaries passed to PlotTopologies.")

        inferred_paleo_boundary_lines = shapelify_feature_lines(
            self.inferred_paleo_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame(
            {"geometry": inferred_paleo_boundary_lines}, geometry="geometry"
        )
        return gdf

    def plot_inferred_paleo_boundaries(self, ax, color="black", **kwargs):
        """Plot inferred paleo boundaries on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with inferred paleo boundaries plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_inferred_paleo_boundaries(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_terrane_boundaries(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed terrane boundary lines.

        Notes
        -----
        The `terrane_boundaries` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `terrane_boundaries` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `terrane_boundaries` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `terrane_boundaries` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No terrane boundaries have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.terrane_boundaries is None:
            raise ValueError("No terrane boundaries passed to PlotTopologies.")

        terrane_boundary_lines = shapelify_feature_lines(
            self.terrane_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame(
            {"geometry": terrane_boundary_lines}, geometry="geometry"
        )
        return gdf

    def plot_terrane_boundaries(self, ax, color="black", **kwargs):
        """Plot terrane boundaries on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with terrane boundaries plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_terrane_boundaries(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_transitional_crusts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed transitional crust lines.

        Notes
        -----
        The `transitional_crusts` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `transitional_crusts` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `transitional_crusts` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `transitional_crusts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No transitional crusts have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.transitional_crusts is None:
            raise ValueError("No transitional crusts passed to PlotTopologies.")

        transitional_crust_lines = shapelify_feature_lines(
            self.transitional_crusts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame(
            {"geometry": transitional_crust_lines}, geometry="geometry"
        )
        return gdf

    def plot_transitional_crusts(self, ax, color="black", **kwargs):
        """Plot transitional crust on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with transitional crust plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_transitional_crusts(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_orogenic_belts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed orogenic belt lines.

        Notes
        -----
        The `orogenic_belts` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `orogenic_belts` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `orogenic_belts` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `orogenic_belts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No orogenic belts have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.orogenic_belts is None:
            raise ValueError("No orogenic belts passed to PlotTopologies.")

        orogenic_belt_lines = shapelify_feature_lines(
            self.orogenic_belts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": orogenic_belt_lines}, geometry="geometry")
        return gdf

    def plot_orogenic_belts(self, ax, color="black", **kwargs):
        """Plot orogenic belts on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with orogenic belts plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_orogenic_belts(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_sutures(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed suture lines.

        Notes
        -----
        The `sutures` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `sutures` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `sutures` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `sutures` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No sutures have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.sutures is None:
            raise ValueError("No sutures passed to PlotTopologies.")

        suture_lines = shapelify_feature_lines(
            self.sutures,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": suture_lines}, geometry="geometry")
        return gdf

    def plot_sutures(self, ax, color="black", **kwargs):
        """Plot sutures on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with sutures plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_sutures(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_continental_crusts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed continental crust lines.

        Notes
        -----
        The `continental_crusts` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `continental_crusts` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `continental_crusts` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `continental_crusts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No continental crust topologies have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.continental_crusts is None:
            raise ValueError(
                "No continental crust topologies passed to PlotTopologies."
            )

        continental_crust_lines = shapelify_feature_lines(
            self.continental_crusts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame(
            {"geometry": continental_crust_lines}, geometry="geometry"
        )
        return gdf

    def plot_continental_crusts(self, ax, color="black", **kwargs):
        """Plot continental crust lines on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with continental crust lines plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_continental_crusts(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_extended_continental_crusts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed extended continental crust lines.

        Notes
        -----
        The `extended_continental_crusts` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `extended_continental_crusts` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `extended_continental_crusts` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `extended_continental_crusts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No extended continental crust topologies have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.extended_continental_crusts is None:
            raise ValueError(
                "No extended continental crust topologies passed to PlotTopologies."
            )

        extended_continental_crust_lines = shapelify_feature_lines(
            self.extended_continental_crusts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame(
            {"geometry": extended_continental_crust_lines}, geometry="geometry"
        )
        return gdf

    def plot_extended_continental_crusts(self, ax, color="black", **kwargs):
        """Plot extended continental crust lines on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with extended continental crust lines plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_extended_continental_crusts(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_passive_continental_boundaries(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed passive continental boundary lines.

        Notes
        -----
        The `passive_continental_boundaries` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `passive_continental_boundaries` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `passive_continental_boundaries` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `passive_continental_boundaries` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No passive continental boundaries have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.passive_continental_boundaries is None:
            raise ValueError(
                "No passive continental boundaries passed to PlotTopologies."
            )

        passive_continental_boundary_lines = shapelify_feature_lines(
            self.passive_continental_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame(
            {"geometry": passive_continental_boundary_lines}, geometry="geometry"
        )
        return gdf

    def plot_passive_continental_boundaries(self, ax, color="black", **kwargs):
        """Plot passive continental boundaries on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with passive continental boundaries plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_passive_continental_boundaries(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_slab_edges(self, central_meridian=0.0, tessellate_degrees=None):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed slab edge lines.

        Notes
        -----
        The `slab_edges` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `slab_edges` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `slab_edges` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `slab_edges` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No slab edges have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.slab_edges is None:
            raise ValueError("No slab edges passed to PlotTopologies.")

        slab_edge_lines = shapelify_feature_lines(
            self.slab_edges,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": slab_edge_lines}, geometry="geometry")
        return gdf

    def plot_slab_edges(self, ax, color="black", **kwargs):
        """Plot slab edges on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with slab edges plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_slab_edges(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_misc_transforms(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed misc transform lines.

        Notes
        -----
        The `misc_transforms` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `misc_transforms` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `misc_transforms` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `misc_transforms` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No miscellaneous transforms have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.misc_transforms is None:
            raise ValueError("No miscellaneous transforms passed to PlotTopologies.")

        misc_transform_lines = shapelify_feature_lines(
            self.misc_transforms,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame({"geometry": misc_transform_lines}, geometry="geometry")
        return gdf

    def plot_misc_transforms(self, ax, color="black", **kwargs):
        """Plot miscellaneous transform boundaries on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with miscellaneous transform boundaries plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_misc_transforms(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_unclassified_features(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed unclassified feature lines.

        Notes
        -----
        The `unclassified_features` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `unclassified_features` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `unclassified_features` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `unclassified_features` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No unclassified features have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.unclassified_features is None:
            raise ValueError("No unclassified features passed to PlotTopologies.")

        unclassified_feature_lines = shapelify_feature_lines(
            self.unclassified_features,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        gdf = gpd.GeoDataFrame(
            {"geometry": unclassified_feature_lines}, geometry="geometry"
        )
        return gdf

    def plot_unclassified_features(self, ax, color="black", **kwargs):
        """Plot GPML unclassified features on a standard map projection.

        Parameters
        ----------
        ax : instance of <cartopy.mpl.geoaxes.GeoAxes> or <cartopy.mpl.geoaxes.GeoAxesSubplot>
            A subclass of `matplotlib.axes.Axes` which represents a map Projection.
            The map should be set at a particular Cartopy projection.

        color : str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs :
            Keyword arguments for parameters such as `alpha`, etc.
            for plotting trench geometries.
            See `Matplotlib` keyword arguments
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map
            with unclassified features plotted onto the chosen map projection.
        """
        if "transform" in kwargs.keys():
            warnings.warn(
                "'transform' keyword argument is ignored by PlotTopologies",
                UserWarning,
            )
            kwargs.pop("transform")
        tessellate_degrees = kwargs.pop("tessellate_degrees", None)
        central_meridian = kwargs.pop("central_meridian", None)
        if central_meridian is None:
            central_meridian = _meridian_from_ax(ax)

        gdf = self.get_unclassified_features(
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )
        if hasattr(ax, "projection"):
            gdf = _clean_polygons(data=gdf, projection=ax.projection)
        else:
            kwargs["transform"] = self.base_projection
        return gdf.plot(ax=ax, facecolor="none", edgecolor=color, **kwargs)

    def get_all_topologies(
        self,
        central_meridian=0.0,
        tessellate_degrees=1,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed unclassified feature lines.

        Notes
        -----
        The `topologies` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `topologies` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `topologies` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `topologies` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No topologies have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.topologies is None:
            raise ValueError("No topologies passed to PlotTopologies.")

        all_topologies = shapelify_features(
            self.topologies,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

        # get plate IDs and feature types to add to geodataframe
        plate_IDs = []
        feature_types = []
        feature_names = []
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

    @append_docstring(PLOT_DOCSTRING)
    def plot_all_topologies(self, ax, color="black", **kwargs):
        """Plot all topologies on a standard map projection."""

        return self._plot_feature(
            ax,
            get_feature_func=self.get_all_topologies,
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    def get_all_topological_sections(
        self,
        central_meridian=0.0,
        tessellate_degrees=1,
    ):
        """Create a geopandas.GeoDataFrame object containing geometries of
        resolved topological sections.

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
            A pandas.DataFrame that has a column with `topologies` geometry.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to
            `PlotTopologies`. This is needed to construct `topologies`
            to the requested `time` and thus populate the GeoDataFrame.

        Notes
        -----
        The `topologies` needed to produce the GeoDataFrame are automatically
        constructed if the optional `time` parameter is passed to the
        `PlotTopologies` object before calling this function. `time` can be
        passed either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 # Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `topologies`
        are reconstructed, they are converted into Shapely lines whose
        coordinates are passed to a geopandas GeoDataFrame.
        """
        if self._time is None:
            raise ValueError(
                "No topologies have been resolved. Set `PlotTopologies.time` to construct them."
            )

        if self.topologies is None:
            raise ValueError("No topologies passed to PlotTopologies.")

        topologies_list = [
            *self.ridge_transforms,
            *self.ridges,
            *self.transforms,
            *self.trenches,
            *self.trench_left,
            *self.trench_right,
            *self.other,
        ]

        # get plate IDs and feature types to add to geodataframe
        geometries = []
        plate_IDs = []
        feature_types = []
        feature_names = []
        for topo in topologies_list:
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

    @append_docstring(PLOT_DOCSTRING)
    def plot_all_topological_sections(self, ax, color="black", **kwargs):
        """Plot all topologies on a standard map projection."""

        return self._plot_feature(
            ax,
            get_feature_func=self.get_all_topological_sections,
            color=color,
            **kwargs,
        )

    def get_ridges(self, central_meridian=0.0, tessellate_degrees=1):
        """Create a geopandas.GeoDataFrame object containing geometries of reconstructed ridge lines.

        Notes
        -----
        The `ridges` needed to produce the GeoDataFrame are automatically constructed if the optional `time`
        parameter is passed to the `PlotTopologies` object before calling this function. `time` can be passed
        either when `PlotTopologies` is first called...

            gplot = gplately.PlotTopologies(..., time=100,...)

        or anytime afterwards, by setting:

            time = 100 #Ma
            gplot.time = time

        ...after which this function can be re-run. Once the `ridges` are reconstructed, they are
        converted into Shapely lines whose coordinates are passed to a geopandas GeoDataFrame.

        Returns
        -------
        gdf : instance of <geopandas.GeoDataFrame>
            A pandas.DataFrame that has a column with `ridges` geometry.
        central_meridian : float
            Central meridian around which to perform wrapping; default: 0.0.
        tessellate_degrees : float or None
            If provided, geometries will be tessellated to this resolution prior
            to wrapping.

        Raises
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `ridges` to the requested `time` and thus populate the GeoDataFrame.

        """
        return self.get_feature(
            self.ridges,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @append_docstring(PLOT_DOCSTRING)
    def plot_ridges(self, ax, color="black", **kwargs):
        """Plot reconstructed ridge polylines onto a standard map Projection.

        Notes
        -----
        The `ridges` for plotting are accessed from the `PlotTopologies` object's
        `ridges` attribute. These `ridges` are reconstructed to the `time`
        passed to the `PlotTopologies` object and converted into Shapely polylines.
        The reconstructed `ridges` are plotted onto the GeoAxes or GeoAxesSubplot map
        `ax` using GeoPandas. Map presentation details (e.g. `facecolor`, `edgecolor`, `alpha`…)
        are permitted as keyword arguments.

        Ridge geometries are wrapped to the dateline using
        pyGPlates' [DateLineWrapper](https://www.gplates.org/docs/pygplates/generated/pygplates.datelinewrapper)
        by splitting a polyline into multiple polylines at the dateline. This is to avoid
        horizontal lines being formed between polylines at longitudes of -180 and 180 degrees.
        Point features near the poles (-89 & 89 degree latitude) are also clipped to ensure
        compatibility with Cartopy.

        """
        if not self.plate_reconstruction.topology_features:
            logger.warn(
                "Plate model does not have topology features. Unable to plot_ridges."
            )
            return ax

        return self._plot_feature(
            ax,
            get_feature_func=self.get_ridges,
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )
