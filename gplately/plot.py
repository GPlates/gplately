#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
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

"""
This sub-module contains tools for reconstructing and plotting geological features and feature data through time.
"""

import logging
import math
import os
import warnings
from functools import partial
from typing import Union

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false

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
from .grids import Raster, read_netcdf_grid
from .mapping.cartopy_plot import DEFAULT_CARTOPY_PROJECTION, CartopyPlotEngine
from .mapping.plot_engine import PlotEngine
from .tools import EARTH_RADIUS
from .utils.feature_utils import shapelify_features as _shapelify_features
from .utils.plot_utils import _meridian_from_ax, PLOT_DOCSTRING, GET_DATE_DOCSTRING
from .utils.plot_utils import plot_subduction_teeth as _plot_subduction_teeth
from .utils import deprecated
from .utils.io_utils import to_geographic_data_array, load_data_array_from_netcdf

logger = logging.getLogger("gplately")


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
    """Read, reconstruct and plot topology features at specific reconstruction times.

    The :class:`PlotTopologies` class:

    * Read features held in GPlates GPML (GPlates Markup Language) files and ESRI shapefiles.
    * Reconstruct the locations of these features as they migrate through geological time.
    * Turn these reconstructed features into Shapely geometries for plotting.

    To create the :class:`PlotTopologies` object, supply:

    * an instance of the :class:`PlateReconstruction` object

    and optionally,

    * a coastline_filename
    * a continent_filename
    * a COB_filename
    * a reconstruction time
    * an anchor_plate_id

    For example:

    .. code-block:: python
        :linenos:

        # create a PlotTopologies object
        gplot = gplately.plot.PlotTopologies(
            plate_reconstruction,
            coastline_filename=None,
            continent_filename=None,
            COB_filename=None,
            time=100,
            anchor_plate_id=None,
        )

        # Setting a new reconstruction time
        gplot.time = 20  # Ma

    The ``coastline_filename``, ``continent_filename`` and ``COB_filename`` can be paths (:class:`str`) to GPML and/or shapefiles,
    as well as instances of `pygplates.FeatureCollection`_.

    Some features for plotting (like plate boundaries) are taken from the :attr:`PlateReconstruction.topology_features`.
    They have already been reconstructed to the given ``time``.
    Simply provide a new reconstruction time by changing the ``time`` attribute, e.g.

    .. code-block:: python
        :linenos:

        gplot.time = 20 # Ma

    which will automatically reconstruct all topologies to the specified time.
    You **MUST** set ``gplot.time`` before plotting anything.

    A variety of geological features can be plotted on maps, including:

    * subduction boundaries & subduction polarity teeth
    * mid-ocean ridge boundaries
    * transform boundaries
    * miscellaneous boundaries
    * coastline polylines
    * continental polygons and
    * continent-ocean boundary polylines
    * topological plate velocity vector fields
    * netCDF4 `MaskedArray`_ or ndarray raster data:
        - seafloor age grids
        - paleo-age grids
        - global relief (topography and bathymetry)
    * assorted reconstructable feature data, for example:
        - seafloor fabric
        - large igneous provinces
        - volcanic provinces

    .. _pygplates.Feature: https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature
    .. _pygplates.ReconstructedFeatureGeometry: https://www.gplates.org/docs/pygplates/generated/pygplates.reconstructedfeaturegeometry
    .. _geopandas.GeoDataFrame: https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html
    .. _MaskedArray: https://numpy.org/doc/stable/reference/maskedarray.generic.html
    """

    def __init__(
        self,
        plate_reconstruction,
        coastlines=None,
        continents=None,
        COBs=None,
        time: float = 140.0,
        anchor_plate_id=None,
        plot_engine: PlotEngine = CartopyPlotEngine(),
    ):
        """Constructor. Create a :class:`PlotTopologies` object.

        Parameters
        ----------
        plate_reconstruction : PlateReconstruction
            The :class:`PlateReconstruction` object will be used to access a plate
            ``rotation_model`` and a set of ``topology_features`` which contains plate boundary
            features like trenches, ridges and transforms.
        coastlines : str, or `pygplates.FeatureCollection`_
            The full string path to a coastline feature file. Coastline features can also
            be passed as a `pygplates.FeatureCollection`_ object (this is
            the case if these features are sourced from the :class:`DataServer` object).
        continents : str, or `pygplates.FeatureCollection`_
            The full string path to a continent feature file. Continent features can also
            be passed as a `pygplates.FeatureCollection`_ object (this is
            the case if these features are sourced from the :class:`DataServer` object).
        COBs : str, or `pygplates.FeatureCollection`_
            The full string path to a COB feature file. COB features can also be passed
            as a `pygplates.FeatureCollection`_ object (this is the case
            if these features are sourced from the :class:`DataServer` object).
        time: float
            The time (Ma) to reconstruct and plot geological features to.
        anchor_plate_id : int, optional
            Anchor plate ID for reconstruction.
            Defaults to the current anchor plate ID of ``plate_reconstruction`` (its ``anchor_plate_id`` attribute).
        plot_engine : :class:`PlotEngine`, default=CartopyPlotEngine()
            Use Cartopy (:class:`CartopyPlotEngine`) or pyGMT (:class:`PygmtPlotEngine`) to plot the map.


        .. _pygplates.FeatureCollection: https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection
        """

        self._plot_engine = plot_engine
        self.plate_reconstruction = plate_reconstruction
        """The :class:`PlateReconstruction` object will be used to access a plate ``rotation_model`` and a set of ``topology_features``
        which contains plate boundary features like trenches, ridges and transforms.

        :type: PlateReconstruction
        """

        if self.plate_reconstruction.topology_features is None:
            self.plate_reconstruction.topology_features = []
            logger.warning("Plate model does not have topology features.")

        self.base_projection = DEFAULT_CARTOPY_PROJECTION
        """Cartopy map projection. Defaults to ``cartopy.crs.PlateCarree``.
        
        .. seealso::
        
            `Cartopy projection list <https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html>`__.
        """

        # Pickle the __init__ arguments used to create coastline/continent/COB FeatureCollection's (see __reduce__ for details).
        self._coastlines_pickle = coastlines
        self._continents_pickle = continents
        self._COBs_pickle = COBs
        # Load the coastline/continent/COB FeatureCollection's.
        #
        # These will get used to create reconstructed geometries when the time is updated.
        self._coastlines = _load_FeatureCollection(coastlines)
        self._continents = _load_FeatureCollection(continents)
        self._COBs = _load_FeatureCollection(COBs)

        self.coastlines = None
        """A list containing coastline features reconstructed to the specified ``time`` attribute.

        :type: iterable/list of `pygplates.ReconstructedFeatureGeometry`_
        """
        self.continents = None
        """A list containing continent features reconstructed to the specified ``time`` attribute.
        
        :type: iterable/list of `pygplates.ReconstructedFeatureGeometry`_
        """
        self.COBs = None
        """A list containing COB features reconstructed to the specified ``time`` attribute.

        :type: iterable/list of `pygplates.ReconstructedFeatureGeometry`_
        """

        self.trenches = []
        """A list of `pygplates.FeatureType.gpml_subduction_zone` features regardless of their subduction polarity.
        Normally, this list equals `trench_left` + `trench_right` but it can be different if there are some trench features 
        with undefined subduction polarity in the model.
        
        :type: iterable/list of `pygplates.Feature`_
        """

        self.trench_left = []
        """A list of `pygplates.FeatureType.gpml_subduction_zone` features whose subduction polarity is "left".
        
        :type: iterable/list of `pygplates.Feature`_
        """

        self.trench_right = []
        """A list of `pygplates.FeatureType.gpml_subduction_zone` features whose subduction polarity is "right". 

        :type: iterable/list of `pygplates.Feature`_
        """

        self.other = []
        """A list of miscellaneous topological plate boundary sections, 
        which are not subduction zones or mid-ocean ridges (ridge/transform). 
        This can include features like unclassified features, extended continental crusts,
        continental rifts, faults, orogenic belts, fracture zones, inferred paleo boundaries, terrane
        boundaries and passive continental boundaries.
        
        :type: iterable/list of `pygplates.Feature`_
        """
        self._topological_plate_boundaries = None
        self._topologies = None
        self._ridges = []
        self._transforms = []

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
        self.unclassified_features = []

        self._ridges_do_not_use_for_now = None
        self._transforms_do_not_use_for_now = None

        self._plot_engine = plot_engine

        if anchor_plate_id is None:
            # Default to the anchor plate of 'self.plate_reconstruction'.
            self._anchor_plate_id = None
        else:
            self._anchor_plate_id = self._check_anchor_plate_id(anchor_plate_id)

        self._time = None
        # Note: This calls '_update_time()' which reconstructs our features and topologies.
        self.time = time

    def __reduce__(self):
        # Arguments for __init__.
        #
        # Note: Instead of pickling the coastline/continent/COB FeatureCollection's we're pickling the __init__ arguments
        #       used to create them. And if they are just filenames then it ends up being faster to rebuild the FeatureCollection's
        #       (from files when __init__ is called during unpickling) than it is to pickle/unpickle the FeatureCollection's.
        init_args = (
            self.plate_reconstruction,
            self._coastlines_pickle,  # pickle the __init__ argument
            self._continents_pickle,  # pickle the __init__ argument
            self._COBs_pickle,  # pickle the __init__ argument
            self._time,  # is always valid (due to 'self.time = time' in __init__)
            self._anchor_plate_id,  # can be None (meaning use PlateReconstruction's anchor plate ID)
            self._plot_engine,
        )

        # The reason we implement __reduce__ is so that __init__ gets called.
        #
        # We do this instead of implementing __getstate__ (where __init__ doesn't get called) since it's easier to keep
        # things in sync, otherwise __setstate__ would need to stay in sync with __init__ whenever __init_'s code gets modified.
        #
        # It also means we don't have to worry about pickling the *reconstructed* features and topologies.
        # These are unpickable items because pygplates reconstructed feature geometries and resolved topological geometries
        # cannot be pickled (only the *unreconstructed* pygplates features and features collections can be pickled).
        #
        # Note: We set the 3rd tuple element (the object state) to None so that __setstate__ doesn't get called since
        #       we currently don't need to set any extra state beyond what __init__ sets up for us.
        return PlotTopologies, init_args, None

    @property
    def topological_plate_boundaries(self):
        """Resolved topologies for rigid boundaries ONLY."""
        return self._topological_plate_boundaries

    @property
    def topologies(self):
        """Resolved topologies for BOTH rigid boundaries and networks.
        A list containing assorted topologies like:

        - pygplates.FeatureType.gpml_topological_network
        - pygplates.FeatureType.gpml_oceanic_crust
        - pygplates.FeatureType.gpml_topological_slab_boundary
        - pygplates.FeatureType.gpml_topological_closed_plate_boundary

        :type: iterable/list of `pygplates.Feature`_
        """
        return self._topologies

    @property
    def ridges(self):
        """Mid-ocean ridge features (all the features which are labelled as ``gpml:MidOceanRidge`` in the model).
        A list containing ridge and transform boundary sections of type `pygplates.FeatureType.gpml_mid_ocean_ridge`.

        :type: iterable/list of `pygplates.Feature`_
        """
        warnings.warn(
            "The 'ridges' property has been changed since GPlately 1.3.0. "
            "You need to check your workflow to make sure the new 'ridges' property still suits your purpose. "
            "In earlier releases of GPlately, we used an algorithm to identify the 'ridges' and 'transforms' within the gpml:MidOceanRidge features. "
            "Unfortunately, the algorithm did not work very well. So we have removed the algorithm and now the 'ridges' property contains all the features "
            "which are labelled as gpml:MidOceanRidge in the reconstruction model.",
            UserWarning,
            stacklevel=2,
        )
        return self._ridges

    @property
    def transforms(self):
        """Transform boundary features (all the features which are labelled as ``gpml:Transform`` in the model).
        A list containing transform boundary sections of type `pygplates.FeatureType.gpml_transforms`.

        :type: iterable/list of `pygplates.Feature`_
        """
        warnings.warn(
            "The 'transforms' property has been changed since GPlately 1.3.0. "
            "You need to check your workflow to make sure the new 'transforms' property still suits your purpose. "
            "In earlier releases of GPlately, we used an algorithm to identify the 'ridges' and 'transforms' within the gpml:MidOceanRidge features. "
            "Unfortunately, the algorithm did not work very well. So we have removed the algorithm and now the 'transforms' property contains all the features "
            "which are labelled as gpml:Transform in the reconstruction model.",
            UserWarning,
            stacklevel=2,
        )
        return self._transforms

    @property
    def time(self):
        """The time (Ma) to reconstruct and plot geological features to.

        :type: float

        .. note::

            You can either set the ``time`` attribute when creating the :class:`PlotTopologies` object or anytime afterwards.

            .. code-block:: python
                :linenos:

                # set the reconstruction time when creating the PlotTopologies object.
                gplot = gplately.PlotTopologies(..., time=100,...)

                # set the reconstruction time after PlotTopologies object is created.
                gplot.time = 100
        """
        return self._time

    @time.setter
    def time(self, new_time: float):
        """Set the new reconstruction time.

        Reconstruct and resolve topologies to this new time.

        Raises
        ------
        ValueError
            If the reconstruction time is invalid or less then 0.
        """
        try:
            new_time_f = float(new_time)
        except:
            raise ValueError(f"The new 'time' ({new_time}) is not a number.")

        if new_time_f < 0:
            raise ValueError(f"The new 'time' ({new_time}) must be greater than 0.")

        if self._time is None or not math.isclose(new_time_f, self._time):
            self._time = new_time_f
            self._update_time()
        else:
            logger.debug(
                "The new reconstruction 'time' is the same with the old one. Do nothing!"
            )

    @property
    def anchor_plate_id(self):
        """Anchor plate ID for reconstruction. Must be an integer >= 0.
        Defaults to :attr:`PlateReconstruction.anchor_plate_id`.

        :type: int
        """
        if self._anchor_plate_id is None:
            # Default to anchor plate of 'self.plate_reconstruction'.
            return self.plate_reconstruction.anchor_plate_id

        return self._anchor_plate_id

    @anchor_plate_id.setter
    def anchor_plate_id(self, anchor_plate):
        if anchor_plate is None:
            # We'll use the anchor plate of 'self.plate_reconstruction'.
            self._anchor_plate_id = None
        else:
            self._anchor_plate_id = self._check_anchor_plate_id(anchor_plate)

        # Reconstructed/resolved geometries depend on the anchor plate.
        self._update_time()

    @staticmethod
    def _check_anchor_plate_id(id):
        id = int(id)
        if id < 0:
            raise ValueError("Invalid anchor plate ID: {}".format(id))
        return id

    def _update_time(self):
        """Rereconstruct features and topologies to the :attr:`gplately.PlotTopologies.time` attribute.

        The following :class:`gplately.PlotTopologies` attributes are updated whenever a reconstruction ``time`` attribute is set:

        - resolved topology features (topological plates and networks)
        - ridge and transform boundary sections (resolved features)
        - ridge boundary sections (resolved features)
        - transform boundary sections (resolved features)
        - subduction boundary sections (resolved features)
        - left subduction boundary sections (resolved features)
        - right subduction boundary sections (resolved features)
        - other boundary sections (resolved features) that are not subduction zones or mid-ocean ridges(ridge/transform)

        Moreover, coastlines, continents and COBs are reconstructed to the new :attr:`gplately.PlotTopologies.time`.
        """

        # Get the topological snapshot (of resolved topologies) for the current time (and our anchor plate ID).
        topological_snapshot = self.plate_reconstruction.topological_snapshot(
            self.time,
            # If our anchor plate is None then this will use the anchor plate of 'self.plate_reconstruction'...
            anchor_plate_id=self._anchor_plate_id,
        )

        #
        # NOTE: If you add a new data member here that's a pygplates reconstructable feature geometry or resolved topological geometry,
        #       then be sure to also include it in __getstate__/()__setstate__()
        #       (basically anything reconstructed or resolved by pygplates since those cannot be pickled).
        #

        # Extract (from the topological snapshot) resolved topologies for BOTH rigid boundaries and networks.
        self._topologies = [
            t.get_resolved_feature()
            for t in topological_snapshot.get_resolved_topologies(
                resolve_topology_types=pygplates.ResolveTopologyType.boundary
                | pygplates.ResolveTopologyType.network
            )
        ]

        (
            self._topological_plate_boundaries,
            self._ridges,
            self._ridges_do_not_use_for_now,  # the algorithm to separate ridges and transforms has not been ready yet
            self._transforms_do_not_use_for_now,
            self.trenches,
            self.trench_left,
            self.trench_right,
            self.other,
        ) = ptt.resolve_topologies.resolve_topological_snapshot_into_features(
            topological_snapshot,
            # use ResolveTopologyType.boundary parameter to resolve rigid plate boundary only
            # because the Mid-ocean ridges(and transforms) should not contain lines from topological networks
            resolve_topology_types=pygplates.ResolveTopologyType.boundary,  # type: ignore
        )

        for topol in self.other:
            if topol.get_feature_type() == pygplates.FeatureType.gpml_continental_rift:  # type: ignore
                self.continental_rifts.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_fault:  # type: ignore
                self.faults.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_fracture_zone:  # type: ignore
                self.fracture_zones.append(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_inferred_paleo_boundary  # type: ignore
            ):
                self.inferred_paleo_boundaries.append(topol)

            elif (
                topol.get_feature_type() == pygplates.FeatureType.gpml_terrane_boundary  # type: ignore
            ):
                self.terrane_boundaries.append(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_transitional_crust  # type: ignore
            ):
                self.transitional_crusts.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_orogenic_belt:  # type: ignore
                self.orogenic_belts.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_suture:  # type: ignore
                self.sutures.append(topol)

            elif (
                topol.get_feature_type() == pygplates.FeatureType.gpml_continental_crust  # type: ignore
            ):
                self.continental_crusts.append(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_extended_continental_crust  # type: ignore
            ):
                self.extended_continental_crusts.append(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_passive_continental_boundary  # type: ignore
            ):
                self.passive_continental_boundaries.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_slab_edge:  # type: ignore
                self.slab_edges.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_transform:  # type: ignore
                self._transforms.append(topol)
                # remove transform features from 'other' since 'other' is supposed to contain features that
                # are not subduction zones or mid-ocean ridges(ridge/transform)
                self.other.remove(topol)

            elif (
                topol.get_feature_type()
                == pygplates.FeatureType.gpml_unclassified_feature  # type: ignore
            ):
                self.unclassified_features.append(topol)

        # reconstruct other important polygons and lines
        if self._coastlines:
            self.coastlines = self.plate_reconstruction.reconstruct(
                self._coastlines,
                self.time,
                from_time=0,
                # If our anchor plate is None then this will use the anchor plate of 'self.plate_reconstruction'...
                anchor_plate_id=self._anchor_plate_id,
            )

        if self._continents:
            self.continents = self.plate_reconstruction.reconstruct(
                self._continents,
                self.time,
                from_time=0,
                # If our anchor plate is None then this will use the anchor plate of 'self.plate_reconstruction'...
                anchor_plate_id=self._anchor_plate_id,
            )

        if self._COBs:
            self.COBs = self.plate_reconstruction.reconstruct(
                self._COBs,
                self.time,
                from_time=0,
                # If our anchor plate is None then this will use the anchor plate of 'self.plate_reconstruction'...
                anchor_plate_id=self._anchor_plate_id,
            )

    @validate_reconstruction_time
    def get_feature(
        self,
        feature,
        central_meridian=0.0,
        tessellate_degrees=None,
        validate_reconstruction_time=True,
    ):
        """Convert feature(s) to a `geopandas.GeoDataFrame`_ object.

        Parameters
        ----------
        feature : `pygplates.Feature`_, `pygplates.ReconstructedFeatureGeometry`_, pygplates.GeometryOnSphere or iterable of the aforementioned three
            Feature object(s).

        Returns
        -------
        gdf : geopandas.GeoDataFrame
            A `geopandas.GeoDataFrame`_ object contaning the feature geometries.


        .. note::

            The feature(s) needed to produce the `geopandas.GeoDataFrame`_ object should already be reconstructed to a ``time``.
            This function converts the feature(s) into a set of Shapely geometries and put them into a `geopandas.GeoDataFrame`_ object.

        """

        if not feature:
            raise ValueError(
                "No feature to convert. Make sure the `feature` parameter contains feature(s)."
            )
        if not validate_reconstruction_time:
            logger.debug(
                "The 'validate_reconstruction_time' parameter is set to False. The reconstruction time is not validated."
            )
        shp = shapelify_features(
            feature,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

        return gpd.GeoDataFrame({"geometry": shp}, geometry="geometry", crs="EPSG:4326")  # type: ignore

    @append_docstring(PLOT_DOCSTRING.format("feature"))
    def plot_feature(
        self, ax, feature, feature_name="", color: Union[str, list] = "black", **kwargs
    ):
        """Plot `pygplates.FeatureCollection`_ or `pygplates.Feature`_ onto a map."""
        if not feature:
            logger.warning(
                f"The feature ({feature_name}) is empty and will not be plotted."
            )
            return ax
        else:
            if "edgecolor" not in kwargs.keys():
                kwargs["edgecolor"] = color
            if "facecolor" not in kwargs.keys():
                kwargs["facecolor"] = "none"
            return self._plot_feature(ax, partial(self.get_feature, feature), **kwargs)

    def _plot_feature(self, ax, get_feature_func, **kwargs) -> None:
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

        if not isinstance(gdf, gpd.GeoDataFrame):
            raise Exception(
                f"Expecting a GeoDataFrame object, but the gdf is {type(gdf)}"
            )

        if len(gdf) == 0:
            logger.debug("No feature found for plotting. Do nothing and return.")
            return ax

        self._plot_engine.plot_geo_data_frame(ax, gdf, **kwargs)

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("coastlines"))
    def get_coastlines(self, central_meridian=0.0, tessellate_degrees=None):
        """Return the reconstructed coastlines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.coastlines,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("coastlines"))
    def plot_coastlines(self, ax, color: Union[str, list] = "black", **kwargs):
        """Plot reconstructed coastlines on a map."""
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
        """Return the reconstructed continental polygons as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.continents,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("continents"))
    def plot_continents(self, ax, color="black", **kwargs):
        """Plot the reconstructed continental polygons on a map."""
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
        """Return the reconstructed continent-ocean boundaries as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.COBs,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("continent ocean boundaries"))
    def plot_continent_ocean_boundaries(
        self, ax, color: Union[str, list] = "black", **kwargs
    ):
        """Plot the reconstructed continent-ocean boundaries (COBs) on a map."""
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
        """Return the reconstructed mid-ocean ridge lines (gpml:MidOceanRidge) as a `geopandas.GeoDataFrame`_ object.

        .. note::

            This function only returns the reconstructed mid-ocean ridges as polylines.
            If you need to plot gap-free topological plate boundaries, see :meth:`plot_all_topological_sections`,
            which allows you to plot all topological sections with different styles.
        """
        warnings.warn(
            "The 'get_ridges' function has been changed since GPlately 1.3.0. "
            "You need to check your workflow to make sure the new 'get_ridges' function still suits your purpose. "
            "In earlier releases of GPlately, we used an algorithm to identify the 'ridges' and 'transforms' within the gpml:MidOceanRidge features. "
            "Unfortunately, the algorithm did not work very well. So we have removed the algorithm and now the 'get_ridges' function returns all the features "
            "which are labelled as gpml:MidOceanRidge in the reconstruction model.",
            UserWarning,
            stacklevel=5,
        )

        return self.get_feature(
            self.ridges,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @validate_topology_availability("ridges")
    @append_docstring(PLOT_DOCSTRING.format("ridges"))
    def plot_ridges(self, ax, color="black", **kwargs):
        """Plot the reconstructed mid-ocean ridge lines(gpml:MidOceanRidge) on a map.

        .. note::

            This function only plots the reconstructed ridges as polylines.
            If you need to plot gap-free topological plate boundaries, see :meth:`plot_all_topological_sections`,
            which allows you to plot all topological sections with different styles.
        """

        warnings.warn(
            "The 'plot_ridges' function has been changed since GPlately 1.3.0. "
            "You need to check your workflow to make sure the new 'plot_ridges' function still suits your purpose. "
            "In earlier releases of GPlately, we used an algorithm to identify the 'ridges' and 'transforms' within the gpml:MidOceanRidge features. "
            "Unfortunately, the algorithm did not work very well. So we have removed the algorithm and now the 'plot_ridges' function plots all the features "
            "which are labelled as gpml:MidOceanRidge in the reconstruction model.",
            UserWarning,
            stacklevel=5,
        )

        return self.plot_feature(
            ax,
            self._ridges,
            feature_name="ridges",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("trenches"))
    def get_trenches(self, central_meridian=0.0, tessellate_degrees=1):
        """Return the reconstructed trench lines as a `geopandas.GeoDataFrame`_ object.

        .. note::

            This function only returns the reconstructed trenches as polylines.
            If you need to plot gap-free topological plate boundaries, see :meth:`plot_all_topological_sections`,
            which allows you to plot all topological sections with different styles.
        """
        return self.get_feature(
            self.trenches,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_topology_availability("trenches")
    @append_docstring(PLOT_DOCSTRING.format("trenches"))
    def plot_trenches(self, ax, color="black", **kwargs):
        """Plot the reconstructed trenches on a map.

        .. note::

            This function only plots the reconstructed trenches as polylines.
            If you need to plot gap-free topological plate boundaries, see :meth:`plot_all_topological_sections`,
            which allows you to plot all topological sections with different styles.
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
        """Return the reconstructed "other" lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.other,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("other"))
    def plot_misc_boundaries(self, ax, color="black", **kwargs):
        """Plot the reconstructed miscellaneous plate boundary sections(polylines) on a map.
        Sections are neither subduction zones nor mid-ocean ridges, nor transforms.
        """
        return self.plot_feature(
            ax,
            self.other,
            feature_name="miscellaneous topological plate boundary sections",
            facecolor="none",
            edgecolor=color,
            **kwargs,
        )

    @validate_reconstruction_time
    def get_subduction_direction(self, central_meridian=0.0, tessellate_degrees=None):
        """Return the :attr:`PlotTopologies.trench_left` and
        :attr:`PlotTopologies.trench_right` as a `geopandas.GeoDataFrame`_ object.

        Parameters
        ----------
        central_meridian : float, default=0.0
            The central meridian of the map. This will affect the dateline wrapping.
        tessellate_degrees : float or None, default=1.0
            If provided, geometries will be tessellated to this resolution prior to wrapping.

        Returns
        -------
        gdf_left :geopandas.GeoDataFrame
            A `geopandas.GeoDataFrame`_ object containing `trench_left` geometries.
        gdf_right : geopandas.GeoDataFrame
            A `geopandas.GeoDataFrame`_ object containing `trench_right` geometries.

        """
        if self.trench_left is None or self.trench_right is None:
            raise Exception(
                "No subduction zone/trench data is found. Make sure the plate model has topology feature."
            )

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

        gdf_left = gpd.GeoDataFrame(
            {"geometry": trench_left_features}, geometry="geometry", crs="EPSG:4326"
        )  # type: ignore
        gdf_right = gpd.GeoDataFrame(
            {"geometry": trench_right_features}, geometry="geometry", crs="EPSG:4326"
        )  # type: ignore

        return gdf_left, gdf_right

    @validate_reconstruction_time
    @validate_topology_availability("Subduction Zones")
    def plot_subduction_teeth(
        self, ax, spacing=0.07, size=None, aspect=None, color="black", **kwargs
    ) -> None:
        """Plot subduction teeth.

        .. note::

            Subduction teeth are tessellated from :attr:`PlotTopologies.trench_left` and
            :attr:`PlotTopologies.trench_right` attributes, and transformed into Shapely polygons for plotting.

        Parameters
        ----------
        ax :
            Cartopy ax or pygmt figure object.

        spacing : float, default=0.07
            The tessellation threshold (in radians). Parametrises subduction tooth density.
            Triangles are generated only along line segments with distances that exceed the given threshold ``spacing``.

        size : float, default=None
            Length of teeth triangle base (in radians). If kept at ``None``, then ``size = 0.5*spacing``.

        aspect : float, default=None
            Aspect ratio of teeth triangles. If kept at ``None``, then ``aspect = 2/3*size``.

        color : str, default='black'
            The colour of the teeth. By default, it is set to black.

        **kwargs :
            Keyword arguments for plotting subduction teeth.
            See Matplotlib keyword arguments `here <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html>`__.
        """
        kwargs["spacing"] = spacing
        kwargs["size"] = size
        kwargs["aspect"] = aspect

        central_meridian = _meridian_from_ax(ax)
        tessellate_degrees = np.rad2deg(spacing)
        gdf_subduction_left, gdf_subduction_right = self.get_subduction_direction(
            tessellate_degrees=tessellate_degrees, central_meridian=central_meridian
        )

        self._plot_engine.plot_subduction_zones(
            ax, gdf_subduction_left, gdf_subduction_right, color=color, **kwargs
        )

    def plot_plate_polygon_by_id(self, ax_or_fig, plate_id, color="black", **kwargs):
        """Plot a plate polygon with the given ``plate_id`` on a map.

        Parameters
        ----------
        ax_or_fig :
            Cartopy ax or PyGMT figure object.

        plate_id : int
            A plate ID that identifies the continental polygon to plot. See the
            `Global EarthByte plate IDs list <https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/SampleData/FeatureCollections/Rotations/Global_EarthByte_PlateIDs_20071218.pdf>`__
            for a full list of plate IDs.

        **kwargs :
            Keyword arguments for plotting.
        """
        features = []
        if self.topologies:
            features = (
                [
                    feature
                    for feature in self.topologies
                    if feature.get_reconstruction_plate_id() == plate_id
                ],
            )
        self.plot_feature(
            ax_or_fig,
            features,
            color=color,
            **kwargs,
        )

    def plot_grid(
        self, ax_or_fig, grid, extent=(-180, 180, -90, 90), shading=None, **kwargs
    ):
        """Plot a grid onto a map.

        Parameters
        ----------
        ax_or_fig :
            Cartopy ax or PyGMT figure object.

        grid : NumPy `MaskedArray`_, GPlately `Raster` or a time-dependent raster name.
            A `MaskedArray`_ with elements that define a grid. The number of rows in the raster
            corresponds to the number of latitudinal coordinates, while the number of raster
            columns corresponds to the number of longitudinal coordinates.
            Alternatively, a GPlately `Raster` object can be provided.
            If a raster name is provided, the raster will be looked up from the time-dependent rasters registered in the Plate Model Manager.
            The :class:`gplately.PlateReconstruction` object must be created with a valid :class:`PlateModel` object.

        extent : tuple, default=(-180, 180, -90, 90)
            A tuple of 4 (min_lon, max_lon, min_lat, max_lat) representing the extent of gird.


        .. note::

            The parameters of this function are different for different plot engines. See :meth:`CartopyPlotEngine.plot_grid`
            and :meth:`PyGMTPlotEngine.plot_grid` for details.

        """

        # TODO: the parameters of this function need to be unified for different plot engines.

        if isinstance(grid, str):  # grid is a raster name
            if not self.plate_reconstruction.plate_model:
                raise Exception(
                    "The 'plate_reconstruction' does not have a valid 'plate_model' object. "
                    "Cannot look up the raster by name. Make sure to create the 'plate_reconstruction' with a valid 'plate_model' object."
                )

            grid = Raster(
                data=self.plate_reconstruction.plate_model.get_raster(grid, self.time),
                plate_reconstruction=self.plate_reconstruction,
                extent=(-180, 180, -90, 90),
            )

        return self._plot_engine.plot_grid(
            ax_or_fig,
            grid,
            extent=extent,
            shading=shading,
            projection=self.base_projection,
            **kwargs,
        )

    def plot_grid_from_netCDF(
        self,
        ax_or_fig,
        filename,
        var_name="",
        shading=None,
        shading_var_name="",
        **kwargs,
    ):
        """Read raster data from a netCDF file and plot it on a map.

        Parameters
        ----------
        ax_or_fig :
            Cartopy ax or PyGMT figure.
        filename : str
            Full path to a netCDF file.
        var_name : str, default=""
            The variable name of the raster data in the netCDF file. If not provided, the first variable in the netCDF file will be used.
        shading : str, default=None
            Full path to a netCDF file containing illumination data, such as topography grid. If not provided, no shading will be applied.
            If you are using PyGMT plot engine, the valid shading parameters for PyGMT can also be used here, such as ``-I+d`` or ``"+a315+ne0.6"``.
        shading_var_name : str, default=""
            The variable name of the shading data in the shading netCDF file. If not provided, the first variable in the netCDF file will be used.
        **kwargs :
            Keyword arguments for Matplotlib's ``imshow()`` or PyGMT's ``grdimage()``.
        """

        raster_xr = load_data_array_from_netcdf(filename, var_name=var_name)
        raster_xr = to_geographic_data_array(raster_xr)

        if (
            shading
            and isinstance(shading, (str, bytes, os.PathLike))
            and os.path.isfile(shading)
        ):
            shading = load_data_array_from_netcdf(shading, var_name=shading_var_name)

        lon_coords = raster_xr.coords["lon"].values
        lat_coords = raster_xr.coords["lat"].values
        extent = [
            float(np.nanmin(lon_coords)),
            float(np.nanmax(lon_coords)),
            float(np.nanmin(lat_coords)),
            float(np.nanmax(lat_coords)),
        ]

        # Override matplotlib default origin ('upper')
        origin = kwargs.pop("origin", "lower")
        if lat_coords[0] < lat_coords[-1]:
            origin = "lower"
        else:
            origin = "upper"

        return self.plot_grid(
            ax_or_fig,
            raster_xr,
            shading=shading,
            extent=extent,
            origin=origin,
            **kwargs,
        )

    def plot_plate_motion_vectors(
        self,
        ax_or_fig,
        color="black",
        lons=[],
        lats=[],
        spacingX=10,
        spacingY=10,
        normalise=False,
        delta_time=5.0,
        **kwargs,
    ):
        """Calculate plate motion velocity vector fields at a particular geological time and plot them on a map.

        .. note::

            The :meth:`plot_plate_motion_vectors` generates a MeshNode domain of point features using
            given spacings in the X and Y directions (``spacingX`` and ``spacingY``). Each point in
            the domain is assigned a plate ID, and these IDs are used to obtain equivalent stage
            rotations of identified tectonic plates over a 5 Ma time interval. Each point and
            its stage rotation are used to calculate plate velocities at a particular geological
            time. Velocities for each domain point are represented in the north-east-down
            coordinate system and plotted on a GeoAxes.

            Vector fields can be optionally normalised by setting ``normalise`` to ``True``. This
            makes vector arrow lengths uniform.

        Parameters
        ----------
        ax_or_fig :
            Cartopy ax or PyGMT figure.
        color : str, default='black'
            The colour of the velocity vectors. By default, it is set to black.
        lons : list of float, default=[]
            A list of longitudes for the velocity domain point features. If not provided, a global mesh grid of points will be generated with the given spacings in the X and Y directions (``spacingX`` and ``spacingY``).
        lats : list of float, default=[]
            A list of latitudes for the velocity domain point features. If not provided, a global mesh grid of points will be generated with the given spacings in the X and Y directions (``spacingX`` and ``spacingY``).
        spacingX : int, default=10
            The spacing in the X direction used to make the velocity domain point feature mesh.
        spacingY : int, default=10
            The spacing in the Y direction used to make the velocity domain point feature mesh.
        normalise : bool, default=False
            Choose whether to normalise the velocity magnitudes so that vector lengths are all equal.
        delta_time : float, default=5.0
            The time interval (in Ma) used to calculate the velocity vectors. The stage rotation between the current time and the time interval is used to calculate the velocity vectors. By default, it is set to 5 Ma.
        **kwargs :
            Keyword arguments for Matplotlib quiver() or PyGMT velo().
        """

        if len(lons) != 0 or len(lats) != 0:
            assert len(lons) == len(
                lats
            ), "The 'lons' and 'lats' parameters must be of the same length if they are provided."

        # if no longitude and latitude arrays are provided, generate a global mesh grid of points with the given spacings in the X and Y directions (spacingX and spacingY)
        if not lons or not lats:
            lon_arr = np.arange(-180, 180 + spacingX, spacingX)
            lat_arr = np.arange(-90, 90 + spacingY, spacingY)
            lons = np.tile(lon_arr, len(lat_arr))
            lats = np.repeat(lat_arr, len(lon_arr))

        velocity_lons, velocity_lats = self.plate_reconstruction.get_point_velocities(
            lons,
            lats,
            self.time,
            delta_time=delta_time,
            # Match previous implementation that used ptt.velocity_tools.get_plate_velocities()...
            velocity_units=pygplates.VelocityUnits.kms_per_my,
            return_east_north_arrays=True,
        )

        if normalise:
            mag = np.hypot(velocity_lons, velocity_lats)
            mag[mag == 0] = 1
            velocity_lons /= mag
            velocity_lats /= mag

        data = gpd.GeoDataFrame(
            {
                "lon": lons,
                "lat": lats,
                "east_vel": velocity_lons,
                "north_vel": velocity_lats,
                "east_sigma": [0] * len(lons),
                "north_sigma": [0] * len(lats),
                "correlation": [0] * len(lons),
            }
        )

        return self._plot_engine.plot_plate_motion_vectors(
            ax_or_fig, data, projection=self.base_projection, color=color, **kwargs
        )

    def plot_pole(self, ax_or_fig, lon, lat, a95, color="green", **kwargs):
        """
        Plot a paleomagnetic pole onto a map.

        Parameters
        ----------
        ax_or_fig :
            Matplotlib axes or PyGMT figure.
        lon : float
            Longitudinal coordinate of the pole.
        lat : float
            Latitudinal coordinate of the pole.
        a95 : float
            The radius(in degrees) of the cone of 95% confidence around a paleomagnetic pole.

        Returns
        -------
        Matplotlib axes or PyGMT figure.
        """
        return self._plot_engine.plot_pole(ax_or_fig, lon, lat, a95, color=color)

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("continental rifts"))
    def get_continental_rifts(
        self,
        central_meridian=0.0,
        tessellate_degrees=None,
    ):
        """Return the reconstructed continental rift lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.continental_rifts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("continental rifts"))
    def plot_continental_rifts(self, ax, color="black", **kwargs):
        """Plot continental rifts on a map."""
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
        """Return the reconstructed fault lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.faults,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("faults"))
    def plot_faults(self, ax, color="black", **kwargs):
        """Plot faults on a map."""
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
        """Return the reconstructed fracture zone lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.fracture_zones,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("fracturezones"))
    def plot_fracture_zones(self, ax, color="black", **kwargs):
        """Plot fracture zones on a map."""
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
        """Return the reconstructed inferred paleo boundary lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.inferred_paleo_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("inferred paleo-boundaries"))
    def plot_inferred_paleo_boundaries(self, ax, color="black", **kwargs):
        """Plot inferred paleo boundaries on a map."""
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
        """Return the reconstructed terrane boundary lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.terrane_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("terrane boundaries"))
    def plot_terrane_boundaries(self, ax, color="black", **kwargs):
        """Plot terrane boundaries on a map."""
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
        """Return the reconstructed transitional crust lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.transitional_crusts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("transitional crusts"))
    def plot_transitional_crusts(self, ax, color="black", **kwargs):
        """Plot transitional crust on a map."""
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
        """Return the reconstructed orogenic belt lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.orogenic_belts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("orogenic belts"))
    def plot_orogenic_belts(self, ax, color="black", **kwargs):
        """Plot orogenic belts on a map."""
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
        """Return the reconstructed suture lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.sutures,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("sutures"))
    def plot_sutures(self, ax, color="black", **kwargs):
        """Plot sutures on a map."""
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
        """Return the reconstructed continental crust lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.continental_crusts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("continental crusts"))
    def plot_continental_crusts(self, ax, color="black", **kwargs):
        """Plot continental crust lines on a map."""
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
        """Return the reconstructed extended continental crust lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.extended_continental_crusts,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("extended continental crusts"))
    def plot_extended_continental_crusts(self, ax, color="black", **kwargs):
        """Plot extended continental crust lines on a map."""
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
        """Return the reconstructed passive continental boundary lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.passive_continental_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("passive continental boundaries"))
    def plot_passive_continental_boundaries(self, ax, color="black", **kwargs):
        """Plot passive continental boundaries on a map."""
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
        """Return the reconstructed slab edge lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.slab_edges,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("slab edges"))
    def plot_slab_edges(self, ax, color="black", **kwargs):
        """Plot slab edges on a map."""
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
        """Return the reconstructed transform lines(gpml:Transform) as a `geopandas.GeoDataFrame`_ object.

        .. note::

            This function only returns the reconstructed transforms as polylines.
            If you need to plot gap-free topological plate boundaries, see :meth:`plot_all_topological_sections`,
            which allows you to plot all topological sections with different styles.
        """
        warnings.warn(
            "The 'get_transforms' function has been changed since GPlately 1.3.0. "
            "You need to check your workflow to make sure the new 'get_transforms' function still suits your purpose. "
            "In earlier releases of GPlately, we used an algorithm to identify the 'ridges' and 'transforms' within the gpml:MidOceanRidge features. "
            "Unfortunately, the algorithm did not work very well. So we have removed the algorithm and now the 'get_transforms' function returns all the features "
            "which are labelled as gpml:Transform in the reconstruction model.",
            UserWarning,
            stacklevel=5,
        )

        return self.get_feature(
            self._transforms,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("transforms"))
    def plot_transforms(self, ax, color="black", **kwargs):
        """Plot transform boundaries(gpml:Transform) on a map.

        .. note::

            This function only plots the reconstructed transforms as polylines.
            If you need to plot gap-free topological plate boundaries, see :meth:`plot_all_topological_sections`,
            which allows you to plot all topological sections with different styles.
        """

        warnings.warn(
            "The 'plot_transforms' function has been changed since GPlately 1.3.0. "
            "You need to check your workflow to make sure the new 'plot_transforms' function still suits your purpose. "
            "In earlier releases of GPlately, we used an algorithm to identify the 'ridges' and 'transforms' within the gpml:MidOceanRidge features. "
            "Unfortunately, the algorithm did not work very well. So we have removed the algorithm and now the 'plot_transforms' function plots all the features "
            "which are labelled as gpml:Transform in the reconstruction model.",
            UserWarning,
            stacklevel=5,
        )

        return self.plot_feature(
            ax,
            self._transforms,
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
        """Return the reconstructed unclassified feature lines as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self.unclassified_features,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_reconstruction_time
    @append_docstring(PLOT_DOCSTRING.format("unclassified features"))
    def plot_unclassified_features(self, ax, color="black", **kwargs):
        """Plot GPML unclassified features on a map."""
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
        """Return the reconstructed topological features listed below as a `geopandas.GeoDataFrame`_ object.

        - pygplates.FeatureType.gpml_topological_network
        - pygplates.FeatureType.gpml_oceanic_crust
        - pygplates.FeatureType.gpml_topological_slab_boundary
        - pygplates.FeatureType.gpml_topological_closed_plate_boundary

        """

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
            crs="EPSG:4326",
        )  # type: ignore
        return gdf

    @validate_topology_availability("all topologies")
    @append_docstring(PLOT_DOCSTRING.format("topologies"))
    def plot_all_topologies(self, ax, color: Union[str, list] = "black", **kwargs):
        """Plot the reconstructed topological features listed below on a map.

        - pygplates.FeatureType.gpml_topological_network
        - pygplates.FeatureType.gpml_oceanic_crust
        - pygplates.FeatureType.gpml_topological_slab_boundary
        - pygplates.FeatureType.gpml_topological_closed_plate_boundary

        """
        if "edgecolor" not in kwargs.keys():
            kwargs["edgecolor"] = color
        if "facecolor" not in kwargs.keys():
            kwargs["facecolor"] = "none"

        return self._plot_feature(
            ax,
            self.get_all_topologies,
            **kwargs,
        )

    @validate_reconstruction_time
    @append_docstring(GET_DATE_DOCSTRING.format("all topological sections"))
    def get_all_topological_sections(
        self,
        central_meridian=0.0,
        tessellate_degrees=1,
    ):
        """Return the reconstructed topological plate boundary sections listed below as a `geopandas.GeoDataFrame`_ object.

        - ridges and transforms
        - subduction zones (trenches)
        - other miscellaneous boundary sections that are neither subduction zones nor mid-ocean ridges (ridge/transform)
        """

        # get plate IDs and feature types to add to geodataframe
        geometries = []
        plate_IDs = []
        feature_types = []
        feature_names = []
        for topo in [
            *self.ridges,
            *self.trenches,
            *self.transforms,
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
                    converted = linemerge(converted, directed=True)
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
            crs="EPSG:4326",
        )  # type: ignore
        return gdf

    @validate_topology_availability("all topological sections")
    def plot_all_topological_sections(
        self,
        ax,
        color="black",
        plot_subduction_teeth=False,
        ridge_kwargs=None,
        transform_kwargs=None,
        trench_kwargs=None,
        other_kwargs=None,
        **kwargs,
    ):
        """Plot the reconstructed topological boundary sections(polylines) on a map.

        This function draws gap-free topological plate boundaries by layering the
        topological sections in the following order:

        1. Miscellaneous topological boundary sections that are neither subduction zones nor
           mid-ocean ridges
        2. Ridges (gpml:MidOceanRidge)
        3. Transforms (gpml:Transform)
        4. Trenches (gpml:SubductionZone)
        5. Optional subduction teeth when ``plot_subduction_teeth=True``

        Parameters
        ----------
        ax :
            Cartopy GeoAxes or PyGMT Figure object to plot on.
        color : str or list, default="black"
            Default colour used for any layer that does not specify its own ``color`` or
            ``edgecolor``.
        plot_subduction_teeth : bool, default=False
            If ``True``, also plot subduction teeth using the trench layer styling.
        ridge_kwargs, transform_kwargs, trench_kwargs, other_kwargs : dict or None, optional
            Layer-specific keyword arguments forwarded to the corresponding plotting call.
            If a layer-specific dict provides ``color`` or ``edgecolor``, it overrides the
            shared ``color`` argument for that layer. Keys in these dictionaries override
            keys passed through ``**kwargs``.
        **kwargs :
            Shared keyword arguments forwarded to each layer. Layer-specific dictionaries
            take precedence over overlapping keys in this shared set.

        Returns
        -------
        ax
            The input axes/figure object, for convenience and chaining.

        Examples
        --------

        .. code-block:: python

            gplot.plot_all_topological_sections(
                ax,
                plot_subduction_teeth=True,
                other_kwargs={"color": "lightgrey", "linewidth": 0.8},
                ridge_kwargs={"color": "red", "linewidth": 1.0},
                transform_kwargs={"color": "green", "linewidth": 1.0},
                trench_kwargs={"color": "blue", "linewidth": 1.0},
            )
        """
        ridge_kwargs = {} if ridge_kwargs is None else dict(ridge_kwargs)
        transform_kwargs = {} if transform_kwargs is None else dict(transform_kwargs)
        trench_kwargs = {} if trench_kwargs is None else dict(trench_kwargs)
        other_kwargs = {} if other_kwargs is None else dict(other_kwargs)
        shared_kwargs = dict(kwargs)

        other_color = other_kwargs.pop("color", None)
        if other_color is None:
            other_color = other_kwargs.pop("edgecolor", None)
        if other_color is None:
            other_color = color
        other_plot_kwargs = dict(shared_kwargs)
        other_plot_kwargs.update(other_kwargs)
        self.plot_feature(
            ax,
            self.other,
            feature_name="other miscellaneous topological boundary sections",
            color=other_color,
            **other_plot_kwargs,
        )

        ridge_color = ridge_kwargs.pop("color", None)
        if ridge_color is None:
            ridge_color = ridge_kwargs.pop("edgecolor", None)
        if ridge_color is None:
            ridge_color = color
        ridge_plot_kwargs = dict(shared_kwargs)
        ridge_plot_kwargs.update(ridge_kwargs)
        self.plot_ridges(ax, color=ridge_color, **ridge_plot_kwargs)

        transform_color = transform_kwargs.pop("color", None)
        if transform_color is None:
            transform_color = transform_kwargs.pop("edgecolor", None)
        if transform_color is None:
            transform_color = color
        transform_plot_kwargs = dict(shared_kwargs)
        transform_plot_kwargs.update(transform_kwargs)
        self.plot_transforms(ax, color=transform_color, **transform_plot_kwargs)

        trench_color = trench_kwargs.pop("color", None)
        if trench_color is None:
            trench_color = trench_kwargs.pop("edgecolor", None)
        if trench_color is None:
            trench_color = color
        trench_plot_kwargs = dict(shared_kwargs)
        trench_plot_kwargs.update(trench_kwargs)
        if plot_subduction_teeth:
            trench_plot_kwargs_without_gmtlabel = trench_plot_kwargs.copy()
            trench_plot_kwargs_without_gmtlabel.pop("gmtlabel", None)
            self.plot_trenches(
                ax, color=trench_color, **trench_plot_kwargs_without_gmtlabel
            )
            trench_teeth_kwargs = dict(shared_kwargs)
            trench_teeth_kwargs.update(trench_kwargs)
            self.plot_subduction_teeth(ax, color=trench_color, **trench_teeth_kwargs)
        else:
            self.plot_trenches(ax, color=trench_color, **trench_plot_kwargs)
        return ax

    def get_topological_plate_boundaries(
        self, central_meridian=0.0, tessellate_degrees=1
    ):
        """Return the reconstructed topological plate boundaries(polygons) as a `geopandas.GeoDataFrame`_ object."""
        return self.get_feature(
            self._topological_plate_boundaries,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
        )

    @validate_topology_availability("topological plate boundaries")
    @append_docstring(PLOT_DOCSTRING.format("topological plate boundaries"))
    def plot_topological_plate_boundaries(self, ax, color="black", **kwargs):
        """Plot the topological plate boundaries(polygons)."""
        return self.plot_feature(
            ax,
            self._topological_plate_boundaries,
            feature_name="topological plate boundaries",
            color=color,
            **kwargs,
        )

    #
    # The functions below are deprecated and will be removed in the future.
    # Please do not use them in new code and check your workflow if you are using them in existing code.
    #
    @property
    @append_docstring(deprecated._ridge_transforms_impl.__doc__ or "")
    def ridge_transforms(self):
        return deprecated._ridge_transforms_impl(self)

    @property
    @append_docstring(deprecated._plot_misc_transforms_impl.__doc__ or "")
    def misc_transforms(self):
        return deprecated._misc_transforms_impl(self)

    @append_docstring(deprecated._plot_misc_transforms_impl.__doc__ or "")
    def plot_misc_transforms(self, ax, color="black", **kwargs):
        return deprecated._plot_misc_transforms_impl(self, ax, color=color, **kwargs)

    @append_docstring(deprecated._get_misc_transforms_impl.__doc__ or "")
    def get_misc_transforms(self, central_meridian=0.0, tessellate_degrees=None): # fmt: skip
        return deprecated._get_misc_transforms_impl(self, central_meridian=central_meridian, tessellate_degrees=tessellate_degrees) # fmt: skip

    @append_docstring(deprecated._plot_ridges_and_transforms_impl.__doc__ or "")
    def plot_ridges_and_transforms(self, ax, color="black", **kwargs):
        return deprecated._plot_ridges_and_transforms_impl(self, ax, color=color, **kwargs) # fmt: skip

    @append_docstring(deprecated._plot_plate_id_impl.__doc__ or "")
    def plot_plate_id(self, *args, **kwargs):
        return deprecated._plot_plate_id_impl(self, *args, **kwargs)

    @append_docstring(deprecated._get_ridges_and_transforms_impl.__doc__ or "")
    def get_ridges_and_transforms(self, central_meridian=0.0, tessellate_degrees=1):
        return deprecated._get_ridges_and_transforms_impl(self, central_meridian=central_meridian, tessellate_degrees=tessellate_degrees )# fmt: skip

    @append_docstring(deprecated._plot_subduction_teeth_deprecated.__doc__ or "")
    def _plot_subduction_teeth_deprecated(self, ax, spacing=0.1, size=2.0, aspect=1, color="black", **kwargs): # fmt: skip
        return deprecated._plot_subduction_teeth_deprecated(self, ax, spacing=spacing, size=size, aspect=aspect, color=color, **kwargs) # fmt: skip

    @append_docstring(deprecated._tessellate_triangles_impl.__doc__ or "")
    def _tessellate_triangles(self, features, tesselation_radians, triangle_base_length, triangle_aspect=1.0): # fmt: skip
        return deprecated._tessellate_triangles_impl(self, features, tesselation_radians, triangle_base_length, triangle_aspect) # fmt: skip
