"""Tools for reconstructing and plotting geological features and feature data through time.

Methods in `plot.py` reconstruct geological features using 
[pyGPlates' `reconstruct` function](https://www.gplates.org/docs/pygplates/generated/pygplates.reconstruct.html),
turns them into plottable Shapely geometries, and plots them onto 
Cartopy GeoAxes using Shapely and GeoPandas.

Classes
-------
PlotTopologies
"""
import re
import warnings

import pygplates
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import ptt
from shapely.geometry import Point, Polygon
from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
from shapely.ops import linemerge

from .pygplates import FeatureCollection as _FeatureCollection
from .pygplates import _is_string
from .reconstruction import PlateReconstruction as _PlateReconstruction
from .geometry import pygplates_to_shapely
import geopandas as gpd
from .io import (
    get_valid_geometries,  # included for backwards compatibility
    get_geometries as _get_geometries,
)
from .tools import EARTH_RADIUS
from .gpml import _load_FeatureCollection


def plot_subduction_teeth(
    geometries,
    width,
    polarity=None,
    height=None,
    spacing=None,
    projection="auto",
    transform=None,
    ax=None,
    **kwargs
):
    """Add subduction teeth to a plot.

    The subduction polarity used for subduction teeth can be specified
    manually or detected automatically if `geometries` is a
    `geopandas.GeoDataFrame` object with a `polarity` column.

    Parameters
    ----------
    geometries : geopandas.GeoDataFrame, sequence of shapely geometries, or str
        If a `geopandas.GeoDataFrame` is given, its geometry attribute
        will be used. If `geometries` is a string, it must be the path to
        a file, which will be loaded with `geopandas.read_file`. Otherwise,
        `geometries` must be a sequence of shapely geometry objects (instances
        of the `shapely.geometry.base.BaseGeometry` class).
    width : float
        The (approximate) width of the subduction teeth. If a projection is
        used, this value will be in projected units.
    polarity : {"left", "l", "right", "r", None}, default None
        The subduction polarity of the geometries. If no polarity is provided,
        and `geometries` is a `geopandas.GeoDataFrame`, this function will
        attempt to find a `polarity` column in the data frame and use the
        values given there. If `polarity` is not manually specified and no
        appropriate column can be found, an error will be raised.
    height : float, default None
        If provided, the height of the subduction teeth. As with `width`,
        this value should be given in projected units. If no value is given,
        the height of the teeth will be equal to 0.6 * `width`.
    spacing : float, default None
        If provided, the spacing between the subduction teeth. As with
        `width` and `height`, this value should be given in projected units.
        If no value is given, `spacing` will default to `width`, producing
        tightly packed subduction teeth.
    projection : cartopy.crs.Transform, "auto", or None, default "auto"
        The projection of the plot. If the plot has no projection, this value
        can be explicitly given as `None`. The default value is "auto", which
        will acquire the projection automatically from the plot axes.
    transform : cartopy.crs.Transform, or None, default None
        If the plot is projected, a `transform` value is usually needed.
        Frequently, the appropriate value is an instance of
        `cartopy.crs.PlateCarree`.
    ax : matplotlib.axes.Axes, or None, default None
        The axes on which the subduction teeth will be drawn. By default,
        the current axes will be acquired using `matplotlib.pyplot.gca`.
    **kwargs
        Any further keyword arguments will be passed to
        `matplotlib.patches.Polygon`.

    Raises
    ------
    ValueError
        If `width` <= 0, or if `polarity` is an invalid value or could not
        be determined.
    """
    if ax is None:
        ax = plt.gca()

    if projection == "auto":
        try:
            projection = ax.projection
        except AttributeError:
            projection = None
    elif isinstance(projection, str):
        raise ValueError("Invalid projection: {}".format(projection))

    if polarity is None:
        polarity_column = _find_polarity_column(geometries.columns.values)
        if polarity_column is None:
            raise ValueError(
                "Could not automatically determine polarity; "
                + "it must be defined manually instead."
            )
        triangles = []
        for p in geometries[polarity_column].unique():
            if p.lower() not in {"left", "l", "right", "r"}:
                continue
            gdf_polarity = geometries[geometries[polarity_column] == p]
            triangles.extend(
                _tessellate_triangles(
                    gdf_polarity,
                    width,
                    p,
                    height,
                    spacing,
                    projection,
                    transform,
                )
            )
    else:
        triangles = _tessellate_triangles(
            geometries,
            width,
            polarity,
            height,
            spacing,
            projection,
            transform,
        )

    if projection is not None:
        domain = projection.domain
        triangles = [domain.intersection(i) for i in triangles]

    if hasattr(ax, "add_geometries") and projection is not None:
        ax.add_geometries(triangles, crs=projection, **kwargs)
    else:
        for triangle in triangles:
            ax.fill(*triangle.exterior.xy, **kwargs)


def _tessellate_triangles(
    geometries,
    width,
    polarity="left",
    height=None,
    spacing=None,
    projection=None,
    transform=None,
):
    """Generate subduction teeth triangles for plotting.

    Forms continuous trench geometries and identifies their subduction polarities.
    Subduction teeth triangles can be customised with a given spacing and width.  
    Their apexes point in the identified polarity directions.

    Parameters
    ----------
    geometries : geopandas.GeoDataFrame, sequence of shapely geometries, or str
        If a `geopandas.GeoDataFrame` is given, its geometry attribute
        will be used. If `geometries` is a string, it must be the path to
        a file, which will be loaded with `geopandas.read_file`. Otherwise,
        `geometries` must be a sequence of shapely geometry objects (instances
        of the `shapely.geometry.base.BaseGeometry` class).
    width : float
        The (approximate) width of the subduction teeth. If a projection is
        used, this value will be in projected units.
    polarity : {"left", "l", "right", "r", None}, default "left"
        The subduction polarity of the geometries. If no polarity is provided,
        and `geometries` is a `geopandas.GeoDataFrame`, this function will
        attempt to find a `polarity` column in the data frame and use the
        values given there. If `polarity` is not manually specified and no
        appropriate column can be found, an error will be raised.
    height : float, default None
        If provided, the height of the subduction teeth. As with `width`,
        this value should be given in projected units. If no value is given,
        the height of the teeth will be equal to 0.6 * `width`.
    spacing : float, default None
        If provided, the spacing between the subduction teeth. As with
        `width` and `height`, this value should be given in projected units.
        If no value is given, `spacing` will default to `width`, producing
        tightly packed subduction teeth.
    projection : cartopy.crs.Transform, "auto", or None, default None
        The projection of the plot. If the plot has no projection, this value
        can be explicitly given as `None`. The default value is "auto", which
        will acquire the projection automatically from the plot axes.
    transform : cartopy.crs.Transform, or None, default None
        If the plot is projected, a `transform` value is usually needed.
        Frequently, the appropriate value is an instance of
        `cartopy.crs.PlateCarree`.

    Returns
    -------
    results : list of shapely Polygon objects
        Subduction teeth generated for the given geometries. 
    """
    if width <= 0.0:
        raise ValueError("Invalid `width` argument: {}".format(width))
    polarity = _parse_polarity(polarity)
    geometries = _parse_geometries(geometries)
    if height is None:
        height = width * 2.0/3.0
    if spacing is None:
        spacing = width

    if projection is not None:
        geometries_new = []
        for i in geometries:
            geometries_new.extend(
                _project_geometry(i, projection, transform)
            )
        geometries = geometries_new
        del geometries_new
    geometries = linemerge(geometries)
    if isinstance(geometries, BaseMultipartGeometry):
        geometries = list(geometries.geoms)
    elif isinstance(geometries, BaseGeometry):
        geometries = [geometries]
    results = _calculate_triangle_vertices(
        geometries,
        width,
        spacing,
        height,
        polarity,
    )
    return results


def _project_geometry(geometry, projection, transform=None):
    """Project shapely geometries onto a certain Cartopy CRS map projection. 

    Uses a coordinate system ("transform"), if given. 

    Parameters
    ----------
    geometry : shapely.geometry.base.BaseGeometry
        An instance of a shapely geometry.
    projection : cartopy.crs.Transform, 
        The projection of the plot. 
    transform : cartopy.crs.Transform or None, default None
        If the plot is projected, a `transform` value is usually needed.
        Frequently, the appropriate value is an instance of
        `cartopy.crs.PlateCarree`.

    Returns
    -------
    projected : list
        The provided shapely geometries projected onto a Cartopy CRS map projection.
    """
    if transform is None:
        transform = ccrs.PlateCarree()
    result = [projection.project_geometry(geometry, transform)]
    projected = []
    for i in result:
        if isinstance(i, BaseMultipartGeometry):
            projected.extend(list(i.geoms))
        else:
            projected.append(i)
    return projected


def _calculate_triangle_vertices(
    geometries,
    width,
    spacing,
    height,
    polarity,
):
    """Generate vertices of subduction teeth triangles.

    Triangle bases are set on shapely BaseGeometry trench instances with their apexes 
    pointing in directions of subduction polarity. Triangle dimensions are set by a 
    specified width, spacing and height (either provided by the user or set as default
    values from _tessellate_triangles). The teeth are returned as shapely polygons.

    Parameters
    ----------
    geometries : list of shapely geometries (instances of the
        shapely.geometry.base.BaseGeometry or shapely.geometry.base.BaseMultipartGeometry
        class)
        Trench geometries projected onto a certain map projection (using a 
        coordinate system if specified), each with identified subduction polarities. 
        Teeth triangles will be generated only on the BaseGeometry instances. 
    width : float
        The (approximate) width of the subduction teeth. If a projection is
        used, this value will be in projected units.
    spacing : float,
        The spacing between the subduction teeth. As with
        `width` and `height`, this value should be given in projected units.
    height : float, default None
        The height of the subduction teeth. This value should also be given in projected
        units.
    polarity : {"left", "right"}
        The subduction polarity of the shapely geometries. 
    
    Returns
    -------
    triangles : list of shapely polygons
        The subduction teeth generated along the supplied trench geometries. 
    """
    if isinstance(geometries, BaseGeometry):
        geometries = [geometries]
    triangles = []
    for geometry in geometries:
        if not isinstance(geometry, BaseGeometry):
            continue
        length = geometry.length
        tessellated_x = []
        tessellated_y = []
        for distance in np.arange(spacing, length, spacing):
            point = Point(geometry.interpolate(distance))
            tessellated_x.append(point.x)
            tessellated_y.append(point.y)
        tessellated_x = np.array(tessellated_x)
        tessellated_y = np.array(tessellated_y)

        for i in range(len(tessellated_x) - 1):
            normal_x = tessellated_y[i] - tessellated_y[i + 1]
            normal_y = tessellated_x[i + 1] - tessellated_x[i]
            normal = np.array((normal_x, normal_y))
            normal_mag = np.sqrt((normal ** 2).sum())
            if normal_mag == 0:
                continue
            normal *= height / normal_mag
            midpoint = np.array((tessellated_x[i], tessellated_y[i]))
            if polarity == "right":
                normal *= -1.0
            apex = midpoint + normal

            next_midpoint = np.array((tessellated_x[i + 1], tessellated_y[i + 1]))
            line_vector = np.array(next_midpoint - midpoint)
            line_vector_mag = np.sqrt((line_vector ** 2).sum())
            line_vector /= line_vector_mag
            triangle_point_a = midpoint + width * 0.5 * line_vector
            triangle_point_b = midpoint - width * 0.5 * line_vector
            triangle_points = np.array(
                (
                    triangle_point_a,
                    triangle_point_b,
                    apex,
                )
            )
            triangles.append(Polygon(triangle_points))
    return triangles


def _parse_polarity(polarity):
    """Ensure subduction polarities have valid strings as labels - either "left", "l", "right" or "r".

    The geometries' subduction polarities are either provided by the user in plot_subduction_teeth
    or found automatically in a geopandas.GeoDataFrame column by _find_polarity_column, if such a
    column exists.
 
    Parameters
    ----------
    polarity : {"left", "l", "right", "r"}
        The subduction polarity of the geometries (either set by the user or found automatically
        from the geometries' data frame). 

    Returns
    -------
    polarity : {"left", "right"}
        Returned if the provided polarity string is one of {"left", "l", "right", "r"}. "l" and "r"
        are classified and returned as "left" and "right" respectively.  

    Raises
    ------
    TypeError
        If the provided polarity is not a string type.
    ValueError
        If the provided polarity is not valid ("left", "l", "right" or "r").
    """
    if not isinstance(polarity, str):
        raise TypeError(
            "Invalid `polarity` argument type: {}".format(type(polarity))
        )
    if polarity.lower() in {"left", "l"}:
        polarity = "left"
    elif polarity.lower() in {"right", "r"}:
        polarity = "right"
    else:
        valid_args = {"left", "l", "right", "r"}
        err_msg = "Invalid `polarity` argument: {}".format(
            polarity
        ) + "\n(must be one of: {})".format(valid_args)
        raise ValueError(err_msg)
    return polarity


def _find_polarity_column(columns):
    """Search for a 'polarity' column in a geopandas.GeoDataFrame to extract subduction
    polarity values.

    Subduction polarities can be used for tessellating subduction teeth.

    Parameters
    ----------
    columns : geopandas.GeoDataFrame.columns.values instance
        A list of geopandas.GeoDataFrame column header strings.

    Returns
    -------
    column : list
        If found, returns a list of all subduction polarities ascribed to the supplied
        geometry data frame.
    None
        if a 'polarity' column was not found in the data frame. In this case, subduction
        polarities will have to be manually provided to plot_subduction_teeth.

    """
    pattern = "polarity"
    for column in columns:
        if re.fullmatch(pattern, column) is not None:
            return column
    return None


def _parse_geometries(geometries):
    """Resolve a geopandas.GeoSeries object into shapely BaseGeometry and/or
    BaseMutipartGeometry instances.

    Parameters
    ----------
    geometries : geopandas.GeoDataFrame, sequence of shapely geometries, or str
        If a `geopandas.GeoDataFrame` is given, its geometry attribute
        will be used. If `geometries` is a string, it must be the path to
        a file, which will be loaded with `geopandas.read_file`. Otherwise,
        `geometries` must be a sequence of shapely geometry objects (instances
        of the `shapely.geometry.base.BaseGeometry` class).

    Returns
    -------
    out : list
        Resolved shapely BaseMutipartGeometry and/or BaseGeometry instances.
    """
    geometries = _get_geometries(geometries)
    if isinstance(geometries, gpd.GeoSeries):
        geometries = list(geometries)

    # Explode multi-part geometries
    # Weirdly the following seems to be faster than
    # the equivalent explode() method from GeoPandas:
    out = []
    for i in geometries:
        if isinstance(i, BaseMultipartGeometry):
            out.extend(list(i.geoms))
        else:
            out.append(i)
    return out


def shapelify_features(features, central_meridian=0.0, tessellate_degrees=None):
    """Generate Shapely `MultiPolygon` or `MultiLineString` geometries
    from reconstructed feature polygons.
    
    Notes
    -----
    Some Shapely polygons generated by `shapelify_features` cut longitudes of 180 
    or -180 degrees. These features may appear unclosed at the dateline, so Shapely 
    "closes" these polygons by connecting any of their open ends with lines. These 
    lines may manifest on GeoAxes plots as horizontal lines that span the entire 
    global extent. To prevent this, `shapelify_features` uses pyGPlates' 
    [DateLineWrapper](https://www.gplates.org/docs/pygplates/generated/pygplates.datelinewrapper)
    to split a feature polygon into multiple closed polygons if it happens to cut the 
    antimeridian.
    Another measure taken to ensure features are valid is to order exterior coordinates 
    of Shapely polygons anti-clockwise. 

    Parameters
    ----------
    features : iterable of <pygplates.Feature>, <ReconstructedFeatureGeometry> or <GeometryOnSphere>
        Iterable containing reconstructed polygon features.
    central_meridian : float
        Central meridian around which to perform wrapping; default: 0.0.
    tessellate_degrees : float or None
        If provided, geometries will be tessellated to this resolution prior
        to wrapping.

    Returns
    -------
    all_geometries : list of `shapely.geometry.BaseGeometry`
        Shapely geometries converted from the given reconstructed features. Any
        geometries at the dateline are split. 

    See Also
    --------
    geometry.pygplates_to_shapely : convert PyGPlates geometry objects to
    Shapely geometries.
    """
    if isinstance(
        features,
        (
            pygplates.Feature,
            pygplates.ReconstructedFeatureGeometry,
            pygplates.GeometryOnSphere,
        ),
    ):
        features = [features]

    geometries = []
    for feature in features:
        if isinstance(feature, pygplates.Feature):
            geometries.extend(feature.get_all_geometries())
        elif isinstance(feature, pygplates.ReconstructedFeatureGeometry):
            geometries.append(feature.get_reconstructed_geometry())
        elif isinstance(feature, (pygplates.GeometryOnSphere, pygplates.LatLonPoint)):
            geometries.append(feature)
        elif isinstance(feature, pygplates.DateLineWrapper.LatLonMultiPoint):
            geometries.append(
                pygplates.MultiPointOnSphere(
                    [i.to_lat_lon() for i in feature.get_points()]
                )
            )
        elif isinstance(feature, pygplates.DateLineWrapper.LatLonPolyline):
            geometries.append(pygplates.PolylineOnSphere(feature.get_points()))
        elif isinstance(feature, pygplates.DateLineWrapper.LatLonPolygon):
            geometries.append(
                pygplates.PolygonOnSphere(
                    [i.to_lat_lon() for i in feature.get_exterior_points()]
                )
            )

    return [
        pygplates_to_shapely(
            i,
            force_ccw=True,
            validate=True,
            central_meridian=central_meridian,
            tessellate_degrees=tessellate_degrees,
            explode=False,
        )
        for i in geometries
    ]


shapelify_feature_lines = shapelify_features
shapelify_feature_polygons = shapelify_features


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
            raise ValueError("Plate model must have topology features.")

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
        if time is not None:
            self.time = time
        else:
            self._time = None

    def __getstate__(self):

        filenames = self.plate_reconstruction.__getstate__()

        # add important variables from Points object
        if self._coastlines:
            filenames["coastlines"] = self._coastlines.filenames
        if self._continents:
            filenames["continents"] = self._continents.filenames
        if self._COBs:
            filenames["COBs"] = self._COBs.filenames
        filenames['time'] = self.time
        filenames['plate_id'] = self._anchor_plate_id

        return filenames

    def __setstate__(self, state):

        self.plate_reconstruction = _PlateReconstruction(state['rotation_model'], state['topology_features'], state['static_polygons'])

        self._coastlines = None
        self._continents = None
        self._COBs = None
        self.coastlines = None
        self.continents = None
        self.COBs = None

        # reinstate unpicklable items
        if 'coastlines' in state:
            self._coastlines = _FeatureCollection()
            for feature in state['coastlines']:
                self._coastlines.add( _FeatureCollection(feature) )

        if 'continents' in state:
            self._continents = _FeatureCollection()
            for feature in state['continents']:
                self._continents.add( _FeatureCollection(feature) )

        if 'COBs' in state:
            self._COBs = _FeatureCollection()
            for feature in state['COBs']:
                self._COBs.add( _FeatureCollection(feature) )


        self._anchor_plate_id = state["plate_id"]
        self.base_projection = ccrs.PlateCarree()
        self._time = None


    @property
    def time(self):
        """ The reconstruction time."""
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
        if var >= 0:
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
            raise ValueError(
                "Invalid anchor plate ID: {}".format(id)
            )
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
            self.time)

        self.topologies, self.ridge_transforms, self.ridges, self.transforms, self.trenches, self.trench_left, self.trench_right, self.other = resolved_topologies

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
                
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_inferred_paleo_boundary:
                self.inferred_paleo_boundaries.append(topol)

            elif topol.get_feature_type() == pygplates.FeatureType.gpml_terrane_boundary:
                self.terrane_boundaries.append(topol)
                
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_transitional_crust:
                self.transitional_crusts.append(topol)
            
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_orogenic_belt:
                self.orogenic_belts.append(topol)
                
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_suture:
                self.sutures.append(topol)
                
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_continental_crust:
                self.continental_crusts.append(topol)
            
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_extended_continental_crust:
                self.extended_continental_crusts.append(topol)
            
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_passive_continental_boundary:
                self.passive_continental_boundaries.append(topol)
            
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_slab_edge:
                self.slab_edges.append(topol)
                
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_transform:
                self.misc_transforms.append(topol)
                
            elif topol.get_feature_type() == pygplates.FeatureType.gpml_unclassified_feature:
                self.unclassified_features.append(topol)

        # reconstruct other important polygons and lines
        if self._coastlines:
            self.coastlines = self.plate_reconstruction.reconstruct(
                self._coastlines, self.time, from_time=0, anchor_plate_id=self.anchor_plate_id)

        if self._continents:
            self.continents = self.plate_reconstruction.reconstruct(
                self._continents, self.time, from_time=0, anchor_plate_id=self.anchor_plate_id)

        if self._COBs:
            self.COBs = self.plate_reconstruction.reconstruct(
                self._COBs, self.time, from_time=0, anchor_plate_id=self.anchor_plate_id)


    # subduction teeth
    def _tessellate_triangles(self, features, tesselation_radians, triangle_base_length, triangle_aspect=1.0):
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
                    pts = np.array([(p.get_longitude(), p.get_latitude()) for p in line.get_points()])

                    for p in range(0, len(pts) - 1):
                        A = pts[p]
                        B = pts[p+1]

                        AB_dist = B - A
                        AB_norm = AB_dist / np.hypot(*AB_dist)
                        cum_distance += np.hypot(*AB_dist)

                        # create a new triangle if cumulative distance is exceeded.
                        if cum_distance >= tesselation_degrees:

                            C = A + triangle_base_length*AB_norm

                            # find normal vector
                            AD_dist = np.array([AB_norm[1], -AB_norm[0]])
                            AD_norm = AD_dist / np.linalg.norm(AD_dist)

                            C0 = A + 0.5*triangle_base_length*AB_norm

                            # project point along normal vector
                            D = C0 + triangle_base_length*triangle_aspect*AD_norm

                            triangle_pointsX.append( [A[0], C[0], D[0]] )
                            triangle_pointsY.append( [A[1], C[1], D[1]] )

                            cum_distance = 0.0

        return np.array(triangle_pointsX), np.array(triangle_pointsY)


    def get_feature(self, feature):
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
        shp = shapelify_features(feature)
        gdf = gpd.GeoDataFrame({'geometry': shp}, geometry='geometry')
        return gdf

    def plot_feature(self, ax, feature, **kwargs):
        """ Plot pygplates.FeatureCollection  or pygplates.Feature onto a map.

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
        gdf = self.get_feature(feature)
        return gdf.plot(ax=ax, transform=self.base_projection, **kwargs)


    def get_coastlines(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `coastlines` to the requested `time` and thus populate the GeoDataFrame.

        """
        if self._time is None:
            raise ValueError("No coastlines have been resolved. Set `PlotTopologies.time` to construct coastlines.")

        if self.coastlines is None:
            raise ValueError("Supply coastlines to PlotTopologies object")

        coastline_polygons = shapelify_feature_polygons(self.coastlines)
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
        gdf = self.get_coastlines()
        return gdf.plot(ax=ax, transform=self.base_projection, **kwargs)


    def get_continents(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `continents` to the requested `time` and thus populate the GeoDataFrame.

        """
        if self._time is None:
            raise ValueError("No continents have been resolved. Set `PlotTopologies.time` to construct continents.")

        if self.continents is None:
            raise ValueError("Supply continents to PlotTopologies object")

        continent_polygons = shapelify_feature_polygons(self.continents)
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
        gdf = self.get_continents()
        return gdf.plot(ax=ax, transform=self.base_projection, **kwargs)


    def get_continent_ocean_boundaries(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `COBs` to the requested `time` and thus populate the GeoDataFrame.

        """
        if self._time is None:
            raise ValueError("No geometries have been resolved. Set `PlotTopologies.time` to construct topologies.")

        if self.COBs is None:
            raise ValueError("Supply COBs to PlotTopologies object")

        COB_lines = shapelify_feature_lines(self.COBs)
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
        gdf = self.get_continent_ocean_boundaries()
        return gdf.plot(ax=ax, transform=self.base_projection, **kwargs)


    def get_ridges(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `ridges` to the requested `time` and thus populate the GeoDataFrame.

        """
        if self._time is None:
            raise ValueError("No ridges have been resolved. Set `PlotTopologies.time` to construct ridges.")

        if self.ridges is None:
            raise ValueError("No ridge topologies passed to PlotTopologies.")

        ridge_lines = shapelify_feature_lines(self.ridges)
        gdf = gpd.GeoDataFrame({"geometry": ridge_lines}, geometry="geometry")
        return gdf


    def plot_ridges(self, ax, color='black', **kwargs):
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
        gdf = self.get_ridges()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_ridges_and_transforms(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `ridge_transforms` to the requested `time` and thus populate the GeoDataFrame.

        """
        if self._time is None:
            raise ValueError("No ridges and transforms have been resolved. Set `PlotTopologies.time` to construct ridges and transforms.")

        if self.ridge_transforms is None:
            raise ValueError("No ridge and transform topologies passed to PlotTopologies.")

        ridge_transform_lines = shapelify_feature_lines(self.ridge_transforms)
        gdf = gpd.GeoDataFrame({"geometry": ridge_transform_lines}, geometry="geometry")
        return gdf


    def plot_ridges_and_transforms(self, ax, color='black', **kwargs):
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
        gdf = self.get_ridges_and_transforms()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_transforms(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `transforms` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No transforms have been resolved. Set `PlotTopologies.time` to construct transforms.")

        if self.transforms is None:
            raise ValueError("No transform topologies passed to PlotTopologies.")

        transform_lines = shapelify_feature_lines(self.transforms)
        gdf = gpd.GeoDataFrame({"geometry": transform_lines}, geometry="geometry")
        return gdf


    def plot_transforms(self, ax, color='black', **kwargs):
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
        gdf = self.get_transforms()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_trenches(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `trenches` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No trenches have been resolved. Set `PlotTopologies.time` to construct trenches.")

        if self.trenches is None:
            raise ValueError("No trenches passed to PlotTopologies.")

        trench_lines = shapelify_feature_lines(self.trenches)
        gdf = gpd.GeoDataFrame({"geometry": trench_lines}, geometry="geometry")
        return gdf


    def plot_trenches(self, ax, color='black', **kwargs):
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
        gdf = self.get_trenches()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_misc_boundaries(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `other` geometries to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No miscellaneous topologies have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.other is None:
            raise ValueError("No miscellaneous topologies passed to PlotTopologies.")

        lines = shapelify_features(self.other)
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
        gdf = self.get_misc_boundaries()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def plot_subduction_teeth_deprecated(self, ax, spacing=0.1, size=2.0, aspect=1, color='black', **kwargs):
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
            triangle_aspect=-aspect)
        subd_xR, subd_yR = self._tessellate_triangles(
            self.trench_right,
            tesselation_radians=spacing,
            triangle_base_length=size,
            triangle_aspect=aspect)
        
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
            raise ValueError("No miscellaneous topologies have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.trench_left is None or self.trench_right is None:
            raise ValueError("No trench_left or trench_right topologies passed to PlotTopologies.")

        trench_left_features  = shapelify_feature_lines(self.trench_left)
        trench_right_features = shapelify_feature_lines(self.trench_right)

        gdf_left  = gpd.GeoDataFrame({"geometry": trench_left_features},  geometry="geometry")
        gdf_right = gpd.GeoDataFrame({"geometry": trench_right_features}, geometry="geometry")

        return gdf_left, gdf_right


    def plot_subduction_teeth(self, ax, spacing=0.07, size=None, aspect=None, color='black', **kwargs):
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

        size : float, default=None
            Length of teeth triangle base (in radians). If kept at `None`, then
            `size = 0.5*spacing`.

        aspect : float, default=None
            Aspect ratio of teeth triangles. If kept at `None`, then `aspect = 2/3*size`.

        color : str, default=’black’
            The colour of the teeth. By default, it is set to black.

        **kwargs : 
            Keyword arguments parameters such as ‘alpha’, etc. 
            for plotting subduction tooth polygons.
            See `Matplotlib` keyword arguments 
            [here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html).

        Returns
        -------
        ax : instance of <geopandas.GeoDataFrame.plot>
            A standard cartopy.mpl.geoaxes.GeoAxes or cartopy.mpl.geoaxes.GeoAxesSubplot map 
            with subduction teeth plotted onto the chosen map projection.
        """
        if self._time is None:
            raise ValueError("No miscellaneous topologies have been resolved. Set `PlotTopologies.time` to construct them.")

        spacing = spacing * EARTH_RADIUS * 1e3

        if aspect is None:
            aspect = 2.0/3.0
        if size is None:
            size = spacing*0.5

        height = size*aspect

        trench_left_features  = shapelify_feature_lines(self.trench_left)
        trench_right_features = shapelify_feature_lines(self.trench_right)

        return(
            plot_subduction_teeth(trench_left_features,  size, 'l', height, spacing, ax=ax, color=color, **kwargs),
            plot_subduction_teeth(trench_right_features,  size, 'r', height, spacing, ax=ax, color=color, **kwargs)
        )


    def plot_plate_id(self, ax, plate_id, **kwargs):
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
        for feature in self.topologies:
            if feature.get_reconstruction_plate_id() == plate_id:
                ft_plate = shapelify_feature_polygons([feature])
                return ax.add_geometries(ft_plate, crs=self.base_projection, **kwargs)


    def plot_grid(self, ax, grid, extent=[-180,180,-90,90], **kwargs):
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


    def plot_plate_motion_vectors(self, ax, spacingX=10, spacingY=10, normalise=False, **kwargs):
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
        
        lons = np.arange(-180, 180+spacingX, spacingX)
        lats = np.arange(-90, 90+spacingY, spacingY)
        lonq, latq = np.meshgrid(lons, lats)

        # create a feature from all the points
        velocity_domain_features = ptt.velocity_tools.make_GPML_velocity_feature(lonq.ravel(), latq.ravel())

        rotation_model = self.plate_reconstruction.rotation_model
        topology_features = self.plate_reconstruction.topology_features

        delta_time = 5.0
        all_velocities = ptt.velocity_tools.get_plate_velocities(
            velocity_domain_features,
            topology_features,
            rotation_model,
            self.time,
            delta_time,
            'vector_comp')

        X, Y, U, V = ptt.velocity_tools.get_x_y_u_v(lons, lats, all_velocities)

        if normalise:
            mag = np.hypot(U, V)
            mag[mag == 0] = 1
            U /= mag
            V /= mag

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            quiver = ax.quiver(
                X, Y,
                U, V,
                transform=self.base_projection,
                **kwargs
            )
        return quiver


    def get_continental_rifts(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `continental_rifts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No continental rifts have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.continental_rifts is None:
            raise ValueError("No continental rifts passed to PlotTopologies.")

        continental_rift_lines = shapelify_feature_lines(self.continental_rifts)
        gdf = gpd.GeoDataFrame({"geometry": continental_rift_lines}, geometry="geometry")
        return gdf


    def plot_continental_rifts(self, ax, color='black', **kwargs):
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
        gdf = self.get_continental_rifts()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_faults(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `faults` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No faults have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.faults is None:
            raise ValueError("No faults passed to PlotTopologies.")

        fault_lines = shapelify_feature_lines(self.faults)
        gdf = gpd.GeoDataFrame({"geometry": fault_lines}, geometry="geometry")
        return gdf


    def plot_faults(self, ax, color='black', **kwargs):
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
        gdf = self.get_faults()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_fracture_zones(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `fracture_zones` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No fracture zones have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.fracture_zones is None:
            raise ValueError("No fracture zones passed to PlotTopologies.")

        fracture_zone_lines = shapelify_feature_lines(self.fracture_zones)
        gdf = gpd.GeoDataFrame({"geometry": fracture_zone_lines}, geometry="geometry")
        return gdf


    def plot_fracture_zones(self, ax, color='black', **kwargs):
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
        gdf = self.get_fracture_zones()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_inferred_paleo_boundaries(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `inferred_paleo_boundaries` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No inferred paleo boundaries have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.inferred_paleo_boundaries is None:
            raise ValueError("No inferred paleo boundaries passed to PlotTopologies.")

        inferred_paleo_boundary_lines = shapelify_feature_lines(self.inferred_paleo_boundaries)
        gdf = gpd.GeoDataFrame({"geometry": inferred_paleo_boundary_lines}, geometry="geometry")
        return gdf


    def plot_inferred_paleo_boundaries(self, ax, color='black', **kwargs):
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
        gdf = get_inferred_paleo_boundaries()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_terrane_boundaries(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `terrane_boundaries` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No terrane boundaries have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.terrane_boundaries is None:
            raise ValueError("No terrane boundaries passed to PlotTopologies.")

        terrane_boundary_lines = shapelify_feature_lines(self.terrane_boundaries)
        gdf = gpd.GeoDataFrame({"geometry": terrane_boundary_lines}, geometry="geometry")
        return gdf


    def plot_terrane_boundaries(self, ax, color='black', **kwargs):
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
        gdf = self.get_terrane_boundaries()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_transitional_crusts(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `transitional_crusts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No transitional crusts have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.transitional_crusts is None:
            raise ValueError("No transitional crusts passed to PlotTopologies.")

        transitional_crust_lines = shapelify_feature_lines(self.transitional_crusts)
        gdf = gpd.GeoDataFrame({"geometry": transitional_crust_lines}, geometry="geometry")
        return gdf 


    def plot_transitional_crusts(self, ax, color='black', **kwargs):
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
        gdf = self.get_transitional_crusts()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_orogenic_belts(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `orogenic_belts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No orogenic belts have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.orogenic_belts is None:
            raise ValueError("No orogenic belts passed to PlotTopologies.")

        orogenic_belt_lines = shapelify_feature_lines(self.orogenic_belts)
        gdf = gpd.GeoDataFrame({"geometry": orogenic_belt_lines}, geometry="geometry")
        return gdf


    def plot_orogenic_belts(self, ax, color='black', **kwargs):
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
        gdf = self.get_orogenic_belts()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_sutures(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `sutures` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No sutures have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.sutures is None:
            raise ValueError("No sutures passed to PlotTopologies.")

        suture_lines = shapelify_feature_lines(self.sutures)
        gdf = gpd.GeoDataFrame({"geometry": suture_lines}, geometry="geometry")
        return gdf


    def plot_sutures(self, ax, color='black', **kwargs):
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
        gdf = self.get_sutures()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_continental_crusts(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `continental_crusts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No continental crust topologies have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.continental_crusts is None:
            raise ValueError("No continental crust topologies passed to PlotTopologies.")

        continental_crust_lines = shapelify_feature_lines(self.continental_crusts)
        gdf = gpd.GeoDataFrame({"geometry": continental_crust_lines}, geometry="geometry")
        return gdf


    def plot_continental_crusts(self, ax, color='black', **kwargs):
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
        gdf = self.get_continental_crusts()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_extended_continental_crusts(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `extended_continental_crusts` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No extended continental crust topologies have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.extended_continental_crusts is None:
            raise ValueError("No extended continental crust topologies passed to PlotTopologies.")

        extended_continental_crust_lines = shapelify_feature_lines(self.extended_continental_crusts)
        gdf = gpd.GeoDataFrame({"geometry": extended_continental_crust_lines}, geometry="geometry")
        return gdf

    def plot_extended_continental_crusts(self, ax, color='black', **kwargs): 
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
        gdf = self.get_extended_continental_crusts()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_passive_continental_boundaries(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `passive_continental_boundaries` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No passive continental boundaries have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.passive_continental_boundaries is None:
            raise ValueError("No passive continental boundaries passed to PlotTopologies.")

        passive_continental_boundary_lines = shapelify_feature_lines(self.passive_continental_boundaries)
        gdf = gpd.GeoDataFrame({"geometry": passive_continental_boundary_lines}, geometry="geometry")
        return gdf


    def plot_passive_continental_boundaries(self, ax, color='black', **kwargs): 
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
        gdf = self.get_passive_continental_boundaries()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_slab_edges(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `slab_edges` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No slab edges have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.slab_edges is None:
            raise ValueError("No slab edges passed to PlotTopologies.")

        slab_edge_lines = shapelify_feature_lines(self.slab_edges)
        gdf = gpd.GeoDataFrame({"geometry": slab_edge_lines}, geometry="geometry")
        return gdf


    def plot_slab_edges(self, ax, color='black', **kwargs): 
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
        gdf = self.get_slab_edges()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_misc_transforms(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `misc_transforms` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No miscellaneous transforms have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.misc_transforms is None:
            raise ValueError("No miscellaneous transforms passed to PlotTopologies.")

        misc_transform_lines = shapelify_feature_lines(self.misc_transforms)
        gdf = gpd.GeoDataFrame({"geometry": misc_transform_lines}, geometry="geometry")
        return gdf


    def plot_misc_transforms(self, ax, color='black', **kwargs): 
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
        gdf = self.get_misc_transforms()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_unclassified_features(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `unclassified_features` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No unclassified features have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.unclassified_features is None:
            raise ValueError("No unclassified features passed to PlotTopologies.")

        unclassified_feature_lines = shapelify_feature_lines(self.unclassified_features)
        gdf = gpd.GeoDataFrame({"geometry": unclassified_feature_lines}, geometry="geometry")
        return gdf


    def plot_unclassified_features(self, ax, color='black', **kwargs): 
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
        gdf = self.get_unclassified_features()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)


    def get_all_topologies(self):
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

        Raises 
        ------
        ValueError
            If the optional `time` parameter has not been passed to `PlotTopologies`. This is needed to construct
            `topologies` to the requested `time` and thus populate the GeoDataFrame.
        """
        if self._time is None:
            raise ValueError("No topologies have been resolved. Set `PlotTopologies.time` to construct them.")

        if self.topologies is None:
            raise ValueError("No topologies passed to PlotTopologies.")

        all_topologies = shapelify_features(self.topologies)

        # get plate IDs and feature types to add to geodataframe
        plate_IDs = []
        feature_types = []
        feature_names = []
        for topo in self.topologies:
            ft_type = topo.get_feature_type()

            plate_IDs.append(topo.get_reconstruction_plate_id())
            feature_types.append(ft_type)
            feature_names.append(ft_type.get_name())

        gdf = gpd.GeoDataFrame({"geometry": all_topologies,
                                "reconstruction_plate_ID": plate_IDs,
                                "feature_type": feature_types,
                                "feature_name": feature_names},
                                geometry="geometry")
        return gdf


    def plot_all_topologies(self, ax, color='black', **kwargs):
        """Plot all topologies on a standard map projection.

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
        gdf = self.get_all_topologies()
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)
