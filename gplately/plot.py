"""The “plot” module offers simple shortcuts to pyGplates functionalities for reconstructing geological features,
transforming them into Shapely geometries and plotting them onto maps.

Methods in this module perform a combination of reconstructing topological features, generating and sorting through 
valid Shapely geometries, and plotting them. Some methods conduct all the steps above, while others only accept 
reconstructed features as parameters. The methods that reconstruct geometries using pyGplates before plotting include:
    * add_coastlines
    * add_continents
    * add_ridges
The features to reconstruct must be defined above where these methods are called under the variable names ‘coastlines’,
‘continents’ and ‘topology_features’ respectively. 

The methods that assume geometries are reconstructed already are:
    * another version of add_ridges
    * add_trenches
For both methods, reconstructed ridges and topologies are assumed to be defined in the specific shapefile paths:
    * "reconstructed_topologies/ridge_transform_boundaries_{reconstruction_time}Ma.shp"
    * "reconstructed_topologies/subduction_boundaries_{reconstruction_time}Ma.shp"
    * "reconstructed_topologies/subduction_boundaries_sL_{reconstruction_time}Ma.shp"
    * "reconstructed_topologies/subduction_boundaries_sR_{reconstruction_time}Ma.shp"
for a specific reconstruction_time. For example, reconstructed ridges at 40Ma should be contained in the path:
reconstructed_topologies/ridge_transform_boundaries_40Ma.shp

add_quiver calculates and plots plate motion velocity vectors, executing Plate Tectonic Tools’ velocity vector
functionalities in a single line. A rotation model and topology feature must be defined above this method as 
“rotation_model” and “topology_features” respectively (if using GPlately, these would be held in the
PlateReconstruction object).

shapelify_feature_lines and shapelify_feature_polygons turn reconstructed features into Shapely geometries for plotting. 

In the PlotTopologies class, 
    * plot_coastlines
    * plot_continents
    * Plot_continent_ocean_boundaries
are plotting methods that use shapefile strings given to the PlotTopologies object. Moreover,
    * plot_ridges
    * plot_ridges_and_transforms
    * plot_transforms
    * plot_trenches
    * plot_subduction_teeth
are plotting methods that use the PlotTopologies feature attributes: self.topologies, self.ridge_transforms,
self.ridges, self.transforms, self.trenches, self.trench_left, self.trench_right and self.other. These features
have been resolved using Plate Tectonics Tools’ ptt.resolve_topologies.resolve_topologies_into_features method. 

The following methods:
    *plot_grid
    *plot_grid_from_netCDF
plot ndarray-like and netCDF grids respectively onto maps. The latter method uses GPlately’s grids module to
extract a grid from a given netCDF file path.

The method ‘plot_plate_motion_vectors’ constructs and plots plate motion velocity vectors using Plate Tectonics Tools functionalities. 

Classes
-------
PlotTopologies
"""
import re

import pygplates
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import ptt
from shapely.geometry import Point, Polygon
from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
from shapely.ops import linemerge

from .geometry import pygplates_to_shapely
import geopandas as gpd
from .io import (
    get_valid_geometries,  # included for backwards compatibility
    get_geometries as _get_geometries,
)
from .tools import EARTH_RADIUS

    
def add_coastlines(ax, reconstruction_time, **kwargs):
    """Reconstructs coastline geometries and plots them onto a standard map. 

    The coastlines for plotting must be defined in a variable "coastlines" above the method, along with a "rotation_model".
    "coastlines" should be a feature collection, or filename, or feature, or sequence of features, 
    or a list or tuple of any combination of those four types. These coastlines are reconstructed, transformed into shapely
    geometries and added onto the chosen map. Map presentation details (e.g. facecolor, edgecolor, alpha…) permitted.

    Parameters
    ----------
    ax : GeoAxis
        A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots. 
        Should be set at a particular Cartopy map projection.

    reconstruction_time : float
        A particular geological time (Ma) at which to reconstruct the feature coastlines. 

    **export_wrap_to_dateline : bool, default=True
        Wrap/clip reconstructed geometries to the dateline.

    **kwargs 
        Keyword arguments that allow control over parameters such as ‘facecolor’, ‘alpha’, etc. for plotting coastline 
        geometries.

    Returns
    -------
    ax : GeoAxis
        The map with coastline features plotted onto the chosen projection. 
    """
    # write shapefile
    reconstructed_coastlines = []
    pygplates.reconstruct(coastlines, rotation_model, reconstructed_coastlines, float(reconstruction_time),
                          export_wrap_to_dateline=True)
    coastlines_geometries = shapelify_feature_polygons(reconstructed_coastlines)
    ax.add_geometries(coastlines_geometries, crs=ccrs.PlateCarree(), **kwargs)
    
def add_continents(ax, reconstruction_time, **kwargs):
    """Reconstructs continental geometries and plots them onto a standard map. 

    The continents for plotting must be defined in a variable "continents" above the method, along with a "rotation_model".
    "continents" should be a feature collection, or filename, or feature, or sequence of features, or a list or tuple 
    of any combination of those four types. These coastlines are reconstructed, transformed into shapely
    geometries and added onto the chosen map. Map presentation details (e.g. facecolor, edgecolor, alpha…) permitted.

    Parameters
    ----------
    ax : GeoAxis
        A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots. 
        Should be set at a particular Cartopy map projection.

    reconstruction_time : float
        A particular geological time (Ma) at which to reconstruct the feature continents. 

    **export_wrap_to_dateline : bool, default=True
        Wrap/clip reconstructed geometries to the dateline.

    **kwargs 
        Keyword arguments that allow control over parameters such as ‘facecolor’, ‘alpha’, etc. for plotting continental
        geometries.

    Returns
    -------
    ax : GeoAxis
        The map with continental features plotted onto the chosen projection. 
    """
    reconstructed_continents = []
    pygplates.reconstruct(continents, rotation_model, reconstructed_continents, float(reconstruction_time),
                          export_wrap_to_dateline=True)
    continent_geometries = shapelify_feature_polygons(reconstructed_continents)
    ax.add_geometries(continent_geometries, crs=ccrs.PlateCarree(), **kwargs)

    
def add_ridges(ax, reconstruction_time, **kwargs):
    """Reconstructs ridge features to a specific geological time, converts them into shapely geometries and plots them 
    onto a standard map.
    
    The ridges for plotting must be defined in a variable "topology_features" above the method, along with a "rotation_model".
    "topology_features" must be a feature collection, or filename, or feature, or sequence of features, 
    or a list or tuple of any combination of those four types. Map presentation details (e.g. facecolor, edgecolor, alpha…)
    permitted.

    Exterior coordinates are ordered anti-clockwise and only valid geometries are passed to ensure compatibility with Cartopy.

    Parameters
    ----------
    ax : GeoAxis
        A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots. 
        Should be set at a particular Cartopy map projection.

    reconstruction_time : float
        A particular geological time (Ma) at which to reconstruct the feature ridges. 

    **kwargs 
        Keyword arguments that allow control over parameters such as ‘facecolor’, ‘alpha’, etc. for plotting ridge geometries.

    Returns
    -------
    ax
        The map with ridge features plotted onto the chosen projection. 

    """
    import shapely
    reconstructed_ridges = get_ridge_transforms(topology_features, rotation_model, float(reconstruction_time))
    all_geometries = []
    for feature in reconstructed_ridges:
        geometry = feature.get_all_geometries()[0].to_lat_lon_array()[::-1,::-1]
        
        # construct shapely geometry
        geom = shapely.geometry.LineString(geometry)

        # we need to make sure the exterior coordinates are ordered anti-clockwise
        # and the geometry is valid otherwise it will screw with cartopy
        if geom.is_valid:
            all_geometries.append(geom)
    
    ax.add_geometries(all_geometries, crs=ccrs.PlateCarree(), **kwargs)
    
def add_ridges(ax, reconstruction_time, **kwargs):
    """Reads a shapefile containing reconstructed ridge features. Plots them onto a standard map.

    Ridge features to be plotted must be reconstructed to a specific geological time already. These are assumed to be 
    held in the filename string
    "reconstructed_topologies/ridge_transform_boundaries_{reconstruction_time}Ma.shp" for each reconstruction time. This
    shape filename is read and its ridge geometries are extracted, turned into shapely geometries, and plotted onto the 
    given map.

    Parameters
    ----------
    ax : GeoAxis
        A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
        Should be set at a particular Cartopy map projection.

    reconstruction_time : float
        A particular geological time (Ma) at which to reconstruct the feature ridges. 

    **kwargs 
        Keyword arguments that allow control over parameters such as ‘facecolor’, ‘alpha’, etc. for plotting ridge 
        geometries.

    Returns
    -------
    ax : GeoAxis
        The map with ridge features plotted onto the chosen projection.
    """
    shp_name = "reconstructed_topologies/ridge_transform_boundaries_{:.2f}Ma.shp".format(reconstruction_time)
    shp_continents = shpreader.Reader(shp_name).geometries()
    ft_continents  = cfeature.ShapelyFeature(shp_continents, ccrs.PlateCarree())
    ax.add_feature(ft_continents, **kwargs)

def add_trenches(ax, reconstruction_time, color='k', linewidth=2, **kwargs):
    """Reads 3 shapefile containing reconstructed subduction boundaries, left subduction trenches and right subduction
    trenches. Generates subduction teeth along these boundaries and plots both subduction zones and teeth onto a standard map.

    The subduction boundary features to be plotted should be reconstructed to a specific geological time already. They are
    assumed to be held in the filename string "reconstructed_topologies/subduction_boundaries_{reconstruction_time}Ma.shp"
    for a specific reconstruction time. 
    
    The left subduction features are assumed to be held in the filename string: 
    "reconstructed_topologies/subduction_boundaries_sL_{reconstruction_time}Ma.shp", while the right subduction features are
    assumed to be held in the filename string: "reconstructed_topologies/subduction_boundaries_sR_{reconstruction_time}Ma.shp".
    These left and right features are tessellated into subduction teeth. 

    Parameters
    ----------
    ax : GeoAxis
        A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
        Should be set at a particular Cartopy map projection.

    reconstruction_time : float
        A particular geological time (Ma) at which to reconstruct the feature ridges.

    color : str, default='k'
        Sets the subduction boundaries to a specific edge color.

    linewidth : float, default=2
        Defines the thickness of the subduction boundaries 

    **kwargs 
        Keyword arguments that allow control over map features like ‘alpha’, etc. for plotting trench geometries.

    Returns
    -------
    ax : GeoAxis
        The map with subduction boundary features & subduction teeth plotted onto the chosen projection.
    """
    shp_name = "reconstructed_topologies/subduction_boundaries_{:.2f}Ma.shp".format(reconstruction_time)
    shp_subd = shpreader.Reader(shp_name).geometries()
    ft_subd  = cfeature.ShapelyFeature(shp_subd, ccrs.PlateCarree())
    ax.add_feature(ft_subd, facecolor='none', edgecolor=color, linewidth=linewidth, zorder=5)
    # add Subduction Teeth
    subd_xL, subd_yL = tesselate_triangles(
        "reconstructed_topologies/subduction_boundaries_sL_{:.2f}Ma.shp".format(reconstruction_time),
        tesselation_radians=0.1, triangle_base_length=2.0, triangle_aspect=-1.0)
    subd_xR, subd_yR = tesselate_triangles(
        "reconstructed_topologies/subduction_boundaries_sR_{:.2f}Ma.shp".format(reconstruction_time),
        tesselation_radians=0.1, triangle_base_length=2.0, triangle_aspect=1.0)
    
    for tX, tY in zip(subd_xL, subd_yL):
        triangle_xy_points = np.c_[tX, tY]
        patch = plt.Polygon(triangle_xy_points, color=color, transform=ccrs.PlateCarree(), zorder=6)
        ax.add_patch(patch)
    for tX, tY in zip(subd_xR, subd_yR):
        triangle_xy_points = np.c_[tX, tY]
        patch = plt.Polygon(triangle_xy_points, color=color, transform=ccrs.PlateCarree(), zorder=6)
        ax.add_patch(patch)
    
    
def add_quiver(ax, reconstruction_time, **kwargs):
    """Calculates velocity data for a specific reconstruction time and plots velocity vectors onto a map. 

    A rotation model and topology feature must be defined above this method as “rotation_model” and “topology_features”
    respectively. These are used to generate velocity domain feature points. Velocities at each domain point are calculated
    for a specific geological time. X and Y nodes (coordinates of domain points) and u and v components of velocity are
    generated to plot velocity vectors onto the given map. Velocity magnitudes are not normalised. 

    Parameters
    -----------
    ax : GeoAxis
        A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
        Should be set at a particular Cartopy map projection.

    reconstruction_time : float
        A particular geological time (Ma) at which to reconstruct the feature ridges.

    **kwargs 
        Keyword arguments that allow control over quiver features like ‘alpha’, etc. for velocity vectors.

    Returns
    -------
    ax : GeoAxis
        The map with velocity vectors plotted onto the chosen projection.
    """
    Xnodes, Ynodes, U, V = ptt.velocity_tools.get_velocity_x_y_u_v(reconstruction_time, rotation_model,
                                                                   topology_features)
    mag = np.hypot(U, V)
#     mag = np.clip(mag, 1.0, 1e99)
#     mag[mag==0] = 1 #to avoid 0 divisor
#     U = U/mag
#     V = V/mag
    
    if mag.any():
        ax.quiver(Xnodes, Ynodes, U, V, transform=ccrs.PlateCarree(), **kwargs)



# subduction teeth
def tesselate_triangles(shapefilename, tesselation_radians, triangle_base_length, triangle_aspect=1.0):
    """Places subduction teeth along subduction boundary line segments within a MultiLineString shapefile. 

    Parameters
    ----------
    shapefilename  : str  
       Path to shapefile containing the reconstructed subduction boundary features.

    tesselation_radians : float
        Parametrises subduction teeth density. Triangles are generated only along line segments with distances that 
        exceed the given threshold tessellation_radians.

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

    import shapefile

    with shapefile.Reader(shapefilename) as shp:
        tesselation_degrees = np.degrees(tesselation_radians)
        triangle_pointsX = []
        triangle_pointsY = []

        for i in range(len(shp)):
            pts = np.array(shp.shape(i).points)

            cum_distance = 0.0
            for p in range(len(pts) - 1):

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
                _tesselate_triangles(
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
        triangles = _tesselate_triangles(
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


def _tesselate_triangles(
    geometries,
    width,
    polarity="left",
    height=None,
    spacing=None,
    projection=None,
    transform=None,
):
    """Generates subduction teeth triangles for plotting.

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
    """Projects shapely geometries onto a certain Cartopy CRS map projection. 

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
    """Generates vertices of subduction teeth triangles.

    Triangle bases are set on shapely BaseGeometry trench instances with their apexes 
    pointing in directions of subduction polarity. Triangle dimensions are set by a 
    specified width, spacing and height (either provided by the user or set as default
    values from _tesselate_triangles). The teeth are returned as shapely polygons.

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
        tesselated_x = []
        tesselated_y = []
        for distance in np.arange(spacing, length, spacing):
            point = Point(geometry.interpolate(distance))
            tesselated_x.append(point.x)
            tesselated_y.append(point.y)
        tesselated_x = np.array(tesselated_x)
        tesselated_y = np.array(tesselated_y)

        for i in range(len(tesselated_x) - 1):
            normal_x = tesselated_y[i] - tesselated_y[i + 1]
            normal_y = tesselated_x[i + 1] - tesselated_x[i]
            normal = np.array((normal_x, normal_y))
            normal_mag = np.sqrt((normal ** 2).sum())
            if normal_mag == 0:
                continue
            normal *= height / normal_mag
            midpoint = np.array((tesselated_x[i], tesselated_y[i]))
            if polarity == "right":
                normal *= -1.0
            apex = midpoint + normal

            next_midpoint = np.array((tesselated_x[i + 1], tesselated_y[i + 1]))
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
    """Ensures subduction polarities are valid strings - either "left", "l", "right" or "r".

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
    """Searches for a 'polarity' column in a geopandas.GeoDataFrame to extract subduction
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
    """Resolves a geopandas.GeoSeries object into shapely BaseGeometry and/or
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
    """Generates shapely MultiPolygon or MultiLineString geometries
    from reconstructed feature polygons.
    
    Wraps geometries around a central meridian by splitting a polygon into
    multiple polygons at the antimeridian. This is to avoid horizontal lines
    being formed between geometries at longitudes of -180 and 180 degrees.
    Exterior coordinates are ordered anti-clockwise and only valid geometries
    are passed to ensure compatibility with Cartopy.

    Parameters
    ----------
    features : iterable of `pygplates.Feature`, `ReconstructedFeatureGeometry`,
    or `GeometryOnSphere`
        Iterable containing reconstructed polygon features.
    central_meridian : float
        Central meridian around which to perform wrapping; default: 0.0.
    tessellate_degrees : float or None
        If provided, geometries will be tessellated to this resolution prior
        to wrapping.

    Returns
    -------
    all_geometries : list of `shapely.geometry.BaseGeometry`
        Shapely geometries converted from the given reconstructed features.

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
    """Provides methods to read, reconstruct and plot geological topology features.

    Accesses PyGplates and Shapely functionalities to reconstruct feature shapefiles to specific geological times. A
    variety of geological features can be plotted on GeoAxes maps as Shapely MultiLineString or MultiPolygon geometries,
    including subduction trenches & subduction teeth, mid-ocean ridge boundaries, transform boundaries, coastlines,
    continents and COBs. Moreover, netCDF4 and ndarray grids can be plotted onto GeoAxes maps. Plate motion velocity
    vectors can be generated and plotted using this object.

    The rotation models and topology features are used with pyGplates to extract resolved topology features, ridge and
    transform boundary sections, ridge boundary sections, transform boundary sections, subduction boundary sections,
    left subduction boundary sections, right subduction boundary sections and other boundary sections that are not
    subduction zones or mid-ocean ridges (ridge/transform). The continents, coastlines and COBs are read from file
    names given to this class.


    Attributes
    ----------
    self.PlateReconstruction_object
    self.base_projection
    self.coastline_filename
    self.continent_filename
    self.COB_filename
    self.time
    self.topologies
    self.ridge_transforms
    self.ridges
    self.transforms
    self.trenches
    self.trench_left
    self.trench_right
    self.other

    Methods
    ---------
    __init__(self, PlateReconstruction_object, time, coastline_filename=None, continent_filename=None, COB_filename=None)
        Constructs all necessary attributes for the PlotTopologies object. 
        
    update_time(self, time)
        
    _tesselate_triangles(self, features, tesselation_radians, triangle_base_length, triangle_aspect=1.0)
        Places subduction teeth along subduction boundary line segments within a MultiLineString shapefile.
        
    plot_coastlines(self, ax, **kwargs)
    plot_continents(self, ax, **kwargs)
    plot_continent_ocean_boundaries(self, ax, **kwargs)
    plot_ridges(self, ax, color='black', **kwargs)
    plot_ridges_and_transforms(self, ax, color='black', **kwargs)
    plot_transforms(self, ax, color='black', **kwargs)
    plot_trenches(self, ax, color='black', **kwargs)
    plot_subduction_teeth(self, ax, spacing=0.1, size=2.0, aspect=1, color='black', **kwargs)
    plot_grid(self, ax, grid, extent=[-180,180,-90,90], **kwargs)
    plot_grid_from_netCDF(self, ax, filename, **kwargs)
    plot_plate_motion_vectors(self, ax, spacingX=10, spacingY=10, normalise=False, **kwargs)
    """
    def __init__(self, PlateReconstruction_object, time, coastline_filename=None, continent_filename=None, COB_filename=None):
        """Constructs all necessary attributes for the PlateTopologies object.
        """
        import ptt
        import cartopy.crs as ccrs

        self.PlateReconstruction_object = PlateReconstruction_object
        self.base_projection = ccrs.PlateCarree()

        self.coastline_filename = coastline_filename
        self.continent_filename = continent_filename
        self.COB_filename = COB_filename

        # store topologies for easy access
        # setting time runs the update_time routine
        self.time = time

    @property
    def time(self):
        """ Reconstruction time """
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


    def update_time(self, time):
        """Reconstructs features and resolves important topologies at the set reconstruction time. Updates all
        reconstructions/resolutions accordingly whenever the time attribute is changed.

        Returns
        -------
        The following class attributes are generated and updated whenever a reconstruction time attribute is set:
        topologies, ridge_transforms, ridges, transforms, trenches, trench_left, trench_right, or other : lists 
            Common topology features resolved to the set geological time. Specifically, these are lists of:
                - resolved topology features (topological plates and networks)
                - ridge and transform boundary sections (resolved features)
                - ridge boundary sections (resolved features)
                - transform boundary sections (resolved features)
                - subduction boundary sections (resolved features)
                - left subduction boundary sections (resolved features)
                - right subduction boundary sections (resolved features)
                - other boundary sections (resolved features) that are not subduction zones or mid-ocean ridges 
                (ridge/transform)

        coastlines, continents, COBs : lists
            Reconstructed features appended to a Python list. 
        """
        self._time = float(time)
        resolved_topologies = ptt.resolve_topologies.resolve_topologies_into_features(
            self.PlateReconstruction_object.rotation_model,
            self.PlateReconstruction_object.topology_features,
            self.time)

        self.topologies, self.ridge_transforms, self.ridges, self.transforms, self.trenches, self.trench_left, self.trench_right, self.other = resolved_topologies

        # reconstruct other important polygons and lines
        if self.coastline_filename:
            self.coastlines = self.PlateReconstruction_object.reconstruct(
                self.coastline_filename, self.time, from_time=0, anchor_plate_id=0)

        if self.continent_filename:
            self.continents = self.PlateReconstruction_object.reconstruct(
                self.continent_filename, self.time, from_time=0, anchor_plate_id=0)

        if self.COB_filename:
            self.COBs = self.PlateReconstruction_object.reconstruct(
                self.COB_filename, self.time, from_time=0, anchor_plate_id=0)


    # subduction teeth
    def _tesselate_triangles(self, features, tesselation_radians, triangle_base_length, triangle_aspect=1.0):
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

    def plot_coastlines(self, ax, **kwargs):
        """Plots reconstructed coastline polygons onto a standard map. 

        The reconstructed coastlines for plotting must be contained in a shape filename. These coastlines are transformed
        into shapely geometries and added onto the chosen map for a specific geological time (supplied to the PlotTopologies
        object). Map presentation details (e.g. facecolor, edgecolor, alpha…) are permitted.

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        **kwargs 
            Keyword arguments that allow control over parameters such as ‘facecolor’, ‘alpha’, etc. for plotting 
            coastline geometries.

        Returns
        -------
        cartopy.mpl.feature_artist.FeatureArtist
            The `FeatureArtist` instance responsible for drawing the geometries.
        """
        if self.coastline_filename is None:
            raise ValueError("Supply coastline_filename to PlotTopologies object")

        coastline_polygons = shapelify_feature_polygons(self.coastlines)
        gdf = gpd.GeoDataFrame({"geometry": coastline_polygons}, geometry="geometry")
        return gdf.plot(ax=ax, transform=self.base_projection, **kwargs)

    def plot_continents(self, ax, **kwargs):
        """Plots reconstructed continental polygons onto a standard map. 

        The reconstructed continents for plotting must be contained in a shape filename. These continents are transformed
        into shapely geometries and added onto the chosen map for a specific geological time (supplied to the PlotTopologies
        object). Map presentation details (e.g. facecolor, edgecolor, alpha…) are permitted.

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        **kwargs 
            Keyword arguments that allow control over parameters such as ‘facecolor’, ‘alpha’, etc. for plotting continental
            geometries.

        Returns
        -------
        cartopy.mpl.feature_artist.FeatureArtist
            The `FeatureArtist` instance responsible for drawing the geometries.
        """
        if self.continent_filename is None:
            raise ValueError("Supply continent_filename to PlotTopologies object")

        continent_polygons = shapelify_feature_polygons(self.continents)
        gdf = gpd.GeoDataFrame({"geometry": continent_polygons}, geometry="geometry")
        return gdf.plot(ax=ax, transform=self.base_projection, **kwargs)

    def plot_continent_ocean_boundaries(self, ax, **kwargs):
        """Plots reconstructed continent-ocean boundary (COB) polygons onto a standard map. 

        The reconstructed COBs for plotting must be contained in a shape filename. These COBs are transformed into shapely
        geometries and added onto the chosen map for a specific geological time (supplied to the PlotTopologies object). Map
        presentation details (e.g. facecolor, edgecolor, alpha…) are permitted.

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        **kwargs 
            Keyword arguments that allow control over parameters such as ‘facecolor’, ‘alpha’, etc. for plotting COBs.

        Returns
        -------
        cartopy.mpl.feature_artist.FeatureArtist
            The `FeatureArtist` instance responsible for drawing the geometries.
        """
        if self.COB_filename is None:
            raise ValueError("Supply COB_filename to PlotTopologies object")

        COB_lines = shapelify_feature_lines(self.COBs)
        gdf = gpd.GeoDataFrame({"geometry": COB_lines}, geometry="geometry")
        return gdf.plot(ax=ax, transform=self.base_projection, **kwargs)

    def plot_ridges(self, ax, color='black', **kwargs):
        """Plots reconstructed ridge polylines onto a standard map. 

        The reconstructed ridges for plotting must be a list held in the ‘ridges’ attribute. These ridges are transformed 
        into shapely geometries and added onto the chosen map for a specific geological time (supplied to the PlotTopologies
        object). Map presentation details (e.g. color, alpha…) are permitted.

        Note: Ridge geometries are wrapped to the dateline by splitting a polyline into multiple polylines at the dateline.
        This is to avoid horizontal lines being formed between polylines at longitudes of -180 and 180 degrees. Lat-lon 
        feature points near the poles (-89 & 89 latitude) are clipped to ensure compatibility with Cartopy. 

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        color : str, default=’black’
            The colour of the ridge lines. By default, it is set to black.

        **kwargs 
            Keyword arguments that allow control over parameters such as ‘alpha’, etc. for plotting ridge geometries.

        Returns
        -------
        cartopy.mpl.feature_artist.FeatureArtist
            The `FeatureArtist` instance responsible for drawing the geometries.
        """
        ridge_lines = shapelify_feature_lines(self.ridges)
        gdf = gpd.GeoDataFrame({"geometry": ridge_lines}, geometry="geometry")
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)

    def plot_ridges_and_transforms(self, ax, color='black', **kwargs):
        """Plots reconstructed ridge & transform boundary polylines onto a standard map. 

        The reconstructed ridge & transform sections for plotting must be a list held in the ‘ridge_transforms’ attribute. 
        These sections are transformed into shapely geometries and added onto the chosen map for a specific geological time
        (supplied to the PlotTopologies object). Map presentation details (e.g. color, alpha…) are permitted.

        Note: Ridge & transform geometries are wrapped to the dateline by splitting a polyline into multiple polylines at 
        the dateline. This is to avoid horizontal lines being formed between polylines at longitudes of -180 and 180 degrees.
        Lat-lon feature points near the poles (-89 & 89 latitude) are clipped to ensure compatibility with Cartopy. 

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        color : str, default=’black’
            The colour of the ridge & transform lines. By default, it is set to black.

        **kwargs 
            Keyword arguments that allow control over parameters such as ‘alpha’, etc. for plotting ridge & transform geometries.

        Returns
        -------
        cartopy.mpl.feature_artist.FeatureArtist
            The `FeatureArtist` instance responsible for drawing the geometries.
        """
        ridge_transform_lines = shapelify_feature_lines(self.ridge_transforms)
        gdf = gpd.GeoDataFrame({"geometry": ridge_transform_lines}, geometry="geometry")
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)

    def plot_transforms(self, ax, color='black', **kwargs):
        """Plots reconstructed transform boundary polylines onto a standard map. 

        The reconstructed transform sections for plotting must be a list held in the ‘transforms’ attribute. These sections
        are transformed into shapely geometries and added onto the chosen map for a specific geological time (supplied to
        the PlotTopologies object). Map presentation details (e.g. color, alpha…) are permitted.

        Note: Transform geometries are wrapped to the dateline by splitting a polyline into multiple polylines at the dateline.
        This is to avoid horizontal lines being formed between polylines at longitudes of -180 and 180 degrees. Lat-lon feature
        points near the poles (-89 & 89 latitude) are clipped to ensure compatibility with Cartopy. 

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        color = str, default=’black’
            The colour of the transform lines. By default, it is set to black.

        **kwargs 
            Keyword arguments that allow control over parameters such as ‘alpha’, etc. for plotting transform geometries.

        Returns
        -------
        cartopy.mpl.feature_artist.FeatureArtist
            The `FeatureArtist` instance responsible for drawing the geometries.
        """
        transform_lines = shapelify_feature_lines(self.transforms)
        gdf = gpd.GeoDataFrame({"geometry": transform_lines}, geometry="geometry")
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)

    def plot_trenches(self, ax, color='black', **kwargs):
        """Plots reconstructed subduction trench polylines onto a standard map. 

        The reconstructed trenches for plotting must be a list held in the ‘trenches’ attribute. These sections are 
        transformed into shapely geometries and added onto the chosen map for a specific geological time (supplied to the
        PlotTopologies object). Map presentation details (e.g. color, alpha…) are permitted.

        Note: Trench geometries are wrapped to the dateline by splitting a polyline into multiple polylines at the dateline.
        This is to avoid horizontal lines being formed between polylines at longitudes of -180 and 180 degrees. Lat-lon 
        feature points near the poles (-89 & 89 latitude) are clipped to ensure compatibility with Cartopy. 

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        color = str, default=’black’
            The colour of the trench lines. By default, it is set to black.

        **kwargs 
            Keyword arguments that allow control over parameters such as ‘alpha’, etc. for plotting trench geometries.

        Returns
        -------
        cartopy.mpl.feature_artist.FeatureArtist
            The `FeatureArtist` instance responsible for drawing the geometries.
        """
        trench_lines = shapelify_feature_lines(self.trenches)
        gdf = gpd.GeoDataFrame({"geometry": trench_lines}, geometry="geometry")
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)

    def plot_misc_boundaries(self, ax, color="black", **kwargs):
        """Plots reconstructed miscellaneous plate boundary polylines onto a standard map.

        The reconstructed boundaries for plotting must be a list held in the `other` attribute. These sections are
        transformed into shapely geometries and added onto the chosen map for a specific geological time (supplied to the
        PlotTopologies object). Map presentation details (e.g. color, alpha…) are permitted.

        Note: Boundary geometries are wrapped to the dateline by splitting a polyline into multiple polylines at the dateline.
        This is to avoid horizontal lines being formed between polylines at longitudes of -180 and 180 degrees. Lat-lon
        feature points near the poles (-89 & 89 latitude) are clipped to ensure compatibility with Cartopy.

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        color = str, default=’black’
            The colour of the boundary lines. By default, it is set to black.

        **kwargs
            Keyword arguments that allow control over parameters such as ‘alpha’, etc. for plotting boundary geometries.

        Returns
        -------
        cartopy.mpl.feature_artist.FeatureArtist
            The `FeatureArtist` instance responsible for drawing the geometries.
        """
        lines = shapelify_features(self.other)
        gdf = gpd.GeoDataFrame({"geometry": lines}, geometry="geometry")
        return gdf.plot(ax=ax, facecolor='none', edgecolor=color, transform=self.base_projection, **kwargs)

    def plot_subduction_teeth_deprecated(self, ax, spacing=0.1, size=2.0, aspect=1, color='black', **kwargs):
        """Plots subduction teeth onto a standard map. 

        To plot subduction teeth, the left and right sides of resolved subduction boundary sections must be accessible
        from the “left_trench” and “right_trench” attributes. Teeth are created and transformed into Shapely polygons for
        plotting. 

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        spacing : float, default=0.1 
            The tessellation threshold (in radians). Parametrises subduction tooth density. Triangles are generated only
            along line segments with distances that exceed the given threshold ‘spacing’.

        size : float, default=2.0
            Length of teeth triangle base.

        aspect : float, default=1
            Aspect ratio of teeth triangles. Ratio is 1.0 by default. 

        color = str, default=’black’
            The colour of the teeth. By default, it is set to black.

        **kwargs 
            Keyword arguments that allow control over parameters such as ‘alpha’, etc. for plotting subduction tooth polygons.

        Returns
        -------
        ax : GeoAxis
            The map with subduction teeth plotted onto the chosen projection (transformed using PlateCarree).
        """
        import shapely

        # add Subduction Teeth
        subd_xL, subd_yL = self._tesselate_triangles(
            self.trench_left,
            tesselation_radians=spacing,
            triangle_base_length=size,
            triangle_aspect=-aspect)
        subd_xR, subd_yR = self._tesselate_triangles(
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


    def plot_subduction_teeth(self, ax, spacing=0.07, size=None, aspect=None, color='black', **kwargs):
        """Plots subduction teeth onto a standard map. 

        To plot subduction teeth, the left and right sides of resolved subduction boundary sections must be accessible
        from the “left_trench” and “right_trench” attributes. Teeth are created and transformed into Shapely polygons for
        plotting. 

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        spacing : float, default=0.1 
            The tessellation threshold (in radians). Parametrises subduction tooth density. Triangles are generated only
            along line segments with distances that exceed the given threshold ‘spacing’.

        size : float, default=None
            Length of teeth triangle base (in radians). If None then size=0.5*spacing

        aspect : float, default=None
            Aspect ratio of teeth triangles. If None then aspect=2/3*size

        color = str, default=’black’
            The colour of the teeth. By default, it is set to black.

        **kwargs 
            Keyword arguments that allow control over parameters such as ‘alpha’, etc. for plotting subduction tooth polygons.
        """

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
        """Plot a plate ID"""
        for feature in self.topologies:
            if feature.get_reconstruction_plate_id() == plate_id:
                ft_plate = shapelify_feature_polygons([feature])
                return ax.add_geometries(ft_plate, crs=self.base_projection, **kwargs)


    def plot_grid(self, ax, grid, extent=[-180,180,-90,90], **kwargs):
        """Plots an ndarray of gridded data onto a standard map. 

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        grid : ndarray
            An array with elements that define a grid. The number of rows corresponds to the number of latitudinal points, 
            while the number of columns corresponds to the number of longitudinal points.

        extent : 1d array, default=[-180,180,-90,90]
            A four-element array to specify the [min lon, max lon, min lat, max lat] with which to constrain the grid image.
            If no extents are supplied, full global extent is assumed. 

        **kwargs
            Keyword arguments that allow control over map presentation parameters such as ‘alpha’, etc. for plotting the grid.

        Returns
        -------
        ax : GeoAxis
            The map with the grid plotted onto the chosen projection (transformed using PlateCarree).
        """
        return ax.imshow(grid, origin='lower', extent=extent, transform=self.base_projection, **kwargs)


    def plot_grid_from_netCDF(self, ax, filename, **kwargs):
        """Reads in a netCDF file, extracts a raster grid and its associated latitude and longitude arrays and plots the
        raster grid onto a standard map. 

        The lat-lon arrays are used to define the extent of the grid image.

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        filename : str
            Path to a netCDF filename.

        **kwargs
            Keyword arguments that allow control over map presentation parameters such as ‘alpha’, etc. for plotting the grid.

        Returns
        -------
        ax : GeoAxis
            The map with the extracted grid plotted onto the chosen projection (transformed using PlateCarree).
        """
        from .grids import read_netcdf_grid

        raster, lon_coords, lat_coords = read_netcdf_grid(filename, return_grids=True)
        extent = [lon_coords.min(), lon_coords.max(), lat_coords.min(), lat_coords.max()]
        return self.plot_grid(ax, raster, extent=extent, **kwargs)


    def plot_plate_motion_vectors(self, ax, spacingX=10, spacingY=10, normalise=False, **kwargs):
        """Calculates plate motion velocity vector fields at a particular geological time and plots them onto a 
        standard map. 
        
        Generates velocity domain feature collections (MeshNode-type features) from given spacing in the X and Y 
        directions. Domain points are extracted from these features and assigned plate IDs, which are used to obtain
        equivalent stage rotations of identified tectonic plates over a 5 Ma time interval. Each point and its stage 
        rotation are used to calculate plate velocities at a particular geological time. Obtained velocities for each
        domain point are represented in the north-east-down coordinate system.
        
        Vector fields can be optionally normalised to have uniform arrow lengths. 

        Parameters
        ----------
        ax : GeoAxis
            A standard map for lat-lon data built on a matplotlib figure. Can be for a single plot or for multiple subplots.
            Should be set at a particular Cartopy map projection.

        spacingX : int, default=10
            The spacing in the X direction used to make the latitude-longitude velocity domain feature mesh. 

        spacingY : int, default=10
            The spacing in the Y direction used to make the latitude-longitude velocity domain feature mesh. 

        normalise : bool, default=False
            Choose whether to normalise the velocity magnitudes so that vector lengths are all equal. 

        **kwargs
            Keyword arguments that allow control over quiver presentation parameters such as ‘alpha’, etc. for plotting
            the vector field.

        Returns
        -------
        ax : GeoAxis
            The map with the plate velocity vector field plotted onto the chosen projection (transformed using PlateCarree).
        """
        
        lons = np.arange(-180, 180+spacingX, spacingX)
        lats = np.arange(-90, 90+spacingY, spacingY)
        lonq, latq = np.meshgrid(lons, lats)

        # create a feature from all the points
        velocity_domain_features = ptt.velocity_tools.make_GPML_velocity_feature(lonq.ravel(), latq.ravel())

        rotation_model = self.PlateReconstruction_object.rotation_model
        topology_features = self.PlateReconstruction_object.topology_features

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

        return ax.quiver(X, Y, U, V, transform=self.base_projection, **kwargs)
