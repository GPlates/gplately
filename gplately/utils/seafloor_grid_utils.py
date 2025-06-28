#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
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
Auxiliary functions for SeafloorGrid
"""

import numpy as np
import pygplates

from .. import ptt


def create_icosahedral_mesh(refinement_levels):
    """Define a global point mesh with Stripy's
    [icosahedral triangulated mesh](https://github.com/underworldcode/stripy/blob/294354c00dd72e085a018e69c345d9353c6fafef/stripy/spherical_meshes.py#L27)
    and turn all mesh domains into pyGPlates MultiPointOnSphere types.

    This global mesh will be masked with a set of continental or COB terrane
    polygons to define the ocean basin at a given reconstruction time.
    The `refinement_levels` integer is proportional to the resolution of the
    mesh and the ocean/continent boundary.

    Parameters
    ----------
    refinement_levels : int
        Refine the number of points in the triangulation. The larger the
        refinement level, the sharper the ocean basin resolution.

    Returns
    -------
    multi_point : instance of <pygplates.MultiPointOnSphere>
        The longitues and latitudes that make up the icosahedral ocean mesh
        collated into a MultiPointOnSphere object.
    icosahedral_global_mesh : instance of <stripy.spherical_meshes.icosahedral_mesh>
        The original global icosahedral triangulated mesh.
    """
    import stripy

    # Create the ocean basin mesh using Stripy's icosahedral spherical mesh
    icosahedral_global_mesh = stripy.spherical_meshes.icosahedral_mesh(
        refinement_levels, include_face_points=False, trisection=False, tree=False
    )
    # Get lons and lats of mesh, and turn them into a MultiPointOnSphere
    lats_arr = np.rad2deg(icosahedral_global_mesh.lats)
    lons_arr = np.rad2deg(icosahedral_global_mesh.lons)
    multi_point = pygplates.MultiPointOnSphere(zip(lats_arr, lons_arr))

    return multi_point, icosahedral_global_mesh


def ensure_polygon_geometry(reconstructed_polygons, rotation_model, time):
    """Ensure COB terrane/continental polygon geometries are polygons
    with reconstruction plate IDs and valid times.

    Notes
    -----
    This step must be done so that the initial set of ocean basin points
    (the Stripy icosahedral mesh) can be partitioned into plates using
    each reconstruction plate ID for the given plate `model`.

    This allows for an oceanic point-in-continental
    polygon query for every identified plate ID. See documentation for
    `point_in_polygon_routine` for more details.

    `ensure_polygon_geometry` works as follows:
    COB terrane/continental polygons are assumed to have been reconstructed
    already in `reconstructed_polygons` (a list of
    type <pygplates.ReconstructedFeatureGeometry>). The list contents are
    turned into a <pygplates.FeatureCollection> to be ascribed a
    `PolygonOnSphere` geometry, a reconstruction plate ID, and a valid time.
    Once finished, this feature collection is turned back into a list of
    instance <pygplates.ReconstructedFeatureGeometry> and returned.

    This revert must be completed for compatibility with the subsequent
    point-in-polygon routine.

    Parameters
    ----------
    reconstructed_polygons : list of instance <pygplates.ReconstructedFeatureGeometry>
        If used in `SeafloorGrid`, these are automatically obtained from the
        `PlotTopologies.continents` attribute (the reconstructed continental
        polygons at the current reconstruction time).

    rotation_model : instance of <pygplates.RotationModel>
        A parameter for turning the <pygplates.FeatureCollection> back into a
        list of instance <pygplates.ReconstructedFeatureGeometry> for
        compatibility with the point-in-polygon routine.

    """
    continent_FeatCol = []
    # self._PlotTopologies_object.continents
    for n in reconstructed_polygons:
        continent_FeatCol.append(n.get_feature())

    polygon_feats = pygplates.FeatureCollection(continent_FeatCol)

    # From GPRM's force_polygon_geometries(); set feature attributes
    # like valid times and plate IDs to each masking polygon
    polygons = []
    for feature in polygon_feats:
        for geom in feature.get_all_geometries():
            polygon = pygplates.Feature(feature.get_feature_type())
            polygon.set_geometry(pygplates.PolygonOnSphere(geom))
            polygon.set_reconstruction_plate_id(feature.get_reconstruction_plate_id())
            # Avoid features in COBTerranes with invalid time
            if feature.get_valid_time()[0] >= feature.get_valid_time()[1]:
                polygon.set_valid_time(
                    feature.get_valid_time()[0], feature.get_valid_time()[1]
                )
                polygons.append(polygon)
    cobter_polygon_features = pygplates.FeatureCollection(polygons)

    # Turn the feature collection back into ReconstructedFeatureGeometry
    # objects otherwise it will not work with PIP
    reconstructed_cobter_polygons = []
    pygplates.reconstruct(
        cobter_polygon_features, rotation_model, reconstructed_cobter_polygons, time
    )
    return reconstructed_cobter_polygons


def point_in_polygon_routine(multi_point, COB_polygons):
    """Perform Plate Tectonic Tools' point in polygon routine to partition
    points in a `multi_point` MultiPointOnSphere feature based on whether
    they are inside or outside the polygons in `COB_polygons`.

    Notes
    -----
    Assuming the `COB_polygons` have passed through `ensure_polygon_geometry`,
    each polygon should have a plate ID assigned to it.

    This PIP routine serves two purposes for `SeafloorGrid`:

    1) It identifies continental regions in the icosahedral global mesh
    MultiPointOnSphere feature and 'erases' in-continent oceanic points
    for the construction of a continental mask at each timestep;

    2) It identifies oceanic points in the icosahedral global mesh.
    These points will be passed to a function that calculates each point's
    proximity to its nearest MOR segment (if any) within the polygonal domain
    of its allocated plate ID. Each distance is divided by half the
    `initial_ocean_mean_spreading_rate` (an attribute of `SeafloorGrids`) to
    determine a simplified seafloor age for each point.

    Number 2) only happens once at the start of the gridding process to
    momentarily fill the gridding region with initial ocean points that have
    set ages (albeit not from a plate model file). After multiple time steps
    of reconstruction, the ocean basin will be filled with new points (with
    plate-model prescribed ages) that emerge from ridge topologies.


    Returns
    -------
    pygplates.MultiPointOnSphere(points_in_arr) : instance <pygplates.MultiPointOnSphere>
        Point features that are within COB terrane polygons.
    pygplates.MultiPointOnSphere(points_out_arr) : instance <pygplates.MultiPointOnSphere>
        Point features that are outside COB terrane polygons.
    zvals : list
        A binary list. If an entry is == 0, its corresponing point in the
        MultiPointOnSphere object is on the ocean. If == 1, the point is
        in the COB terrane polygon.
    """
    # Convert MultiPointOnSphere to array of PointOnSphere
    multi_point = np.array(multi_point.get_points(), dtype="object")

    # Collect reconstructed geometries of continental polygons
    polygons = np.empty(len(COB_polygons), dtype="object")
    for ind, i in enumerate(COB_polygons):
        if isinstance(i, pygplates.ReconstructedFeatureGeometry):
            geom = i.get_reconstructed_geometry()
        elif isinstance(i, pygplates.GeometryOnSphere):
            geom = i
        else:  # e.g. ndarray of coordinates
            geom = pygplates.PolygonOnSphere(i)
        polygons[ind] = geom
    proxies = np.ones(polygons.size)

    pip_result = ptt.utils.points_in_polygons.find_polygons(
        multi_point, polygons, proxies, all_polygons=False
    )  # 1 for points in polygons, None for points outside
    zvals = np.array(
        pip_result,
        dtype="float",
    ).ravel()
    zvals[np.isnan(zvals)] = 0.0
    zvals = zvals.astype("int")
    points_in_arr = multi_point[zvals == 1]
    points_out_arr = multi_point[zvals != 1]

    return (
        pygplates.MultiPointOnSphere(points_in_arr),
        pygplates.MultiPointOnSphere(points_out_arr),
        zvals,
    )
