#
#     Copyright (C) 2016-2025 The University of Sydney, Australia
#
#     This program is free software; you can redistribute it and/or modify it under
#     the terms of the GNU General Public License, version 2, as published by
#     the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to Free Software Foundation, Inc.,
#     51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

import numpy as np
import pygplates


def plot_velocities_uv(x, y, u, v, ax):
    """draw the velocity vectors in a map
    Some arrows are long and some are very short.
    To make the plot clearer, we nomalize the velocity magnitude and use color to denote the different speed.

    Parameters
    ----------
    x, y: the domain points coordinates

    u, v: the u and v component of the velocities

    ax: the matplotlib axes object

    Returns
    -------
    A colour bar object
    """
    import cartopy.crs as ccrs

    u = np.array(u)
    v = np.array(v)
    mag = np.sqrt(u * u + v * v)
    mag[mag == 0] = 1  # to avoid 0 divisor
    u = u / mag
    v = v / mag
    cb = ax.quiver(
        x, y, u, v, mag, transform=ccrs.PlateCarree(), cmap="jet", zorder=100
    )
    return cb


def plot_velocities(x, y, velocities, ax):
    """draw the velocity vectors in a map
    Some arrows are long and some are very short.
    To make the plot clearer, we nomalize the velocity magnitude and use color to denote the different speed.

    Parameters
    ----------
    x, y: the domain points coordinates

    velocities: the velocity data

    ax: the matplotlib axes object

    Returns
    -------
    A colour bar object
    """
    x_, y_, u, v = get_x_y_u_v(x, y, velocities)
    return plot_velocities_uv(x_, y_, u, v, ax)


# function to make a velocity mesh nodes at an arbitrary set of points defined in Lat
# Long and Lat are assumed to be 1d arrays.
def make_GPML_velocity_feature(Long, Lat):
    # Add points to a multipoint geometry
    multi_point = pygplates.MultiPointOnSphere(
        [(float(lat), float(lon)) for lat, lon in zip(Lat, Long)]
    )

    # Create a feature containing the multipoint feature, and defined as MeshNode type
    meshnode_feature = pygplates.Feature(
        pygplates.FeatureType.create_from_qualified_string("gpml:MeshNode")
    )
    meshnode_feature.set_geometry(multi_point)
    meshnode_feature.set_name("Velocity Mesh Nodes from pygplates")

    output_feature_collection = pygplates.FeatureCollection(meshnode_feature)

    # NB: at this point, the feature could be written to a file using
    # output_feature_collection.write('myfilename.gpmlz')

    # for use within the notebook, the velocity domain feature is returned from the function
    return output_feature_collection


def get_plate_velocities(
    velocity_domain_features,
    topology_features,
    rotation_model,
    time,
    delta_time,
    rep="vector_comp",
):
    """function to get velocites via pygplates"""

    # Extract points from velocity domain features.
    points = []
    for velocity_domain_feature in velocity_domain_features:
        # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
        # Iterate over them all.
        for velocity_domain_geometry in velocity_domain_feature.get_geometries():
            for velocity_domain_point in velocity_domain_geometry.get_points():
                points.append(velocity_domain_point)

    # Topological snapshot containing resolved topologies.
    topological_snapshot = pygplates.TopologicalSnapshot(
        topology_features,
        rotation_model,
        time,
    )

    # Calculate point velocities using topological snapshot.
    point_velocities = topological_snapshot.get_point_velocities(
        points,
        velocity_delta_time=delta_time,
        velocity_units=pygplates.VelocityUnits.kms_per_my,
        earth_radius_in_kms=pygplates.Earth.mean_radius_in_kms,
    )

    # Replace any missing velocities with zero velocity.
    #
    # If a point does not intersect a topological plate (or network) then its velocity is None.
    for point_index in range(len(points)):
        if point_velocities[point_index] is None:
            point_velocities[point_index] = pygplates.Vector3D.zero

    if rep == "mag_azim":
        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
        return pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
            points, point_velocities
        )

    elif rep == "vector_comp":
        # Convert global 3D velocity vectors to local (North, East, Down) vectors (one per point).
        return pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
            points, point_velocities
        )


def get_velocities(
    time, rotation_model, topology_filenames, delta_time=5.0, Xnodes=[], Ynodes=[]
):
    """get velocity data

    Parameters
    ----------
    time: number
        the reconstruction time

    rotation_model
        the rotation model

    topology_filenames
        the topology files

    delta_time: number
        the time increment

    Xnodes, Ynodes: the coordinates of domain points

    Returns
    -------
    the coordinates of domain points and velocity data

    """
    if len(Xnodes) == 0 or len(Ynodes) == 0:
        Xnodes = np.arange(-180, 180, 10)
        Ynodes = np.arange(-90, 90, 10)
    Xg, Yg = np.meshgrid(Xnodes, Ynodes)
    Xg = Xg.flatten()
    Yg = Yg.flatten()
    velocity_domain_features = make_GPML_velocity_feature(Xg, Yg)

    # Load the topological plate polygon features.
    topology_features = []
    for fname in topology_filenames:
        for f in pygplates.FeatureCollection(fname):
            topology_features.append(f)

    # Call the function we created above to get the velocities
    return (
        Xnodes,
        Ynodes,
        get_plate_velocities(
            velocity_domain_features,
            topology_features,
            rotation_model,
            time,
            delta_time,
            "vector_comp",
        ),
    )


def get_x_y_u_v(Xnodes, Ynodes, all_velocities):
    """get velocity u,v components from velocity data

    Parameters
    ----------
    all_velocities: the velocity data
    Xnodes, Ynodes: the coordinates of domain points

    Returns
    -------
    the coordinates of domain points and velocity u, v components
    """
    uu = []
    vv = []
    for vel in all_velocities:
        if not hasattr(vel, "get_y"):
            uu.append(vel[1])
            vv.append(vel[0])
        else:
            uu.append(vel.get_y())
            vv.append(vel.get_x())
    u = np.asarray([uu]).reshape((Ynodes.shape[0], Xnodes.shape[0]))
    v = np.asarray([vv]).reshape((Ynodes.shape[0], Xnodes.shape[0]))

    return Xnodes, Ynodes, u, v


def get_velocity_x_y_u_v(
    time, rotation_model, topology_filenames, delta_time=5.0, Xnodes=[], Ynodes=[]
):
    """get velocity data in x ,y, u, v form

    Parameters
    ----------
    time: number
        the reconstruction time

    rotation_model
        the rotation model

    topology_filenames
        the topology files

    delta_time: number
        the time increment

    Xnodes, Ynodes: the coordinates of domain points

    Returns
    -------
    the coordinates of domain points and velocity u, v components
    """
    Xnodes, Ynodes, all_velocities = get_velocities(
        time,
        rotation_model,
        topology_filenames,
        delta_time=5.0,
        Xnodes=Xnodes,
        Ynodes=Ynodes,
    )
    return get_x_y_u_v(Xnodes, Ynodes, all_velocities)
