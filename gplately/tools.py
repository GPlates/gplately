"""A module that offers tools for executing common geological calculations, mathematical conversions and numpy conversions.
"""
from itertools import repeat
from multiprocessing import cpu_count, Pool, set_start_method
from platform import system

import numpy as np
import pygplates
from scipy.optimize import root_scalar

EARTH_RADIUS = pygplates.Earth.mean_radius_in_kms

_DEFAULT_PLATE_THICKNESS = 125.0e3
_DEFAULT_T_MANTLE = 1350.0
_DEFAULT_KAPPA = 8.04e-7
_SEC_PER_MYR = 3.15576e13


def plate_temp(
    age,
    z,
    plate_thickness,
    kappa=_DEFAULT_KAPPA,
    t_mantle=_DEFAULT_T_MANTLE,
    t_surface=0.0,
):
    """Computes the temperature in a cooling plate for a given age and plate thickness at a depth = z. 

    Assumes a mantle temperature of 1350 degrees, and a 0 degree surface temperature. Kappa is 0.804e-6.

    Parameters
    ----------
    age : float
        The geological time (Ma) at which to calculate plate temperature.

    z : float
        The plate depth (m) at which to calculate temperature.

    plate_thickness : float
        The thickness (m) of the plate in consideration.

    Returns
    -------
    float
        Plate temperature at given age and depth.
    """

    age = age * _SEC_PER_MYR
    if age <= 0.0:
        return t_mantle

    sine_arg = np.pi * z / plate_thickness
    exp_arg = -kappa * (np.pi ** 2) * age / (plate_thickness ** 2)
    k = np.arange(1, 20).reshape(-1, 1)
    cumsum = (np.sin(k * sine_arg) * np.exp((k ** 2) * exp_arg) / k).sum(axis=0)

    result = (
        t_surface
        + 2.0 * cumsum * (t_mantle - t_surface) / np.pi
        + (t_mantle - t_surface) * z / plate_thickness
    )
    return result


def plate_isotherm_depth(
    ages,
    temp=_DEFAULT_T_MANTLE,
    plate_thickness=_DEFAULT_PLATE_THICKNESS,
    kappa=_DEFAULT_KAPPA,
    t_mantle=_DEFAULT_T_MANTLE,
    t_surface=0.0,
    n_jobs=1,
    **kwargs
):
    """Compute the depth to the temp - isotherm in a cooling plate model.

    This function uses `scipy.optimize.root_scalar` to efficiently solve
    `gplately.tools.plate_temp` to find the appropriate depths for the
    given ages.

    Parameters
    ----------
    ages : array_like
        Geological ages (Ma) at which to compute depths.
    temp : float, default: 1350.0
        Temperature for which to calculate the isotherm depth.
    plate_thickness : float, default: 125.0e3
        Thickness of the plate, in metres. Defaults to 125 km, as in
        Parsons/Sclater.
    kappa : float, default: 8.04e-7
        Thermal diffusivity coefficient, in metres per second.
    t_mantle : float, default: 1350.0
        Mantle temperature, in degrees Celsius or kelvin.
    t_surface : float, default: 0.0
        Crust surface temperature, in degrees Celsius or kelvin.
    n_jobs : int, default: 1
        Number of processes to use; -1 corresponds to all processes.
        Uses `multiprocessing` module.
    **kwargs
        Any further keyword arguments are passed to
        `scipy.optimize.root_scalar`.

    Returns
    -------
    depths : ndarray
        The calculated isotherm depths. This is a scalar if `ages` is a scalar.
    """
    ages = np.atleast_1d(ages)
    if ages.size == 1:
        age = np.atleast_1d(ages.squeeze())[0]
        return _plate_isotherm_depth_single(
            age,
            temp=temp,
            plate_thickness=plate_thickness,
            kappa=kappa,
            t_mantle=t_mantle,
            t_surface=t_surface,
            **kwargs
        )

    if n_jobs == -1:
        n_jobs = cpu_count()
    elif n_jobs == 0:
        n_jobs = 1

    if n_jobs == 1:
        return np.array(
            [
                _plate_isotherm_depth_single(
                    age=age,
                    temp=temp,
                    plate_thickness=plate_thickness,
                    kappa=kappa,
                    t_mantle=t_mantle,
                    t_surface=t_surface,
                    **kwargs
                )
                for age in ages
            ]
        )
    if system() == "Darwin":
        # Avoid multiprocessing errors on MacOS
        try:
            set_start_method("spawn")
        except RuntimeError:
            pass

    n = len(ages)
    chunk_size = max([int(n / n_jobs / 2), 1])
    fn_args = []
    fn_kwargs = []
    for age in ages:
        fn_args.append(
            (
                age,
                temp,
                plate_thickness,
                kappa,
                t_mantle,
                t_surface,
            )
        )
        fn_kwargs.append(kwargs)

    with Pool(n_jobs) as pool:
        out = _starstarmap(
            pool,
            _plate_isotherm_depth_single,
            fn_args,
            fn_kwargs,
            chunksize=chunk_size,
        )
        pool.close()
        pool.join()
    out = np.array(out).reshape(ages.shape)
    return out


def _plate_isotherm_depth_single(
    age,
    temp=_DEFAULT_T_MANTLE,
    plate_thickness=_DEFAULT_PLATE_THICKNESS,
    kappa=_DEFAULT_KAPPA,
    t_mantle=_DEFAULT_T_MANTLE,
    t_surface=0.0,
    **kwargs
):
    if age == 0.0:
        return 0.0

    f = lambda z: (
        plate_temp(
            age=age,
            z=z,
            plate_thickness=plate_thickness,
            kappa=kappa,
            t_mantle=t_mantle,
            t_surface=t_surface,
        )
        - temp
    )
    solution = root_scalar(f, bracket=[0.0, plate_thickness], **kwargs)
    if not solution.converged:
        err = "Solution did not converge for age {} Ma".format(age)
        err += "\n\t{}".format(solution.flag)
        raise RuntimeError(err)
    return solution.root


def points_to_features(lons, lats, plate_ID=None):
    """Creates point features represented on a unit length sphere in 3D cartesian coordinates from a latitude and 
    longitude list.

    Parameters
    ----------
    lons : list
        The longitudes of needed point features.
        
    lats : list
        The latitudes of needed point features.
        
    plate_ID : int, default=None
        The reconstruction plate ID to assign to the needed point features.

    Returns
    -------
    point_features : list
        Topological point features resolved from the given lists of lat-lon point coordinates.
    """
    try:
        len(lons)
    except TypeError:
        lons = [lons]
    try:
        len(lats)
    except TypeError:
        lats = [lats]
    if len(lons) != len(lats):
        raise ValueError(
            "'lons' and 'lats' must be of equal length"
            + " ({} != {})".format(len(lons), len(lats))
        )
    if len(lons) == 0:
        raise ValueError("'lons' and 'lats' must be provided")

    try:
        len(plate_ID)
    except TypeError:
        plate_ID = [plate_ID] * len(lons)

    # create point features
    point_features = []
    for lon, lat, id in zip(lons, lats, plate_ID):
        point_feature = pygplates.Feature()
        point_feature.set_geometry(pygplates.PointOnSphere(float(lat), float(lon)))
        if id is not None:
            point_feature.set_reconstruction_plate_id(id)
        point_features.append(point_feature)

    if len(point_features) == 1:
        return point_features[0]
    return point_features

def extract_feature_lonlat(features):
    """Extracts the latitudes and longitudes of topological feature points.

    Parameters
    ----------
    features : list
        A list of topological point features to extract lat-lon coordinates from.

    Returns
    -------
    rlon, rlat : lists
        Extracted longitudes and latitudes of feature points. Each index corresponds to different points. Same length
        as the input feature array.
    """
    rlon = np.zeros(len(features))
    rlat = np.zeros(len(features))

    for i, feature in enumerate(features):
        geometry = feature.get_reconstructed_geometry()
        rlat[i], rlon[i] = geometry.to_lat_lon()

    return rlon, rlat


def lonlat2xyz(lon, lat, degrees=True):
    """Convert lon / lat (radians) for spherical triangulation into Cartesian (x,y,z) coordinates on the unit sphere.

    Parameters
    ----------
    lon, lat : lists
        Longitudes and latitudes of feature points in radians.

    Returns
    -------
    xs, ys, zs : lists
        Cartesian coordinates of each feature point in all 3 dimensions.
    """
    if degrees:
        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)
    cosphi = np.cos(lat)
    xs = cosphi * np.cos(lon)
    ys = cosphi * np.sin(lon)
    zs = np.sin(lat)
    return xs, ys, zs

def xyz2lonlat(x, y, z, validate=False, degrees=True):
    """Converts Cartesian (x,y,z) representation of points (on the unit sphere) for spherical triangulation into 
    lon / lat (radians).

    Note: No check is made here that (x,y,z) are unit vectors - it is assumed.

    Parameters
    ----------
    x, y, z : lists
        Cartesian coordinates of each feature point in all 3 dimensions.

    Returns
    -------
    lon, lat : lists
        Longitudes and latitudes of feature points in radians.

    Notes
    -----
    No check is made here that (x,y,z) are unit vectors, unless validate=True is specified
    """
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    if validate:
        mags = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        ones = np.full_like(mags, 1)
        if not np.all(np.equal(mags, ones)):
            raise ValueError("All (x, y, z) must be unit vectors")
    lons = np.arctan2(y, x)
    lats = np.arcsin(z)
    if degrees:
        lons = np.rad2deg(lons)
        lats = np.rad2deg(lats)
    if lons.size == 1:
        lons = np.atleast_1d(np.squeeze(lons))[0]
    if lats.size == 1:
        lats = np.atleast_1d(np.squeeze(lats))[0]
    return lons, lats


def haversine_distance(lon1, lon2, lat1, lat2, degrees=True):
    """Computes the Haversine distance (the shortest distance on the surface of an ideal spherical Earth) between two 
    points given their latitudes and longitudes.

    Sources
    -------
    https://en.wikipedia.org/wiki/Haversine_formula
    https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters

    Parameters
    ----------
    lon1, lon2 : float
        Longitudes of both points

    lat1, lat2 : float
        Latitudes of both points

    Returns
    -------
    d : float
        The Haversine distance in metres.

    Notes
    -----
    Default behaviour assumes values in degrees; for radians specify degrees=False
    """
    if degrees:
        dLat = np.deg2rad(lat2) - np.deg2rad(lat1)
        dLon = np.deg2rad(lon2) - np.deg2rad(lon1)
    else:
        dLat = lat2 - lat1
        dLon = lon2 - lon1
    a = (
        np.sin(dLat / 2) ** 2
        + np.cos(lat1 * np.pi / 180)
        * np.cos(lat2 * np.pi / 180)
        * np.sin(dLon / 2) ** 2
    )
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    d = EARTH_RADIUS * c
    return d * 1000


def geocentric_radius(lat, degrees=True):
    """ Calculates the latitude-dependent radius of an ellipsoid Earth. 

    Parameters
    ----------
    lat : float
        The geodetic latitude at which to calculate the Earth's radius
    degrees : bool, default=True
        Specify whether the given latitude is in degrees.

    Returns
    -------
    earth_radius : float
        The Earth's geocentric radius (in metres) at the given geodetic latitude.
    """
    if degrees:
        rlat = np.radians(lat)
    else:
        rlat = lat

    coslat = np.cos(rlat)
    sinlat = np.sin(rlat)
    r1 = 6384.4e3
    r2 = 6352.8e3
    num = (r1**2*coslat)**2 + (r2**2*sinlat)**2
    den = (r1*coslat)**2 + (r2*sinlat)**2
    earth_radius = np.sqrt(num/den)
    return earth_radius


def _apply_fn_args_kwargs(fn, args, kwargs):
    """Needed for _starstarmap function."""
    return fn(*args, **kwargs)


def _starstarmap(pool, fn, fn_args, fn_kwargs, **starmap_kwargs):
    """Use kwargs with Pool.starmap.

    https://stackoverflow.com/q/45718523"""
    starmap_args = zip(repeat(fn), fn_args, fn_kwargs)
    return pool.starmap(_apply_fn_args_kwargs, starmap_args, **starmap_kwargs)
