"""A module that offers tools for executing common geological calculations, mathematical conversions and numpy conversions.
"""
import numpy as np
import pygplates
import scipy

EARTH_RADIUS = pygplates.Earth.mean_radius_in_kms

_DEFAULT_PLATE_THICKNESS = 125.0e3
_DEFAULT_T_MANTLE = 1350.0
_DEFAULT_KAPPA = 8.04e-7
_SEC_PER_MYR = 3.15576e13


def plate_temp(
    age,
    z,
    plate_thickness=_DEFAULT_PLATE_THICKNESS,
    kappa=_DEFAULT_KAPPA,
    t_mantle=_DEFAULT_T_MANTLE,
    t_surface=0.0,
):
    """Compute the temperature in a cooling plate for a given age, plate
    thickness, and depth.

    By default, assumes a mantle temperature of 1350 degrees, and a 0 degree
    surface temperature. Kappa defaults to 0.804e-6.

    Parameters
    ----------
    age : float
        The geological time (Ma) at which to calculate plate temperature.
    z : array_like
        The plate depth(s) (m) at which to calculate temperature.
    plate_thickness : float, default: 125.0e3
        The thickness (m) of the plate in consideration.
    kappa : float, default: 0.804e-6
    t_mantle : float, default: 1350.0
        Mantle temperature at the Moho, in degrees Celsius.
    t_surface : float, default: 0.0
        Temperature at the surface, in degrees Celsius.

    Returns
    -------
    ndarray
        The plate temperature at the given age and depth. This is a scalar
        if `z` is a scalar.
    """
    aged = age * _SEC_PER_MYR

    z = np.atleast_1d(z)

    sine_arg = np.pi * z / plate_thickness
    exp_arg = -kappa * (np.pi ** 2) * aged / (plate_thickness ** 2)
    k = np.arange(1, 20).reshape(-1, 1)
    cumsum = (np.sin(k * sine_arg) * np.exp((k ** 2) * exp_arg) / k).sum(axis=0)

    result = (
        t_surface
        + 2.0 * cumsum * (t_mantle - t_surface) / np.pi
        + (t_mantle - t_surface) * z / plate_thickness
    )
    result = result.reshape(z.shape)
    if result.size == 1:  # input was a scalar
        return result.flatten()[0]
    return result


def plate_isotherm_depth(
    age,
    temp=1150.,
    plate_thickness=_DEFAULT_PLATE_THICKNESS,
    maxiter=50,
    tol=0.001,
    require_convergence=False,
    **kwargs
):
    """Computes the depth to the temp - isotherm in a cooling plate mode. Solution by iteration. 

    By default the plate thickness is 125 km as in Parsons/Sclater.

    Parameters
    ----------
    age : array_like
        Geological ages (Ma) at which to compute depths.
    temp : float, default: 1150.0
        The temperature of a temp-isotherm for which to calculate depth,
        in degrees Celsius.
    plate_thickness : float, default: 125.0e3
        Thickness of the plate, in metres.
    maxiter : int, default: 50
        Maximum number of iterations.
    tol: float, default: 0.001
        Tolerance for convergence of `plate_temp`.
    require_convergence: bool, default: False
        If True, raise a `RuntimeError` if convergence is not reached within
        `maxiter` iterations; if False, a `RuntimeWarning` is issued instead.
    **kwargs
        Any further keyword arguments are passed on to `plate_temp`.

    Returns
    -------
    z : ndarray
        Array of depths to the chosen temperature isotherm at the given ages.
        This is a scalar if `age` is a scalar.

    Raises
    ------
    RuntimeError
        If `require_convergence` is True and convergence is not reached within
        `maxiter` iterations.
    """

    # If `rtol` is passed, use that as `tol`
    # `rtol` actually means 'relative tolerance' in scipy,
    # while `xtol` is used for absolute tolerance
    tol = kwargs.pop("rtol", tol)

    # For backwards compatibility, also allow `n` as a keyword argument
    maxiter = kwargs.pop("n", maxiter)

    maxiter = int(maxiter)
    if maxiter <= 0:
        raise ValueError(
            "`maxiter` must be greater than zero ({})".format(maxiter)
        )
    age = np.atleast_1d(age)

    non_zero_ages = age[age > 0.0]

    zi = np.zeros_like(non_zero_ages, dtype=float)  # starting depth is 0

    z_too_small = np.zeros_like(non_zero_ages, dtype=float)
    z_too_big = np.full_like(non_zero_ages, plate_thickness, dtype=float)

    t_diff = np.full_like(non_zero_ages, np.nan)
    for _ in range(maxiter):
        zi = 0.5 * (z_too_small + z_too_big)
        ti = plate_temp(non_zero_ages, zi, plate_thickness, **kwargs)
        t_diff = temp - ti
        z_too_big[t_diff < -tol] = zi[t_diff < -tol]
        z_too_small[t_diff > tol] = zi[t_diff > tol]

        if (np.abs(t_diff) < tol).all():
            break

    # convergence warning
    if (np.abs(t_diff) > tol).any():
        failed_ages = non_zero_ages[np.abs(t_diff) > tol]
        message = (
            "Solution did not converge below tol={}".format(tol)
            + " within maxiter={} iterations".format(maxiter)
            + " for the following ages: {}".format(failed_ages)
        )
        if require_convergence:
            raise RuntimeError(message)
        import warnings

        warnings.warn(message, category=RuntimeWarning)

    # protect against negative ages
    out = np.zeros_like(age)
    out[age > 0.0] = zi
    out = out.reshape(age.shape)
    if out.size == 1:  # input age was a scalar
        return out.flatten()[0]
    return out


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


def plate_partitioner_for_point(lat_lon_tuple, topology_features, rotation_model):
    """ Determine the present-day plate ID of a (lat, lon) coordinate pair if 
    it is not specified.
    """
    plate_partitioner = pygplates.PlatePartitioner(
        pygplates.FeatureCollection(topology_features), 
        pygplates.RotationModel(rotation_model), 
        reconstruction_time=float(0)
    )
    partitioning_plate = plate_partitioner.partition_point(
        pygplates.PointOnSphere(
            (float(lat_lon_tuple[0]), float(lat_lon_tuple[1]))
        )
    )
    plate_id_at_present_day = partitioning_plate.get_feature().get_reconstruction_plate_id()
    return(plate_id_at_present_day)


# Auxiliary functions for the Muller et al. 2022 paper "Evolution of Earth’s tectonic carbon conveyor belt"
def surface_area_oblate_spheroid(r1, r2):
    e = np.sqrt(1.0 - r2**2/r1**2)
    return 2.0*np.pi*r1**2*(1.0 + (1.0-e**2)/e * np.arctanh(e))


def lat_area_function(latitude_one, latitude_two, longitude_resolution):
    '''
    Calculates the point area of an evenly gridded lat/lon mesh
    Longitude resolution is lon2 - lon1
    '''
    dlat = np.sin(np.radians(latitude_two)) - np.sin(np.radians(latitude_one))
    lat_area = 2 * np.pi * 6371.009e3**2 * np.abs(dlat)/longitude_resolution
    return lat_area


def smooth_1D(array, sigma=3.0, axis=0):
    """ Gaussian filter with standard deviation """
    return scipy.ndimage.gaussian_filter1d(array, sigma, axis=axis)


def My2s(Ma):
    return Ma*3.1536e13


def update_progress(progress):
    from IPython.display import clear_output
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1
    
    bar_length = 20
    block = int(round(bar_length * progress))

    clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    print(text)


def read_csv(filename, readcols):
    """ read csv and reorder from 0 - 250 Ma """
    Ma = np.loadtxt(filename, delimiter=',', usecols=(0,), skiprows=1, dtype=int, unpack=True)
    data = np.loadtxt(filename, delimiter=',', usecols=readcols, skiprows=1, unpack=False)
    return data[Ma]


