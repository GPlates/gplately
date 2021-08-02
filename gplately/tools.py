import numpy as np
import pygplates

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
    """
    Computes the temperature in a cooling plate for age = t
    and at a depth = z.

    Notes:
        - age is given in Myr
        - z and plate_thickness are given in metres
        - t_mantle and t_surface are given in degrees Celsius or Kelvin
    """

    age *= _SEC_PER_MYR

    sine_arg = np.pi * z / plate_thickness
    exp_arg = -kappa * (np.pi ** 2) * age / (plate_thickness ** 2)
    k = np.arange(1, 20).reshape(-1, 1)
    cumsum = (np.sin(k * sine_arg) * np.exp((k ** 2) * exp_arg) / k).sum(axis=0)

    result = (
        t_surface
        + 2.0 * cumsum * (t_mantle - t_surface) / np.pi
        + (t_mantle - t_surface) * z / plate_thickness
    )
    if result.size == 1:
        return result[0]
    return result


def plate_isotherm_depth(
    age,
    temp=_DEFAULT_T_MANTLE,
    plate_thickness=_DEFAULT_PLATE_THICKNESS,
    n=20,
    rtol=0.001,
):
    """
    Computes the depth to the temp - isotherm in a cooling plate mode.
    Solution by iteration. By default the plate thickness is 125 km as
    in Parsons/Sclater.
    """

    n = int(n)
    if n <= 0:
        raise ValueError("n must be greater than zero (n = {})".format(n))
    age = np.atleast_1d(age)

    zi = np.atleast_1d(np.zeros_like(age, dtype=float))  # starting depth is 0

    z_too_small = np.atleast_1d(np.zeros_like(age, dtype=float))
    z_too_big = np.atleast_1d(np.full_like(age, plate_thickness, dtype=float))

    func = np.vectorize(plate_temp)
    for _ in range(n):
        zi = 0.5 * (z_too_small + z_too_big)
        ti = func(age, zi, plate_thickness)
        t_diff = temp - ti
        z_too_big[t_diff < -rtol] = zi[t_diff < -rtol]
        z_too_small[t_diff > rtol] = zi[t_diff > rtol]

        if (np.abs(t_diff) < rtol).all():
            break

    # convergence warning
    if np.abs(t_diff) > rtol:
        import warnings
        warnings.warn("Iterations did not converge below rtol={}".format(rtol))

    # protect against negative ages
    zi[age <= 0] = 0
    zi = np.atleast_1d(np.squeeze(zi))
    if zi.size == 1:
        return zi[0]
    return zi


def points_to_features(lons, lats, plate_ID=None):
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
        point_feature.set_geometry(pygplates.PointOnSphere(lat, lon))
        if id is not None:
            point_feature.set_reconstruction_plate_id(id)
        point_features.append(point_feature)

    if len(point_features) == 1:
        return point_features[0]
    return point_features


def extract_feature_lonlat(features):

    rlon = np.zeros(len(features))
    rlat = np.zeros(len(features))

    for i, feature in enumerate(features):
        geometry = feature.get_reconstructed_geometry()
        rlat[i], rlon[i] = geometry.to_lat_lon()

    return rlon, rlat


def lonlat2xyz(lon, lat, degrees=True):
    """
    Convert lon / lat (radians) for the spherical triangulation into x, y, z
    on the unit sphere
    """
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    if degrees:
        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)
    cosphi = np.cos(lat)
    xs = cosphi * np.cos(lon)
    ys = cosphi * np.sin(lon)
    zs = np.sin(lat)
    if xs.size == 1:
        xs = np.atleast_1d(np.squeeze(xs))[0]
    if ys.size == 1:
        ys = np.atleast_1d(np.squeeze(ys))[0]
    if zs.size == 1:
        zs = np.atleast_1d(np.squeeze(zs))[0]
    return xs, ys, zs


def xyz2lonlat(x, y, z, validate=False, degrees=True):
    """
    Convert x,y,z representation of points *on the unit sphere* of the
    spherical triangulation to lon / lat (in radians, by default).

    Notes:
        - no check is made here that (x,y,z) are unit vectors, unless
            validate=True is specified
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
    """
    From  https://en.wikipedia.org/wiki/Haversine_formula and
    https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters

    Notes:
        - Default behaviour assumes values in degrees;
            for radians specify degrees=False
        - Value returned is in metres
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
