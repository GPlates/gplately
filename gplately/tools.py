"""A module that offers tools for executing common geological calculations, mathematical conversions and numpy conversions.
"""
import numpy as np
import pygplates

def plate_temp(age, z, PLATE_THICKNESS) :
    """Computes the temperature in a cooling plate for a given age and plate thickness at a depth = z. 

    Assumes a mantle temperature of 1380 degrees, and a 0 degree surface temperature. Kappa is 0.804e-6. 

    Parameters
    ----------
    age : float
        The geological time (Ma) at which to calculate plate temperature.

    z : float
        The plate depth (m) at which to calculate temperature.

    PLATE_THICKNESS : float
        The thickness (m) of the plate in consideration.

    Returns
    -------
    list
        A list enclosing ONE floating-point number equal to plate temperature (e.g. [1367.33962383]).
    """

    KAPPA = 0.804e-6
    T_MANTLE = 1350.0
    T_SURFACE = 0.0

    sine_arg = np.pi * z / PLATE_THICKNESS
    exp_arg = -KAPPA * np.pi * np.pi * age / (PLATE_THICKNESS * PLATE_THICKNESS)
    k = np.ones_like(age)*np.arange(1, 20).reshape(-1,1)
    cumsum = ( np.sin(k * sine_arg) * np.exp(k*k*exp_arg)/k ).sum(axis=0)

    return T_SURFACE + 2.0 * cumsum * (T_MANTLE - T_SURFACE)/np.pi + (T_MANTLE - T_SURFACE) * z/PLATE_THICKNESS

def plate_isotherm_depth(age, temp=1350.0):
    """Computes the depth to the temp - isotherm in a cooling plate mode. Solution by iteration. 

    By default the plate thickness is 125 km as in Parsons/Sclater.

    Parameters
    ----------
    age : ndarray
        An array of geological ages (Ma) at which to compute depths. 

    temp : float, default=1350.0
        The temperature of a temp-isotherm to calculate the depth to. Defaults to 1350 degrees.

    Returns
    -------
    zi : ndarray
        An array of depths to the chosen temperature isotherm. Each entry corresponds to each unique ‘age’ given. 
    """
    PLATE_THICKNESS = 125e3
    
    z = 0.0 # starting depth is 0
    rtol = 0.001 # error tolerance
    
    z_too_small = np.atleast_1d(np.zeros_like(age, dtype=np.float))
    z_too_big = np.atleast_1d(np.full_like(age, PLATE_THICKNESS, dtype=np.float))
    
    for i in range(20):
        zi = 0.5 * (z_too_small + z_too_big)
        ti = plate_temp (age, zi, PLATE_THICKNESS)
        t_diff = temp - ti
        z_too_big[t_diff < -rtol] = zi[t_diff < -rtol]
        z_too_small[t_diff > rtol] = zi[t_diff > rtol]
        
        if (np.abs(t_diff) < rtol).all():
            break
            
    # protect against negative ages
    zi[age <= 0] = 0
    return np.squeeze(zi)

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
    # create point features
    point_features = []
    for lon, lat in zip(lons, lats):
        point_feature = pygplates.Feature()
        point_feature.set_geometry(pygplates.PointOnSphere(lat, lon))
        if plate_ID:
            point_feature.set_reconstruction_plate_id(plate_ID)
        point_features.append(point_feature)

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


def lonlat2xyz(lon, lat):
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
    cosphi = np.cos(lat)
    xs = cosphi*np.cos(lon)
    ys = cosphi*np.sin(lon)
    zs = np.sin(lat)
    return xs, ys, zs

def xyz2lonlat(x,y,z):
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
    """
    lons = np.arctan2(ys, xs)
    lats = np.arcsin(zs)
    return lons, lats


def haversine_distance(lon1, lon2, lat1, lat2):
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
    """
    R = 6378.137 # radius of earth in km
    dLat = lat2*np.pi/180 - lat1*np.pi/180
    dLon = lon2*np.pi/180 - lon1*np.pi/180
    a = np.sin(dLat/2)**2 + np.cos(lat1*np.pi/180)*np.cos(lat2*np.pi/180) * np.sin(dLon/2)**2
    c = 2.0*np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    d = R*c
    return d*1000
