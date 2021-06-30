import numpy as np
import pygplates

def plate_temp(age, z, PLATE_THICKNESS) :
    "Computes the temperature in a cooling plate for age = t\
    and at a depth = z."

    KAPPA = 0.804e-6
    T_MANTLE = 1350.0
    T_SURFACE = 0.0

    sine_arg = np.pi * z / PLATE_THICKNESS
    exp_arg = -KAPPA * np.pi * np.pi * age / (PLATE_THICKNESS * PLATE_THICKNESS)
    k = np.ones_like(age)*np.arange(1, 20).reshape(-1,1)
    cumsum = ( np.sin(k * sine_arg) * np.exp(k*k*exp_arg)/k ).sum(axis=0)

    return T_SURFACE + 2.0 * cumsum * (T_MANTLE - T_SURFACE)/np.pi + (T_MANTLE - T_SURFACE) * z/PLATE_THICKNESS

def plate_isotherm_depth(age, temp=1350.0):
    "Computes the depth to the temp - isotherm in a cooling plate mode.\
    Solution by iteration. By default the plate thickness is 125 km as\
    in Parsons/Sclater."

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

    rlon = np.zeros(len(features))
    rlat = np.zeros(len(features))

    for i, feature in enumerate(features):
        geometry = feature.get_reconstructed_geometry()
        rlat[i], rlon[i] = geometry.to_lat_lon()

    return rlon, rlat


def lonlat2xyz(lon, lat):
    """
    Convert lon / lat (radians) for the spherical triangulation into x,y,z
    on the unit sphere
    """
    cosphi = np.cos(lat)
    xs = cosphi*np.cos(lon)
    ys = cosphi*np.sin(lon)
    zs = np.sin(lat)
    return xs, ys, zs

def xyz2lonlat(x,y,z):
    """
    Convert x,y,z representation of points *on the unit sphere* of the
    spherical triangulation to lon / lat (radians).

    Notes:
        no check is made here that (x,y,z) are unit vectors
    """
    lons = np.arctan2(ys, xs)
    lats = np.arcsin(zs)
    return lons, lats


def haversine_distance(lon1, lon2, lat1, lat2):
    """
    from  https://en.wikipedia.org/wiki/Haversine_formula
    https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters
    """
    R = 6378.137 # radius of earth in km
    dLat = lat2*np.pi/180 - lat1*np.pi/180
    dLon = lon2*np.pi/180 - lon1*np.pi/180
    a = np.sin(dLat/2)**2 + np.cos(lat1*np.pi/180)*np.cos(lat2*np.pi/180) * np.sin(dLon/2)**2
    c = 2.0*np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    d = R*c
    return d*1000
