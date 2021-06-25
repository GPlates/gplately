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