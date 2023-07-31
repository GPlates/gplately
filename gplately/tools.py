"""A module that offers tools for executing common geological calculations, mathematical conversions and numpy conversions.
"""
import numpy as np
import pygplates
import pandas as pd
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
    for lon, lat, pid in zip(lons, lats, plate_ID):
        point_feature = pygplates.Feature()
        point_feature.set_geometry(pygplates.PointOnSphere(float(lat), float(lon)))
        if pid is not None:
            point_feature.set_reconstruction_plate_id(pid)
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
    if not isinstance(rotation_model, pygplates.RotationModel):
        rotation_model = pygplates.RotationModel(rotation_model)
    plate_partitioner = pygplates.PlatePartitioner(
        pygplates.FeatureCollection(topology_features), 
        rotation_model,
        reconstruction_time=float(0)
    )
    partitioning_plate = plate_partitioner.partition_point(
        pygplates.PointOnSphere(
            (float(lat_lon_tuple[0]), float(lat_lon_tuple[1]))
        )
    )
    plate_id_at_present_day = partitioning_plate.get_feature().get_reconstruction_plate_id()
    return(plate_id_at_present_day)


def read_rotation_file_pandas(rotation_file_paths):
    """ Written by Nicky M Wright. Extract data from one rotation file, and write 
    it to a pandas dataframe.
    """
    rotation_file = pd.read_csv(
        rotation_file_paths, 
        names = ['reconstruction_plate_id', 'age', 'lat', 'lon', 'angle', 'anchor_plate_id', 'comment'], 
        delim_whitespace=True, 
        comment='!'
    )
    with open(rotation_file_paths, 'r') as f:
        lines = f.readlines()
        output = []

        comment = '!'
        for line in lines:
            head, sep, tail = line.partition(comment)
            tail = tail.strip('\n')
            output.append(tail)
    
    rotation_file['comment'] = output
    
    return rotation_file


def correct_longitudes_for_dateline(lons):
    lons[lons < 0] += 360 # correct for dateline
    return lons

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


def smooth_1D_gaussian(
    input_data, time_window, 
    axis=-1, output=None, mode="reflect", truncate=4.0):
    """Smooth every data element in the 1D `input_data` array over a 
    specified `time_window`, using a one-dimensional, zeroth order 
    Gaussian filter. 
    
    The `time_window`, or Gaussian kernel diameter, is used to 
    calculate a sigma for the Gaussian kernel.
    
    For example, if a `time_window` of 20 Myr is supplied, an array with 21
    elements (10 Myr on each side of a central point, including the central 
    point) is produced. Each element is filled with the Gaussian evaluated 
    at that point. 
    
    This Gaussian is correlated to each data point in the input_array to
    smooth the data.
    
    Parameters
    ----------
    input_data : 1d array
        A one-dimensional array of input data to smooth using a Gaussian filter.
    time_window : float, default = 5 (Myr)
        A float or integer to specify the full width of the Gaussian filter with
        which to smooth each data point in `input_data`. 1 pixel : 1 Myr. 
    axis : int, default = -1
        The axis of `input_data` along which to smooth. Default is -1.
    output : array or dtype, optional
        The array in which to place the output, or the dtype of the returned array. 
        By default an array of the same dtype as input will be created.
    mode : str, default "reflect"
        The way the input_array is extended beyond its bounds to ensure all data
        points can be convolved with the created Gaussian. See 
        [scipy's docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html)
        to see the full list of extension methods. By default, `"reflect"` extends
        `input_array` such that [a b c d] becomes (d c b a | a b c d | d c b a),
        reflecting about the edge of the last pixel.
    truncate : float or int, default = 4.0
        A multiplicative factor that truncates the number of standard deviations, or
        sigmas, that the Gaussian distribution spans in either direction from the
        centre. This impacts the number of pixels that the Gaussian kernel spans. 
        By default, this is set to 4 standard deviations.
    """
    # Sigma, the magnitude of standard deviation in pixel units, is truncated
    # by a specified number of permissible standard deviations; 4 by default
    # half the time window defines the width of the gaussian kernel
    radius = time_window/2
    sigma = int(radius-0.5)/int(truncate)
    
    # Produce the gaussian kernel given a sigma and kernel full-width/diameter 
    # (time_window). The coefficient 1/sqrt(2*pi)*sigma is omitted as it 
    # just scales the kernel distribution and does not impact smoothing weights
    time_window_array = np.arange(-radius, radius+1)
    kernel = np.exp(-0.5 / sigma**2 * time_window_array ** 2)
    
    # Ensure sum of kernel is normalised to unity
    kernel = kernel / kernel.sum()
    
    # Correlate the input data to the Gaussian kernel
    smoothed_data = scipy.ndimage.correlate1d(
        input_data, 
        kernel, 
        axis, output, mode,
        origin=0, # 0 centers the filter over the pixel
        cval=0.0 # Only applicable if mode is 'constant', extends filter by 0s everywhere. 
    )
    return smoothed_data


# From Simon Williams' GPRM
def find_distance_to_nearest_ridge(resolved_topologies,shared_boundary_sections,
                                   point_features,fill_value=5000.):

    all_point_distance_to_ridge = []
    all_point_lats = []
    all_point_lons = []

    for topology in resolved_topologies:
        plate_id = topology.get_resolved_feature().get_reconstruction_plate_id()

        # Section to isolate the mid-ocean ridge segments that bound the current plate
        mid_ocean_ridges_on_plate = []
        for shared_boundary_section in shared_boundary_sections:

            if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('MidOceanRidge'):
                for shared_subsegment in shared_boundary_section.get_shared_sub_segments():
                    sharing_resolved_topologies = shared_subsegment.get_sharing_resolved_topologies()
                    for resolved_polygon in sharing_resolved_topologies:
                        if resolved_polygon.get_feature().get_reconstruction_plate_id() == plate_id:
                            mid_ocean_ridges_on_plate.append(shared_subsegment.get_resolved_geometry())

        point_distance_to_ridge = []
        point_lats = []
        point_lons = []

        for point_feature in point_features:

            for points in point_feature.get_geometries():
                for point in points:

                    if topology.get_resolved_geometry().is_point_in_polygon(point):

                        if len(mid_ocean_ridges_on_plate)>0:

                            min_distance_to_ridge = None

                            for ridge in mid_ocean_ridges_on_plate:
                                distance_to_ridge = pygplates.GeometryOnSphere.distance(point,ridge,min_distance_to_ridge)

                                if distance_to_ridge is not None:
                                    min_distance_to_ridge = distance_to_ridge

                            point_distance_to_ridge.append(min_distance_to_ridge*pygplates.Earth.mean_radius_in_kms)
                            point_lats.append(point.to_lat_lon()[0])
                            point_lons.append(point.to_lat_lon()[1])

                        else:

                            # Originally, give points in plate IDs without MORs a fill value
                            point_distance_to_ridge.append(fill_value)
                            point_lats.append(point.to_lat_lon()[0])
                            point_lons.append(point.to_lat_lon()[1])

                            # Try allocating a NoneType to points instead
                            # point_lats.append(None)
                            # point_lons.append(None)

                            # Try skipping the point (this causes the workflow to use nearest neighbour interpolation to fill these regions)
                            #continue

        all_point_distance_to_ridge.extend(point_distance_to_ridge)
        all_point_lats.extend(point_lats)
        all_point_lons.extend(point_lons)


    return all_point_lons,all_point_lats,all_point_distance_to_ridge


def calculate_spreading_rates(
    time,
    lons,
    lats,
    left_plates,
    right_plates,
    rotation_model,
    delta_time=1.0,
    units=pygplates.VelocityUnits.kms_per_my,
    earth_radius=pygplates.Earth.mean_radius_in_kms,
):

    VELOCITY_UNITS = {
    pygplates.VelocityUnits.kms_per_my,
    pygplates.VelocityUnits.cms_per_yr,
    }
    SPREADING_FEATURES = {
        "gpml:MidOceanRidge",
        "gpml:ContinentalRift",
    }

    if units not in VELOCITY_UNITS:
        raise ValueError("Invalid `units` argument: {}".format(units))

    time = float(time)
    #delta_time = float(time)
    delta_time = float(delta_time)

    if len(set(len(i) for i in (lons, lats, left_plates, right_plates))) > 1:
        err_msg = (
            "Length mismatch: len(lons) = {}".format(len(lons))
            + ", len(lats) = {}".format(len(lats))
            + ", len(left_plates) = {}".format(len(left_plates))
            + ", len(right_plates) = {}".format(len(right_plates))
        )
        raise ValueError(err_msg)

    plate_pairs = {
        pair: {
            "lons": [],
            "lats": [],
            "rotation": _get_rotation(
                pair,
                rotation_model,
                time,
                delta_time,
                use_identity_for_missing_plate_ids=False,
            ),
            "indices": [],
        }
        for pair in set(zip(left_plates, right_plates))
    }
    out = {}
    for index, (lon, lat, left_plate, right_plate) in enumerate(
        zip(
            lons,
            lats,
            left_plates,
            right_plates,
        )
    ):
        pair = (left_plate, right_plate)
        plate_pairs[pair]["lons"].append(lon)
        plate_pairs[pair]["lats"].append(lat)
        plate_pairs[pair]["indices"].append(index)

    for pair in plate_pairs.keys():
        rotation = plate_pairs[pair]["rotation"]
        pair_lons = plate_pairs[pair]["lons"]
        pair_lats = plate_pairs[pair]["lats"]
        pair_points = list(zip(pair_lats, pair_lons))
        indices = plate_pairs[pair]["indices"]

        if rotation is not None:
            velocities = pygplates.calculate_velocities(
                domain_points=pair_points,
                finite_rotation=rotation,
                time_interval_in_my=delta_time,
                velocity_units=units,
                earth_radius_in_kms=earth_radius,
            )
            velocities = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                local_origins=pair_points,
                vectors=velocities,
            )
            velocities = np.abs(np.array([i[0] for i in velocities]))
        else:
            velocities = np.full(len(indices), np.nan)

        for index, velocity in zip(indices, velocities):
            out[index] = velocity

    #print(out)
    #print(sorted(out))
    #print("\n")

    return [out[i] for i in sorted(out)], indices #[indices[i] for i in sorted(indices)]


def _get_rotation(
    plate_pair,
    rotation_model,
    time,
    delta_time=1.0,
    **kwargs
):
    to_time = time - delta_time * 0.5
    from_time = time + delta_time * 0.5

    if to_time <= 0.0:
        to_time = 0.0
        from_time = delta_time

    rotation = rotation_model.get_rotation(
        to_time=to_time,
        from_time=from_time,
        moving_plate_id=plate_pair[0],
        anchor_plate_id=plate_pair[1],
        **kwargs
    )

    return rotation


def _deg2pixels(deg_res, deg_min, deg_max):
    return int(np.floor((deg_max - deg_min) / deg_res)) + 1

def _pixels2deg(spacing_pixel, deg_min, deg_max):
    return (deg_max - deg_min) / np.floor(int(spacing_pixel - 1))
    

