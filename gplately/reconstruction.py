"""The “reconstruction” module offers simple shortcuts to pyGplates and Plate Tectonic Tools functionalities for reconstructing features, working with point data, and calculating plate velocities at specific geological times. 

Classes
-------
PlateReconstruction
Points
"""
import pygplates
import numpy as np
import ptt
import warnings

from . import tools as _tools



class PlateReconstruction(object):
    """The PlateReconstruction class contains methods to reconstruct topology features at a specific geological time using
    a given rotation model, a feature/feature collection and a set of static polygons. Velocity data for specific times can
    also be extracted from these reconstructed features. 

    Attributes
    ----------
    rotation_model : str or list
        A pygplates rotation model
    topology_features : list
        A feature collection list containing pygplates features (aggregated using pygplates.FeatureCollection())
    static_polygons : str 
        Path to static polygon file
        
    Methods
    -------
    __init__(self, rotation_model=None, topology_features=None, static_polygons=None)
        Constructs all necessary attributes for the plate reconstruction object.

    tesselate_subduction_zones(self, time, tessellation_threshold_radians=0.001, anchor_plate_id=0) 
        Samples points along subduction zone trenches and obtains both convergence and absolute velocities at a
        particular geological time.

    tesselate_mid_ocean_ridges(self, time, tessellation_threshold_radians=0.001, anchor_plate_id=0)
        Samples points along resolved spreading features (e.g. mid-ocean ridges) and calculates spreading rate and length
        of ridge segments at a particular geological time.

    reconstruct(self, feature, to_time, from_time=0, anchor_plate_id=0, **kwargs)
        Reconstructs regular geological features, motion paths or flowlines to a specific geological time.

    get_point_velocities(self, lons, lats, time, delta_time=1.0)
        Generates a velocity domain feature collection, resolves them into points, and calculates the north and east 
        components of the velocity vector for each point in the domain at a particular geological time. 
    """
    
    def __init__(self, rotation_model=None, topology_features=None, static_polygons=None):
        """Constructs all necessary attributes for the plate reconstruction object.

        Parameters
        ----------
        rotation_model : str, or :class:`FeatureCollection`, or :class:`Feature`, or sequence of :class:`Feature`, 
        or sequence of any combination of those four types, default=None 
            Can be provided as a rotation filename, or rotation feature collection, or rotation feature, or sequence of 
            rotation features, or a sequence (eg, a list or tuple) of any combination of those four types.

        topology_features : str, or a sequence (eg, ``list`` or ``tuple``) of :class:`Feature`, or a single :class:`Feature`,
        default=None
            Can be provided as an optional topology-feature filename, or sequence of features, or a single feature. 
            Note: since a :class:`FeatureCollection` is an iterable sequence of features it can be used in the 
            *features* argument.

        static_polygons : :class:`FeatureCollection`, or str, or :class:`Feature`, or sequence of :class:`Feature`, 
        or a sequence of any combination of those four types, default=None
            Can be provided as a static polygon feature collection, or optional filename, or a single feature, or a sequence of
            features.


        Raises
        ------
        OpenFileForReadingError 
            if any file is not readable (when filenames specified)

        FileFormatNotSupportedError 
            if any file format (identified by the filename extensions) does not support reading (when filenames specified)
        """

        rotation_model = pygplates.RotationModel(rotation_model)

        default_topology_features = pygplates.FeatureCollection()
        for topology in topology_features:
            default_topology_features.add( pygplates.FeatureCollection(topology) )

        self.rotation_model = rotation_model
        self.topology_features = default_topology_features
        self.static_polygons = static_polygons


    def tesselate_subduction_zones(self, time, tessellation_threshold_radians=0.001, ignore_warnings=False, **kwargs):
        """Samples points along subduction zone trenches and obtains both convergence and absolute velocities at a particular
        geological time.
        
        Resolves topologies at 'time', tessellates all resolved subducting features to within 'tessellation_threshold_radians'
        radians and obtains the following information for each sampled point along a trench:
    
        0 - longitude of sampled point
        1 - latitude of sampled point
        2 - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
        3 - subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
        4 - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
        5 - trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
        6 - length of arc segment (in degrees) that current point is on
        7 - trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
        8 - subducting plate ID
        9 - trench plate ID

        The obliquity angles are in the range (-180 180). The range (0, 180) goes clockwise (when viewed from above the Earth)
        from the trench normal direction to the velocity vector. The range (0, -180) goes counter-clockwise.
        You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.
        The trench normal is perpendicular to the trench and pointing toward the overriding plate.
    
        Note that the convergence velocity magnitude is negative if the plates are diverging (if convergence obliquity angle 
        is greater than 90 or less than -90). And note that the absolute velocity magnitude is negative if the trench 
        (subduction zone) is moving towards the overriding plate (if absolute obliquity angle is less than 90 or greater 
        than -90) - note that this ignores the kinematics of the subducting plate.
        
        The delta time interval used for velocity calculations is, by default, assumed to be 1Ma.

        Parameters
        ----------
        time : float
            The reconstruction time (Ma) at which to query subduction convergence.

        tessellation_threshold_radians : float, default=0.001 
            The threshold sampling distance along the subducting trench (in radians).

        Returns
        -------
        subduction_data : a list of vertically-stacked tuples
            The results for all tessellated points sampled along the trench.
            The size of the returned list is equal to the number of tessellated points.
            Each tuple in the list corresponds to a tessellated point and has the following tuple items:

            * longitude of sampled point
            * latitude of sampled point
            * subducting convergence (relative to trench) velocity magnitude (in cm/yr)
            * subducting convergence velocity obliquity angle (angle between trench normal vector and convergence 
            velocity vector)
            * trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
            * trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute 
            velocity vector)
            * length of arc segment (in degrees) that current point is on
            * trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
            * subducting plate ID
            * trench plate ID


        Notes
        -----
        Each point in the output is the midpoint of a great circle arc between two adjacent points in the trench polyline.
        The trench normal vector used in the obliquity calculations is perpendicular to the great circle arc of each point 
        (arc midpoint) and pointing towards the overriding plate (rather than away from it).

        Each trench is sampled at approximately uniform intervals along its length (specified via a threshold sampling distance).
        The sampling along the entire length of a trench is not exactly uniform. Each segment along a trench is sampled
        such that the samples have a uniform spacing that is less than or equal to the threshold sampling distance. However each
        segment in a trench might have a slightly different spacing distance (since segment lengths are not integer multiples of
        the threshold sampling distance).

        The trench normal (at each arc segment mid-point) always points *towards* the overriding plate.
        """
        if ignore_warnings:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                subduction_data = ptt.subduction_convergence.subduction_convergence(
                    self.rotation_model,
                    self.topology_features,
                    tessellation_threshold_radians,
                    float(time),
                    **kwargs)

        else:
            subduction_data = ptt.subduction_convergence.subduction_convergence(
                self.rotation_model,
                self.topology_features,
                tessellation_threshold_radians,
                float(time),
                **kwargs)

        subduction_data = np.vstack(subduction_data)
        return subduction_data


    def tesselate_mid_ocean_ridges(self, time, tessellation_threshold_radians=0.001, ignore_warnings=False, **kwargs):
        """Samples points along resolved spreading features (e.g. mid-ocean ridges) and calculates spreading rate and 
        length of ridge segments at a particular geological time.
         
        Resolves topologies at 'time', tessellates all resolved spreading features to within 'tessellation_threshold_radians'
        radians and obtains the following data:
    
        0 - longitude of sampled point
        1 - latitude of sampled point
        2 - spreading velocity magnitude (in cm/yr)
        3 - length of arc segment (in degrees) that current point is on
        
        All spreading feature types are considered. The transform segments of spreading features are ignored. 
        Note: by default, the function assumes that a segment can deviate 45 degrees from the stage pole before it is 
        considered a transform segment.

        Parameters
        ----------
        time : float
            The reconstruction time (Ma) at which to query subduction convergence.

        tessellation_threshold_radians : float, default=0.001 
            The threshold sampling distance along the subducting trench (in radians).

        anchor_plate_id : int, default=0
            The anchor plate of the reconstruction model.

        Returns
        -------
        ridge_data : a list of vertically-stacked tuples
            The results for all tessellated points sampled along the trench.
            The size of the returned list is equal to the number of tessellated points.
            Each tuple in the list corresponds to a tessellated point and has the following tuple items:
            
            * longitude of sampled point
            * latitude of sampled point
            * spreading velocity magnitude (in cm/yr)
            * length of arc segment (in degrees) that current point is on
        """
        if ignore_warnings:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                ridge_data = ptt.ridge_spreading_rate.spreading_rates(
                    self.rotation_model,
                    self.topology_features,
                    float(time),
                    tessellation_threshold_radians,
                    **kwargs)

        else:
            ridge_data = ptt.ridge_spreading_rate.spreading_rates(
                self.rotation_model,
                self.topology_features,
                float(time),
                tessellation_threshold_radians,
                **kwargs)

        ridge_data = np.vstack(ridge_data)
        return ridge_data

    def reconstruct(self, feature, to_time, from_time=0, anchor_plate_id=0, **kwargs):
        """Reconstructs regular geological features, motion paths or flowlines to a specific geological time.
        
        Parameters
        ----------
        feature : :class:`FeatureCollection`, or string, or :class:`Feature`, or sequence of :class:`Feature`, or sequence 
        of any combination of those four types.
            The features to reconstruct. Can be provided as a feature collection, or filename, or feature, or sequence of
            features, 
            or a sequence (eg, a list or tuple) of any combination of those four types.

        to_time : float, or :class:`GeoTimeInstant`
            The specific geological time to reconstruct to.

        from_time : float, default=0
            The specific geological time to reconstruct from. By default, this is set to present day. Raises 
            NotImplementedError if from_time not equal to 0.0.

        anchor_plate_id : int, default=0
            The anchor plate of the reconstruction model.

        **reconstruct_type : ReconstructType, default=ReconstructType.feature_geometry
            The specific reconstruction type to generate based on input feature geometry type. Can be provided as 
            pygplates.ReconstructType.feature_geometry to only reconstruct regular feature geometries, or
            pygplates.ReconstructType.motion_path to only reconstruct motion path features, or 
            pygplates.ReconstructType.flowline to only reconstruct flowline features. 
            Generates :class:`reconstructed feature geometries<ReconstructedFeatureGeometry>’, or :class:`reconstructed 
            motion paths<ReconstructedMotionPath>’, or :class:`reconstructed flowlines<ReconstructedFlowline>’ respectively.

        **group_with_feature : bool, default=False
            Used to group reconstructed geometries with their features. This can be useful when a feature has more than one
            geometry and hence more than one reconstructed geometry. The output *reconstructed_geometries* then becomes a 
            list of tuples where each tuple contains a :class:`feature<Feature>` and a ``list`` of reconstructed geometries. 
            Note: this keyword argument only applies when *reconstructed_geometries* is a list because exported files are 
            always grouped with their features. This is applicable to all ReconstructType features.
        
        **export_wrap_to_dateline : bool, default=True
            Wrap/clip reconstructed geometries to the dateline (currently ignored).

        Returns
        -------
        reconstructed_features : list
            Reconstructed geometrical features (generated by the reconstruction) are appended to the Python list.
            The reconstructed geometries are output in the same order as that of their respective input features (in the 
            parameter “features”). The order across input feature collections is also retained. This happens regardless 
            of whether *features* and *reconstructed_features* include files or not. Note: if keyword argument 
            group_with_feature=True then the list contains tuples that group each :class:`feature<Feature>` with a list 
            of its reconstructed geometries.

        Raises
        ------
        NotImplementedError
            if the starting time for reconstruction “from_time” not equal to 0.0
        """
        from_time, to_time = float(from_time), float(to_time)
        if from_time != 0.0:
            raise NotImplementedError("Soon...")

        reconstructed_features = []
        pygplates.reconstruct(feature, self.rotation_model, reconstructed_features, to_time,\
            anchor_plate_id=anchor_plate_id, **kwargs)
        return reconstructed_features


    def get_point_velocities(self, lons, lats, time, delta_time=1.0):
        """Generates a velocity domain feature collection, resolves them into points, and calculates the north and east 
        components of the velocity vector for each point in the domain at a particular geological time. 
        
        Velocity domain feature collections are MeshNode-type features. These are produced from lat-lon points represented as
        multi-point geometries (projections onto the surface of the unit length sphere). These features are resolved into 
        domain points and assigned plate IDs, which are used to obtain the equivalent stage rotations of identified tectonic
        plates over a time interval. Each velocity domain point and its stage rotation are used to calculate its velocity at 
        a particular geological time. Obtained velocities for each domain point are represented in the north-east-down 
        coordinate system. 

        Parameters
        ----------
        lons : array
            A 1D array of longitude points.

        lats : array
            A 1D array of latitude points.

        time : float
            The specific geological time (Ma) at which to calculate plate velocities.

        delta_time : float, default=1.0
            The time increment used for generating partitioning plate stage rotations. 1.0Ma by default.

        Returns
        -------
        all_velocities : 1D list of tuples
            For each velocity domain feature point, a tuple of (north, east, down) velocity components is generated and 
            appended to a list of velocity data. The length of all_velocities is equivalent to the number of domain points
            resolved from the lat-lon array parameters.
        """
        # Add points to a multipoint geometry

        time = float(time)

        multi_point = pygplates.MultiPointOnSphere([(float(lat),float(lon)) for lat, lon in zip(lats,lons)])

        # Create a feature containing the multipoint feature, and defined as MeshNode type
        meshnode_feature = pygplates.Feature(pygplates.FeatureType.create_from_qualified_string('gpml:MeshNode'))
        meshnode_feature.set_geometry(multi_point)
        meshnode_feature.set_name('Velocity Mesh Nodes from pygplates')

        velocity_domain_features = pygplates.FeatureCollection(meshnode_feature)
        
        # NB: at this point, the feature could be written to a file using
        # output_feature_collection.write('myfilename.gpmlz')
        
        
        # All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
        all_domain_points = []
        all_velocities = []

        # Partition our velocity domain features into our topological plate polygons at the current 'time'.
        plate_partitioner = pygplates.PlatePartitioner(self.topology_features, self.rotation_model, time)

        for velocity_domain_feature in velocity_domain_features:
            # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
            # Iterate over them all.
            for velocity_domain_geometry in velocity_domain_feature.get_geometries():

                for velocity_domain_point in velocity_domain_geometry.get_points():

                    all_domain_points.append(velocity_domain_point)

                    partitioning_plate = plate_partitioner.partition_point(velocity_domain_point)
                    if partitioning_plate:

                        # We need the newly assigned plate ID
                        # to get the equivalent stage rotation of that tectonic plate.
                        partitioning_plate_id = partitioning_plate.get_feature().get_reconstruction_plate_id()

                        # Get the stage rotation of partitioning plate from 'time + delta_time' to 'time'.
                        equivalent_stage_rotation = self.rotation_model.get_rotation(time,
                                                                                     partitioning_plate_id,
                                                                                     time + delta_time)

                        # Calculate velocity at the velocity domain point.
                        # This is from 'time + delta_time' to 'time' on the partitioning plate.
                        velocity_vectors = pygplates.calculate_velocities(
                            [velocity_domain_point],
                            equivalent_stage_rotation,
                            delta_time)

                        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples
                        # (one tuple per point).
                        velocities =pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                                [velocity_domain_point],
                                velocity_vectors)
                        all_velocities.append((velocities[0].get_x(), velocities[0].get_y()))

                    else:
                        all_velocities.append((0,0))
                        
        return np.array(all_velocities)



class Points(object):
    """The Points class offers simple methods to work with point data. It reconstructs geological features and can extract their
    plate velocities at a specific geological time. Note that the Points class must be used with the PlateReconstruction class
    since the rotation model and static polygons needed for Point object methods are sourced from the PlateReconstruction
    class.

    Attributes
    ----------
    PlateReconstruction_object : object pointer
    lons, lats : 1d array
    time : float
    plate_id : int
    x, y, z : 1d array
    lonlat, xyz : list
    rotation_model : str, or list (accessed using PlateReconstruction_object.rotation_model)
    static_polygons : str, or list (accessed using PlateReconstruction_object.static_polygons)
    features : list
    FeatureCollection : list 

    Methods
    -------
     __init__(self, PlateReconstruction_object, lons, lats, time=0, plate_id=None)
        Constructs all necessary attributes for the points object.
        
    reconstruct(self, time, anchor_plate_id=0, **kwargs)
        Reconstructs regular geological features, motion paths or flowlines to a specific geological time and extracts
        the latitudinal and longitudinal points of these features.
        
    plate_velocity(self, time, delta_time=1)
        Calculates the x and y components of tectonic plate velocities at a particular geological time.
        
    save(self, filename)
        Saves the feature collection used in the Points object under a given filename to the current directory. 
    """
    def __init__(self, PlateReconstruction_object, lons, lats, time=0, plate_id=None):
        """Constructs all necessary attributes for the points object.

        Parameters
        ----------
        PlateReconstruction_object : object pointer
            Allows for the accessibility of PlateReconstruction object attributes. Namely, PlateReconstruction object 
            attributes rotation_model, topology_featues and static_polygons can be used in the points object if called using
            “self.PlateReconstruction_object.X”, where X is the attribute.

        lons : float, or 1D array
            A single point, or a 1D array of longitude points.

        lats : float 1D array
            A single point, or a 1D array of latitude points.

        time : float, default=0
            The specific geological time (Ma) at which to start reconstructing. Default time is present day (0).

        plate_id : int, default=None
            The plate ID of a particular tectonic plate. Defaults to none.

        Returns
        -------
        An extension of accessible points object attributes, such as:

        x, y, z : float, or 1D array
            Cartesian coordinate equivalents of supplied lat, lon points scaled to mean Earth radius in km.

        lonlat, xyz : list
            Concatenated arrays of [lat, lon] points or their [x, y, z] Cartesian coordinate equivalents. Each list element 
            can be a float or 1D array.

        rotation_model, static_polygons : :class:`FeatureCollection`, or str, or :class:`Feature`, or sequence of
        :class:`Feature`, or a sequence of any combination of those four types, default=None
            Can be provided as a rotation model / static polygon feature collection, or optional filename, or a single feature, 
            or a sequence of features. Accessible with the points object after calling 
            “PlateReconstruction_object.rotation_model”, or “PlateReconstruction_object.static_polygons”.

        features : a sequence (eg, list or tuple) of :class:`Feature`, or a single :class:`Feature`
            A single, or list of point features with spherical geometry generated from a given tectonic plate ID 
            (default=None) elsewhere (eg. in another method) and given a set of lat-lon point(s). 
            If a plate ID is not supplied, "features" is instead a single set or list of features partitioned using 
            static polygons. 

        FeatureCollection : list
            A list of a set of features aggregated into a feature collection. 
        """
        self.lons = lons
        self.lats = lats
        self.time = time

        # get Cartesian coordinates
        self.x, self.y, self.z = _tools.lonlat2xyz(lons, lats, degrees=False)

        # scale by average radius of the Earth
        self.x *= _tools.EARTH_RADIUS
        self.y *= _tools.EARTH_RADIUS
        self.z *= _tools.EARTH_RADIUS

        # store concatenated arrays
        self.lonlat = np.c_[self.lons, self.lats]
        self.xyz = np.c_[self.x, self.y, self.z]


        rotation_model = PlateReconstruction_object.rotation_model
        static_polygons = PlateReconstruction_object.static_polygons
        self.PlateReconstruction_object = PlateReconstruction_object

        features = _tools.points_to_features(lons, lats, plate_id)

        if plate_id:
            self.features = features
        else:
            # partition using static polygons
            # being careful to observe 'from time'
            partitioned_features = pygplates.partition_into_plates(static_polygons, rotation_model, features)
            self.features = partitioned_features

            plate_id = np.empty(len(self.lons), dtype=int)
            for i, feature in enumerate(partitioned_features):
                plate_id[i] = feature.get_reconstruction_plate_id()

        self.plate_id = plate_id
        self.FeatureCollection = pygplates.FeatureCollection(self.features)


    def reconstruct(self, time, anchor_plate_id=0, **kwargs):
        """Reconstructs regular geological features, motion paths or flowlines to a specific geological time and extracts 
        the latitudinal and longitudinal points of these features.

        Note: this method accesses and uses the rotation model attribute from the PointReconstruction object, and reconstructs 
        the feature lat-lon point attributes of the Points object.

        Parameters
        ----------
        time : float
            The specific geological time (Ma) to reconstruct features to.

        anchor_plate_id : int, default=0
            The anchor plate of the reconstruction model.

        **reconstruct_type : ReconstructType, default=ReconstructType.feature_geometry
            The specific reconstruction type to generate based on input feature geometry type. Can be provided as
            ReconstructType.feature_geometry to only reconstruct regular feature geometries, or ReconstructType.MotionPath to
            only reconstruct motion path features, or ReconstructType.Flowline to only reconstruct flowline features. Generates
            :class:`reconstructed feature geometries<ReconstructedFeatureGeometry>’, or :class:`reconstructed motion
            paths<ReconstructedMotionPath>’, or :class:`reconstructed flowlines<ReconstructedFlowline>’ respectively.

        **group_with_feature : bool, default=False
            Used to group reconstructed geometries with their features. This can be useful when a feature has more than one
            geometry and hence more than one reconstructed geometry. The output *reconstructed_geometries* then becomes a 
            list of tuples where each tuple contains a :class:`feature<Feature>` and a ``list`` of reconstructed geometries. 
            Note: this keyword argument only applies when *reconstructed_geometries* is a list because exported files are 
            always grouped with their features. This is applicable to all ReconstructType features.

        **export_wrap_to_dateline : bool, default=True
            Wrap/clip reconstructed geometries to the dateline (currently ignored).

        Returns
        -------
        rlons, rlats : lists
            Two 1D numpy arrays enclosing all reconstructed feature points transformed into lat-lon points. 

        Raises
        ------
        NotImplementedError
            if the starting time for reconstruction “from_time” not equal to 0.0
        """
        from_time = self.time
        to_time = time
        reconstructed_features = self.PlateReconstruction_object.reconstruct(
            self.features, to_time, from_time, anchor_plate_id=anchor_plate_id, **kwargs)

        rlons, rlats = _tools.extract_feature_lonlat(reconstructed_features)
        return rlons, rlats


    def reconstruct_to_birth_age(self, ages, anchor_plate_id=0, **kwargs):
        from_time = self.time
        ages = np.array(ages)

        if len(ages) != len(self.features):
            raise ValueError("Number of points and ages must be identical")

        unique_ages = np.unique(ages)
        rlons = np.zeros(ages.shape)
        rlats = np.zeros(ages.shape)

        for age in unique_ages:
            mask_age = ages == age

            reconstructed_features = self.PlateReconstruction_object.reconstruct(
                self.features, age, from_time, anchor_plate_id=anchor_plate_id, **kwargs)

            lons, lats = _tools.extract_feature_lonlat(reconstructed_features)

            rlons[mask_age] = lons[mask_age]
            rlats[mask_age] = lats[mask_age]


        return rlons, rlats

    def plate_velocity(self, time, delta_time=1):
        """Calculates the x and y components of tectonic plate velocities at a particular geological time.

        This method accesses and uses the rotation_model attribute from the PointReconstruction object, and uses the features
        attribute of this Points object. Feature points are extracted and assigned plate IDs that are used to obtain the
        equivalent stage rotations of identified tectonic plates over a time interval. Each feature point and its stage rotation
        are used to calculate its plate velocity at a particular geological time. Obtained velocities for each domain point are
        represented in the north-east-down coordinate system, and their x,y Cartesian coordinate components are extracted. 

        Parameters
        ----------
        time : float
            The specific geological time (Ma) at which to calculate plate velocities.

        delta_time : float, default=1.0
            The time increment used for generating partitioning plate stage rotations. 1.0Ma by default.
            

        Returns
        -------
        all_velocities.T : 2D numpy list
            A transposed 2D numpy list with two rows and a number of columns equal to the number of x,y Cartesian velocity 
            components obtained (and thus the number of feature points extracted from a supplied feature). Each list column 
            stores one point’s x,y, velocity components along its two rows.
        """
        rotation_model = self.PlateReconstruction_object.rotation_model
        all_velocities = np.empty((len(self.features), 2))

        for i, feature in enumerate(self.features):
            geometry = feature.get_geometry()
            partitioning_plate_id = feature.get_reconstruction_plate_id()
            equivalent_stage_rotation = rotation_model.get_rotation(time, partitioning_plate_id, time+delta_time)
            
            velocity_vectors = pygplates.calculate_velocities(
                [geometry],
                equivalent_stage_rotation,
                delta_time,
                pygplates.VelocityUnits.cms_per_yr)
            
            velocities = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                [geometry],
                velocity_vectors)

            all_velocities[i] = velocities[0].get_x(), velocities[0].get_y()

        return list(all_velocities.T)


    def save(self, filename):
        """Saves the feature collection used in the Points object under a given filename to the current directory. 

        The needed file format to save to is determined from the filename extension. 

        Parameters
        ----------
        filename : string
            Can be provided as a string including the filename and the file format needed.

        Returns
        -------
        Feature collection saved under given filename to current directory.
        """
        self.FeatureCollection.save(filename)
