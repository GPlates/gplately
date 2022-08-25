"""Tools that wrap up pyGplates and Plate Tectonic Tools functionalities for reconstructing features,
 working with point data, and calculating plate velocities at specific geological times. 
"""
import pygplates
import numpy as np
import ptt
import warnings

from . import tools as _tools



class PlateReconstruction(object):
    """The `PlateReconstruction` class contains methods to reconstruct topology features to specific 
    geological times given a `rotation_model`, a set of `topology_features` and a set of 
    `static_polygons`. Topological plate velocity data at specific geological times can also be 
    calculated from these reconstructed features. 

    Attributes
    ----------
    rotation_model : str, or instance of <pygplates.FeatureCollection>, or <pygplates.Feature>, or sequence of <pygplates.Feature>, or instance of <pygplates.RotationModel>, default None
        A rotation model to query equivalent and/or relative topological plate rotations
        from a time in the past relative to another time in the past or to present day. Can be 
        provided as a rotation filename, or rotation feature collection, or rotation feature, or 
        sequence of rotation features, or a sequence (eg, a list or tuple) of any combination of 
        those four types.
    topology_features : str, or a sequence (eg, `list` or `tuple`) of instances of <pygplates.Feature>, or a single instance of <pygplates.Feature>, or an instance of <pygplates.FeatureCollection>, default None
        Reconstructable topological features like trenches, ridges and transforms. Can be provided 
        as an optional topology-feature filename, or sequence of features, or a single feature. 
    static_polygons : str, or instance of <pygplates.Feature>, or sequence of <pygplates.Feature>,or an instance of <pygplates.FeatureCollection>, default None
        Present-day polygons whose shapes do not change through geological time. They are
        used to cookie-cut dynamic polygons into identifiable topological plates (assigned 
        an ID) according to their present-day locations. Can be provided as a static polygon feature 
        collection, or optional filename, or a single feature, or a sequence of
        features.
    """
    
    def __init__(self, rotation_model=None, topology_features=None, static_polygons=None):
        rotation_model = pygplates.RotationModel(rotation_model)

        default_topology_features = pygplates.FeatureCollection()
        for topology in topology_features:
            default_topology_features.add( pygplates.FeatureCollection(topology) )

        self.rotation_model = rotation_model
        self.topology_features = default_topology_features
        self.static_polygons = static_polygons


    def tesselate_subduction_zones(self, time, tessellation_threshold_radians=0.001, ignore_warnings=False, **kwargs):
        """Samples points along subduction zone trenches and obtains subduction data at a particular
        geological time.
        
        Resolves topologies at `time`, tessellates all resolved subducting features to within 'tessellation_threshold_radians'
        radians and obtains the following information for each sampled point along a trench:
    
        `tesselate_subduction_zones` returns a list of 10 vertically-stacked tuples with the following data per sampled trench point:

        * Col. 0 - longitude of sampled trench point
        * Col. 1 - latitude of sampled trench point
        * Col. 2 - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
        * Col. 3 - subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
        * Col. 4 - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
        * Col. 5 - trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
        * Col. 6 - length of arc segment (in degrees) that current point is on
        * Col. 7 - trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
        * Col. 8 - subducting plate ID
        * Col. 9 - trench plate ID


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

            * Col. 0 - longitude of sampled trench point
            * Col. 1 - latitude of sampled trench point
            * Col. 2 - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
            * Col. 3 - subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
            * Col. 4 - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
            * Col. 5 - trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
            * Col. 6 - length of arc segment (in degrees) that current point is on
            * Col. 7 - trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
            * Col. 8 - subducting plate ID
            * Col. 9 - trench plate ID

        Notes
        -----
        Each sampled point in the output is the midpoint of a great circle arc between two adjacent points in the trench polyline.
        The trench normal vector used in the obliquity calculations is perpendicular to the great circle arc of each point 
        (arc midpoint) and pointing towards the overriding plate (rather than away from it).

        Each trench is sampled at approximately uniform intervals along its length (specified via a threshold sampling distance).
        The sampling along the entire length of a trench is not exactly uniform. Each segment along a trench is sampled
        such that the samples have a uniform spacing that is less than or equal to the threshold sampling distance. However each
        segment in a trench might have a slightly different spacing distance (since segment lengths are not integer multiples of
        the threshold sampling distance).

        The trench normal (at each arc segment mid-point) always points *towards* the overriding plate.
        The obliquity angles are in the range (-180 180). The range (0, 180) goes clockwise (when viewed from above the Earth)
        from the trench normal direction to the velocity vector. The range (0, -180) goes counter-clockwise.
        You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.
        The trench normal is perpendicular to the trench and pointing toward the overriding plate.
    
        Note that the convergence velocity magnitude is negative if the plates are diverging (if convergence obliquity angle 
        is greater than 90 or less than -90). And note that the absolute velocity magnitude is negative if the trench 
        (subduction zone) is moving towards the overriding plate (if absolute obliquity angle is less than 90 or greater 
        than -90) - note that this ignores the kinematics of the subducting plate.
        
        The delta time interval used for velocity calculations is, by default, assumed to be 1Ma.
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


    def total_subduction_zone_length(self, time, use_ptt=False, ignore_warnings=False):
        """Calculates the total length of all mid-ocean ridges (km) at the specified geological time (Ma).

        if `use_ptt` is True
        
        Uses Plate Tectonic Tools' `subduction_convergence` module to calculate trench segment lengths on a unit sphere. 
        The aggregated total subduction zone length is scaled to kilometres using the geocentric radius.

        Otherwise

        Resolves topology features ascribed to the `PlateReconstruction` model and extracts their shared boundary sections.
        The lengths of each trench boundary section are appended to the total subduction zone length.
        The total length is scaled to kilometres using a latitude-dependent (geocentric) Earth radius.


        Parameters
        ----------
        time : int
            The geological time at which to calculate total mid-ocean ridge lengths.
        use_ptt : bool, default=False
            If set to `True`, the PTT method is used.
        ignore_warnings : bool, default=False
            Choose whether to ignore warning messages from PTT's `subduction_convergence` workflow. These warnings alert the user 
            when certain subduction sub-segments are ignored - this happens when the trench segments have unidentifiable subduction 
            polarities and/or subducting plates. 

        Raises
        ------
        ValueError
            If neither `use_pygplates` or `use_ptt` have been set to `True`.

        Returns
        -------
        total_subduction_zone_length_kms : float
            The total subduction zone length (in km) at the specified `time`.

        """
        if use_ptt:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                subduction_data = self.tesselate_subduction_zones(time, ignore_warnings=ignore_warnings)

            trench_arcseg = subduction_data[:,6]
            trench_pt_lat = subduction_data[:,1]
            
            total_subduction_zone_length_kms = 0
            for i, segment in enumerate(trench_arcseg):
                earth_radius = _tools.geocentric_radius(trench_pt_lat[i])/1e3
                total_subduction_zone_length_kms += np.deg2rad(segment)*earth_radius 
                
            return total_subduction_zone_length_kms

        else:
            resolved_topologies = []
            shared_boundary_sections = []
            pygplates.resolve_topologies(self.topology_features, self.rotation_model, resolved_topologies, time, shared_boundary_sections)

            total_subduction_zone_length_kms = 0.0
            for shared_boundary_section in shared_boundary_sections:
                if shared_boundary_section.get_feature().get_feature_type() != pygplates.FeatureType.gpml_subduction_zone:
                    continue
                for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                    clat, clon = shared_sub_segment.get_resolved_geometry().get_centroid().to_lat_lon()
                    earth_radius = _tools.geocentric_radius(clat) / 1e3
                    total_subduction_zone_length_kms += shared_sub_segment.get_resolved_geometry().get_arc_length()*earth_radius

            return total_subduction_zone_length_kms


    def total_continental_arc_length(self, time, continental_grid=None, trench_arc_distance=0.0, ignore_warnings=True):
        """Calculates the total length of all global continental arcs (km) at the specified geological time (Ma).

        Uses Plate Tectonic Tools' `subduction_convergence` workflow to sample a given plate model's trench features into 
        point features and obtain their subduction polarities. The resolved points are projected out by the `trench_arc_distance`
        and their new locations are linearly interpolated onto the supplied `continental_grid`. If the projected trench 
        points lie in the grid, they are considered continental arc points, and their arc segment lengths are appended 
        to the total continental arc length for the specified `time`. The total length is scaled to kilometres using the geocentric 
        Earth radius. 

        Parameters
        ----------
        time : int
            The geological time at which to calculate total continental arc lengths.
        continental_grid: str or MaskedArray, default=None
            A MaskedArray or full string path to a continental grid with which to interpolate projected trench points on 
            (thereby identifying continental arc points). 
        trench_arc_distance : float, default=0.0
            The trench-to-arc distance (in kilometres) to project sampled trench points out by in the direction of their 
            subduction polarities. 
        ignore_warnings : bool, default=False
            Choose whether to ignore warning messages from PTT's subduction_convergence workflow that alerts the user of 
            subduction sub-segments that are ignored due to unidentified polarities and/or subducting plates. 

        Raises
        ------
        ValueError
            * If a continental grid directory is not supplied.
            * If the trench_arc_distance is not supplied or kept at 0.0km.

        Returns
        -------
        total_continental_arc_length_kms : float
            The continental arc length (in km) at the specified time.
        """
        from . import grids as _grids
        # Right now, raster reconstruction is not supported.
        if continental_grid is None:
            raise ValueError("Please provide a directory to a continental grid or a masked array for the current time.")
        
        elif isinstance(continental_grid, np.ma.MaskedArray):
            # Process the masked continental grid
            graster = _grids.Raster(PlateReconstruction, array=continental_grid, extent=[-180,180,-90,90])

        elif isinstance(continental_grid, str):
            # Process the continental grid directory
            graster = _grids.Raster(PlateReconstruction, continental_grid, extent=[-180,180,-90,90])

        # Obtain trench data with Plate Tectonic Tools
        trench_data = self.tesselate_subduction_zones(time, ignore_warnings=ignore_warnings)
        
        # Extract trench data
        trench_normal_azimuthal_angle = trench_data[:,7]
        trench_arcseg = trench_data[:,6]
        trench_pt_lon = trench_data[:,0]
        trench_pt_lat = trench_data[:,1]
        
        # Modify the trench-arc distance using the geocentric radius
        arc_distance = trench_arc_distance / (_tools.geocentric_radius(trench_pt_lat)/1000)
        
        # Project trench points out along trench-arc distance, and obtain their new lat-lon coordinates
        dlon = arc_distance*np.sin(np.radians(trench_normal_azimuthal_angle))
        dlat = arc_distance*np.cos(np.radians(trench_normal_azimuthal_angle))
        ilon = trench_pt_lon + np.degrees(dlon)
        ilat = trench_pt_lat + np.degrees(dlat)
        
        # Linearly interpolate projected points onto continental grids, and collect the indices of points that lie
        # within the grids.
        sampled_points = graster.interpolate(ilon, ilat, method='linear', return_indices=True, return_distances=False)     
        in_raster = [i for i, point in enumerate(sampled_points[0]) if point > 0]
        
        # Define arrays + total arc length
        lat_in = []
        lon_in = []
        total_continental_arc_length_kms = 0
        subd_we_count_lat = []
        subd_we_count_lon = []
        
        # Loop through all successful in-raster indices
        for index in in_raster:
            
            # Get the lat-lon coordinate of the in-raster point, and the corresponding trench point
            lat_in.append(ilat[index])
            lon_in.append(ilon[index])
            subd_we_count_lat.append(trench_pt_lat[index])
            subd_we_count_lon.append(trench_pt_lon[index])

            # Append the continental trench segment to the total arc length
            earth_radius = _tools.geocentric_radius(trench_pt_lat[index])/1000
            total_continental_arc_length_kms += np.deg2rad(trench_data[:,6][index])*earth_radius
            
        return total_continental_arc_length_kms


    def tesselate_mid_ocean_ridges(self, time, tessellation_threshold_radians=0.001, ignore_warnings=False, **kwargs):
        """Samples points along resolved spreading features (e.g. mid-ocean ridges) and calculates spreading rates and 
        lengths of ridge segments at a particular geological time.
         
        Resolves topologies at `time`, tessellates all resolved spreading features to within 'tessellation_threshold_radians'
        radians. Returns a 4-column vertically stacked tuple with the following data.

        * Col. 0 - longitude of sampled ridge point
        * Col. 1 - latitude of sampled ridge point
        * Col. 2 - spreading velocity magnitude (in cm/yr)
        * Col. 3 - length of arc segment (in degrees) that current point is on
        
        All spreading feature types are considered. The transform segments of spreading features are ignored. 
        Note: by default, the function assumes that a segment can deviate 45 degrees from the stage pole before it is 
        considered a transform segment.

        Parameters
        ----------
        time : float
            The reconstruction time (Ma) at which to query subduction convergence.

        tessellation_threshold_radians : float, default=0.001 
            The threshold sampling distance along the subducting trench (in radians).

        ignore_warnings : bool, default=False
            Choose to ignore warnings from Plate Tectonic Tools' ridge_spreading_rate workflow. 

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
                spreading_feature_types = [pygplates.FeatureType.gpml_mid_ocean_ridge]
                ridge_data = ptt.ridge_spreading_rate.spreading_rates(
                    self.rotation_model,
                    self.topology_features,
                    float(time),
                    tessellation_threshold_radians,
                    spreading_feature_types,
                    **kwargs)

        else:
            spreading_feature_types = [pygplates.FeatureType.gpml_mid_ocean_ridge]
            ridge_data = ptt.ridge_spreading_rate.spreading_rates(
                self.rotation_model,
                self.topology_features,
                float(time),
                tessellation_threshold_radians,
                spreading_feature_types,
                **kwargs)

        ridge_data = np.vstack(ridge_data)
        return ridge_data


    def total_ridge_length(self, time, use_ptt=False, ignore_warnings=False):
        """Calculates the total length of all mid-ocean ridges (km) at the specified geological time (Ma).

        if `use_ptt` is True
        
        Uses Plate Tectonic Tools' `ridge_spreading_rate` workflow to calculate ridge segment lengths. Scales lengths to
        kilometres using the geocentric radius.

        Otherwise

        Resolves topology features of the PlateReconstruction model and extracts their shared boundary sections.
        The lengths of each GPML mid-ocean ridge shared boundary section are appended to the total ridge length.
        Scales lengths to kilometres using the geocentric radius.


        Parameters
        ----------
        time : int
            The geological time at which to calculate total mid-ocean ridge lengths.
        use_ptt : bool, default=False
            If set to `True`, the PTT method is used. 
        ignore_warnings : bool, default=False
            Choose whether to ignore warning messages from PTT's `ridge_spreading_rate` workflow.

        Raises
        ------
        ValueError
            If neither `use_pygplates` or `use_ptt` have been set to `True`.

        Returns
        -------
        total_ridge_length_kms : float
            The total length of global mid-ocean ridges (in kilometres) at the specified time.
        """

        if use_ptt is True:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                ridge_data = self.tesselate_mid_ocean_ridges(time)

            ridge_arcseg = ridge_data[:,3]
            ridge_pt_lat = ridge_data[:,1]

            total_ridge_length_kms = 0
            for i, segment in enumerate(ridge_arcseg):
                earth_radius = _tools.geocentric_radius(ridge_pt_lat[i])/1e3
                total_ridge_length_kms += np.deg2rad(segment)*earth_radius 

            return total_ridge_length_kms

        else:
            resolved_topologies = []
            shared_boundary_sections = []
            pygplates.resolve_topologies(self.topology_features, self.rotation_model, resolved_topologies, time, shared_boundary_sections)

            total_ridge_length_kms = 0.0
            for shared_boundary_section in shared_boundary_sections:
                if shared_boundary_section.get_feature().get_feature_type() != pygplates.FeatureType.gpml_mid_ocean_ridge:
                    continue
                for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                    clat, clon = shared_sub_segment.get_resolved_geometry().get_centroid().to_lat_lon()
                    earth_radius = _tools.geocentric_radius(clat) / 1e3
                    total_ridge_length_kms += shared_sub_segment.get_resolved_geometry().get_arc_length()*earth_radius

            return total_ridge_length_kms


    def reconstruct(self, feature, to_time, from_time=0, anchor_plate_id=0, **kwargs):
        """Reconstructs regular geological features, motion paths or flowlines to a specific geological time.
        
        Parameters
        ----------
        feature : str, or instance of <pygplates.FeatureCollection>, or <pygplates.Feature>, or sequence of <pygplates.Feature>
            The geological features to reconstruct. Can be provided as a feature collection, or filename, 
            or feature, or sequence of features, or a sequence (eg, a list or tuple) of any combination of 
            those four types.

        to_time : float, or :class:`GeoTimeInstant`
            The specific geological time to reconstruct to.

        from_time : float, default=0
            The specific geological time to reconstruct from. By default, this is set to present day. Raises 
            `NotImplementedError` if `from_time` is not set to 0.0 Ma (present day).

        anchor_plate_id : int, default=0
            Reconstruct features with respect to a certain anchor plate. By default, reconstructions are made 
            with respect to the absolute reference frame, like a stationary geological element (e.g. a mantle
            plume). This frame is given the plate ID of 0. 

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
            Wrap/clip reconstructed geometries to the dateline.

        Returns
        -------
        reconstructed_features : list
            Reconstructed geological features (generated by the reconstruction) are appended to the list.
            The reconstructed geometries are output in the same order as that of their respective input features (in the 
            parameter “features”). The order across input feature collections is also retained. This happens regardless 
            of whether *features* and *reconstructed_features* include files or not. Note: if keyword argument 
            group_with_feature=True then the list contains tuples that group each :class:`feature<Feature>` with a list 
            of its reconstructed geometries.

        Raises
        ------
        NotImplementedError
            if the starting time for reconstruction `from_time` is not equal to 0.0.
        """
        from_time, to_time = float(from_time), float(to_time)

        reconstructed_features = []
        pygplates.reconstruct(feature, self.rotation_model, reconstructed_features, to_time,\
            anchor_plate_id=anchor_plate_id, **kwargs)
        return reconstructed_features


    def get_point_velocities(self, lons, lats, time, delta_time=1.0):
        """Generates a velocity domain feature collection, resolves them into points, and calculates the north and east 
        components of the velocity vector for each point in the domain at a particular geological `time`. 
        
        Notes
        -----
        Velocity domain feature collections are MeshNode-type features. These are produced from `lons` and `lats` pairs
        represented as multi-point geometries (projections onto the surface of the unit length sphere). These features are 
        resolved into  domain points and assigned plate IDs which are used to obtain the equivalent stage rotations of 
        identified tectonic plates over a time interval (`delta_time`). Each velocity domain point and its stage rotation 
        are used to calculate the point's plate velocity at a particular `time`. Obtained velocities for each domain point 
        are represented in the north-east-down coordinate system. 

        Parameters
        ----------
        lons : array
            A 1D array of point data's longitudes.

        lats : array
            A 1D array of point data's latitudes.

        time : float
            The specific geological time (Ma) at which to calculate plate velocities.

        delta_time : float, default=1.0
            The time increment used for generating partitioning plate stage rotations. 1.0Ma by default.

        Returns
        -------
        all_velocities : 1D list of tuples
            For each velocity domain feature point, a tuple of (north, east, down) velocity components is generated and 
            appended to a list of velocity data. The length of `all_velocities` is equivalent to the number of domain points
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


    def create_motion_path(self, lons, lats, time_array, plate_id=None, anchor_plate_id=0, return_rate_of_motion=False):
        """ Create a path of points to mark the trajectory of a plate's motion 
        through geological time.
        
        Parameters
        ----------
        lons : arr
            An array containing the longitudes of seed points on a plate in motion.
        lats : arr
            An array containing the latitudes of seed points on a plate in motion.
        time_array : arr
            An array of reconstruction times at which to determine the trajectory 
            of a point on a plate. For example:
                
                import numpy as np
                min_time = 30
                max_time = 100
                time_step = 2.5
                time_array = np.arange(min_time, max_time + time_step, time_step)

        plate_id : int, default=None
            The ID of the moving plate. If this is not passed, the plate ID of the 
            seed points are ascertained using pygplates' `PlatePartitioner`.
        anchor_plate_id : int, default=0
            The ID of the anchor plate.
        return_rate_of_motion : bool, default=False
            Choose whether to return the rate of plate motion through time for each
            
        Returns
        -------
        rlons : ndarray
            An n-dimensional array with columns containing the longitudes of 
            the seed points at each timestep in `time_array`. There are n 
            columns for n seed points. 
        rlats : ndarray
            An n-dimensional array with columns containing the latitudes of 
            the seed points at each timestep in `time_array`. There are n 
            columns for n seed points. 
        StepTimes
        StepRates
            
        Examples
        --------
        To access the latitudes and longitudes of each seed point's motion path:
        
            for i in np.arange(0,len(seed_points)):
                current_lons = lon[:,i]
                current_lats = lat[:,i]
        """
        lons = np.atleast_1d(lons)
        lats = np.atleast_1d(lats)
        time_array = np.atleast_1d(time_array)

        # ndarrays to fill with reconstructed points and 
        # rates of motion (if requested)
        rlons = np.empty((len(time_array), len(lons)))
        rlats = np.empty((len(time_array), len(lons)))
        StepTimes = np.empty(((len(time_array)-1)*2, len(lons)))
        StepRates = np.empty(((len(time_array)-1)*2, len(lons)))

        seed_points = list(zip(lats,lons))
        for i, lat_lon in enumerate(seed_points):

            seed_points_at_digitisation_time = pygplates.PointOnSphere(
                pygplates.LatLonPoint(float(lat_lon[0]), float(lat_lon[1]))
            )
            # Allocate the present-day plate ID to the PointOnSphere if 
            # it was not given.
            if plate_id is None:
                plate_id = _tools.plate_partitioner_for_point(
                    lat_lon,
                    self.topology_features,
                    self.rotation_model
                )
            # Create the motion path feature. enforce float and int for C++ signature.
            motion_path_feature = pygplates.Feature.create_motion_path(
                seed_points_at_digitisation_time, 
                time_array, 
                valid_time=(time_array.max(), time_array.min()),
                relative_plate=int(anchor_plate_id),
                reconstruction_plate_id=int(plate_id))

            reconstructed_motion_paths = self.reconstruct(
                motion_path_feature, 
                to_time=0, 
                reconstruct_type=pygplates.ReconstructType.motion_path,
                anchor_plate_id=int(anchor_plate_id),
            )
            # Turn motion paths in to lat-lon coordinates
            for reconstructed_motion_path in reconstructed_motion_paths:
                trail = reconstructed_motion_path.get_motion_path().to_lat_lon_array()

            lon, lat = np.flipud(trail[:, 1]), np.flipud(trail[:, 0])

            rlons[:,i] = lon
            rlats[:,i] = lat

            # Obtain step-plot coordinates for rate of motion
            if return_rate_of_motion is True:
                
                # Get timestep
                TimeStep = []
                for j in range(len(time_array) - 1):
                    diff = time_array[j + 1] - time_array[j]
                    TimeStep.append(diff)

                # Iterate over each segment in the reconstructed motion path, get the distance travelled by the moving
                # plate relative to the fixed plate in each time step
                Dist = []
                for reconstructed_motion_path in reconstructed_motion_paths:
                    for segment in reconstructed_motion_path.get_motion_path().get_segments():
                        Dist.append(segment.get_arc_length() * _tools.geocentric_radius(segment.get_start_point().to_lat_lon()[0]) / 1e3)

                # Note that the motion path coordinates come out starting with the oldest time and working forwards
                # So, to match our 'times' array, we flip the order
                Dist = np.flipud(Dist)

                # Get rate of motion as distance per Myr
                Rate = np.asarray(Dist)/TimeStep
                
                # Manipulate arrays to get a step plot
                StepRate = np.zeros(len(Rate)*2)
                StepRate[::2] = Rate
                StepRate[1::2] = Rate

                StepTime = np.zeros(len(Rate)*2)
                StepTime[::2] = time_array[:-1]
                StepTime[1::2] = time_array[1:]
                
                # Append the nth point's step time and step rate coordinates to the ndarray
                StepTimes[:,i] = StepTime
                StepRates[:,i] = StepRate

        if return_rate_of_motion is True:
            return np.squeeze(rlons), np.squeeze(rlats), np.squeeze(StepTimes), np.squeeze(StepRates)
        else:
            return np.squeeze(rlons), np.squeeze(rlats)


    def create_flowline(self, lons, lats, time_array, left_plate_ID, right_plate_ID, return_rate_of_motion=False):
        """ Create a path of points to track plate motion away from 
        spreading ridges over time using half-stage rotations.

        Parameters
        ----------
        lons : arr
            An array of longitudes of points along spreading ridges.
        lats : arr
            An array of latitudes of points along spreading ridges.
        time_array : arr
            A list of times to obtain seed point locations at.
        left_plate_ID : int
            The plate ID of the polygon to the left of the spreading
            ridge.
        right_plate_ID : int
            The plate ID of the polygon to the right of the spreading
            ridge.
        return_rate_of_motion : bool, default False
            Choose whether to return a step time and step rate array
            for a step plot of motion. 

        Returns
        -------
        left_lon : ndarray
            The longitudes of the __left__ flowline for n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in `time_array`. 
        left_lat : ndarray
            The latitudes of the __left__ flowline of n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in `time_array`.
        right_lon : ndarray
            The longitudes of the __right__ flowline of n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in `time_array`.
        right_lat : ndarray
            The latitudes of the __right__ flowline of n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in `time_array`.

        Examples
        --------
        To access the ith seed point's left and right latitudes and
        longitudes:

            for i in np.arange(0,len(seed_points)):
                left_flowline_longitudes = left_lon[:,i] 
                left_flowline_latitudes = left_lat[:,i] 
                right_flowline_longitudes = right_lon[:,i]
                right_flowline_latitudes = right_lat[:,i]
        """
        lats = np.atleast_1d(lats)
        lons = np.atleast_1d(lons)
        time_array = np.atleast_1d(time_array)
        
        seed_points = list(zip(lats, lons))
        multi_point = pygplates.MultiPointOnSphere(seed_points)
        
        start = 0
        if time_array[0] != 0:
            start = 1
            time_array = np.hstack([[0], time_array])

        # Create the flowline feature
        flowline_feature = pygplates.Feature.create_flowline(
            multi_point,
            time_array.tolist(),
            valid_time=(time_array.max(), time_array.min()),
            left_plate=left_plate_ID,
            right_plate=right_plate_ID)

        # reconstruct the flowline in present-day coordinates
        reconstructed_flowlines = self.reconstruct(
            flowline_feature, 
            to_time=0,
            reconstruct_type=pygplates.ReconstructType.flowline)

        # Wrap things to the dateline, to avoid plotting artefacts.
        date_line_wrapper = pygplates.DateLineWrapper()

        # Create lat-lon ndarrays to store the left and right lats and lons of flowlines
        left_lon  = np.empty((len(time_array), len(lons)))
        left_lat  = np.empty((len(time_array), len(lons)))
        right_lon = np.empty((len(time_array), len(lons)))
        right_lat = np.empty((len(time_array), len(lons)))
        StepTimes = np.empty(((len(time_array)-1)*2, len(lons)))
        StepRates = np.empty(((len(time_array)-1)*2, len(lons)))


        # Iterate over the reconstructed flowlines. Each seed point results in a 'left' and 'right' flowline 
        for i, reconstructed_flowline in enumerate(reconstructed_flowlines):

            # Get the points for the left flowline only
            left_latlon = reconstructed_flowline.get_left_flowline().to_lat_lon_array()
            left_lon[:,i] = left_latlon[:,1]
            left_lat[:,i] = left_latlon[:,0]

            # Repeat for the right flowline points
            right_latlon = reconstructed_flowline.get_right_flowline().to_lat_lon_array()
            right_lon[:,i] = right_latlon[:,1]
            right_lat[:,i] = right_latlon[:,0]

        if return_rate_of_motion:
            for i, reconstructed_motion_path in enumerate(reconstructed_flowlines):
                distance = []
                for segment in reconstructed_motion_path.get_left_flowline().get_segments():
                    distance.append(segment.get_arc_length() * _tools.geocentric_radius(segment.get_start_point().to_lat_lon()[0]) / 1e3)
                    
                # Get rate of motion as distance per Myr
                # Need to multiply rate by 2, since flowlines give us half-spreading rate
                time_step = time_array[1]-time_array[0]
                Rate = (np.asarray(distance)/time_step) * 2 # since we created the flowline at X increment

                # Manipulate arrays to get a step plot
                StepRate = np.zeros(len(Rate)*2)
                StepRate[::2] = Rate
                StepRate[1::2] = Rate

                StepTime = np.zeros(len(Rate)*2)
                StepTime[::2] = time_array[:-1]
                StepTime[1::2] = time_array[1:]

                # Append the nth point's step time and step rate coordinates to the ndarray
                StepTimes[:,i] = StepTime
                StepRates[:,i] = StepRate

            return left_lon[start:], left_lat[start:], right_lon[start:], right_lat[start:], StepTimes, StepRates

        else:
            return left_lon[start:], left_lat[start:], right_lon[start:], right_lat[start:]



class Points(object):
    """`Points` contains methods to reconstruct and work with with geological point data. For example, the 
    locations and plate velocities of point data can be calculated at a specific geological `time`. The `Points` 
    object requires the `PlateReconstruction` object to work because it holds the `rotation_model` and `static_polygons` 
    needed to classify topological plates and quantify feature rotations through time. 

    Attributes
    ----------
    PlateReconstruction_object : object pointer
        Allows for the accessibility of `PlateReconstruction` object attributes: `rotation_model`, `topology_featues` 
        and `static_polygons` for use in the `Points` object if called using “self.PlateReconstruction_object.X”, 
        where X is the attribute.

    lons : float, or 1D array
        A single float, or a 1D array containing the longitudes of point data.

    lats : float 1D array
        A single float, or a 1D array containing the latitudes of point data.

    time : float, default=0
        The specific geological time (Ma) at which to reconstruct the point data. By default, it is set to
        the present day (0 Ma).

    plate_id : int, default=None
        The plate ID of a particular tectonic plate on which point data lies, if known. This is obtained in `init` 
        if not provided.

    """
    def __init__(self, PlateReconstruction_object, lons, lats, time=0, plate_id=None):
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
            plate_id = np.atleast_1d(plate_id)
            self.features = features
        else:
            # partition using static polygons
            # being careful to observe 'from time'
            partitioned_features = pygplates.partition_into_plates(
                static_polygons,
                rotation_model,
                features,
                reconstruction_time=time)
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
            Reconstruct features with respect to a certain anchor plate. By default, reconstructions are made 
            with respect to the absolute reference frame, like a stationary geological element (e.g. a mantle
            plume). This frame is given the plate ID of 0.

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
        rlons : list of float
            A 1D numpy array enclosing all reconstructed point features' longitudes. 

        rlats : list of float
            A 1D numpy array enclosing all reconstructed point features' latitudes. 

        Raises
        ------
        NotImplementedError
            if the starting time for reconstruction `from_time` is not equal to 0.0
        """
        from_time = self.time
        to_time = time
        reconstructed_features = self.PlateReconstruction_object.reconstruct(
            self.features, to_time, from_time, anchor_plate_id=anchor_plate_id, **kwargs)

        rlons, rlats = _tools.extract_feature_lonlat(reconstructed_features)
        return rlons, rlats


    def reconstruct_to_birth_age(self, ages, anchor_plate_id=0, **kwargs):
        """ Reconstructs point features supplied to the `Points` object from the supplied initial time (`self.time`)
        to a range of times. The number of supplied times must equal the number of point features supplied to the Points object. 

        Attributes
        ----------
        ages : array
            Geological times to reconstruct features to. Must have the same length as the `Points `object's `self.features` attribute 
            (which holds all point features represented on a unit length sphere in 3D Cartesian coordinates).
        anchor_plate_id : int, default=0
            Reconstruct features with respect to a certain anchor plate. By default, reconstructions are made 
            with respect to the absolute reference frame, like a stationary geological element (e.g. a mantle
            plume). This frame is given the plate ID of 0.
        **kwargs 
            Additional keyword arguments for the `gplately.PlateReconstruction.reconstruct` method.

        Raises
        ------
        ValueError
            If the number of ages and number of point features supplied to the Points object are not identical.

        Returns
        -------
        rlons, rlats : float
            The longitude and latitude coordinate lists of all point features reconstructed to all specified ages.

        Examples
        --------
        To reconstruct n seed points' locations to B Ma (for this example n=2, with (lon,lat) = (78,30) and (56,22) at time=0 Ma,
        and we reconstruct to B=10 Ma):

            # Longitude and latitude of n=2 seed points
            pt_lon = np.array([78., 56])
            pt_lat = np.array([30., 22])

            # Call the Points object!
            gpts = gplately.Points(model, pt_lon, pt_lat)
            print(gpts.features[0].get_all_geometries())   # Confirms we have features represented as points on a sphere

            ages = numpy.linspace(10,10, len(pt_lon))
            rlons, rlats = gpts.reconstruct_to_birth_age(ages)

        """
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

        This method accesses and uses the `rotation_model` attribute from the `PlateReconstruction` object and uses the `Points` 
        object's `self.features` attribute. Feature points are extracted and assigned plate IDs. These IDs are used to obtain the
        equivalent stage rotations of identified tectonic plates over a time interval `delta_time`. Each feature point and its stage
        rotation are used to calculate the point's plate velocity at a particular geological time. Obtained velocities for each domain
        point are represented in the north-east-down coordinate system, and their x,y Cartesian coordinate components are extracted. 

        Parameters
        ----------
        time : float
            The specific geological time (Ma) at which to calculate plate velocities.

        delta_time : float, default=1.0
            The time increment used for generating partitioning plate stage rotations. 1.0 Ma by default.
            

        Returns
        -------
        all_velocities.T : 2D numpy list
            A transposed 2D numpy list with two rows and a number of columns equal to the number of x,y Cartesian velocity 
            components obtained (and thus the number of feature points extracted from a supplied feature). Each list column 
            stores one point’s x,y, velocity components along its two rows.
        """
        time = float(time)

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

            all_velocities[i] = velocities[0].get_y(), velocities[0].get_x()

        return list(all_velocities.T)


    def motion_path(self, time_array, anchor_plate_id=0, return_rate_of_motion=False):
        """ Create a path of points to mark the trajectory of a plate's motion 
        through geological time.
        
        Parameters
        ----------
        time_array : arr
            An array of reconstruction times at which to determine the trajectory 
            of a point on a plate. For example:
                
                import numpy as np
                min_time = 30
                max_time = 100
                time_step = 2.5
                time_array = np.arange(min_time, max_time + time_step, time_step)

        anchor_plate_id : int, default=0
            The ID of the anchor plate.
        return_rate_of_motion : bool, default=False
            Choose whether to return the rate of plate motion through time for each
            
        Returns
        -------
        rlons : ndarray
            An n-dimensional array with columns containing the longitudes of 
            the seed points at each timestep in `time_array`. There are n 
            columns for n seed points. 
        rlats : ndarray
            An n-dimensional array with columns containing the latitudes of 
            the seed points at each timestep in `time_array`. There are n 
            columns for n seed points. 
        """
        time_array = np.atleast_1d(time_array)
        
        # ndarrays to fill with reconstructed points and 
        # rates of motion (if requested)
        rlons = np.empty((len(time_array),   len(self.lons)))
        rlats = np.empty((len(time_array),   len(self.lons)))
        StepTimes = np.empty(((len(time_array)-1)*2, len(self.lons)))
        StepRates = np.empty(((len(time_array)-1)*2, len(self.lons)))

        seed_points = list(zip(self.lats, self.lons))
        for i, lat_lon in enumerate(seed_points):
            
            seed_points_at_digitisation_time = pygplates.PointOnSphere(
                pygplates.LatLonPoint(float(lat_lon[0]), float(lat_lon[1]))
            )

            # Create the motion path feature
            motion_path_feature = pygplates.Feature.create_motion_path(
                seed_points_at_digitisation_time, 
                time_array.tolist(),
                valid_time=(time_array.max(), time_array.min()),
                #relative_plate=int(self.plate_id[i]),
                #reconstruction_plate_id=int(anchor_plate_id))
                relative_plate=int(anchor_plate_id),
                reconstruction_plate_id=int(self.plate_id[i]))

            reconstructed_motion_paths = self.PlateReconstruction_object.reconstruct(
                motion_path_feature, 
                to_time=0, 
                #from_time=0, 
                reconstruct_type=pygplates.ReconstructType.motion_path)

            # Turn motion paths in to lat-lon coordinates
            for reconstructed_motion_path in reconstructed_motion_paths:
                trail = reconstructed_motion_path.get_motion_path().to_lat_lon_array()

            lon, lat = np.flipud(trail[:, 1]), np.flipud(trail[:, 0])

            rlons[:,i] = lon
            rlats[:,i] = lat

            # Obtain step-plot coordinates for rate of motion
            if return_rate_of_motion is True:
                
                # Get timestep
                TimeStep = []
                for j in range(len(time_array) - 1):
                    diff = time_array[j + 1] - time_array[j]
                    TimeStep.append(diff)

                # Iterate over each segment in the reconstructed motion path, get the distance travelled by the moving
                # plate relative to the fixed plate in each time step
                Dist = []
                for reconstructed_motion_path in reconstructed_motion_paths:
                    for segment in reconstructed_motion_path.get_motion_path().get_segments():
                        Dist.append(segment.get_arc_length() * _tools.geocentric_radius(segment.get_start_point().to_lat_lon()[0]) / 1e3)

                # Note that the motion path coordinates come out starting with the oldest time and working forwards
                # So, to match our 'times' array, we flip the order
                Dist = np.flipud(Dist)

                # Get rate of motion as distance per Myr
                Rate = np.asarray(Dist)/TimeStep
                
                # Manipulate arrays to get a step plot
                StepRate = np.zeros(len(Rate)*2)
                StepRate[::2] = Rate
                StepRate[1::2] = Rate

                StepTime = np.zeros(len(Rate)*2)
                StepTime[::2] = time_array[:-1]
                StepTime[1::2] = time_array[1:]
                
                # Append the nth point's step time and step rate coordinates to the ndarray
                StepTimes[:,i] = StepTime
                StepRates[:,i] = StepRate
        
        if return_rate_of_motion is True:
            return np.squeeze(rlons), np.squeeze(rlats), np.squeeze(StepTimes), np.squeeze(StepRates)
        else:
            return np.squeeze(rlons), np.squeeze(rlats)


    def flowline(self, time_array, left_plate_ID, right_plate_ID, return_rate_of_motion=False):
        """ Create a path of points to track plate motion away from 
        spreading ridges over time using half-stage rotations.
        
        Parameters
        ----------
        lons : arr
            An array of longitudes of points along spreading ridges.
        lats : arr
            An array of latitudes of points along spreading ridges.
        time_array : arr
            A list of times to obtain seed point locations at.
        left_plate_ID : int
            The plate ID of the polygon to the left of the spreading
            ridge.
        right_plate_ID : int
            The plate ID of the polygon to the right of the spreading
            ridge.
        return_rate_of_motion : bool, default False
            Choose whether to return a step time and step rate array for
            a step-plot of flowline motion.
            
        Returns
        -------
        left_lon : ndarray
            The longitudes of the __left__ flowline for n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in `time_array`. 
        left_lat : ndarray
            The latitudes of the __left__ flowline of n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in `time_array`.
        right_lon : ndarray
            The longitudes of the __right__ flowline of n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in `time_array`.
        right_lat : ndarray
            The latitudes of the __right__ flowline of n seed points.
            There are n columns for n seed points, and m rows
            for m time steps in `time_array`.
        
        Examples
        --------
        To access the ith seed point's left and right latitudes and
        longitudes:
        
            for i in np.arange(0,len(seed_points)):
                left_flowline_longitudes = left_lon[:,i] 
                left_flowline_latitudes = left_lat[:,i] 
                right_flowline_longitudes = right_lon[:,i]
                right_flowline_latitudes = right_lat[:,i]
        """
        model = self.PlateReconstruction_object
        return model.create_flowline(self.lons, self.lats, time_array, left_plate_ID, right_plate_ID, return_rate_of_motion)


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
        filename = str(filename)

        numpy_types = {'.csv': ',', '.txt': ' ', '.dat': ' '}

        if filename.endswith(tuple(numpy_types.keys())):
            data = np.c_[self.lons, self.lats, self.plate_id]

            # find appropriate delimiter
            delimiter = numpy_types[filename[filename.index('.'):]]

            header = "Longitude{0}Latitude{0}Plate_ID".format(delimiter)
            np.savetxt(filename, data, header=header, delimiter=delimiter, comments='')

        elif filename.endswith('.gpml'):
            self.FeatureCollection.write(filename)

        else:
            raise ValueError("Cannot save to specified file type, use csv or gpml file extension.")