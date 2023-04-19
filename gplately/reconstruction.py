"""Tools that wrap up pyGplates and Plate Tectonic Tools functionalities for reconstructing features,
 working with point data, and calculating plate velocities at specific geological times. 
"""
import math
import os
import warnings

import numpy as np
import ptt
import pygplates

from . import tools as _tools
from .gpml import _load_FeatureCollection
from .pygplates import RotationModel as _RotationModel
from .pygplates import FeatureCollection as _FeatureCollection


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
    
    def __init__(self, rotation_model, topology_features=None, static_polygons=None):

        if hasattr(rotation_model, "reconstruction_identifier"):
            self.name = rotation_model.reconstruction_identifier
        else:
            self.name = None

        self.rotation_model = _RotationModel(rotation_model)
        self.topology_features = _load_FeatureCollection(topology_features)
        self.static_polygons = _load_FeatureCollection(static_polygons)

    def __getstate__(self):

        filenames = {"rotation_model": self.rotation_model.filenames}

        if self.topology_features:
            filenames['topology_features'] = self.topology_features.filenames
        if self.static_polygons:
            filenames['static_polygons'] = self.static_polygons.filenames

        # # remove unpicklable items
        # del self.rotation_model, self.topology_features, self.static_polygons
        
        # # really make sure they're gone
        # self.rotation_model = None
        # self.topology_features = None
        # self.static_polygons = None

        return filenames


    def __setstate__(self, state):

        # reinstate unpicklable items
        self.rotation_model = _RotationModel(state['rotation_model'])

        self.topology_features = None
        self.static_polygons = None

        if 'topology_features' in state:
            self.topology_features = _FeatureCollection()
            for topology in state['topology_features']:
                self.topology_features.add( _FeatureCollection(topology) )

        if 'static_polygons' in state:
            self.static_polygons = _FeatureCollection()
            for polygon in state['static_polygons']:
                self.static_polygons.add( _FeatureCollection(polygon) )


    def tessellate_subduction_zones(self, time, tessellation_threshold_radians=0.001, ignore_warnings=False, return_geodataframe=False, **kwargs):
        """Samples points along subduction zone trenches and obtains subduction data at a particular
        geological time.
        
        Resolves topologies at `time`, tessellates all resolved subducting features to within 'tessellation_threshold_radians'
        radians and obtains the following information for each sampled point along a trench:
    
        `tessellate_subduction_zones` returns a list of 10 vertically-stacked tuples with the following data per sampled trench point:

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

        if return_geodataframe:
            from shapely import geometry
            import geopandas as gpd
            coords = [geometry.Point(lon,lat) for lon,lat in zip(subduction_data[:,0], subduction_data[:,1])]
            d = {"geometry": coords}

            labels = ['convergence velocity (cm/yr)',
                      'convergence obliquity angle (degrees)',
                      'trench velocity (cm/yr)',
                      'trench obliquity angle (degrees)',
                      'length (degrees)',
                      'trench normal angle (degrees)',
                      'subducting plate ID',
                      'overriding plate ID']

            for i, label in enumerate(labels):
                index = 2+i
                d[label] = subduction_data[:,index]

            gdf = gpd.GeoDataFrame(d, geometry="geometry")
            return gdf

        else:
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
                subduction_data = self.tessellate_subduction_zones(time, ignore_warnings=ignore_warnings)

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


    def total_continental_arc_length(
        self,
        time,
        continental_grid,
        trench_arc_distance,
        ignore_warnings=True,
    ):
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
        continental_grid: Raster, array_like, or str
            The continental grid used to identify continental arc points. Must
            be convertible to `Raster`. For an array, a global extent is
            assumed [-180,180,-90,90]. For a filename, the extent is obtained
            from the file.
        trench_arc_distance : float
            The trench-to-arc distance (in kilometres) to project sampled trench points out by in the direction of their 
            subduction polarities. 
        ignore_warnings : bool, default=True
            Choose whether to ignore warning messages from PTT's subduction_convergence workflow that alerts the user of 
            subduction sub-segments that are ignored due to unidentified polarities and/or subducting plates. 

        Returns
        -------
        total_continental_arc_length_kms : float
            The continental arc length (in km) at the specified time.
        """
        from . import grids as _grids

        if isinstance(continental_grid, _grids.Raster):
            graster = continental_grid
        elif isinstance(continental_grid, str):
            # Process the continental grid directory
            graster = _grids.Raster(
                data=continental_grid,
                realign=True,
                time=float(time),
            )
        else:
            # Process the masked continental grid
            try:
                continental_grid = np.array(continental_grid)
                graster = _grids.Raster(
                    data=continental_grid,
                    extent=[-180,180,-90,90],
                    time=float(time),
                )
            except Exception as e:
                raise TypeError(
                    "Invalid type for `continental_grid` (must be Raster,"
                    + " str, or array_like)"
                ) from e
        if (time != graster.time) and (not ignore_warnings):
            raise RuntimeWarning(
                "`continental_grid.time` ({}) ".format(graster.time)
                + "does not match `time` ({})".format(time)
            )

        # Obtain trench data with Plate Tectonic Tools
        trench_data = self.tessellate_subduction_zones(time, ignore_warnings=ignore_warnings)

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
        sampled_points = graster.interpolate(
            ilon,
            ilat,
            method='linear',
            return_indices=False,
        )
        continental_indices = np.where(sampled_points > 0)
        point_lats = ilat[continental_indices]
        point_radii = _tools.geocentric_radius(point_lats) * 1.0e-3  # km
        segment_arclens = np.deg2rad(trench_arcseg[continental_indices])
        segment_lengths = point_radii * segment_arclens
        return np.sum(segment_lengths)


    def tessellate_mid_ocean_ridges(self, time, tessellation_threshold_radians=0.001, ignore_warnings=False, return_geodataframe=False, **kwargs):
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

        if return_geodataframe:
            from shapely import geometry
            import geopandas as gpd

            points = [geometry.Point(lon,lat) for lon,lat in zip(ridge_data[:,0], ridge_data[:,1])]
            gdf = gpd.GeoDataFrame({"geometry": points,
                                    "velocity (cm/yr)": ridge_data[:,2],
                                    "length (degrees)": ridge_data[:,3]},
                                    geometry="geometry")
            return gdf

        else:
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
                ridge_data = self.tessellate_mid_ocean_ridges(time)

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
                StepRates[:,i] = StepRate*0.1 # cm/yr

                # Obseleted by Lauren's changes above (though it is more efficient)
                # multiply arc length of the motion path segment by a latitude-dependent Earth radius
                # use latitude of the segment start point
                # distance.append( segment.get_arc_length() * _tools.geocentric_radius(segment.get_start_point().to_lat_lon()[0]) / 1e3)
                # rate = np.asarray(distance)/np.diff(time_array)
                # rates[:,i] = np.flipud(rate)
                # rates *= 0.1 # cm/yr
        
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
                StepRates[:,i] = StepRate*0.1 # cm/yr

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
    plate_reconstruction : object pointer
        Allows for the accessibility of `PlateReconstruction` object attributes: `rotation_model`, `topology_featues` 
        and `static_polygons` for use in the `Points` object if called using “self.plate_reconstruction.X”, 
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
    def __init__(self, plate_reconstruction, lons, lats, time=0, plate_id=None):
        self.lons = lons
        self.lats = lats
        self.time = time
        self.attributes = dict()

        self.plate_reconstruction = plate_reconstruction

        self._update(lons, lats, time, plate_id)


    def _update(self, lons, lats, time=0, plate_id=None):

        # get Cartesian coordinates
        self.x, self.y, self.z = _tools.lonlat2xyz(lons, lats, degrees=False)

        # scale by average radius of the Earth
        self.x *= _tools.EARTH_RADIUS
        self.y *= _tools.EARTH_RADIUS
        self.z *= _tools.EARTH_RADIUS

        # store concatenated arrays
        self.lonlat = np.c_[self.lons, self.lats]
        self.xyz = np.c_[self.x, self.y, self.z]


        rotation_model = self.plate_reconstruction.rotation_model
        static_polygons = self.plate_reconstruction.static_polygons
        

        features = _tools.points_to_features(lons, lats, plate_id)

        if plate_id is not None:
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

    @property
    def size(self):
        """ Number of points
        """
        return len(self.lons)


    def __getstate__(self):

        filenames = self.plate_reconstruction.__getstate__()

        # add important variables from Points object
        filenames["lons"] = self.lons
        filenames["lats"] = self.lats
        filenames['time'] = self.time
        filenames['plate_id'] = self.plate_id
        for key in self.attributes:
            filenames[key] = self.attributes[key]

        return filenames

    def __setstate__(self, state):

        self.plate_reconstruction = PlateReconstruction(state['rotation_model'], state['topology_features'], state['static_polygons'])

        # reinstate unpicklable items
        self.lons = state['lons']
        self.lats = state['lats']
        self.time = state['time']
        self.plate_id = state['plate_id']
        self.attributes = dict()

        self._update(self.lons, self.lats, self.time, self.plate_id)

        for key in state:
            if key not in ['lons', 'lats', 'time', 'plate_id']:
                self.attributes[key] = state[key]

    def copy(self):
        """ Returns a copy of the Points object

        Returns
        -------
        Points
            A copy of the current Points object
        """
        gpts = Points(self.plate_reconstruction, self.lons.copy(), self.lats.copy(), self.time, self.plate_id.copy())
        gpts.add_attributes(**self.attributes.copy())

    def add_attributes(self, **kwargs):
        """ Adds the value of a feature attribute associated with a key.

        Example
        -------

            # Define latitudes and longitudes to set up a Points object
            pt_lons = np.array([140., 150., 160.])
            pt_lats = np.array([-30., -40., -50.])

            gpts = gplately.Points(model, pt_lons, pt_lats)

            # Add the attributes a, b and c to the points in the Points object
            gpts.add_attributes(
                a=[10,2,2], 
                b=[2,3,3], 
                c=[30,0,0],
            )

            print(gpts.attributes)

        The output would be:

            {'a': [10, 2, 2], 'b': [2, 3, 3], 'c': [30, 0, 0]}

        Parameters
        ----------
        **kwargs : sequence of key=item/s
            A single key=value pair, or a sequence of key=value pairs denoting the name and
            value of an attribute.


        Notes
        -----
        * An **assertion** is raised if the number of points in the Points object is not equal
        to the number of values associated with an attribute key. For example, consider an instance
        of the Points object with 3 points. If the points are ascribed an attribute `temperature`,
        there must be one `temperature` value per point, i.e. `temperature = [20, 15, 17.5]`.

        """
        keys = kwargs.keys()

        for key in kwargs:
            attribute = kwargs[key]

            # make sure attribute is the same size as self.lons
            if type(attribute) is int or type(attribute) is float:
                array = np.full(self.lons.size, attribute)
                attribute = array
            elif isinstance(attribute, np.ndarray):
                if attribute.size == 1:
                    array = np.full(self.lons.size, attribute, dtype=attribute.dtype)
                    attribute = array

            assert len(attribute) == self.lons.size, "Size mismatch, ensure attributes have the same number of entries as Points"
            self.attributes[key] = attribute

        if any(kwargs):
            # add these to the FeatureCollection
            for f, feature in enumerate(self.FeatureCollection):
                for key in keys:
                    # extract value for each row in attribute
                    val = self.attributes[key][f]

                    # set this attribute on the feature
                    feature.set_shapefile_attribute(key, val)

    def get_geopandas_dataframe(self):
        """ Adds a shapely point `geometry` attribute to each point in the `gplately.Points` object.
        pandas.DataFrame that has a column with geometry
        Any existing point attributes are kept.

        Returns
        -------
        GeoDataFrame : instance of `geopandas.GeoDataFrame`
            A pandas.DataFrame with rows equal to the number of points in the `gplately.Points` object,
            and an additional column containing a shapely `geometry` attribute.

        Example
        -------

            pt_lons = np.array([140., 150., 160.])
            pt_lats = np.array([-30., -40., -50.])

            gpts = gplately.Points(model, pt_lons, pt_lats)

            # Add sample attributes a, b and c to the points in the Points object
            gpts.add_attributes(
                a=[10,2,2], 
                b=[2,3,3], 
                c=[30,0,0],
            )

            gpts.get_geopandas_dataframe()

        ...has the output:

                a  b   c                     geometry
            0  10  2  30  POINT (140.00000 -30.00000)
            1   2  3   0  POINT (150.00000 -40.00000)
            2   2  3   0  POINT (160.00000 -50.00000)


        """
        import geopandas as gpd
        from shapely import geometry

        # create shapely points
        points = []
        for lon, lat in zip(self.lons, self.lats):
            points.append( geometry.Point(lon, lat) )

        attributes = self.attributes.copy()
        attributes['geometry'] = points

        return gpd.GeoDataFrame(attributes, geometry='geometry')

    def get_geodataframe(self):
        """ Returns the output of `Points.get_geopandas_dataframe()`.

        Adds a shapely point `geometry` attribute to each point in the `gplately.Points` object.
        pandas.DataFrame that has a column with geometry
        Any existing point attributes are kept.

        Returns
        -------
        GeoDataFrame : instance of `geopandas.GeoDataFrame`
            A pandas.DataFrame with rows equal to the number of points in the `gplately.Points` object,
            and an additional column containing a shapely `geometry` attribute.

        Example
        -------

            pt_lons = np.array([140., 150., 160.])
            pt_lats = np.array([-30., -40., -50.])

            gpts = gplately.Points(model, pt_lons, pt_lats)

            # Add sample attributes a, b and c to the points in the Points object
            gpts.add_attributes(
                a=[10,2,2], 
                b=[2,3,3], 
                c=[30,0,0],
            )

            gpts.get_geopandas_dataframe()

        ...has the output:

                a  b   c                     geometry
            0  10  2  30  POINT (140.00000 -30.00000)
            1   2  3   0  POINT (150.00000 -40.00000)
            2   2  3   0  POINT (160.00000 -50.00000)


        """
        return self.get_geopandas_dataframe()

    def reconstruct(self, time, anchor_plate_id=0, return_array=False, **kwargs):
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

        return_array : bool, default False
            Return a `numpy.ndarray`, rather than a `Points` object.

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
        reconstructed_features = self.plate_reconstruction.reconstruct(
            self.features, to_time, from_time, anchor_plate_id=anchor_plate_id, **kwargs)

        rlons, rlats = _tools.extract_feature_lonlat(reconstructed_features)
        if return_array:
            return rlons, rlats
        else:
            gpts = Points(self.plate_reconstruction, rlons, rlats, time=to_time, plate_id=self.plate_id)
            gpts.add_attributes(**self.attributes.copy())
            return gpts


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

            reconstructed_features = self.plate_reconstruction.reconstruct(
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

        rotation_model = self.plate_reconstruction.rotation_model
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

            reconstructed_motion_paths = self.plate_reconstruction.reconstruct(
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
                StepRates[:,i] = StepRate*0.1 # cm/yr
        
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
        model = self.plate_reconstruction
        return model.create_flowline(self.lons, self.lats, time_array, left_plate_ID, right_plate_ID, return_rate_of_motion)


    def _get_dataframe(self):
        import geopandas as gpd

        data = dict()
        data["Longitude"] = self.lons
        data["Latitude"]  = self.lats
        data["Plate_ID"]  = self.plate_id
        for key in self.attributes:
            data[key] = self.attributes[key]

        return gpd.GeoDataFrame(data)

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

        if filename.endswith(('.csv', '.txt', '.dat')):
            df = self._get_dataframe()
            df.to_csv(filename, index=False)

        elif filename.endswith(('.xls', '.xlsx')):
            df = self._get_dataframe()
            df.to_excel(filename, index=False)

        elif filename.endswith('xml'):
            df = self._get_dataframe()
            df.to_xml(filename, index=False)

        elif filename.endswith('.gpml') or filename.endswith('.gpmlz'):
            self.FeatureCollection.write(filename)

        else:
            raise ValueError("Cannot save to specified file type. Use csv, gpml, or xls file extension.")


# FROM RECONSTRUCT_BY_TOPOLOGIES.PY
class _DefaultCollision(object):
    """
    Default collision detection function class (the function is the '__call__' method).
    """

    DEFAULT_GLOBAL_COLLISION_PARAMETERS = (7.0, 10.0)
    """
    Default collision parameters for all feature types.

    This is a 2-tuple of (threshold velocity delta in kms/my, threshold distance to boundary per My in kms/my):
    Here we default to the same constants used internally in GPlates 2.0 (ie, 7.0 and 10.0).
    """

    def __init__(
            self,
            global_collision_parameters = DEFAULT_GLOBAL_COLLISION_PARAMETERS,
            feature_specific_collision_parameters = None):
        """
        global_collision_parameters: The collision parameters to use for any feature type not specified in 'feature_specific_collision_parameters'.
                                     Should be a 2-tuple of (threshold velocity delta in kms/my, threshold distance to boundary per My in kms/my).
                                     The first threshold parameter means:
                                        A point that transitions from one plate to another can disappear if the change in velocity exceeds this threshold.
                                     The second threshold parameter means:
                                        Only those transitioning points exceeding the threshold velocity delta and that are close enough to a plate boundary can disappear.
                                        The distance is proportional to the relative velocity (change in velocity), plus a constant offset based on the threshold distance to boundary
                                        to account for plate boundaries that change shape significantly from one time step to the next
                                        (note that some boundaries are meant to do this and others are a result of digitisation).
                                        The actual distance threshold used is (threshold_distance_to_boundary + relative_velocity) * time_interval
                                     Defaults to parameters used in GPlates 2.0, if not specified.
        
        feature_specific_collision_parameters: Optional sequence of collision parameters specific to feature types.
                                               If specified then should be a sequence of 2-tuples, with each 2-tuple specifying (feature_type, collision_parameters).
                                               And where each 'collision_parameters' is a 2-tuple of (threshold velocity delta in kms/my, threshold distance to boundary per My in kms/my).
                                                   See 'global_collision_parameters' for details on these thresholds.
                                               Any feature type not specified here defaults to using 'global_collision_parameters'.
        """

        # Convert list of (feature_type, collision_parameters) tuples to a dictionary.
        if feature_specific_collision_parameters:
            self.feature_specific_collision_parameters = dict(feature_specific_collision_parameters)
        else:
            self.feature_specific_collision_parameters = dict()
        # Fallback for any feature type not specified in the optional feature-specific list.
        self.global_collision_parameters = global_collision_parameters

        # Used to improve performance by caching velocity stage rotations in a dict (for a specific reconstruction time).
        self.velocity_stage_rotation_dict = {}
        self.velocity_stage_rotation_time = None


    def __call__(
            self,
            rotation_model,
            time,
            reconstruction_time_interval,
            prev_point,
            curr_point,
            prev_topology_plate_id,
            prev_resolved_plate_boundary,
            curr_topology_plate_id,
            curr_resolved_plate_boundary):
        """
        Returns True if a collision was detected.

        If transitioning from a rigid plate to another rigid plate with a different plate ID then
        calculate the difference in velocities and continue testing as follows
        (otherwise, if there's no transition, then the point is still active)...
        
        If the velocity difference is below a threshold then we assume the previous plate was split,
        or two plates joined. In this case the point has not subducted (forward in time) or
        been consumed by a mid-ocean (backward in time) and hence is still active.
        
        If the velocity difference is large enough then we see if the distance of the *previous* position
        to the polygon boundary (of rigid plate containing it) exceeds a threshold.
        If the distance exceeds the threshold then the point is far enough away from the boundary that it
        cannot be subducted or consumed by it and hence the point is still active.
        However if the point is close enough then we assume the point was subducted/consumed
        (remember that the point switched plate IDs).
        Also note that the threshold distance increases according to the velocity difference to account for fast
        moving points (that would otherwise tunnel through the boundary and accrete onto the other plate).
        The reason for testing the distance from the *previous* point, and not from the *current* point, is:
        
          (i)  A topological boundary may *appear* near the current point (such as a plate split at the current time)
               and we don't want that split to consume the current point regardless of the velocity difference.
               It won't get consumed because the *previous* point was not near a boundary (because before split happened).
               If the velocity difference is large enough then it might cause the current point to transition to the
               adjacent split plate in the *next* time step (and that's when it should get consumed, not in the current time step).
               An example of this is a mid-ocean ridge suddenly appearing (going forward in time).
        
          (ii) A topological boundary may *disappear* near the current point (such as a plate merge at the current time)
               and we want that merge to consume the current point if the velocity difference is large enough.
               In this case the *previous* point is near a boundary (because before plate merged) and hence can be
               consumed (provided velocity difference is large enough). And since the boundary existed in the previous
               time step, it will affect position of the current point (and whether it gets consumed or not).
               An example of this is a mid-ocean ridge suddenly disappearing (going backward in time).
        
        ...note that items (i) and (ii) above apply both going forward and backward in time.
        """

        # See if a collision occurred.
        if (curr_topology_plate_id != prev_topology_plate_id and
            prev_topology_plate_id is not None and
            curr_topology_plate_id is not None):

            #
            # Speed up by caching velocity stage rotations in a dict.
            #
            if time != self.velocity_stage_rotation_time:
                # We've just switched to a new time so clear the cache.
                #
                # We only cache stage rotations for a specific time.
                # We only really need to cache different plate IDs at the same 'time', so this avoids caching for all times
                # (which would also require including 'time' in the key) and using memory unnecessarily.
                self.velocity_stage_rotation_dict.clear()
                self.velocity_stage_rotation_time = time
            prev_location_velocity_stage_rotation = self.velocity_stage_rotation_dict.get(prev_topology_plate_id)
            if not prev_location_velocity_stage_rotation:
                prev_location_velocity_stage_rotation = rotation_model.get_rotation(time + 1, prev_topology_plate_id, time)
                self.velocity_stage_rotation_dict[prev_topology_plate_id] = prev_location_velocity_stage_rotation
            curr_location_velocity_stage_rotation = self.velocity_stage_rotation_dict.get(curr_topology_plate_id)
            if not curr_location_velocity_stage_rotation:
                curr_location_velocity_stage_rotation = rotation_model.get_rotation(time + 1, curr_topology_plate_id, time)
                self.velocity_stage_rotation_dict[curr_topology_plate_id] = curr_location_velocity_stage_rotation

            # Note that even though the current point is not inside the previous boundary (because different plate ID), we can still
            # calculate a velocity using its plate ID (because we really should use the same point in our velocity comparison).
            prev_location_velocity = pygplates.calculate_velocities(
                    (curr_point,), prev_location_velocity_stage_rotation, 1, pygplates.VelocityUnits.kms_per_my)[0]
            curr_location_velocity = pygplates.calculate_velocities(
                    (curr_point,), curr_location_velocity_stage_rotation, 1, pygplates.VelocityUnits.kms_per_my)[0]

            delta_velocity = curr_location_velocity - prev_location_velocity
            delta_velocity_magnitude = delta_velocity.get_magnitude()
            
            # If we have feature-specific collision parameters then iterate over the boundary sub-segments of the *previous* topological boundary
            # and test proximity to each sub-segment individually (with sub-segment feature type specific collision parameters).
            # Otherwise just test proximity to the entire boundary polygon using the global collision parameters.
            if self.feature_specific_collision_parameters:
                for prev_boundary_sub_segment in prev_resolved_plate_boundary.get_boundary_sub_segments():
                    # Use feature-specific collision parameters if found (falling back to global collision parameters).
                    threshold_velocity_delta, threshold_distance_to_boundary_per_my = self.feature_specific_collision_parameters.get(
                            prev_boundary_sub_segment.get_feature().get_feature_type(),
                            # Default to global collision parameters if no collision parameters specified for sub-segment's feature type...
                            self.global_collision_parameters)
                    
                    # Since each feature type could use different collision parameters we must use the current boundary sub-segment instead of the boundary polygon.
                    if self._detect_collision_using_collision_parameters(
                            reconstruction_time_interval,
                            delta_velocity_magnitude,
                            prev_point,
                            prev_boundary_sub_segment.get_resolved_geometry(),
                            threshold_velocity_delta,
                            threshold_distance_to_boundary_per_my):
                        # Detected a collision.
                        return True
            else:
                # No feature-specific collision parameters so use global fallback.
                threshold_velocity_delta, threshold_distance_to_boundary_per_my = self.global_collision_parameters
                
                # Since all feature types use the same collision parameters we can use the boundary polygon instead of iterating over its sub-segments.
                if self._detect_collision_using_collision_parameters(
                        reconstruction_time_interval,
                        delta_velocity_magnitude,
                        prev_point,
                        prev_resolved_plate_boundary.get_resolved_boundary(),
                        threshold_velocity_delta,
                        threshold_distance_to_boundary_per_my):
                    # Detected a collision.
                    return True
        
        return False
        
        
    def _detect_collision_using_collision_parameters(
            self,
            reconstruction_time_interval,
            delta_velocity_magnitude,
            prev_point,
            prev_boundary_geometry,
            threshold_velocity_delta,
            threshold_distance_to_boundary_per_my):
        
        if delta_velocity_magnitude > threshold_velocity_delta:
            # Add the minimum distance threshold to the delta velocity threshold.
            # The delta velocity threshold only allows those points that are close enough to the boundary to reach
            # it given their current relative velocity.
            # The minimum distance threshold accounts for sudden changes in the shape of a plate boundary
            # which are no supposed to represent a new or shifted boundary but are just a result of the topology
            # builder/user digitising a new boundary line that differs noticeably from that of the previous time period.
            distance_threshold_radians = ((threshold_distance_to_boundary_per_my + delta_velocity_magnitude)
                * reconstruction_time_interval / pygplates.Earth.equatorial_radius_in_kms)

            if pygplates.GeometryOnSphere.distance(prev_point, prev_boundary_geometry, distance_threshold_radians) is not None:
                # Detected a collision.
                return True
        
        return False


_DEFAULT_COLLISION = _DefaultCollision()


class _ContinentCollision(object):
    """
    Continental collision detection function class (the function is the '__call__' method).
    """

    def __init__(
        self,
        grd_output_dir,
        chain_collision_detection=_DEFAULT_COLLISION,
        verbose=False,
    ):
        """
        grd_output_dir: The directory containing the continental grids.

        chain_collision_detection: Another collision detection class/function to reference if we find no collision.
                                   If None then no collision detection is chained. Defaults to the default collision detection.

        verbose: Print progress messages
        """
        
        self.grd_output_dir = grd_output_dir
        self.chain_collision_detection = chain_collision_detection
        self.verbose = verbose
        
        # Load a new grid each time the reconstruction time changes.
        self.grid_time = None

    @property
    def grid_time(self):
        return self._grid_time

    @grid_time.setter
    def grid_time(self, time):
        from .grids import read_netcdf_grid

        if time is None:
            self._grid_time = time
        else:
            filename = '{:s}'.format(self.grd_output_dir.format(time))
            if self.verbose:
                print(
                    'Points masked against grid: {0}'.format(
                        os.path.basename(filename)
                    )
                )
            gridZ, gridX, gridY = read_netcdf_grid(filename, return_grids=True)
            self.gridZ = gridZ
            self.ni, self.nj = gridZ.shape
            self.xmin = np.nanmin(gridX)
            self.xmax = np.nanmax(gridX)
            self.ymin = np.nanmin(gridY)
            self.ymax = np.nanmax(gridY)
            self._grid_time = float(time)


    def __call__(
            self,
            rotation_model,
            time,
            reconstruction_time_interval,
            prev_point,
            curr_point,
            prev_topology_plate_id,
            prev_resolved_plate_boundary,
            curr_topology_plate_id,
            curr_resolved_plate_boundary):
        """
        Returns True if a collision with a continent was detected, or returns result of
        chained collision detection if 'self.chain_collision_detection' is not None.
        """
        # Load the grid for the current time if encountering a new time.
        if time != self.grid_time:
            self.grid_time = time
            self.continent_deletion_count = 0

        # Sample mask grid, which is one over continents and zero over oceans.
        point_lat, point_lon = curr_point.to_lat_lon()
        point_i = (self.ni - 1) * (
            (point_lat - self.ymin)
            / (self.ymax - self.ymin)
        )
        point_j = (self.nj - 1) * (
            (point_lon - self.xmin)
            / (self.xmax - self.xmin)
        )
        point_i_uint = np.rint(point_i).astype(np.uint)
        point_j_uint = np.rint(point_j).astype(np.uint)
        try:
            mask_value = self.gridZ[point_i_uint, point_j_uint]
        except IndexError:
            point_i = np.clip(np.rint(point_i), 0, self.ni - 1).astype(np.int_)
            point_j = np.clip(np.rint(point_j), 0, self.nj - 1).astype(np.int_)
            mask_value = self.gridZ[point_i, point_j]
        if mask_value >= 0.5:
            # Detected a collision.
            self.continent_deletion_count += 1
            return True

        # We didn't find a collision, so ask the chained collision detection if it did (if we have anything chained).
        if self.chain_collision_detection:
            return self.chain_collision_detection(
                    rotation_model,
                    time,
                    reconstruction_time_interval,
                    prev_point,
                    curr_point,
                    prev_topology_plate_id,
                    prev_resolved_plate_boundary,
                    curr_topology_plate_id,
                    curr_resolved_plate_boundary)

        return False


def  reconstruct_points(
            rotation_features_or_model,
            topology_features,
            reconstruction_begin_time,
            reconstruction_end_time,
            reconstruction_time_interval,
            points,
            point_begin_times = None,
            point_end_times = None,
            point_plate_ids = None,
            detect_collisions = _DEFAULT_COLLISION):
    """
    Function to reconstruct points using the ReconstructByTopologies class below.
    
    For description of parameters see the ReconstructByTopologies class below.
    """
    
    topology_reconstruction = _ReconstructByTopologies(
            rotation_features_or_model,
            topology_features,
            reconstruction_begin_time,
            reconstruction_end_time,
            reconstruction_time_interval,
            points,
            point_begin_times,
            point_end_times,
            point_plate_ids,
            detect_collisions)
    
    return topology_reconstruction.reconstruct()


class _ReconstructByTopologies(object):
    """
    Class to reconstruct geometries using topologies.
    
    Currently only points are supported.
    
    use_plate_partitioner: If True then use pygplates.PlatePartitioner to partition points,
                           otherwise use faster points_in_polygons.find_polygons().
    """
    use_plate_partitioner = False
    
    
    def __init__(
            self,
            rotation_features_or_model,
            topology_features,
            reconstruction_begin_time,
            reconstruction_end_time,
            reconstruction_time_interval,
            points,
            point_begin_times = None,
            point_end_times = None,
            point_plate_ids = None,
            detect_collisions = _DEFAULT_COLLISION):
        """
        rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).
        
        topology_features: Topology feature collection(s), or list of features, or filename(s) or any combination of those.
        
        detect_collisions: Collision detection function, or None. Defaults to _DEFAULT_COLLISION.
        """
        
        # Turn rotation data into a RotationModel (if not already).
        if not isinstance(rotation_features_or_model, pygplates.RotationModel):
            rotation_model = pygplates.RotationModel(rotation_features_or_model)
        else:
            rotation_model = rotation_features_or_model
        self.rotation_model = rotation_model
        
        # Turn topology data into a list of features (if not already).
        self.topology_features = pygplates.FeaturesFunctionArgument(topology_features).get_features()
        
        # Set up an array of reconstruction times covering the reconstruction time span.
        self.reconstruction_begin_time = reconstruction_begin_time
        self.reconstruction_end_time = reconstruction_end_time
        if reconstruction_time_interval <= 0.0:
            raise ValueError("'reconstruction_time_interval' must be positive.")
        # Reconstruction can go forward or backward in time.
        if self.reconstruction_begin_time > self.reconstruction_end_time:
            self.reconstruction_time_step = -reconstruction_time_interval
        else:
            self.reconstruction_time_step = reconstruction_time_interval
        # Get number of times including end time if time span is a multiple of time step.
        # The '1' is because, for example, 2 time intervals is 3 times.
        # The '1e-6' deals with limited floating-point precision, eg, we want (3.0 - 0.0) / 1.0 to be 3.0 and not 2.999999 (which gets truncated to 2).
        self.num_times = 1 + int(math.floor(1e-6 + float(self.reconstruction_end_time - self.reconstruction_begin_time) / self.reconstruction_time_step))
        # It's possible the time step is larger than the time span, in which case we change it to equal the time span.
        # This guarantees there'll be at least one time step (which has two times; one at either end of interval).
        if self.num_times == 1:
            self.num_times = 2
            self.reconstruction_time_step = self.reconstruction_end_time - self.reconstruction_begin_time
        self.reconstruction_time_interval = math.fabs(self.reconstruction_time_step)
            
        self.last_time_index = self.num_times - 1
        
        self.points = points
        self.num_points = len(points)
        
        # Use the specified point begin times if provided (otherwise use 'inf').
        self.point_begin_times = point_begin_times
        if self.point_begin_times is None:
            self.point_begin_times = [float('inf')] * self.num_points
        elif len(self.point_begin_times) != self.num_points:
            raise ValueError("Length of 'point_begin_times' must match length of 'points'.")
        
        # Use the specified point end times if provided (otherwise use '-inf').
        self.point_end_times = point_end_times
        if self.point_end_times is None:
            self.point_end_times = [float('-inf')] * self.num_points
        elif len(self.point_end_times) != self.num_points:
            raise ValueError("Length of 'point_end_times' must match length of 'points'.")
        
        # Use the specified point plate IDs if provided (otherwise use '0').
        # These plate IDs are only used when a point falls outside all resolved topologies during a time step.
        self.point_plate_ids = point_plate_ids
        if self.point_plate_ids is None:
            self.point_plate_ids = [0] * self.num_points
        elif len(self.point_plate_ids) != self.num_points:
            raise ValueError("Length of 'point_plate_ids' must match length of 'points'.")
        
        self.detect_collisions = detect_collisions
    
    
    def reconstruct(self):
        
        # Initialise the reconstruction.
        self.begin_reconstruction()
        
        # Loop over the reconstruction times until reached end of the reconstruction time span, or
        # all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        while self.reconstruct_to_next_time():
            pass
        
        return self.get_active_current_points()
    
    
    def begin_reconstruction(self):
        
        self.current_time_index = 0
        
        # Set up point arrays.
        # Store active and inactive points here (inactive points have None in corresponding entries).
        self.prev_points = [None] * self.num_points
        self.curr_points = [None] * self.num_points
        self.next_points = [None] * self.num_points
        
        # Each point can only get activated once (after deactivation it cannot be reactivated).
        self.point_has_been_activated = [False] * self.num_points
        self.num_activated_points = 0
        
        # Set up topology arrays (corresponding to active/inactive points at same indices).
        self.prev_topology_plate_ids = [None] * self.num_points
        self.curr_topology_plate_ids = [None] * self.num_points
        self.prev_resolved_plate_boundaries = [None] * self.num_points
        self.curr_resolved_plate_boundaries = [None] * self.num_points

        # Array to store indices of points found in continents
        self.in_continent_indices = [None] * self.num_points
        self.in_continent_points = [None] * self.num_points

        self.deletedpoints = []
        
        self._activate_deactivate_points()
        self._find_resolved_topologies_containing_points()
    
    
    def get_current_time(self):
        
        return self.reconstruction_begin_time + self.current_time_index * self.reconstruction_time_step
    
    
    def get_all_current_points(self):
        
        return self.curr_points
    
    
    def get_active_current_points(self):
        
        # Return only the active points (the ones that are not None).
        return [point for point in self.get_all_current_points() if point is not None]


    def get_in_continent_indices(self):

        return self.in_continent_points, self.in_continent_indices
    
    
    def reconstruct_to_next_time(self):
        
        # If we're at the last time then there is no next time to reconstruct to.
        if self.current_time_index == self.last_time_index:
            return False
        
        # If all points have been previously activated, but none are currently active then we're finished.
        # This means all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        if (self.num_activated_points == self.num_points and
            not any(self.curr_points)):
            return False
        
        # Cache stage rotations by plate ID.
        # Use different dicts since using different rotation models and time steps, etc.
        reconstruct_stage_rotation_dict = {}
        
        current_time = self.get_current_time()
        
        # Iterate over all points to reconstruct them to the next time step.
        for point_index in range(self.num_points):
            
            curr_point = self.curr_points[point_index]
            if curr_point is None:
                # Current point is not currently active.
                # So we cannot reconstruct to next time.
                self.next_points[point_index] = None
                continue
            
            # Get plate ID of resolved topology containing current point
            # (this was determined in last call to '_find_resolved_topologies_containing_points()').
            curr_plate_id = self.curr_topology_plate_ids[point_index]
            if curr_plate_id is None:
                # Current point is currently active but it fell outside all resolved polygons.
                # So instead we just reconstruct using its plate ID (that was manually assigned by the user/caller).
                curr_plate_id = self.point_plate_ids[point_index]
                continue
            
            # Get the stage rotation that will move the point from where it is at the current time to its
            # location at the next time step, based on the plate id that contains the point at the current time.
            
            # Speed up by caching stage rotations in a dict.
            stage_rotation = reconstruct_stage_rotation_dict.get(curr_plate_id)
            if not stage_rotation:
                stage_rotation = self.rotation_model.get_rotation(
                        # Positive/negative time step means reconstructing backward/forward in time.
                        current_time + self.reconstruction_time_step,
                        curr_plate_id,
                        current_time)
                reconstruct_stage_rotation_dict[curr_plate_id] = stage_rotation
            
            # Use the stage rotation to reconstruct the tracked point from position at current time 
            # to position at the next time step.
            self.next_points[point_index] = stage_rotation * curr_point
        
        #
        # Set up for next loop iteration.
        #
        # Rotate previous, current and next point arrays.
        # The new previous will be the old current.
        # The new current will be the old next.
        # The new next will be the old previous (but values are ignored and overridden in next time step; just re-using its memory).
        self.prev_points, self.curr_points, self.next_points = self.curr_points, self.next_points, self.prev_points
        # Swap previous and current topology arrays.
        # The new previous will be the old current.
        # The new current will be the old previous (but values are ignored and overridden in next time step; just re-using its memory).
        self.prev_topology_plate_ids, self.curr_topology_plate_ids = self.curr_topology_plate_ids, self.prev_topology_plate_ids
        self.prev_resolved_plate_boundaries, self.curr_resolved_plate_boundaries = self.curr_resolved_plate_boundaries, self.prev_resolved_plate_boundaries
        
        # Move the current time to the next time.
        self.current_time_index += 1
        current_time = self.get_current_time()
        
        self._activate_deactivate_points()
        self._find_resolved_topologies_containing_points()
        
        # Iterate over all points to detect collisions.
        if self.detect_collisions:
            for point_index in range(self.num_points):
                
                prev_point = self.prev_points[point_index]
                curr_point = self.curr_points[point_index]
                if prev_point is None or curr_point is None:
                    # If current point is not currently active then no need to detect a collision for it (to deactivate it).
                    # Also previous point might just have been activated now, at end of current time step, and hence
                    # not active at beginning of time step.
                    continue
                
                # Get plate IDs of resolved topology containing previous and current point
                # (this was determined in last call to '_find_resolved_topologies_containing_points()').
                #
                # Note that could be None, so the collision detection needs to handle that.
                prev_plate_id = self.prev_topology_plate_ids[point_index]
                curr_plate_id = self.curr_topology_plate_ids[point_index]
                
                # Detect collisions at the end of the current time step since we need previous, and current, points and topologies.
                # De-activate point (in 'curr_points') if subducted (forward in time) or consumed back into MOR (backward in time).
                if self.detect_collisions(
                        self.rotation_model,
                        current_time,
                        self.reconstruction_time_interval,
                        prev_point, curr_point,
                        prev_plate_id, self.prev_resolved_plate_boundaries[point_index],
                        curr_plate_id, self.curr_resolved_plate_boundaries[point_index]):
                    # An inactive point in 'curr_points' becomes None.
                    # It may have been reconstructed from the previous time step to a valid position
                    # but now we override that result as inactive.
                    self.curr_points[point_index] = None
                    #self.curr_points.remove(self.curr_points[point_index])
                    self.deletedpoints.append(point_index)
                
        # We successfully reconstructed to the next time.
        return True
    
    
    def _activate_deactivate_points(self):
        
        current_time = self.get_current_time()
        
        # Iterate over all points and activate/deactivate as necessary depending on each point's valid time range.
        for point_index in range(self.num_points):
            if self.curr_points[point_index] is None:
                if not self.point_has_been_activated[point_index]:
                    # Point is not active and has never been activated, so see if can activate it.
                    if (current_time <= self.point_begin_times[point_index] and
                        current_time >= self.point_end_times[point_index]):
                        # The initial point is assumed to be the position at the current time
                        # which is typically the point's begin time (approximately).
                        # But it could be the beginning of the reconstruction time span (specified in constructor)
                        # if that falls in the middle of the point's valid time range - in this case the
                        # initial point position is assumed to be in a position that is some time *after*
                        # it appeared (at its begin time) - and this can happen, for example, if you have a
                        # uniform grids of points at some intermediate time and want to see how they
                        # reconstruct to either a younger or older time (remembering that points can
                        # be subducted forward in time and consumed back into a mid-ocean ridge going
                        # backward in time).
                        self.curr_points[point_index] = self.points[point_index]
                        self.point_has_been_activated[point_index] = True
                        self.num_activated_points += 1
            else:
                # Point is active, so see if can deactivate it.
                if not (current_time <= self.point_begin_times[point_index] and
                    current_time >= self.point_end_times[point_index]):
                    self.curr_points[point_index] = None
    
    
    def _find_resolved_topologies_containing_points(self):
        
        current_time = self.get_current_time()
        
        # Resolve the plate polygons for the current time.
        resolved_topologies = []
        pygplates.resolve_topologies(self.topology_features, self.rotation_model, resolved_topologies, current_time)
        
        if _ReconstructByTopologies.use_plate_partitioner:
            # Create a plate partitioner from the resolved polygons.
            plate_partitioner = pygplates.PlatePartitioner(resolved_topologies, self.rotation_model)
        else:
            # Some of 'curr_points' will be None so 'curr_valid_points' contains only the valid (not None)
            # points, and 'curr_valid_points_indices' is the same length as 'curr_points' but indexes into
            # 'curr_valid_points' so we can quickly find which point (and hence which resolved topology)
            # in 'curr_valid_points' is associated with the a particular point in 'curr_points'.
            curr_valid_points = []
            curr_valid_points_indices = [None] * self.num_points
            for point_index, curr_point in enumerate(self.curr_points):
                if curr_point is not None:
                    curr_valid_points_indices[point_index] = len(curr_valid_points)
                    curr_valid_points.append(curr_point)
            # For each valid current point find the resolved topology containing it.
            resolved_topologies_containing_curr_valid_points = ptt.utils.points_in_polygons.find_polygons(
                    curr_valid_points,
                    [resolved_topology.get_resolved_boundary() for resolved_topology in resolved_topologies],
                    resolved_topologies)
        
        # Iterate over all points.
        for point_index, curr_point in enumerate(self.curr_points):
            
            if curr_point is None:
                # Current point is not currently active - so skip it.
                self.curr_topology_plate_ids[point_index] = None
                self.curr_resolved_plate_boundaries[point_index] = None
                continue
            
            # Find the plate id of the polygon that contains 'curr_point'.
            if _ReconstructByTopologies.use_plate_partitioner:
                curr_polygon = plate_partitioner.partition_point(curr_point)
            else:
                curr_polygon = resolved_topologies_containing_curr_valid_points[
                        # Index back into 'curr_valid_points' and hence also into
                        # 'resolved_topologies_containing_curr_valid_points'.
                        curr_valid_points_indices[point_index]]
            self.curr_resolved_plate_boundaries[point_index] = curr_polygon
            
            # If the polygon is None, that means (presumably) that it fell into a crack between 
            # topologies. So it will be skipped and thrown away from future iterations.
            if curr_polygon is None:
                self.curr_topology_plate_ids[point_index] = None
                continue
            
            # Set the plate ID of resolved topology containing current point.
            self.curr_topology_plate_ids[point_index] = curr_polygon.get_feature().get_reconstruction_plate_id()
