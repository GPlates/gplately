import ptt
import numpy as np
import matplotlib.pyplot as plt
import pygplates

def get_valid_geometries(shape_filename):
    """ only return valid geometries """
    import cartopy.io.shapereader as shpreader
    shp_geom = shpreader.Reader(shape_filename).geometries()
    geometries = []
    for record in shp_geom:
        geometries.append(record.buffer(0.0))
    return geometries
    
def add_coastlines(ax, reconstruction_time, **kwargs):
    # write shapefile
    reconstructed_coastlines = []
    pygplates.reconstruct(coastlines, rotation_model, reconstructed_coastlines, float(reconstruction_time),
                          export_wrap_to_dateline=True)
    coastlines_geometries = shapelify_feature_polygons(reconstructed_coastlines)
    
    ax.add_geometries(coastlines_geometries, crs=ccrs.PlateCarree(), **kwargs)
    
def add_continents(ax, reconstruction_time, **kwargs):
    reconstructed_continents = []
    pygplates.reconstruct(continents, rotation_model, reconstructed_continents, float(reconstruction_time),
                          export_wrap_to_dateline=True)
    continent_geometries = shapelify_feature_polygons(reconstructed_continents)
    
    ax.add_geometries(continent_geometries, crs=ccrs.PlateCarree(), **kwargs)
    
def add_ridges(ax, reconstruction_time, **kwargs):
    import shapely
    reconstructed_ridges = get_ridge_transforms(topology_features, rotation_model, float(reconstruction_time))
    all_geometries = []
    for feature in reconstructed_ridges:
        geometry = feature.get_all_geometries()[0].to_lat_lon_array()[::-1,::-1]
        
        # construct shapely geometry
        geom = shapely.geometry.LineString(geometry)

        # we need to make sure the exterior coordinates are ordered anti-clockwise
        # and the geometry is valid otherwise it will screw with cartopy
        if geom.is_valid:
            all_geometries.append(geom)
    
    ax.add_geometries(all_geometries, crs=ccrs.PlateCarree(), **kwargs)
    
def add_ridges(ax, reconstruction_time, **kwargs):
    shp_name = "reconstructed_topologies/ridge_transform_boundaries_{:.2f}Ma.shp".format(reconstruction_time)
    shp_continents = shpreader.Reader(shp_name).geometries()
    ft_continents  = cfeature.ShapelyFeature(shp_continents, ccrs.PlateCarree())
    ax.add_feature(ft_continents, **kwargs)

def add_trenches(ax, reconstruction_time, color='k', linewidth=2, **kwargs):
    shp_name = "reconstructed_topologies/subduction_boundaries_{:.2f}Ma.shp".format(reconstruction_time)
    shp_subd = shpreader.Reader(shp_name).geometries()
    ft_subd  = cfeature.ShapelyFeature(shp_subd, ccrs.PlateCarree())
    ax.add_feature(ft_subd, facecolor='none', edgecolor=color, linewidth=linewidth, zorder=5)
    # add Subduction Teeth
    subd_xL, subd_yL = tesselate_triangles(
        "reconstructed_topologies/subduction_boundaries_sL_{:.2f}Ma.shp".format(reconstruction_time),
        tesselation_radians=0.1, triangle_base_length=2.0, triangle_aspect=-1.0)
    subd_xR, subd_yR = tesselate_triangles(
        "reconstructed_topologies/subduction_boundaries_sR_{:.2f}Ma.shp".format(reconstruction_time),
        tesselation_radians=0.1, triangle_base_length=2.0, triangle_aspect=1.0)
    
    for tX, tY in zip(subd_xL, subd_yL):
        triangle_xy_points = np.c_[tX, tY]
        patch = plt.Polygon(triangle_xy_points, color=color, transform=ccrs.PlateCarree(), zorder=6)
        ax.add_patch(patch)
    for tX, tY in zip(subd_xR, subd_yR):
        triangle_xy_points = np.c_[tX, tY]
        patch = plt.Polygon(triangle_xy_points, color=color, transform=ccrs.PlateCarree(), zorder=6)
        ax.add_patch(patch)
    
    
def add_quiver(ax, reconstruction_time, **kwargs):
    Xnodes, Ynodes, U, V = ptt.velocity_tools.get_velocity_x_y_u_v(reconstruction_time, rotation_model,
                                                                   topology_features)
    mag = np.hypot(U, V)
#     mag = np.clip(mag, 1.0, 1e99)
#     mag[mag==0] = 1 #to avoid 0 divisor
#     U = U/mag
#     V = V/mag
    
    if mag.any():
        ax.quiver(Xnodes, Ynodes, U, V, transform=ccrs.PlateCarree(), **kwargs)



# subduction teeth
def tesselate_triangles(shapefilename, tesselation_radians, triangle_base_length, triangle_aspect=1.0):
    """
    Place subduction teeth along line segments within a MultiLineString shapefile
    
    Parameters
    ----------
        shapefilename  : str  path to shapefile
        tesselation_radians : float
        triangle_base_length : float  length of base
        triangle_aspect : float  aspect ratio
        
    Returns
    -------
        X_points : (n,3) array of triangle x points
        Y_points : (n,3) array of triangle y points
    """

    import shapefile
    shp = shapefile.Reader(shapefilename)

    tesselation_degrees = np.degrees(tesselation_radians)
    triangle_pointsX = []
    triangle_pointsY = []

    for i in range(len(shp)):
        pts = np.array(shp.shape(i).points)

        cum_distance = 0.0
        for p in range(len(pts) - 1):

            A = pts[p]
            B = pts[p+1]

            AB_dist = B - A
            AB_norm = AB_dist / np.hypot(*AB_dist)
            cum_distance += np.hypot(*AB_dist)

            # create a new triangle if cumulative distance is exceeded.
            if cum_distance >= tesselation_degrees:

                C = A + triangle_base_length*AB_norm

                # find normal vector
                AD_dist = np.array([AB_norm[1], -AB_norm[0]])
                AD_norm = AD_dist / np.linalg.norm(AD_dist)

                C0 = A + 0.5*triangle_base_length*AB_norm

                # project point along normal vector
                D = C0 + triangle_base_length*triangle_aspect*AD_norm

                triangle_pointsX.append( [A[0], C[0], D[0]] )
                triangle_pointsY.append( [A[1], C[1], D[1]] )

                cum_distance = 0.0

    shp.close()
    return np.array(triangle_pointsX), np.array(triangle_pointsY)

def shapelify_feature_polygons(features):
    import shapely

    date_line_wrapper = pygplates.DateLineWrapper()

    all_geometries = []
    for feature in features:

        rings = []
        wrapped_polygons = date_line_wrapper.wrap(feature.get_reconstructed_geometry())
        for poly in wrapped_polygons:
            ring = np.array([(p.get_longitude(), p.get_latitude()) for p in poly.get_exterior_points()])
            ring[:,1] = np.clip(ring[:,1], -89, 89) # anything approaching the poles creates artefacts
            ring_polygon = shapely.geometry.Polygon(ring)

            # we need to make sure the exterior coordinates are ordered anti-clockwise
            # and the geometry is valid otherwise it will screw with cartopy
            if not ring_polygon.exterior.is_ccw:
                ring_polygon.exterior.coords = list(ring[::-1])

            rings.append(ring_polygon)

        geom = shapely.geometry.MultiPolygon(rings)

        # we need to make sure the exterior coordinates are ordered anti-clockwise
        # and the geometry is valid otherwise it will screw with cartopy
        all_geometries.append(geom.buffer(0.0)) # add 0.0 buffer to deal with artefacts

    return all_geometries

def shapelify_feature_lines(features):
    import shapely

    date_line_wrapper = pygplates.DateLineWrapper()

    all_geometries = []
    for feature in features:

        rings = []
        wrapped_lines = date_line_wrapper.wrap(feature.get_reconstructed_geometry())
        for line in wrapped_lines:
            ring = np.array([(p.get_longitude(), p.get_latitude()) for p in line.get_points()])
            ring[:,1] = np.clip(ring[:,1], -89, 89) # anything approaching the poles creates artefacts
            ring_linestring = shapely.geometry.LineString(ring)

            rings.append(ring_linestring)

        # construct shapely geometry
        geom = shapely.geometry.MultiLineString(rings)

        if geom.is_valid:
            all_geometries.append(geom)

    return all_geometries



class PlotTopologies(object):

    def __init__(self, PlateReconstruction_object, time, coastline_filename=None, continent_filename=None, COB_filename=None):

        import ptt
        import cartopy.crs as ccrs

        self.PlateReconstruction_object = PlateReconstruction_object
        self.base_projection = ccrs.PlateCarree()

        self.coastline_filename = coastline_filename
        self.continent_filename = continent_filename
        self.COB_filename = COB_filename

        # store topologies for easy access
        # setting time runs the update_time routine
        self.time = time

    @property
    def time(self):
        """ Reconstruction time """
        return self._time

    @time.setter
    def time(self, var):
        if var >= 0:
            self.update_time(var)
        else:
            raise ValueError("Enter a valid time >= 0")


    def update_time(self, time):
        self._time = float(time)
        resolved_topologies = ptt.resolve_topologies.resolve_topologies_into_features(
            self.PlateReconstruction_object.rotation_model,
            self.PlateReconstruction_object.topology_features,
            self.time)

        self.topologies, self.ridge_transforms, self.ridges, self.transforms, self.trenches, self.trench_left, self.trench_right, self.other = resolved_topologies

        # reconstruct other important polygons and lines
        if self.coastline_filename:
            self.coastlines = self.PlateReconstruction_object.reconstruct(
                self.coastline_filename, self.time, from_time=0, anchor_plate_id=0)

        if self.continent_filename:
            self.continents = self.PlateReconstruction_object.reconstruct(
                self.continent_filename, self.time, from_time=0, anchor_plate_id=0)

        if self.COB_filename:
            self.COBs = self.PlateReconstruction_object.reconstruct(
                self.COB_filename, self.time, from_time=0, anchor_plate_id=0)

    def _get_feature_lines(self, features):
        import shapely
        
        date_line_wrapper = pygplates.DateLineWrapper()

        all_geometries = []
        for feature in features:

            # get geometry in lon lat order
            rings = []
            for geometry in feature.get_geometries():
                wrapped_lines = date_line_wrapper.wrap(geometry)
                for line in wrapped_lines:
                    ring = np.array([(p.get_longitude(), p.get_latitude()) for p in line.get_points()])
                    ring[:,1] = np.clip(ring[:,1], -89, 89) # anything approaching the poles creates artefacts
                    ring_linestring = shapely.geometry.LineString(ring)
                    
                    rings.append(ring_linestring)
                
            # construct shapely geometry
            geom = shapely.geometry.MultiLineString(rings)

            if geom.is_valid:
                all_geometries.append(geom)

        return all_geometries

    # subduction teeth
    def _tesselate_triangles(self, features, tesselation_radians, triangle_base_length, triangle_aspect=1.0):
        """
        Place subduction teeth along line segments within a MultiLineString shapefile
        
        Parameters
        ----------
            shapefilename  : str  path to shapefile
            tesselation_radians : float
            triangle_base_length : float  length of base
            triangle_aspect : float  aspect ratio
            
        Returns
        -------
            X_points : (n,3) array of triangle x points
            Y_points : (n,3) array of triangle y points
        """

        tesselation_degrees = np.degrees(tesselation_radians)
        triangle_pointsX = []
        triangle_pointsY = []

        date_line_wrapper = pygplates.DateLineWrapper()


        for feature in features:

            cum_distance = 0.0

            for geometry in feature.get_geometries():
                wrapped_lines = date_line_wrapper.wrap(geometry)
                for line in wrapped_lines:
                    pts = np.array([(p.get_longitude(), p.get_latitude()) for p in line.get_points()])

                    for p in range(0, len(pts) - 1):
                        A = pts[p]
                        B = pts[p+1]

                        AB_dist = B - A
                        AB_norm = AB_dist / np.hypot(*AB_dist)
                        cum_distance += np.hypot(*AB_dist)

                        # create a new triangle if cumulative distance is exceeded.
                        if cum_distance >= tesselation_degrees:

                            C = A + triangle_base_length*AB_norm

                            # find normal vector
                            AD_dist = np.array([AB_norm[1], -AB_norm[0]])
                            AD_norm = AD_dist / np.linalg.norm(AD_dist)

                            C0 = A + 0.5*triangle_base_length*AB_norm

                            # project point along normal vector
                            D = C0 + triangle_base_length*triangle_aspect*AD_norm

                            triangle_pointsX.append( [A[0], C[0], D[0]] )
                            triangle_pointsY.append( [A[1], C[1], D[1]] )

                            cum_distance = 0.0

        return np.array(triangle_pointsX), np.array(triangle_pointsY)

    def plot_coastlines(self, ax, **kwargs):
        if self.coastline_filename is None:
            raise ValueError("Supply coastline_filename to PlotTopologies object")

        coastline_polygons = shapelify_feature_polygons(self.coastlines)
        ax.add_geometries(coastline_polygons, crs=self.base_projection, **kwargs)

    def plot_continents(self, ax, **kwargs):
        if self.continent_filename is None:
            raise ValueError("Supply continent_filename to PlotTopologies object")

        continent_polygons = shapelify_feature_polygons(self.continents)
        ax.add_geometries(continent_polygons, crs=self.base_projection, **kwargs)

    def plot_continent_ocean_boundaries(self, ax, **kwargs):
        if self.COB_filename is None:
            raise ValueError("Supply COB_filename to PlotTopologies object")

        COB_lines = shapelify_feature_lines(self.COBs)
        ax.add_geometries(COB_lines, crs=self.base_projection, facecolor='none', **kwargs)

    def plot_ridges(self, ax, color='black', **kwargs):
        ridge_lines = self._get_feature_lines(self.ridges)
        ax.add_geometries(ridge_lines, crs=self.base_projection, facecolor='none', edgecolor=color, **kwargs)

    def plot_ridges_and_transforms(self, ax, color='black', **kwargs):
        ridge_transform_lines = self._get_feature_lines(self.ridge_transforms)
        ax.add_geometries(ridge_transform_lines, crs=self.base_projection, facecolor='none', edgecolor=color, **kwargs)

    def plot_transforms(self, ax, color='black', **kwargs):
        transform_lines = self._get_feature_lines(self.transforms)
        ax.add_geometries(transform_lines, crs=self.base_projection, facecolor='none', edgecolor=color, **kwargs)

    def plot_trenches(self, ax, color='black', **kwargs):
        trench_lines = self._get_feature_lines(self.trenches)
        ax.add_geometries(trench_lines, crs=self.base_projection, facecolor='none', edgecolor=color, **kwargs)

    def plot_subduction_teeth(self, ax, spacing=0.1, size=2.0, aspect=1, color='black', **kwargs):
        import shapely

        # add Subduction Teeth
        subd_xL, subd_yL = self._tesselate_triangles(
            self.trench_left,
            tesselation_radians=spacing,
            triangle_base_length=size,
            triangle_aspect=-aspect)
        subd_xR, subd_yR = self._tesselate_triangles(
            self.trench_right,
            tesselation_radians=spacing,
            triangle_base_length=size,
            triangle_aspect=aspect)
        
        teeth = []
        for tX, tY in zip(subd_xL, subd_yL):
            triangle_xy_points = np.c_[tX, tY]
            shp = shapely.geometry.Polygon(triangle_xy_points)
            teeth.append(shp)

        for tX, tY in zip(subd_xR, subd_yR):
            triangle_xy_points = np.c_[tX, tY]
            shp = shapely.geometry.Polygon(triangle_xy_points)
            teeth.append(shp)

        ax.add_geometries(teeth, crs=self.base_projection, color=color, **kwargs)

    def plot_grid(self, ax, grid, extent=[-180,180,-90,90], **kwargs):
        ax.imshow(grid, origin='lower', extent=extent, transform=self.base_projection, **kwargs)


    def plot_grid_from_netCDF(self, ax, filename, **kwargs):
        from .grids import read_netcdf_grid

        raster, lon_coords, lat_coords = read_netcdf_grid(filename, return_grids=True)
        extent = [lon_coords.min(), lon_coords.max(), lat_coords.min(), lat_coords.max()]
        self.plot_grid(ax, raster, extent=extent, **kwargs)


    def plot_plate_motion_vectors(self, ax, spacingX, spacingY, **kwargs):
        from .tools import get_point_velocities

        lons = np.arange(-180,180+spacingX,spacingX)
        lats = np.arange(-90,90+spacingY,spacingY)
        lonq, latq = np.meshgrid(lons, lats)

        rotation_model = self.PlateReconstruction_object.rotation_model
        topology_features = self.PlateReconstruction_object.topology_features

        velocities = get_point_velocities(lonq.ravel(), latq.ravel(), topology_features, rotation_model, self.time)
        U = velocities[:,0].reshape(lonq.shape)
        V = velocities[:,1].reshape(latq.shape)

        mag = np.sqrt(U**2 + V**2)
        mag[mag == 0] = 1

        ax.quiver(lonq, latq, U, V, transform=self.base_projection, **kwargs)
