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




class PlotTopologies(object):

    def __init__(self, time, PlateReconstruction_object):

        import ptt

        self.PlateReconstruction_object = PlateReconstruction_object

        # store topologies for easy access




    def update_time(self, time):
        self.time = float(time)
        resolved_topologies = ptt.resolve_topologies.resolve_topologies_into_features(
            self.PlateReconstruction_object.rotation_model,
            self.PlateReconstruction_object.topology_features,
            self.time)

        self.topologies, self.transforms, self.ridges, self.trenches, self.trench_left, self.trench_right, self.other = resolve_topologies


    def _get_reconstructed_lines(feature):
        import shapely

        all_geometries = []
        for feature in reconstructed_feature:

            # get geometry in lon lat order
            geometry = feature.get_reconstructed_geometry().to_lat_lon_array()[::-1,::-1]

            # construct shapely geometry
            geom = shapely.geometry.LineString(geometry)

            # we need to make sure the exterior coordinates are ordered anti-clockwise
            # and the geometry is valid otherwise it will screw with cartopy
            if not geom.exterior.is_ccw:
                geom.exterior.coords = list(geometry[::-1])
            if geom.is_valid:
                all_geometries.append(geom)

        return all_geometries

    def plot_coastlines(self):
        pass

    def plot_ridges(self, ax, **kwargs):
        ax.add_geometries(all_geometries, crs=ccrs.PlateCarree(), **kwargs)

    def plot_trenches(self, ax, **kwargs):
        ax.add_geometries(all_geometries, crs=ccrs.PlateCarree(), **kwargs)

    def plot_subduction_teeth(self, ax, **kwargs):

        # add Subduction Teeth
        subd_xL, subd_yL = tesselate_triangles(
            "reconstructed_topologies/subduction_boundaries_sL_{:.2f}Ma.shp".format(reconstruction_time),
            tesselation_radians=0.1, triangle_base_length=2.0, triangle_aspect=-1.0)
        subd_xR, subd_yR = tesselate_triangles(
            "reconstructed_topologies/subduction_boundaries_sR_{:.2f}Ma.shp".format(reconstruction_time),
            tesselation_radians=0.1, triangle_base_length=2.0, triangle_aspect=1.0)
        
        for tX, tY in zip(subd_xL, subd_yL):
            triangle_xy_points = np.c_[tX, tY]
            patch = plt.Polygon(triangle_xy_points, transform=ccrs.PlateCarree(), **kwargs)
            ax.add_patch(patch)
        for tX, tY in zip(subd_xR, subd_yR):
            triangle_xy_points = np.c_[tX, tY]
            patch = plt.Polygon(triangle_xy_points, transform=ccrs.PlateCarree(), **kwargs)
            ax.add_patch(patch)

    def plot_raster(self, ax, grid, extent=[-180,180,-90,90], **kwargs):
        ax.imshow(grid, origin='lower', extent=extent, transform=ccrs.PlateCarree(), **kwargs)


    def plot_raster_from_netCDF4_file(self, ax, filename, **kwargs):
        from .grids import read_netcdf_grid

        raster, lon_coords, lat_coords = read_netcdf_grid(filename, return_grids=True)
        extent = [lon_coords.min(), lon_coords.max(), lat_coords.min(), lat_coords.max()]
        self.plot_raster(ax, raster, extent=extent, **kwargs)
