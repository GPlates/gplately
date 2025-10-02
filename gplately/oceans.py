#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

import logging
import math
import multiprocessing
import os
import warnings
from functools import partial
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
import pygplates

from . import grids, tools
from .lib.reconstruct_by_topologies import (
    _ContinentCollision,
    _DefaultCollision,
    _ReconstructByTopologiesImpl,
    _ReconstructByTopologicalModelImpl,
)
from .ptt import continent_contours, separate_ridge_transform_segments
from .tools import _deg2pixels, _pixels2deg
from .utils import seafloor_grid_utils
from .utils.log_utils import get_debug_level

logger = logging.getLogger("gplately")


# JUST FOR PYDOC TO GENERATE THE SAME DOC AS BEFORE
def create_icosahedral_mesh(refinement_levels):
    return seafloor_grid_utils.create_icosahedral_mesh(refinement_levels)


def ensure_polygon_geometry(reconstructed_polygons, rotation_model, time):
    return seafloor_grid_utils.ensure_polygon_geometry(
        reconstructed_polygons, rotation_model, time
    )


def point_in_polygon_routine(multi_point, COB_polygons):
    return seafloor_grid_utils.point_in_polygon_routine(multi_point, COB_polygons)


create_icosahedral_mesh.__doc__ = seafloor_grid_utils.create_icosahedral_mesh.__doc__
ensure_polygon_geometry.__doc__ = seafloor_grid_utils.ensure_polygon_geometry.__doc__
point_in_polygon_routine.__doc__ = seafloor_grid_utils.point_in_polygon_routine.__doc__


class SeafloorGrid(object):
    """Generate grids that track data atop global ocean basin points (which emerge from mid-ocean ridges) through geological time.

    This generates grids of seafloor age, seafloor spreading rate and other oceanic data from the :class:`gplately.PlateReconstruction`
    and :class:`gplately.PlotTopologies` objects.

    Gridding methods in this class have been adapted from Simon Williams' development repository for an
    [auto-age-gridding workflow](https://github.com/siwill22/agegrid-0.1).

    The sample jupyter notebook [10-SeafloorGrid](https://github.com/GPlates/gplately/blob/master/Notebooks/10-SeafloorGrids.ipynb)
    demonstrates how the functionalities within this class work. Below you can find documentation for each of these functionalities.

    Methodology
    -----------
    There are two main steps that this class follows to generate grids:

    1. Preparation for reconstruction by topologies
    2. Reconstruction by topologies

    The preparation step involves building a:

    * global domain of initial points that populate the seafloor at ``max_time``,
    * continental mask that separates ocean points from continent regions per timestep, and
    * set of points that emerge to the left and right of mid-ocean ridge segments per timestep, as well as the z-value to allocate to these points.

    First, the global domain of initial points is created using [stripy's](https://github.com/underworldcode/stripy/blob/master/stripy/spherical_meshes.py#L27)
    icosahedral triangulated mesh. The number of points in this mesh can be controlled using a ``refinement_levels`` integer (the larger this integer,
    the more resolved the continent masks will be).

    ![RefinementLevels](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/seafloorgrid_refinement.png)

    These points are spatially partitioned by plate ID so they can be passed into a point-in-polygon routine.
    This identifies points that lie within continental polygon boundaries and those that are in the ocean. From this, continental masks are built
    per timestep, and the initial seed points are allocated ages at the first reconstruction timestep ``max_time``. Each point's initial age
    is calculated by dividing its proximity to the nearest MOR segment by half its assumed spreading rate. This spreading rate
    (``initial_ocean_mean_spreading_rate``) is assumed to be uniform for all points.

    These initial points momentarily fill the global ocean basin, and all have uniform spreading rates.
    Thus, the spreading rate grid at ``max_time`` will be uniformly populated with the ``initial_ocean_mean_spreading_rate`` (mm/yr).
    The age grid at ``max_time`` will look like a series of smooth, linear age gradients clearly partitioned by tectonic plates with unique plate IDs:

    ![MaxTimeGrids](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/max_time_grids.png)

    Ridge "line" topologies are resolved at each reconstruction time step and partitioned into segments with a valid stage rotation.
    Each segment is further divided into points at a specified ridge sampling spacing (``ridge_sampling``).
    Each point is ascribed a latitude, longitude, spreading rate and age (from plate reconstruction model files, as opposed to ages of the initial ocean mesh points),
    a point index and the general z-value that will be gridded onto it.

    ![NewRidgePoints](https://raw.githubusercontent.com/GPlates/gplately/master/Notebooks/NotebookFiles/pdoc_Files/new_ridge_points.png)

    Reconstruction by topologies involves determining which points are active and inactive (collided with a continent or subducted at a trench)
    for each reconstruction time step. This is done using a hidden object in :class:`PlateReconstruction` called ``ReconstructByTopologies``.

    If an ocean point with a certain velocity on one plate ID transitions into another rigid plate ID at another timestep (with another velocity),
    the velocity difference between both plates is calculated. The point may have subducted/collided with a continent if this velocity difference
    is higher than a specified velocity threshold (which can be controlled with ``subduction_collision_parameters``). To ascertain whether the point
    should be deactivated, a displacement test is conducted. If the proximity of the point's previous time position to the polygon boundary it is
    approaching is higher than a set distance threshold, then the point is far enough away from the boundary that it cannot be subducted or consumed by it,
    and hence the point is still active. Otherwise, it is deemed inactive and deleted from the ocean basin mesh.

    With each reconstruction time step, points from mid-ocean ridges (which have more accurate spreading rates and attributed valid times) will spread across
    the ocean floor. Eventually, points will be pushed into continental boundaries or subduction zones, where they are deleted. Ideally, all initial ocean points
    (from the Stripy icosahedral mesh) should be deleted over time. However, not all will be deleted - such points typically reside near continental boundaries.
    This happens if the emerged ridge points do not spread far enough to "phase out" these points at collision regions - likely due to insufficient reconstruction detail.
    These undeleted points form artefacts of anomalously high seafloor age that append over the reconstruction time range.

    Once reconstruction by topologies determines the ocean basin snapshot per timestep, a data frame of all longitudes, latitudes, seafloor ages, spreading rates and
    any other attributed z values will be written to a gridding input file per timestep.

    Each active longitude, latitude and chosen z value is gridded using nearest-neighbour interpolation and written to a netCDF4 format.
    """

    # Keys (column headers) of the gridded data.
    CURRENT_LONGITUDES_KEY = "CURRENT_LONGITUDES"
    CURRENT_LATITUDES_KEY = "CURRENT_LATITUDES"
    SEAFLOOR_AGE_KEY = "SEAFLOOR_AGE"
    SPREADING_RATE_KEY = "SPREADING_RATE"
    BIRTH_LAT_SNAPSHOT_KEY = "BIRTH_LAT_SNAPSHOT"
    POINT_ID_SNAPSHOT_KEY = "POINT_ID_SNAPSHOT"

    # Mid-ocean ridge pickle filename format (accepting time).
    MOR_PKL_FILE_NAME = "MOR_{:0.2f}_Ma.pkl"

    def __init__(
        self,
        PlateReconstruction_object,
        PlotTopologies_object,
        max_time: Union[float, int],
        min_time: Union[float, int],
        ridge_time_step: Union[float, int],
        *,
        save_directory: Union[str, Path] = "seafloor-grid-output",
        file_collection: str = "",
        refinement_levels: int = 5,
        ridge_sampling: float = 0.5,
        extent: Tuple = (-180, 180, -90, 90),
        grid_spacing: float = 0.1,
        subduction_collision_parameters=(5.0, 10.0),
        initial_ocean_mean_spreading_rate: float = 75.0,
        resume_from_checkpoints=False,
        continent_mask_filename=None,
        use_continent_contouring=False,
        nprocs=-2,
        **kwargs,
    ):
        """Constructor. Create a :class:`SeafloorGrid` object.

        Parameters
        ----------
        PlateReconstruction_object : PlateReconstruction
            A :class:`PlateReconstruction` object with a `pygplates.RotationModel`_ and
            a `pygplates.FeatureCollection`_ containing topology features.
        PlotTopologies_object : PlotTopologies
            A :class:`PlotTopologies` object with a continental polygon or COB terrane polygon file to mask grids with.
        max_time : float
            The maximum time for age gridding.
        min_time : float
            The minimum time for age gridding.
        ridge_time_step : float
            The delta time for resolving ridges (and thus age gridding).
        save_directory : str, default=None
            The top-level directory to save all outputs to.
        file_collection : str, default=""
            A string to identify the plate model used (will be automated later).
        refinement_levels : int, default=5
            Control the number of points in the icosahedral mesh (higher integer means higher resolution of continent masks).
        ridge_sampling : float, default=0.5
            Spatial resolution (in degrees) at which points that emerge from ridges are tessellated.
        extent : tuple of 4, default=(-180.,180.,-90.,90.)
            A tuple containing the mininum longitude, maximum longitude, minimum latitude and
            maximum latitude extents for all masking and final grids.
        grid_spacing : float, default=0.1
            The degree spacing/interval with which to space grid points across all masking and
            final grids. If ``grid_spacing`` is provided, all grids will use it. If not, ``grid_spacing`` defaults to 0.1.
        subduction_collision_parameters : len-2 tuple of float, default=(5.0, 10.0)
            A 2-tuple of (threshold velocity delta in kms/my, threshold distance to boundary per My in kms/my)
        initial_ocean_mean_spreading_rate : float, default=75.
            A spreading rate to uniformly allocate to points that define the initial ocean
            basin. These points will have inaccurate ages, but most of them will be phased
            out after points with plate-model prescribed ages emerge from ridges and spread
            to push them towards collision boundaries (where they are deleted).
        resume_from_checkpoints : bool, default=False
            If set to ``True``, and the gridding preparation stage (continental masking and/or
            ridge seed building) is interrupted, SeafloorGrids will resume gridding preparation
            from the last successful preparation time.
            If set to ``False``, SeafloorGrids will automatically overwrite all files in the
            ``save_directory`` if re-run after interruption, or normally re-run, thus beginning
            gridding preparation from scratch. ``False`` will be useful if data allocated to the
            MOR seed points need to be augmented.
        continent_mask_filename : str, optional
            An optional parameter pointing to the full path to a continental mask for each timestep.
            Assuming the time is in the filename, i.e. ``/path/to/continent_mask_0Ma.nc``, it should be
            passed as ``/path/to/continent_mask_{}Ma.nc`` with curly brackets. Include decimal formatting if needed.
        use_continent_contouring: bool, default=False
            Note that this is ignored if ``continent_mask_filename`` is specified.
            If ``True`` then builds the continent mask for a given time using ptt's 'continent contouring' method
            (for more information about 'Continent Contouring', visit https://github.com/EarthByte/continent-contouring).
            If ``False`` then builds the continent masks using the continents of ``PlotTopologies_object``.
        nprocs : int, default=-2
            The number of CPUs to use for parts of the code that are parallelized.
            Must be an integer or convertible to an integer (eg, float is rounded towards zero).
            If positive then uses that many CPUs.
            If ``1`` then executes in serial (ie, is not parallelized).
            If ``0`` then a ``ValueError`` is raised.
            If ``-1`` then all available CPUs are used.
            If ``-2`` then all available CPUs except one are used, etc.
            Defaults to ``-2`` (ie, uses all available CPUs except one to keep system responsive).
        **kwargs
            Handle deprecated arguments such as ``zval_names``.



        .. _pygplates.RotationModel: https://www.gplates.org/docs/pygplates/generated/pygplates.rotationmodel
        .. _pygplates.Feature: https://www.gplates.org/docs/pygplates/generated/pygplates.feature#pygplates.Feature
        .. _pygplates.FeatureCollection: https://www.gplates.org/docs/pygplates/generated/pygplates.featurecollection#pygplates.FeatureCollection
        """

        # Handle deprecated arguments.
        if kwargs.pop("zval_names", None):
            warnings.warn(
                "`zval_names` keyword argument has been deprecated, it is no longer used",
                DeprecationWarning,
            )
        for key in kwargs.keys():
            raise TypeError(
                "SeafloorGrid.__init__() got an unexpected keyword argument '{}'".format(
                    key
                )
            )

        # Determine number of CPUs to use.
        self.num_cpus = _get_num_cpus(nprocs)

        self.plate_reconstruction = PlateReconstruction_object
        self.plot_topologies = PlotTopologies_object

        self.file_collection = file_collection

        if continent_mask_filename:
            # Filename for continental masks that the user can provide instead of building it here
            self.continent_mask_filepath = continent_mask_filename
            self.continent_mask_is_provided = True
        else:
            self.continent_mask_is_provided = False

        self.use_continent_contouring = use_continent_contouring

        self._setup_output_paths(save_directory)

        # Topological parameters
        self.refinement_levels = refinement_levels
        self.ridge_sampling = ridge_sampling
        self.subduction_collision_parameters = subduction_collision_parameters
        self.initial_ocean_mean_spreading_rate = initial_ocean_mean_spreading_rate

        # Gridding parameters
        self.extent = extent

        self._set_grid_resolution(grid_spacing)

        self.resume_from_checkpoints = resume_from_checkpoints

        # Temporal parameters
        self._max_time = max_time
        self._min_time = min_time
        self._ridge_time_step = ridge_time_step
        self._times = np.arange(
            self._max_time, self._min_time - 1e-4, -self._ridge_time_step
        )

        # Get the reconstructed continents at the max time.
        # But first need to set the max time on 'self.plot_topologies'.
        self.plot_topologies.time = self._max_time
        reconstructed_continents = self.plot_topologies.continents

        # Essential features and meshes for the SeafloorGrid
        self.continental_polygons = ensure_polygon_geometry(
            reconstructed_continents,
            self.plate_reconstruction.rotation_model,
            self._max_time,
        )

        self.icosahedral_multi_point = create_icosahedral_mesh(self.refinement_levels)

    def _map_res_to_node_percentage(self, continent_mask_filename):
        """Determine which percentage to use to scale the continent mask resolution at max time."""
        maskY, maskX = grids.read_netcdf_grid(
            continent_mask_filename.format(self._max_time)
        ).shape

        mask_deg = _pixels2deg(maskX, self.extent[0], self.extent[1])

        if mask_deg <= 0.1:
            percentage = 0.1
        elif mask_deg <= 0.25:
            percentage = 0.3
        elif mask_deg <= 0.5:
            percentage = 0.5
        elif mask_deg < 0.75:
            percentage = 0.6
        elif mask_deg >= 1:
            percentage = 0.75
        else:
            raise Exception(
                "Unable to determine percentage. The code here is suspicious. I am unsure about it."
            )
        return mask_deg, percentage

    def _setup_output_paths(self, save_directory):
        """Create various folders for output files."""
        self.save_directory = Path(save_directory)

        # middle ocean ridge files
        self.mid_ocean_ridges_dir = os.path.join(
            self.save_directory, "middle_ocean_ridges"
        )
        Path(self.mid_ocean_ridges_dir).mkdir(parents=True, exist_ok=True)
        if self.file_collection:
            self.mid_ocean_ridges_file_path = os.path.join(
                self.mid_ocean_ridges_dir,
                self.file_collection + "_" + self.MOR_PKL_FILE_NAME,
            )
        else:
            self.mid_ocean_ridges_file_path = os.path.join(
                self.mid_ocean_ridges_dir, self.MOR_PKL_FILE_NAME
            )

        # continent mask files
        # only generate continent mask files if user does not provide them
        if not self.continent_mask_is_provided:
            self.continent_mask_directory = os.path.join(
                self.save_directory, "continent_mask"
            )
            Path(self.continent_mask_directory).mkdir(parents=True, exist_ok=True)

            if self.use_continent_contouring:
                continent_mask_file_basename = (
                    "continent_mask_ptt_continent_contouring_{:0.2f}Ma.nc"
                )
            else:
                continent_mask_file_basename = "continent_mask_{:0.2f}Ma.nc"

            if self.file_collection:
                continent_mask_file_basename = (
                    self.file_collection + "_" + continent_mask_file_basename
                )

            self.continent_mask_filepath = os.path.join(
                self.continent_mask_directory, continent_mask_file_basename
            )

        # gridding input files
        self.gridding_input_directory = os.path.join(
            self.save_directory, "gridding_input"
        )
        Path(self.gridding_input_directory).mkdir(parents=True, exist_ok=True)
        gridding_input_basename = "gridding_input_{:0.2f}Ma.npz"
        if self.file_collection:
            gridding_input_basename = (
                self.file_collection + "_" + gridding_input_basename
            )
        self.gridding_input_filepath = os.path.join(
            self.gridding_input_directory, gridding_input_basename
        )

    def _set_grid_resolution(self, grid_spacing=0.1):
        """Determine the output grid resolution."""
        if not grid_spacing:
            grid_spacing = 0.1
        # A list of degree spacings that allow an even division of the global lat-lon extent.
        divisible_degree_spacings = [0.1, 0.25, 0.5, 0.75, 1.0]

        # If the provided degree spacing is in the list of permissible spacings, use it
        # and prepare the number of pixels in x and y (spacingX and spacingY)
        if grid_spacing in divisible_degree_spacings:
            self.grid_spacing = grid_spacing
            self.spacingX = _deg2pixels(grid_spacing, self.extent[0], self.extent[1])
            self.spacingY = _deg2pixels(grid_spacing, self.extent[2], self.extent[3])

        # If the provided spacing is >>1 degree, use 1 degree
        elif grid_spacing >= divisible_degree_spacings[-1]:
            self.grid_spacing = divisible_degree_spacings[-1]
            self.spacingX = _deg2pixels(
                divisible_degree_spacings[-1], self.extent[0], self.extent[1]
            )
            self.spacingY = _deg2pixels(
                divisible_degree_spacings[-1], self.extent[2], self.extent[3]
            )

            with warnings.catch_warnings():
                warnings.simplefilter("always")
                warnings.warn(
                    f"The provided grid_spacing of {grid_spacing} is quite large. To preserve the grid resolution, a {self.grid_spacing} degree spacing has been employed instead."
                )

        # If the provided degree spacing is not in the list of permissible spacings, but below
        # a degree, find the closest permissible degree spacing. Use this and find
        # spacingX and spacingY.
        else:
            for divisible_degree_spacing in divisible_degree_spacings:
                # The tolerance is half the difference between consecutive divisible spacings.
                # Max is 1 degree for now - other integers work but may provide too little of a
                # grid resolution.
                if abs(grid_spacing - divisible_degree_spacing) <= 0.125:
                    new_deg_res = divisible_degree_spacing
                    self.grid_spacing = new_deg_res
                    self.spacingX = _deg2pixels(
                        new_deg_res, self.extent[0], self.extent[1]
                    )
                    self.spacingY = _deg2pixels(
                        new_deg_res, self.extent[2], self.extent[3]
                    )

            with warnings.catch_warnings():
                warnings.simplefilter("always")
                warnings.warn(
                    f"The provided grid_spacing of {grid_spacing} does not cleanly divide into the global extent. A degree spacing of {self.grid_spacing} has been employed instead."
                )

    # Allow SeafloorGrid time to be updated, and to update the internally-used
    # PlotTopologies' time attribute too. If PlotTopologies is used outside the
    # object, its `time` attribute is not updated.
    @property
    def max_time(self):
        """The reconstruction time.

        :type: float
        """
        return self._max_time

    @property
    def PlotTopologiesTime(self):
        """The :attr:`PlotTopologies.time` attribute.

        :type: float
        """
        return self.plot_topologies.time

    @max_time.setter
    def max_time(self, var):
        if var >= 0:
            self.update_time(var)
        else:
            raise ValueError("Enter a valid time >= 0")

    def update_time(self, max_time: float):
        """Set the new reconstruction time.

        Parameters
        ----------
        max_time: float
            The new reconstruction time.
        """
        self._max_time = float(max_time)
        self.plot_topologies.time = float(max_time)

    def _generate_ocean_points(self):
        """Generate ocean points by using the icosahedral mesh."""
        # Get the reconstructed continents at the max time.
        # But first need to set the max time on 'self.plot_topologies'.
        self.plot_topologies.time = self._max_time
        reconstructed_continents = self.plot_topologies.continents

        # Ensure COB terranes at max time have reconstruction IDs and valid times
        COB_polygons = ensure_polygon_geometry(
            reconstructed_continents,
            self.plate_reconstruction.rotation_model,
            self._max_time,
        )

        # zval is a binary array encoding whether a point
        # coordinate is within a COB terrane polygon or not.
        # Use the icosahedral mesh MultiPointOnSphere attribute
        _, ocean_basin_point_mesh, zvals = point_in_polygon_routine(
            self.icosahedral_multi_point, COB_polygons
        )

        ocean_pt_feature = pygplates.Feature()
        ocean_pt_feature.set_geometry(
            pygplates.MultiPointOnSphere(ocean_basin_point_mesh)
        )
        return [ocean_pt_feature]

    def _get_ocean_points_from_continent_mask(self):
        """Get the ocean points from continent mask grid."""
        max_time_cont_mask = grids.Raster(
            self.continent_mask_filepath.format(self._max_time)
        )
        # If the user provides a continental mask filename, we need to downsize the mask
        # resolution for when we create the initial ocean mesh. The mesh does not need to be high-res.
        # If the input grid is at 0.5 degree uniform spacing, then the input
        # grid is 7x more populated than a 6-level stripy icosahedral mesh and
        # using this resolution for the initial ocean mesh will dramatically slow down reconstruction by topologies.
        # Scale down the resolution based on the input mask resolution
        _, percentage = self._map_res_to_node_percentage(self.continent_mask_filepath)
        max_time_cont_mask.resize(
            int(max_time_cont_mask.shape[0] * percentage),
            int(max_time_cont_mask.shape[1] * percentage),
            inplace=True,
        )

        lat = np.linspace(-90, 90, max_time_cont_mask.shape[0])
        lon = np.linspace(-180, 180, max_time_cont_mask.shape[1])

        llon, llat = np.meshgrid(lon, lat)

        mask_inds = np.where(max_time_cont_mask.data.flatten() == 0)
        mask_vals = max_time_cont_mask.data.flatten()
        mask_lon = llon.flatten()[mask_inds]
        mask_lat = llat.flatten()[mask_inds]

        ocean_pt_feature = pygplates.Feature()
        ocean_pt_feature.set_geometry(
            pygplates.MultiPointOnSphere(zip(mask_lat, mask_lon))
        )
        return [ocean_pt_feature]

    def create_initial_ocean_seed_points(self):
        """Create the initial ocean basin seed point domain (at ``max_time`` only)
        using Stripy's icosahedral triangulation with the specified ``self.refinement_levels``.

        The ocean mesh starts off as a global-spanning Stripy icosahedral mesh.
        ``create_initial_ocean_seed_points`` passes the automatically-resolved-to-current-time
        continental polygons from the :attr:`PlotTopologies.continents` attribute
        (which can be from a COB terrane file or a continental polygon file) into
        Plate Tectonic Tools' point-in-polygon routine. It identifies ocean basin points that lie:

        * outside the polygons (for the ocean basin point domain)
        * inside the polygons (for the continental mask)

        Points from the mesh outside the continental polygons make up the ocean basin seed
        point mesh.

        .. note::

            This point mesh represents ocean basin seafloor that was produced
            before :attr:`SeafloorGrid.max_time`, and thus has unknown properties like valid
            time and spreading rate. As time passes, the plate reconstruction model sees
            points emerging from MORs. These new points spread to occupy the ocean basins,
            moving the initial filler points closer to subduction zones and continental
            polygons with which they can collide. If a collision is detected, these points are deleted.

            Ideally, if a reconstruction tree spans a large time range, **all** initial mesh
            points would collide with a continent or be subducted, leaving behind a mesh of
            well-defined MOR-emerged ocean basin points that data can be attributed to.
            However, some of these initial points situated close to continental boundaries are
            retained through time - these form point artefacts with anomalously high ages. Even
            deep-time plate models (e.g. 1 Ga) will have these artefacts - removing them would
            require more detail to be added to the reconstruction model.

        Returns
        -------
        lons, lats, ages, spreading_rates : ndarray
            The initial ocean seed point locations, ages and spreading rates.
        """

        if (
            os.path.isfile(self.continent_mask_filepath.format(self._max_time))
            and self.continent_mask_is_provided
        ):
            # If a set of continent masks was passed, we can use max_time's continental
            # mask to build the initial profile of seafloor age.
            ocean_points = self._get_ocean_points_from_continent_mask()
        else:
            ocean_points = self._generate_ocean_points()

        # Now that we have ocean points...
        # Determine age of ocean basin points using their proximity to MOR features
        # and an assumed globally-uniform ocean basin mean spreading rate.
        # We need resolved topologies at the `max_time` to pass into the proximity function
        resolved_topologies = self.plate_reconstruction.topological_snapshot(
            self._max_time
        ).get_resolved_topologies()
        lons, lats, distances_to_ridge = tools.find_distance_to_nearest_ridge(
            resolved_topologies,
            ocean_points,
        )

        lons, lats = np.array(lons), np.array(lats)

        # Divide spreading rate by 2 to use half the mean spreading rate
        ages = np.array(distances_to_ridge) / (
            self.initial_ocean_mean_spreading_rate / 2.0
        )

        # Constant initial spreading rate.
        spreading_rates = np.full(len(lons), self.initial_ocean_mean_spreading_rate)

        # the code below is for debug purpose only
        if get_debug_level() > 100:
            #
            # Save initial spreading ages and (constant) spreading rates as a coverage to GPML.
            #

            initial_ocean_coverage_feature = pygplates.Feature()
            initial_ocean_coverage_geometry = pygplates.MultiPointOnSphere(
                zip(lats, lons)
            )
            initial_ocean_coverage_scalars = {
                pygplates.ScalarType.create_gpml("SeafloorAge"): ages,
                pygplates.ScalarType.create_gpml("SpreadingRate"): spreading_rates,
            }
            initial_ocean_coverage_feature.set_geometry(
                (initial_ocean_coverage_geometry, initial_ocean_coverage_scalars)
            )
            # Only visible at initial time (max time).
            initial_ocean_coverage_feature.set_valid_time(
                self._max_time + 0.5 * self._ridge_time_step,
                self._max_time - 0.5 * self._ridge_time_step,
            )  # type: ignore

            basename = "ocean_basin_seed_points_{}_RLs_{}Ma.gpmlz".format(
                self.refinement_levels,
                self._max_time,
            )
            if self.file_collection:
                basename = "{}_{}".format(self.file_collection, basename)

            pygplates.FeatureCollection(initial_ocean_coverage_feature).write(
                os.path.join(self.save_directory, basename)
            )

        return lons, lats, ages, spreading_rates

    def build_all_MOR_seedpoints(self):
        """Resolve mid-ocean ridges for all times between ``min_time`` and ``max_time``, divide them
        into points that make up their shared sub-segments. Rotate these points to the left
        and right of the ridge using their stage rotation so that they spread from the ridge.

        Z-value allocation to each point is done here. In future, a function (like
        the spreading rate function) to calculate general z-data will be an input parameter.

        .. note::

            If MOR seed point building is interrupted, progress is safeguarded as long as
            ``resume_from_checkpoints`` is set to ``True``.

            This assumes that points spread from ridges symmetrically, with the exception of
            large ridge jumps at successive timesteps. Therefore, z-values allocated to ridge-emerging
            points will appear symmetrical until changes in spreading ridge geometries create
            asymmetries.

            In future, this will have a checkpoint save feature so that execution
            (which occurs during preparation for ``ReconstructByTopologies`` and can take several hours)
            can be safeguarded against run interruptions.

        .. seealso::

            `Get tessellated points along a mid ocean ridge <https://github.com/siwill22/agegrid-0.1/blob/master/automatic_age_grid_seeding.py#L117>`__.
        """
        overwrite = True
        if self.resume_from_checkpoints:
            overwrite = False

        if self.num_cpus > 1:
            with multiprocessing.Pool(self.num_cpus) as pool:
                #
                # Temporary hack to avoid a large slowdown due to pickling pygplates.RotationModel and pygplates.FeatureCollection.
                # Instead of pickling pygplates.RotationModel and pygplates.FeatureCollection we instead pickle 'self.plate_reconstruction'
                # which already attempts to avoid pickling those pygplates objects.
                #
                # Note: It only works when 'self.plate_reconstruction` was created using rotation *filenames* and topology *filenames*
                #       (which is the case for the plate models downloaded by the PlateModelManager).
                #       When they're not filenames then pygplates.RotationModel and pygplates.FeatureCollection get pickled anyway.
                #
                # TODO: When pyGPlates can pickle pygplates.RotationModel and pygplates.FeatureCollection noticeably faster (than it currently does in pyGPlates 1.0), then
                #       remove '_generate_mid_ocean_ridge_points_parallel()' and instead use '_generate_mid_ocean_ridge_points()'.
                pool.map(
                    partial(
                        _generate_mid_ocean_ridge_points_parallel,
                        delta_time=self._ridge_time_step,
                        mid_ocean_ridges_file_path=self.mid_ocean_ridges_file_path,
                        plate_reconstruction=self.plate_reconstruction,
                        ridge_sampling=self.ridge_sampling,
                        overwrite=overwrite,
                    ),
                    self._times[1:],  # exclude initial (max) time
                )
        else:
            for time in self._times[1:]:  # exclude initial (max) time
                _generate_mid_ocean_ridge_points(
                    time,
                    delta_time=self._ridge_time_step,
                    mid_ocean_ridges_file_path=self.mid_ocean_ridges_file_path,
                    rotation_model=self.plate_reconstruction.rotation_model,
                    topology_features=self.plate_reconstruction.topology_features,
                    ridge_sampling=self.ridge_sampling,
                    overwrite=overwrite,
                )

    def build_all_continental_masks(self):
        """Create a continental mask to define the ocean basin for all times between ``min_time`` and ``max_time``.

        .. note::

            Continental masking progress is safeguarded if ever masking is interrupted,
            provided that ``resume_from_checkpoints`` is set to ``True``.

            The continental masks will be saved to ``continent_mask_{time}Ma.nc`` as compressed netCDF4 files.
        """
        if not self.continent_mask_is_provided:
            overwrite = True
            if self.resume_from_checkpoints:
                overwrite = False

            if self.num_cpus > 1:
                #
                # Temporary hack to avoid a large slowdown due to pickling pygplates.RotationModel and pygplates.FeatureCollection.
                # Instead of pickling pygplates.RotationModel and pygplates.FeatureCollection we instead pickle 'self.plate_reconstruction'
                # and 'self.plot_topologies' which both already attempt to avoid pickling those pygplates objects.
                #
                # Note: This only works when 'self.plate_reconstruction` and 'self.plot_topologies' were created using rotation *filenames* and topology *filenames*
                #       (which is the case for the plate models downloaded by the PlateModelManager).
                #       When they're not filenames then pygplates.RotationModel and pygplates.FeatureCollection get pickled anyway.
                #
                # TODO: When pyGPlates can pickle pygplates.RotationModel and pygplates.FeatureCollection noticeably faster (than it currently does in pyGPlates 1.0), then
                #       remove the '*_parallel()' versions of the function (and instead use the serial versions, but still in parallel).
                # Note: '_build_continental_mask()', unlike '_build_continental_mask_with_contouring()', passes *reconstructed* continents.
                #       However, *reconstructed* feature geometries cannot be pickled (by pygplates) - only regular features can.
                #       This means '_build_continental_mask_parallel()' may still be required in order to pass 'self.plot_topologies' which is pickled and
                #       then it reconstructs the continents after it is unpickled). Either that or pass the unreconstructed continent features and explicitly
                #       reconstruct them after unpickling.
                #
                if self.use_continent_contouring:
                    with multiprocessing.Pool(self.num_cpus) as pool:
                        pool.map(
                            partial(
                                _build_continental_mask_with_contouring_parallel,
                                continent_mask_filepath=self.continent_mask_filepath,
                                plate_reconstruction=self.plate_reconstruction,
                                plot_topologies=self.plot_topologies,
                                overwrite=overwrite,
                            ),
                            self._times,
                        )
                else:
                    with multiprocessing.Pool(self.num_cpus) as pool:
                        pool.map(
                            partial(
                                _build_continental_mask_parallel,
                                continent_mask_filepath=self.continent_mask_filepath,
                                plate_reconstruction=self.plate_reconstruction,
                                # Note: We can't pickle 'self.plot_topologies.continents' because they are *reconstructed* pygplates objects.
                                #       Instead we pickly 'self.plot_topologies' and then ask it to reconstruct the continents after it's unpickled...
                                plot_topologies=self.plot_topologies,
                                spacingY=self.spacingY,
                                spacingX=self.spacingX,
                                extent=self.extent,
                                overwrite=overwrite,
                            ),
                            self._times,
                        )
            else:
                if self.use_continent_contouring:
                    for time in self._times:
                        _build_continental_mask_with_contouring(
                            time,
                            continent_mask_filepath=self.continent_mask_filepath,
                            rotation_model=self.plate_reconstruction.rotation_model,
                            continent_features=self.plot_topologies._continents,
                            overwrite=overwrite,
                        )
                else:
                    for time in self._times:
                        # Get the reconstructed continents at 'time''.
                        # But first need to set 'time' on 'self.plot_topologies'.
                        self.plot_topologies.time = time
                        reconstructed_continents = self.plot_topologies.continents
                        _build_continental_mask(
                            time,
                            continent_mask_filepath=self.continent_mask_filepath,
                            rotation_model=self.plate_reconstruction.rotation_model,
                            reconstructed_continents=reconstructed_continents,
                            spacingY=self.spacingY,
                            spacingX=self.spacingX,
                            extent=self.extent,
                            overwrite=overwrite,
                        )

    def prepare_for_reconstruction_by_topologies(self):
        """Prepare three main input data sets for seafloor data gridding:

        * Initial ocean seed points (at ``max_time``)
        * Continental masks (from ``max_time`` to ``min_time``)
        * MOR points (from ``max_time`` to ``min_time``)

        Return lists of all attributes for the initial ocean point mesh and
        all ridge points for all times in the reconstruction time array.
        """

        # Create initial ocean seed points -------------------------------------------------

        initial_lons, initial_lats, initial_ages, initial_spreading_rates = (
            self.create_initial_ocean_seed_points()
        )
        logger.info("Finished building initial_ocean_seed_points!")

        # Create continental masks and mid-ocean ridge seed points -------------------------

        # If 'self.resume_from_checkpoints' is True then only those continental mask files and
        # MOR seed point files that don't exist will be generated. This enables a previously
        # interrupted workflow to be continued.

        self.build_all_continental_masks()

        self.build_all_MOR_seedpoints()

        # Load the initial ocean seed points (at initial time) -----------------------------

        all_points = [
            pygplates.PointOnSphere(lat, lon)
            for lon, lat in zip(initial_lons, initial_lats)
        ]
        # Appearance time of initial points is their initial ages plus the initial time (max time).
        all_appearance_times = self._max_time + initial_ages
        all_birth_lats = initial_lats
        all_spreading_rates = initial_spreading_rates

        # Load the MOR seed points (at each time step) -------------------------------------

        for time in self._times[1:]:  # exclude initial (max) time
            # load MOR points for each time step
            mor_df = pd.read_pickle(self.mid_ocean_ridges_file_path.format(time))
            mor_lons = mor_df[self.CURRENT_LONGITUDES_KEY].to_numpy()
            mor_lats = mor_df[self.CURRENT_LATITUDES_KEY].to_numpy()
            mor_spreading_rates = mor_df[self.SPREADING_RATE_KEY].to_numpy()

            all_points.extend(
                pygplates.PointOnSphere(lat, lon)
                for lon, lat in zip(mor_lons, mor_lats)
            )
            all_appearance_times = np.concatenate(
                (all_appearance_times, np.full(len(mor_lons), time))
            )
            all_birth_lats = np.concatenate((all_birth_lats, mor_lats))
            all_spreading_rates = np.concatenate(
                (all_spreading_rates, mor_spreading_rates)
            )

        return all_points, all_appearance_times, all_birth_lats, all_spreading_rates

    def reconstruct_by_topologies(self):
        """Obtain all active ocean seed points which are points that have not been consumed at subduction zones
        or have not collided with continental polygons. Active points' latitudes, longitues, seafloor ages, spreading rates and all
        other general z-values are saved to a gridding input file (.npz).
        """
        self._reconstruct_by_topologies_impl(use_topological_model=False)

    def reconstruct_by_topological_model(self):
        """Use `pygplates.TopologicalModel`_ class to reconstruct seed points.
        This method is an alternative to :meth:`reconstruct_by_topologies()` which uses Python code to do the reconstruction.

        .. _pygplates.TopologicalModel: https://www.gplates.org/docs/pygplates/generated/pygplates.topologicalmodel
        """
        self._reconstruct_by_topologies_impl(use_topological_model=True)

    def _reconstruct_by_topologies_impl(self, use_topological_model):
        logger.info("Preparing all initial files...")

        # Obtain all info from the ocean seed points and all MOR points through time, store in arrays.
        (
            all_points,
            all_appearance_times,
            all_birth_lats,
            all_spreading_rates,
        ) = self.prepare_for_reconstruction_by_topologies()

        # Indices for all points that have existed from `max_time` to `min_time`.
        all_point_ids = np.arange(len(all_points))

        ####  Begin reconstruction by topology process:

        # Specify the default collision detection region as subduction zones
        default_collision = _DefaultCollision(
            feature_specific_collision_parameters=[
                (
                    pygplates.FeatureType.gpml_subduction_zone,
                    self.subduction_collision_parameters,
                )
            ]
        )
        # In addition to the default subduction detection, also detect continental collisions
        collision_spec = _ContinentCollision(
            # This filename string should not have a time formatted into it - this is
            # taken care of later.
            self.continent_mask_filepath,
            default_collision,
            verbose=False,
        )

        # Call the reconstruct by topologies object
        if use_topological_model:
            topology_reconstruction = _ReconstructByTopologicalModelImpl(
                self.plate_reconstruction.rotation_model,
                self.plate_reconstruction.topology_features,
                self._max_time,
                self._min_time,
                self._ridge_time_step,
                all_points,
                point_begin_times=all_appearance_times,
                detect_collisions=collision_spec,
            )
        else:
            topology_reconstruction = _ReconstructByTopologiesImpl(
                self.plate_reconstruction.rotation_model,
                self.plate_reconstruction.topology_features,
                self._max_time,
                self._min_time,
                self._ridge_time_step,
                all_points,
                point_begin_times=all_appearance_times,
                detect_collisions=collision_spec,
            )

        # Initialise the reconstruction.
        topology_reconstruction.begin_reconstruction()

        # Loop over the reconstruction times until the end of the reconstruction time span, or until
        # all points have entered their valid time range *and* either exited their time range or
        # have been deactivated (subducted forward in time or consumed by MOR backward in time).
        while True:
            logger.info(
                f"Reconstruct by topologies: working on time {topology_reconstruction.get_current_time():0.2f} Ma"
            )

            # Get *all* the points at the current time step.
            #
            # This is the same size as 'all_points', which is all initial ocean points at `max_time` and all MOR seed points
            # that have ever emerged from spreading ridges through `max_time` to `min_time`.
            # Here 'all_current_points', despite being the same size, differs in that any *inactive* points at the current time
            # are None and the *active* points are at their reconstructed position at the current time.
            all_current_points = topology_reconstruction.get_all_current_points()

            # Get the indices of the currently *active* points into all the points.
            #
            # Note: Points that are None represent currently *inactive* points.
            #       These are points that have either:
            #         1) not been activated yet because the current time is older then their appearance time, or
            #         2) have been deactivated because the current time is younger than their dissappearance time, or
            #         3) have been deactivated through collision detection (ie, collided with a continent or trench).
            #       Although the points currently don't have a hardwired disappearance time (it's implicitly '-inf').
            #       Instead we rely on collision detection to deactivate points.
            active_point_indices = np.fromiter(
                (
                    point_index
                    for point_index, point in enumerate(all_current_points)
                    if point is not None
                ),
                dtype=int,
            )

            logger.debug(f"The total number of points is :{len(all_current_points)}")
            logger.debug(f"The number of active points is :{len(active_point_indices)}")

            if len(active_point_indices) > 0:

                # Latitudes and longitudes of the active points.
                active_latitudes = np.empty(len(active_point_indices), dtype=float)
                active_longitudes = np.empty(len(active_point_indices), dtype=float)
                for index, point_index in enumerate(active_point_indices):
                    active_latitudes[index], active_longitudes[index] = (
                        all_current_points[point_index].to_lat_lon()
                    )

                active_seafloor_ages = (
                    all_appearance_times[active_point_indices]
                    - topology_reconstruction.get_current_time()
                )

                active_spreading_rates = all_spreading_rates[active_point_indices]

                active_birth_lats = all_birth_lats[active_point_indices]

                active_point_ids = all_point_ids[active_point_indices]

                # Map the key names of the gridding data to their arrays.
                gridding_input_dict = {
                    self.CURRENT_LONGITUDES_KEY: active_longitudes,
                    self.CURRENT_LATITUDES_KEY: active_latitudes,
                    self.SEAFLOOR_AGE_KEY: active_seafloor_ages,
                    self.SPREADING_RATE_KEY: active_spreading_rates,
                    self.BIRTH_LAT_SNAPSHOT_KEY: active_birth_lats,
                    self.POINT_ID_SNAPSHOT_KEY: active_point_ids,
                }

                gridding_input_filename = self.gridding_input_filepath.format(
                    topology_reconstruction.get_current_time()
                )
                # Saving arrays using *keyword* arguments stores the key as the name of the array (rather than just "arr_0", "arr_1", etc).
                # This way when it's loaded (with 'np.load') the arrays will be indexed by these same keys.
                np.savez_compressed(gridding_input_filename, **gridding_input_dict)

                # save debug file
                if get_debug_level() > 100:
                    logger.debug(
                        f"The max and min values of seafloor age are: {np.max(active_seafloor_ages)} - {np.min(active_seafloor_ages)} ({topology_reconstruction.get_current_time()}Ma)"
                    )

                    #
                    # Save currently active sample points and their seafloor ages and spreading rates as a coverage to GPML.
                    #

                    sample_point_coverage_feature = pygplates.Feature()
                    sample_point_coverage_geometry = pygplates.MultiPointOnSphere(
                        zip(active_latitudes, active_longitudes)
                    )
                    sample_point_coverage_scalars = {
                        pygplates.ScalarType.create_gpml(
                            "SeafloorAge"
                        ): active_seafloor_ages,
                        pygplates.ScalarType.create_gpml(
                            "SpreadingRate"
                        ): active_spreading_rates,
                    }
                    sample_point_coverage_feature.set_geometry(
                        (sample_point_coverage_geometry, sample_point_coverage_scalars)
                    )
                    # Only visible at curent time.
                    sample_point_coverage_feature.set_valid_time(
                        topology_reconstruction.get_current_time()
                        + 0.5 * self._ridge_time_step,
                        topology_reconstruction.get_current_time()
                        - 0.5 * self._ridge_time_step,
                    )  # type: ignore

                    # Write gridded data at `time` to gpmlz.
                    sample_point_basename, _ = os.path.splitext(
                        gridding_input_filename
                    )  # remove extension
                    pygplates.FeatureCollection(sample_point_coverage_feature).write(
                        sample_point_basename + ".gpmlz"
                    )

            if not topology_reconstruction.reconstruct_to_next_time():
                break

            logger.info(
                f"Reconstruction has been done for {topology_reconstruction.get_current_time()}!"
            )

    def lat_lon_z_to_netCDF(
        self,
        zval_name,
        time_arr=None,
        unmasked=False,
        nprocs=None,
    ):
        """Produce a netCDF4 grid of a z-value identified by its ``zval_name`` for a given time range in ``time_arr``.

        Seafloor age can be gridded by passing ``zval_name`` as ``SEAFLOOR_AGE``, and spreading
        rate can be gridded with ``SPREADING_RATE``.

        Saves all grids to compressed netCDF format in the attributed directory. Grids
        can be read into ndarray format using :func:`gplately.read_netcdf_grid()`.

        Parameters
        ----------
        zval_name : str
            A string identifier for the z-value.
        time_arr : list of float, default=None
            A time range to turn lons, lats and z-values into netCDF4 grids. If not provided,
            ``time_arr`` defaults to the full ``time_array`` provided to :class:`SeafloorGrid`.
        unmasked : bool, default=False
            Save unmasked grids, in addition to masked versions.
        nprocs : int, optional
            Optional number of CPUs to use for producing grids.
            If not specified then defaults to the number of CPUs specified in :py:meth:`__init__`.
            If specified then must be an integer or convertible to an integer (eg, float is rounded towards zero).
            If positive then uses that many CPUs.
            If ``1`` then executes in serial (ie, is not parallelized).
            If ``0`` then a ``ValueError`` is raised.
            If ``-1`` then all available CPUs are used.
            If ``-2`` then all available CPUs except one are used, etc.
        """

        # User can put any time array within SeafloorGrid bounds, but if none
        # is provided, it defaults to the attributed time array
        if time_arr is None:
            time_arr = self._times

        if nprocs is None:
            # Use default number of CPUs (specified in '__init__()').
            num_cpus = self.num_cpus
        else:
            # Determine number of CPUs to use.
            num_cpus = _get_num_cpus(nprocs)

        if num_cpus > 1:
            with multiprocessing.Pool(num_cpus) as pool:
                pool.map(
                    partial(
                        _lat_lon_z_to_netCDF_time,
                        zval_name=zval_name,
                        file_collection=self.file_collection,
                        save_directory=self.save_directory,
                        extent=self.extent,
                        resX=self.spacingX,
                        resY=self.spacingY,
                        unmasked=unmasked,
                        continent_mask_filename=self.continent_mask_filepath,
                        gridding_input_filename=self.gridding_input_filepath,
                    ),
                    time_arr,
                )
        else:
            for time in time_arr:
                _lat_lon_z_to_netCDF_time(
                    time=time,
                    zval_name=zval_name,
                    file_collection=self.file_collection,
                    save_directory=self.save_directory,
                    extent=self.extent,
                    resX=self.spacingX,
                    resY=self.spacingY,
                    unmasked=unmasked,
                    continent_mask_filename=self.continent_mask_filepath,
                    gridding_input_filename=self.gridding_input_filepath,
                )


def _lat_lon_z_to_netCDF_time(
    time,
    zval_name,
    file_collection,
    save_directory: Union[str, Path],
    extent,
    resX,
    resY,
    continent_mask_filename: str,
    gridding_input_filename: str,
    unmasked=False,
):
    # Read the gridding input made by ReconstructByTopologies:
    # Use pandas to load in lons, lats and z values from npz files
    with np.load(gridding_input_filename.format(time)) as npz_loaded:
        curr_data = pd.DataFrame.from_dict(
            # The key names become the column headers...
            {column: npz_loaded[column] for column in npz_loaded.files},
            orient="columns",
        )

    # drop invalid data
    curr_data = curr_data.replace([np.inf, -np.inf], np.nan)
    curr_data = curr_data.dropna(
        subset=[
            SeafloorGrid.CURRENT_LONGITUDES_KEY,
            SeafloorGrid.CURRENT_LATITUDES_KEY,
            zval_name,
        ]
    )
    # Drop points with age greater than 350.
    if SeafloorGrid.SEAFLOOR_AGE_KEY == zval_name:
        curr_data = curr_data.drop(curr_data[curr_data.SEAFLOOR_AGE > 350].index)

    # Drop duplicate latitudes and longitudes
    unique_data = curr_data.drop_duplicates(
        subset=[SeafloorGrid.CURRENT_LONGITUDES_KEY, SeafloorGrid.CURRENT_LATITUDES_KEY]
    )

    # Acquire lons, lats and zvalues for each time
    lons = unique_data[SeafloorGrid.CURRENT_LONGITUDES_KEY].to_list()
    lats = unique_data[SeafloorGrid.CURRENT_LATITUDES_KEY].to_list()
    zdata = np.array(unique_data[zval_name].to_list())

    # zdata = np.where(zdata > 375, float("nan"), zdata), to deal with vmax in the future
    # zdata = np.nan_to_num(zdata)

    # Create a regular grid on which to interpolate lats, lons and zdata
    extent_globe = extent
    grid_lon = np.linspace(extent_globe[0], extent_globe[1], resX)
    grid_lat = np.linspace(extent_globe[2], extent_globe[3], resY)
    X, Y = np.meshgrid(grid_lon, grid_lat)

    # Interpolate lons, lats and zvals over a regular grid using nearest neighbour interpolation
    Z = tools.griddata_sphere((lons, lats), zdata, (X, Y), method="nearest")
    Z = Z.astype("f4")  # float32 precision

    unmasked_basename = f"{zval_name}_grid_unmasked_{time:0.2f}Ma.nc"
    grid_basename = f"{zval_name}_grid_{time:0.2f}Ma.nc"
    if file_collection:
        unmasked_basename = f"{file_collection}_{unmasked_basename}"
        grid_basename = f"{file_collection}_{grid_basename}"
    output_dir = os.path.join(save_directory, zval_name)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    grid_output_unmasked = os.path.join(output_dir, unmasked_basename)
    grid_output = os.path.join(output_dir, grid_basename)

    if unmasked:
        grids.write_netcdf_grid(
            grid_output_unmasked,
            Z,
            extent=extent,
            significant_digits=2,
            fill_value=None,
        )

    # Identify regions in the grid in the continental mask
    # We need the continental mask to match the number of nodes
    # in the uniform grid defined above. This is important if we
    # pass our own continental mask to SeafloorGrid
    cont_mask = grids.read_netcdf_grid(
        continent_mask_filename.format(time), resize=(resX, resY)
    )

    # Whether to compress the masked grid, or not (None).
    #
    # Note: Previously this was 2 (enabling compression).
    #       However we want reasonable accuracy for the seafloor age grid
    #       (for example) and setting significant digits to 3 - since ages
    #       can be 3 digits (before the decimal point) - doesn't reduce the
    #       grid file sizes by large amounts. For a 0.1 degree age grid it
    #       reduces from ~800KB to ~480KB - which is still significant but not
    #       by a factor of 2 or 3 or 4 or more as reported in
    #       https://github.com/GPlates/gplately/pull/125 (that was likely for
    #       data, such as topography, that is less smooth than an age grid).
    masked_significant_digits = None
    if masked_significant_digits is None:
        masked_fill_value = np.nan
    else:
        # When compressing, we can't seem to use NaN as a fill value because
        # reading the written grid does not seem to mask out the NaN regions.
        # It was reported in https://github.com/GPlates/gplately/pull/125
        # that 2 significant digits was enough to preserve NaN masks, but
        # it doesn't seem to work for me (using netCDF4 1.7.2 on Windows).
        # Even 7 significant digits doesn't work.
        #
        # Instead we just use the default '_FillValue' used by netCDF4.
        masked_fill_value = 9.969209968386869e36

    # Use the continental mask to mask out continents
    Z[cont_mask.astype(bool)] = masked_fill_value

    grids.write_netcdf_grid(
        grid_output,
        Z,
        extent=extent,
        significant_digits=masked_significant_digits,
        fill_value=masked_fill_value,
    )
    logger.info(f"{zval_name} netCDF grids for {time:0.2f} Ma complete!")


def _build_continental_mask(
    time: float,
    continent_mask_filepath,
    rotation_model,
    reconstructed_continents,
    spacingY,
    spacingX,
    extent,
    overwrite=False,
):
    """Create a continental mask for a given time."""
    mask_fn = continent_mask_filepath.format(time)
    if os.path.isfile(mask_fn) and not overwrite:
        logger.info(
            f"Continent mask file exists and will not create again -- {mask_fn}"
        )
        return

    final_grid = grids.rasterise(
        reconstructed_continents,
        rotation_model,
        key=1.0,
        shape=(spacingY, spacingX),
        extent=extent,
        origin="lower",
    )
    final_grid[np.isnan(final_grid)] = 0.0

    grids.write_netcdf_grid(
        continent_mask_filepath.format(time),
        final_grid.astype("i1"),
        extent=(-180, 180, -90, 90),
        fill_value=None,
    )
    logger.info(f"Finished building a continental mask at {time} Ma!")


def _build_continental_mask_parallel(
    time: float,
    continent_mask_filepath,
    plate_reconstruction,
    plot_topologies,
    spacingY,
    spacingX,
    extent,
    overwrite,
):
    # Get the reconstructed continents at 'time''.
    # But first need to set 'time' on our unpickled 'plot_topologies'.
    #
    # Note: This only affects our unpickled copy of the original 'self.plot_topologies'.
    #       The original 'self.plot_topologies' will likely still have a different 'time'.
    plot_topologies.time = time
    reconstructed_continents = plot_topologies.continents

    _build_continental_mask(
        time=time,
        continent_mask_filepath=continent_mask_filepath,
        rotation_model=plate_reconstruction.rotation_model,
        reconstructed_continents=reconstructed_continents,
        spacingY=spacingY,
        spacingX=spacingX,
        extent=extent,
        overwrite=overwrite,
    )


def _build_continental_mask_with_contouring(
    time: float,
    continent_mask_filepath,
    rotation_model,
    continent_features,
    overwrite=False,
):
    """Build the continent mask for a given time using ptt's 'continent contouring' method.
    For more information about 'Continent Contouring', visit https://github.com/EarthByte/continent-contouring.
    """
    mask_fn = continent_mask_filepath.format(time)
    if os.path.isfile(mask_fn) and not overwrite:
        logger.info(
            f"Continent mask file exists and will not create again -- {mask_fn}"
        )
        return

    continent_contouring_point_spacing_degrees = 0.25
    continent_contouring_area_threshold_square_kms = 0
    continent_contouring_area_threshold_steradians = (
        continent_contouring_area_threshold_square_kms
        / (pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms)
    )
    continent_exclusion_area_threshold_square_kms = 800000
    continent_exclusion_area_threshold_steradians = (
        continent_exclusion_area_threshold_square_kms
        / (pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms)
    )
    continent_separation_distance_threshold_radians = (
        continent_contours.DEFAULT_CONTINENT_SEPARATION_DISTANCE_THRESHOLD_RADIANS
    )

    def continent_contouring_buffer_and_gap_distance_radians(time, contoured_continent):
        # One distance for time interval [1000, 300] and another for time interval [250, 0].
        # And linearly interpolate between them over the time interval [300, 250].
        pre_pangea_distance_radians = math.radians(2.5)  # convert degrees to radians
        post_pangea_distance_radians = math.radians(0.0)  # convert degrees to radians
        if time > 300:
            buffer_and_gap_distance_radians = pre_pangea_distance_radians
        elif time < 250:
            buffer_and_gap_distance_radians = post_pangea_distance_radians
        else:
            # Linearly interpolate between 250 and 300 Ma.
            interp = float(time - 250) / (300 - 250)
            buffer_and_gap_distance_radians = (
                interp * pre_pangea_distance_radians
                + (1 - interp) * post_pangea_distance_radians
            )

        # Area of the contoured continent.
        area_steradians = contoured_continent.get_area()

        # Linearly reduce the buffer/gap distance for contoured continents with area smaller than 1 million km^2.
        area_threshold_square_kms = 500000
        area_threshold_steradians = area_threshold_square_kms / (
            pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms
        )
        if area_steradians < area_threshold_steradians:
            buffer_and_gap_distance_radians *= (
                area_steradians / area_threshold_steradians
            )

        return buffer_and_gap_distance_radians

    continent_contouring = continent_contours.ContinentContouring(
        rotation_model,
        continent_features,
        continent_contouring_point_spacing_degrees,
        continent_contouring_area_threshold_steradians,
        continent_contouring_buffer_and_gap_distance_radians,
        continent_exclusion_area_threshold_steradians,
        continent_separation_distance_threshold_radians,
    )
    continent_mask, _ = (
        continent_contouring.get_continent_mask_and_contoured_continents(time)
    )
    grids.write_netcdf_grid(
        continent_mask_filepath.format(time),
        continent_mask.astype("i1"),
        extent=(-180, 180, -90, 90),
        fill_value=None,
    )
    logger.warning(
        f"Finished building a continental mask at {time} Ma using ptt's 'Continent Contouring'!"
        + " For more information about 'Continent Contouring', visit https://github.com/EarthByte/continent-contouring."
    )


# TODO: When pyGPlates can pickle pygplates.RotationModel and pygplates.FeatureCollection noticeably faster (than it currently does in pyGPlates 1.0),
#       then remove this function definition.
def _build_continental_mask_with_contouring_parallel(
    time: float,
    continent_mask_filepath,
    plate_reconstruction,
    plot_topologies,
    overwrite,
):
    _build_continental_mask_with_contouring(
        time=time,
        continent_mask_filepath=continent_mask_filepath,
        rotation_model=plate_reconstruction.rotation_model,
        continent_features=plot_topologies._continents,
        overwrite=overwrite,
    )


def _generate_mid_ocean_ridge_points(
    time: float,
    delta_time: float,
    mid_ocean_ridges_file_path: str,
    rotation_model,
    topology_features,
    ridge_sampling,
    overwrite=True,
):
    """generate middle ocean ridge seed points at a given time"""
    mor_fn = mid_ocean_ridges_file_path.format(time)
    if os.path.isfile(mor_fn) and not overwrite:
        logger.info(
            f"Middle ocean ridge file exists and will not create again.\n{mor_fn}"
        )
        return

    # Points and their z values that emerge from MORs at this time.
    shifted_mor_points = []
    point_spreading_rates = []

    # Resolve topologies to the current time.
    topological_snapshot = pygplates.TopologicalSnapshot(
        topology_features, rotation_model, time
    )
    shared_boundary_sections = topological_snapshot.get_resolved_topological_sections()

    # pygplates.ResolvedTopologicalSection objects.
    for shared_boundary_section in shared_boundary_sections:
        if (
            shared_boundary_section.get_feature().get_feature_type()
            == pygplates.FeatureType.gpml_mid_ocean_ridge
        ):
            spreading_feature = shared_boundary_section.get_feature()

            # Find the stage rotation of the spreading feature in the
            # frame of reference of its geometry at the current
            # reconstruction time (the MOR is currently actively spreading).
            # The stage pole can then be directly geometrically compared
            # to the *reconstructed* spreading geometry.
            stage_rotation = separate_ridge_transform_segments.get_stage_rotation_for_reconstructed_geometry(
                spreading_feature, rotation_model, time
            )
            if not stage_rotation:
                # Skip current feature - it's not a spreading feature.
                continue

            # Get the stage pole of the stage rotation.
            # Note that the stage rotation is already in frame of
            # reference of the *reconstructed* geometry at the spreading time.
            stage_pole, _ = stage_rotation.get_euler_pole_and_angle()

            # One way rotates left and the other right, but don't know
            # which - doesn't matter in our example though.
            rotate_slightly_off_mor_one_way = pygplates.FiniteRotation(
                stage_pole, np.radians(0.01)
            )
            rotate_slightly_off_mor_opposite_way = (
                rotate_slightly_off_mor_one_way.get_inverse()
            )

            # Iterate over the shared sub-segments.
            for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                # Tessellate MOR section.
                mor_points = pygplates.MultiPointOnSphere(
                    shared_sub_segment.get_resolved_geometry().to_tessellated(
                        np.radians(ridge_sampling)
                    )
                )

                coords = mor_points.to_lat_lon_list()
                lats = [i[0] for i in coords]
                lons = [i[1] for i in coords]
                boundary_feature = shared_boundary_section.get_feature()
                left_plate = boundary_feature.get_left_plate(None)
                right_plate = boundary_feature.get_right_plate(None)
                if left_plate is not None and right_plate is not None:
                    # Get the spreading rates for all points in this sub segment
                    (
                        spreading_rates,
                        _,
                    ) = tools.calculate_spreading_rates(
                        time=time,
                        lons=lons,
                        lats=lats,
                        left_plates=[left_plate] * len(lons),
                        right_plates=[right_plate] * len(lons),
                        rotation_model=rotation_model,
                        delta_time=delta_time,
                    )

                else:
                    spreading_rates = [np.nan] * len(lons)

                # Loop through all but the 1st and last points in the current sub segment
                for point, rate in zip(
                    mor_points.get_points()[1:-1],
                    spreading_rates[1:-1],
                ):
                    # Add the point "twice" to the main shifted_mor_points list; once for a L-side
                    # spread, another for a R-side spread. Then add the same spreading rate twice
                    # to the list - this therefore assumes spreading rate is symmetric.
                    shifted_mor_points.append(rotate_slightly_off_mor_one_way * point)
                    shifted_mor_points.append(
                        rotate_slightly_off_mor_opposite_way * point
                    )
                    point_spreading_rates.extend([rate] * 2)

    # save the middle ocean ridges points to .pkl file
    lats_lons = [point.to_lat_lon() for point in shifted_mor_points]
    df = pd.DataFrame(
        lats_lons,
        columns=[
            SeafloorGrid.CURRENT_LATITUDES_KEY,
            SeafloorGrid.CURRENT_LONGITUDES_KEY,
        ],
    )
    df[SeafloorGrid.SPREADING_RATE_KEY] = point_spreading_rates
    df.to_pickle(mor_fn)

    if get_debug_level() > 100:
        #
        # Save newly created mid-ocean ridge points and their spreading rates as a coverage to GPML.
        #

        mor_coverage_feature = pygplates.Feature()
        mor_coverage_geometry = pygplates.MultiPointOnSphere(shifted_mor_points)
        mor_coverage_scalars = {
            pygplates.ScalarType.create_gpml("SpreadingRate"): point_spreading_rates,
        }
        mor_coverage_feature.set_geometry((mor_coverage_geometry, mor_coverage_scalars))
        # Only visible at curent time.
        mor_coverage_feature.set_valid_time(
            time + 0.5 * delta_time,
            time - 0.5 * delta_time,
        )  # type: ignore

        # Write MOR points at `time` to gpmlz
        mor_basename, _ = os.path.splitext(mor_fn)  # remove extension
        pygplates.FeatureCollection(mor_coverage_feature).write(mor_basename + ".gpmlz")

    logger.info(f"Finished building MOR seedpoints at {time} Ma!")


# TODO: When pyGPlates can pickle pygplates.RotationModel and pygplates.FeatureCollection noticeably faster (than it currently does in pyGPlates 1.0),
#       then remove this function definition.
def _generate_mid_ocean_ridge_points_parallel(
    time: float,
    delta_time: float,
    mid_ocean_ridges_file_path: str,
    plate_reconstruction,
    ridge_sampling,
    overwrite,
):
    _generate_mid_ocean_ridge_points(
        time=time,
        delta_time=delta_time,
        mid_ocean_ridges_file_path=mid_ocean_ridges_file_path,
        rotation_model=plate_reconstruction.rotation_model,
        topology_features=plate_reconstruction.topology_features,
        ridge_sampling=ridge_sampling,
        overwrite=overwrite,
    )


def _get_num_cpus(nprocs):
    """Return number of CPUs to use.

    Parameters
    ----------
    nprocs : int
        The number of CPUs to use for parts of the code that are parallelized.
        Must be an integer or convertible to an integer (eg, float is rounded towards zero).
        If positive then uses that many CPUs.
        If ``1`` then executes in serial (ie, is not parallelized).
        If ``0`` then a ``ValueError`` is raised.
        If ``-1`` then all available CPUs are used.
        If ``-2`` then all available CPUs except one are used, etc.
    """

    try:
        nprocs = int(nprocs)
    except ValueError:
        raise TypeError('"nprocs" should be an integer, or convertible to integer')

    if nprocs == 0:
        raise ValueError('"nprocs" should not be zero')

    if nprocs > 0:
        # A positive integer specifying the number of CPUs to use.
        num_cpus = nprocs
    else:  # nprocs < 0
        #
        # A negative integer specifying the number of CPUs to NOT use.
        # '-1' means use all CPUs. '-2' means use all CPUs but one. Etc.
        try:
            num_cpus = multiprocessing.cpu_count() + 1 + nprocs
            # If specified more CPUs to NOT use than there are CPUs available.
            if num_cpus < 1:
                num_cpus = 1
        except NotImplementedError:
            num_cpus = 1

    return num_cpus
