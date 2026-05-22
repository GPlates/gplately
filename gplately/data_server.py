#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
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

"""DataServer class for downloading plate-model assets and rasters."""

import logging
from pathlib import Path
from typing import Union

import numpy as np

from .grids import Raster

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false
import pygplates
from matplotlib import image
from plate_model_manager import PlateModelManager, PresentDayRasterManager
from plate_model_manager.utils.download import FileDownloader
import pooch

logger = logging.getLogger("gplately")


class DataServer(object):
    """
    Download the plate reconstruction models from the `EarthByte server <https://repo.gplates.org/webdav/pmm/>`__.

    The :class:`DataServer` object downloads the model files to the ``GPlately cache folder``.
    If the same model is requested again, a new :class:`DataServer` instance will retrieve the files from the cache --
    provided they haven't been moved or deleted.

    .. seealso::

        - `This table <../plate_models.html>`__ provides a list of available plate reconstruction models.
        - Visit this `EarthByte web page <https://www.earthbyte.org/category/resources/data-models/global-regional-plate-motion-models/>`__ for more information about these plate models.
        - Call :meth:`gplately.auxiliary.get_data_server_cache_path` to see the path to the ``GPlately cache folder``.
    """

    @staticmethod
    def _path_to_cache():
        return pooch.utils.cache_location(
            pooch.os_cache("gplately"), env=None, version=None
        )

    def __init__(self, file_collection, data_dir=None, verbose=True):
        """Constructor. Create a :class:`DataServer` object.

        Example
        -------
        .. code-block:: python

            # create a DataServer object for the Cao2024 model (https://zenodo.org/records/11536686)
            data_server = gplately.DataServer("Cao2024")

        Parameters
        ----------
        file_collection: str
            The model name of interest.

        verbose: bool, default=True
            Toggle print messages regarding server/internet connection status, file availability, etc.
        """

        if not data_dir:
            _data_dir = self._path_to_cache()
        else:
            _data_dir = data_dir

        self.file_collection = file_collection.capitalize()
        self.pmm = PlateModelManager().get_model(
            self.file_collection, data_dir=str(_data_dir)
        )
        if not self.pmm:
            raise Exception(
                f"Unable to get plate model {self.file_collection}. Check if the model name is correct."
            )
        self._available_layers = self.pmm.get_avail_layers()
        self.verbose = verbose

        # initialise empty attributes
        self._rotation_model = None
        self._topology_features = None
        self._static_polygons = None
        self._coastlines = None
        self._continents = None
        self._COBs = None

    def _create_feature_collection(self, file_list):
        feature_collection = pygplates.FeatureCollection()
        for feature in file_list:
            feature_collection.add(pygplates.FeatureCollection(feature))
        return feature_collection

    @staticmethod
    def cache_path():
        """The location of DataServer cache on your computer."""
        return DataServer._path_to_cache()

    @property
    def rotation_model(self):
        """A pygplates.RotationModel object for the plate reconstruction model."""
        if self._rotation_model is None and self.pmm:
            self._rotation_model = pygplates.RotationModel(
                self.pmm.get_rotation_model()
            )
            self._rotation_model.reconstruction_identifier = self.file_collection
            # Setting an attribute on a pyGPlates object produces the following error in version 1.0 of pyGPlates:
            #   RuntimeError: Incomplete pickle support (__getstate_manages_dict__ not set)
            #
            # This is fixed in pyGPlates 1.1 (which implements __getstate__ to copy __dict__ just to be sure),
            # but until that's released we can just set __getstate_manages_dict__ to True.
            #
            # This is because it turns out that Boost.Python (used in pyGPlates) copies the __dict__ in its default __getstate__
            # so we can just manually set __getstate_manages_dict__ to True (I'm not sure why Boost.Python doesn't set it).
            self._rotation_model.__getstate_manages_dict__ = True
        return self._rotation_model

    @property
    def topology_features(self):
        """A pygplates.FeatureCollection object containing topology features."""
        if self._topology_features is None and self.pmm:
            if "Topologies" in self._available_layers:
                self._topology_features = self._create_feature_collection(
                    self.pmm.get_topologies()
                )
            else:
                self._topology_features = pygplates.FeatureCollection()
        return self._topology_features

    @property
    def static_polygons(self):
        """A pygplates.FeatureCollection object containing static polygons."""
        if self._static_polygons is None and self.pmm:
            if "StaticPolygons" in self._available_layers:
                self._static_polygons = self._create_feature_collection(
                    self.pmm.get_static_polygons()
                )
            else:
                self._static_polygons = pygplates.FeatureCollection()
        return self._static_polygons

    @property
    def coastlines(self):
        """A pygplates.FeatureCollection object containing coastlines."""
        if self._coastlines is None and self.pmm:
            if "Coastlines" in self._available_layers:
                self._coastlines = self._create_feature_collection(
                    self.pmm.get_coastlines()
                )
            else:
                self._coastlines = pygplates.FeatureCollection()
        return self._coastlines

    @property
    def continents(self):
        """A pygplates.FeatureCollection object containing continental polygons."""
        if self._continents is None and self.pmm:
            if "ContinentalPolygons" in self._available_layers:
                self._continents = self._create_feature_collection(
                    self.pmm.get_continental_polygons()
                )
            else:
                self._continents = pygplates.FeatureCollection()
        return self._continents

    @property
    def COBs(self):
        """A pygplates.FeatureCollection object containing continent-ocean boundaries."""
        if self._COBs is None and self.pmm:
            if "COBs" in self._available_layers:
                self._COBs = self._create_feature_collection(self.pmm.get_COBs())
            else:
                self._COBs = pygplates.FeatureCollection()
        return self._COBs

    @property
    def from_age(self) -> float:
        """The max age/time of the plate model."""
        if self.pmm:
            return self.pmm.get_big_time()
        else:
            raise Exception(
                "Unable to get max reconstruction age/time. Check the PlateModel object."
            )

    @property
    def to_age(self) -> float:
        """The min age/time of the plate model."""
        if self.pmm:
            return self.pmm.get_small_time()
        else:
            raise Exception(
                "Unable to get min reconstruction age/time. Check the PlateModel object."
            )

    @property
    def time_range(self):
        """Deprecated!!! Use :attr:`DataServer.valid_time` instead.
        Keep consistent with `GML naming <https://www.gplates.org/docs/gpgim/#gml:validTime>`__.
        """
        return self.from_age, self.to_age

    @property
    def valid_times(self):
        """Deprecated!!! Use :attr:`DataServer.valid_time` instead.
        Keep consistent with `GML naming <https://www.gplates.org/docs/gpgim/#gml:validTime>`__.
        """
        return self.from_age, self.to_age

    @property
    def valid_time(self):
        """The period of time the plate model are valid. Return a tuple of (max time, min time)."""
        return self.from_age, self.to_age

    def get_plate_reconstruction_files(self):
        """Download and return a tuple of **rotation_model**, **topology_features** and **static_polygons**.
        These objects can then be used to create :class:`gplately.PlateReconstruction` object.

        Returns
        -------
        rotation_model : pygplates.RotationModel
            A rotation model to query equivalent and/or relative topological plate rotations
            from a time in the past relative to another time in the past or to present day.
        topology_features : pygplates.FeatureCollection
            Topological features including ridges, transforms, subduction zones, etc.
            These features can be used to build topological plate boundaries and networks.
        static_polygons : pygplates.FeatureCollection
            Static polygons which can be used to assign plate IDs for other geometries.
            The plate IDs are essential to tectonic plate reconstruction.


        .. note::

            The example code below downloads ``rotation model``, ``topology features`` and ``static polygons`` files from the
            Muller et al. (2019) plate reconstruction model and create a :class:`gplately.PlateReconstruction` object.

            .. code-block:: python
                :linenos:

                import gplately

                data_server = gplately.DataServer("Muller2019")
                rotation_model, topology_features, static_polygons = (
                    data_server.get_plate_reconstruction_files()
                )

                # create a PlateReconstruction object using the returned objects
                model = gplately.PlateReconstruction(
                    rotation_model, topology_features, static_polygons
                )

            If the requested plate model does not have certain file(s), warning messages will alert user of the missing file(s).
        """

        return self.rotation_model, self.topology_features, self.static_polygons

    def get_topology_geometries(self):
        """Download and return coastlines, continental polygons and COBs (continent-ocean boundary).
        These feature collections can be used to create :class:`gplately.PlotTopologies` object and plot paleomaps.

        Returns
        -------
        coastlines : pygplates.FeatureCollection
            Global coastlines. These coastlines have been assigned plate IDs using static polygons and are ready to
            be reconstructed to a particular geological time.

        continents : pygplates.FeatureCollection
            Continental polygons containing continental crust and volcanically-modified oceanic crust (including island arcs).

        COBs : pygplates.FeatureCollection
            Continent-ocean boundary. The COBs are represented as lines along passive margins and does not include data from active margins.


        .. note::

            The example code below will attempt to download ``coastlines``, ``continents`` and ``COBs`` from the Muller
            et al. (2019) plate reconstruction model and create a :class:`gplately.PlotTopologies` object.

            .. code-block:: python
                :linenos:
                :emphasize-lines: 9, 12

                data_server = gplately.download.DataServer("Muller2019")
                rotation_model, topology_features, static_polygons = (
                    data_server.get_plate_reconstruction_files()
                )
                model = gplately.PlateReconstruction(
                    rotation_model, topology_features, static_polygons
                )

                coastlines, continents, COBs = data_server.get_topology_geometries()

                # create a gplately.PlotTopologies object at 100Ma
                gPlot = gplately.PlotTopologies(model, 100, continents, coastlines, COBs)

            If the requested plate model does not have certain geometries, warning messages will be printed to alert the user.
        """

        return self.coastlines, self.continents, self.COBs

    def get_age_grid(self, times: Union[int, list[int]]):
        """Download the seafloor age grids for the plate model. Save the grids in the ``GPlately cache folder``.

        .. seealso::

            :meth:`DataServer.cache_path`

        The available seafloor age grids are listed below.

        * Muller et al. 2019

            * ``file_collection`` = ``Muller2019``
            * Time range: 0-250 Ma
            * Seafloor age grids in netCDF format.

        * Muller et al. 2016

            * ``file_collection`` = ``Muller2016``
            * Time range: 0-240 Ma
            * Seafloor age grids in netCDF format.

        * Seton et al. 2012

            * ``file_collection`` = ``Seton2012``
            * Time range: 0-200 Ma
            * Seafloor age grids in netCDF format.


        Parameters
        ----------
        times : int, or a list of int
            A reconstruction time or a list of reconstruction times.

        Returns
        -------
        :class:`gplately.Raster` or a list of :class:`gplately.Raster`


        .. note::

            The age grid data can be accessed as a numpy ndarray or MaskedArray via the :attr:`gplately.Raster.data` attribute.

            For example:

            .. code-block:: python
                :linenos:

                data_server = gplately.DataServer("Muller2019")
                graster = data_server.get_age_grid(100)
                graster_data = graster.data

            where ``graster_data`` is a numpy ndarray.

        Raises
        ------
        ValueError
            If the ``times`` parameter contains invalid reconstruction time.


        Example
        -------
        To download  Muller et al. (2019) seafloor age grids for 0Ma, 1Ma and 100 Ma:

            .. code-block:: python
                :linenos:

                data_server = gplately.DataServer("Muller2019")
                age_grids = data_server.get_age_grid([0, 1, 100])

        .. seealso::

            - :meth:`PlateModel.get_raster()`
            - :meth:`PlateModel.get_rasters()`

            .. code-block:: python
                :linenos:
                :emphasize-lines: 4,5

                from gplately import PlateModelManager

                model = PlateModelManager().get_model("Muller2019")
                print(model.get_rasters("AgeGrids", times=[10, 20, 30]))
                print(model.get_raster("AgeGrids", time=100))

        """
        if not self.pmm:
            raise Exception("The plate model object is None. Unable to get age grid.")

        if "AgeGrids" not in self.pmm.get_cfg()["TimeDepRasters"]:
            raise ValueError(
                f"The time-dependent seafloor age grids are not currently available for {self.file_collection}."
            )

        age_grids = []
        for time in np.atleast_1d(times):
            try:
                time_i = int(time)
            except Exception:
                raise ValueError(
                    f"Invalid time {time}. Reconstruction time must be a number."
                )
            if time_i < self.to_age or time_i > self.from_age:
                raise ValueError(
                    f"Invalid time {time}. Reconstruction time must be between {self.time_range}."
                )
            age_grids.append(Raster(data=self.pmm.get_raster("AgeGrids", time_i)))

        if not age_grids:
            raise Exception(f"Unable to get the seafloor age grids for times: {times}")

        if len(age_grids) == 1:
            return age_grids[0]
        else:
            return age_grids

    def get_spreading_rate_grid(self, times):
        """Download seafloor spreading rate grids from the plate reconstruction model and save the grids in the ``GPlately cache folder``.

        .. seealso::

            :meth:`DataServer.cache_path`

        The available seafloor spreading rate grids are listed below.

        * Clennett et al. 2020

            * `file_collection` = `Clennett2020`
            * Time range: 0-250 Ma
            * Seafloor spreading rate grids in netCDF format.


        Parameters
        ----------
        time : int, or list of int
            Request spreading grid(s) for one (an integer) or multiple reconstruction times (a list of integers).

        Returns
        -------
        :class:`gplately.Raster` or a list of :class:`gplately.Raster`
        """

        if not self.pmm:
            raise Exception(
                "The plate model object is None. Unable to get spreading rate grids."
            )

        if "SpreadingRateGrids" not in self.pmm.get_cfg()["TimeDepRasters"]:
            raise ValueError(
                "The time-dependant SpreadingRateGrids are not currently available for {}".format(
                    self.file_collection
                )
            )

        spread_grids = []
        for time in np.atleast_1d(times):
            try:
                time_i = int(time)
            except Exception:
                raise ValueError(
                    f"Invalid time {time}. Reconstruction time must be a number."
                )
            if time_i < self.to_age or time_i > self.from_age:
                raise ValueError(
                    f"Invalid time {time}. Reconstruction time must be between {self.time_range}."
                )
            spread_grids.append(
                Raster(data=self.pmm.get_raster("SpreadingRateGrids", time_i))
            )

        if not spread_grids:
            raise Exception(
                f"Unable to get the seafloor spreading rate grids for times: {times}"
            )

        if len(spread_grids) == 1:
            return spread_grids[0]
        else:
            return spread_grids

    def get_valid_times(self):
        """Deprecated!!! Use :attr:`DataServer.valid_times` instead.
        Return a tuple (max_time, min_time) representing the valid time range of the plate model.
        """
        return self.from_age, self.to_age

    @staticmethod
    def get_raster(raster_name: str):
        """Download rasters that are not associated with any plate reconstruction models. Store the rasters in the ``GPlately cache``.

        The available present-day rasters are listed below.

        * `ETOPO1 <https://www.ngdc.noaa.gov/mgg/global/>`__
            * Filetypes available : TIF, netCDF (GRD)
            * `raster_name` = ``ETOPO1_grd``, ``ETOPO1_tif`` (depending on the requested format)
            * A 1-arc minute global relief model combining lang topography and ocean bathymetry.
            * Citation: doi:10.7289/V5C8276M


        Parameters
        ----------
        raster_name : :class:`str`
            The raster name of interest.

        Returns
        -------
        :class:`gplately.Raster`
            A :class:`gplately.Raster` object containing the raster data which can be accessed as
            a ``numpy ndarray`` or ``MaskedArray`` via :attr:`gplately.Raster.data` attribute.

            For example:

            .. code-block:: python
                :linenos:

                graster = gplately.DataServer.get_raster("ETOPO1_tif")
                graster_data = graster.data

            where ``graster_data`` is a ``numpy ndarray``. This array can be visualised using ``matplotlib.pyplot.imshow`` (see example below).

        Raises
        ------
        Exception
            Raise ``Exception`` when ``raster_name`` is invalid.


        Example
        -------
        Download ETOPO1 and plot it on a map with Mollweide projection.

        .. code-block:: python
            :linenos:

            import cartopy.crs as ccrs
            import matplotlib.pyplot as plt

            import gplately

            etopo1 = gplately.DataServer.get_raster("ETOPO1_tif")
            fig = plt.figure(figsize=(18, 14), dpi=300)
            ax = fig.add_subplot(111, projection=ccrs.Mollweide(central_longitude=-150))
            ax.imshow(etopo1.data, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        """
        if raster_name:
            raster_path = PresentDayRasterManager(
                data_dir=str(DataServer._path_to_cache())
            ).get_raster(raster_name)
            if raster_path.endswith(".grd") or raster_path.endswith(".nc"):
                raster = Raster(data=raster_path)
            else:
                # Otherwise, the raster is an image; use imread to process
                raster = Raster(data=image.imread(raster_path))

            if raster_name.lower() == "etopo1_tif":
                raster.lats = raster.lats[::-1]
            if raster_name.lower() == "etopo1_grd":
                raster._data = raster._data.astype(float)  # type: ignore
            return raster
        else:
            raise Exception("The 'raster_name' parameter is required!")

    @staticmethod
    def get_feature_data(feature_data_id_string: str):
        """Download and load supported geological feature datasets.

        This method fetches a dataset from the EarthByte GPlates sample-data
        server, caches it under the GPlately cache directory, and returns one
        or more ``pygplates.FeatureCollection`` objects from matching files.

        Supported dataset identifiers:

        - ``SeafloorFabric``: Seafloor tectonic fabric data (Matthews et al., 2011).
        - ``Johansson2018``: Large igneous provinces (Johansson et al., 2018).
        - ``Whittaker2015``: Igneous province products interpreted as plume products (Whittaker et al., 2015).
        - ``Hotspots``: Present-day hotspot/plume locations (Whittaker et al., 2013).

        Detailed descriptions can be found below:

        - Large igneous provinces from Johansson et al. (2018)

            Information:
                * Formats: .gpmlz
                * `feature_data_id_string` = `Johansson2018`

            Citations:
                Johansson, L., Zahirovic, S., and Müller, R. D., In Prep, The
                interplay between the eruption and weathering of Large Igneous Provinces and
                the deep-time carbon cycle: Geophysical Research Letters.


        - Large igneous province products interpreted as plume products from Whittaker et al. (2015)

            Information:
                * Formats: .gpmlz, .shp
                * `feature_data_id_string` = `Whittaker2015`

            Citations:
                Whittaker, J. M., Afonso, J. C., Masterton, S., Müller, R. D.,
                Wessel, P., Williams, S. E., & Seton, M. (2015). Long-term interaction between
                mid-ocean ridges and mantle plumes. Nature Geoscience, 8(6), 479-483.
                doi:10.1038/ngeo2437.


        - Seafloor tectonic fabric (fracture zones, discordant zones, V-shaped structures, unclassified V-anomalies, propagating ridge lineations and extinct ridges) from Matthews et al. (2011)

            Information:
                * Formats: .gpmlz
                * `feature_data_id_string` = `SeafloorFabric`

            Citations:
                Matthews, K.J., Müller, R.D., Wessel, P. and Whittaker, J.M., 2011. The
                tectonic fabric of the ocean basins. Journal of Geophysical Research, 116(B12):
                B12109, DOI: 10.1029/2011JB008413.


        - Present day surface hotspot/plume locations from Whittaker et al. (2013)

            Information:
                * Formats: .gpmlz
                * `feature_data_id_string` = `Hotspots`

            Citation:
                Whittaker, J., Afonso, J., Masterton, S., Müller, R., Wessel, P.,
                Williams, S., and Seton, M., 2015, Long-term interaction between mid-ocean ridges and
                mantle plumes: Nature Geoscience, v. 8, no. 6, p. 479-483, doi:10.1038/ngeo2437.

        For examples of plotting downloaded feature data, see
        `notebook 05 <https://github.com/GPlates/gplately/blob/master/Notebooks/05-WorkingWithFeatureGeometries.ipynb>`__.

        Parameters
        ----------
        feature_data_id_string : str
            Dataset identifier to download and load.

        Returns
        -------
        pygplates.FeatureCollection or list[pygplates.FeatureCollection]
            A single feature collection if exactly one matching file is found,
            otherwise a list of feature collections (one per matching file).

        Raises
        ------
        ValueError
            If ``feature_data_id_string`` is empty or not supported.
        Exception
            If no files matching the dataset pattern are found after download.


        .. note::

            Files are only downloaded when local cached data are missing or stale.
        """

        database = {
            "SeafloorFabric": (
                "https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/SeafloorFabric.zip",
                "*.gpmlz",
            ),
            "Johansson2018": (
                "https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/IgneousProvinces.zip",
                "Johansson*.gpmlz",
            ),
            "Whittaker2015": (
                "https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/IgneousProvinces.zip",
                "Whittaker*.gpmlz",
            ),
            "Hotspots": (
                "https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/Hotspots.zip",
                "Hotspots*.gpmlz",
            ),
        }

        if not feature_data_id_string:
            raise ValueError(
                "The 'feature_data_id_string' parameter is required to specify which feature data to download."
            )
        if feature_data_id_string not in database:
            raise ValueError(
                f"Invalid 'feature_data_id_string' parameter. Supported values are: {list(database.keys())}"
            )
        file_url, file_pattern = database[feature_data_id_string]
        file_path = f"{DataServer._path_to_cache()}/{feature_data_id_string}"
        logger.info(file_path)
        downloader = FileDownloader(file_url, f"{file_path}/.metadata.json", file_path)

        # only re-download when necessary
        if downloader.check_if_file_need_update():
            downloader.download_file_and_update_metadata()
        else:
            if downloader.check_if_expire_date_need_update():
                downloader.update_metadata()
            else:
                logger.debug(
                    f"The local files in {file_path} are still good. Will not download again at this moment."
                )

        feature_files = list(Path(file_path).rglob(file_pattern))
        if len(feature_files) == 0:
            raise Exception(
                f"No files matching the pattern '{file_pattern}' were found in the downloaded data for '{feature_data_id_string}'."
            )
        elif len(feature_files) == 1:
            return pygplates.FeatureCollection(str(feature_files[0]))
        else:
            fcs = []
            for f in feature_files:
                fcs.append(pygplates.FeatureCollection(str(f)))
            return fcs
