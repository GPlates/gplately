import pooch
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
from pooch import Decompress as _Decompress
import glob, os
import pygplates
import re
import numpy as np

def fetch_from_web(url):
    """Download file(s) in a given url to the 'gplately' cache folder. Processes
    compressed files using either Pooch's Unzip (if .zip) or Decompress (if .gz, 
    .xz or .bz2)."""
    def pooch_retrieve(url, processor):
        """Downloads file(s) from a URL using Pooch."""
        fnames = _retrieve(
            url=url,
            known_hash=None,  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gplately'),
            processor=processor)
        return fnames

    archive_formats = tuple([".gz", ".xz", ".bz2"])
    if url.endswith(".zip"):
        fnames = pooch_retrieve(url, processor=_Unzip())
    elif url.endswith(archive_formats):
        fnames = pooch_retrieve(url, processor=_Decompress())
    else:
        fnames = pooch_retrieve(url, processor=None)
    return fnames


def collect_file_extension(fnames, file_extension):
    """Searches cached directory for filenames with a specified extension(s)."""
    sorted_fnames = []
    file_extension=tuple(file_extension)
    for file in fnames:
        if file.endswith(file_extension):
            sorted_fnames.append(file)
    return sorted_fnames


def str_in_folder(fnames, strings_to_include=None, strings_to_ignore=None):
    """Filter though files with/without """
    sorted_fnames = []
    for i, fname in enumerate(fnames):
        parent_directory = fname.split("/")[:-1]
        if strings_to_ignore is not None:
            check = [s for s in strings_to_ignore if s in parent_directory]
            if check:
                continue
        if strings_to_include is not None:
            if any(x.lower() in fname.lower() for x in strings_to_include):
                sorted_fnames.append(fname)
        else:
            sorted_fnames.append(fname)
    if sorted_fnames:
        return sorted_fnames
    else:
        return fnames
    

def str_in_filename(fnames, strings_to_include=None, strings_to_ignore=None):
    sorted_fnames = []
    for i, fname in enumerate(fnames):
        filename = fname.split("/")[-1]
        if strings_to_ignore is not None:
            check = [s for s in strings_to_ignore if s in fname]
            if check:
                continue
        if strings_to_include is not None:
            if any(x.lower() in filename.lower() for x in strings_to_include):
                sorted_fnames.append(fname)
        else:
            sorted_fnames.append(fname)
    if sorted_fnames:
        return sorted_fnames
    else:
        return fnames


def check_gpml_or_shp(fnames):
    """For topology features, returns GPML by default. Searches for ESRI Shapefiles 
    instead if GPML files not found."""
    sorted_fnames = []
    for file in fnames:
        if file.endswith(".gpml"):
            sorted_fnames.append(file)
    if not sorted_fnames:
        for file in fnames:
            if file.endswith(".shp"):
                sorted_fnames.append(file)
    return sorted_fnames


def _remove_hash(fname):
    """Removes hashes (32 character file IDs) from cached filenames."""
    split_paths = fname.split("-")
    cache_path = split_paths[0][:-32]
    new_path = cache_path + "-".join(split_paths[1:])
    return new_path


def order_filenames_by_time(fnames):
    """Orders filenames in a list from present day to deeper geological time if they
    are labelled by time."""
    # Collect all digits in each filename.
    filepath_digits=[]
    for i, file in enumerate(fnames):
        digits = []
        for element in re.split('([0-9]+)', _remove_hash(file)):
            if element.isdigit():
                digits.append(int(str(element)))
        filepath_digits.append(digits)

    # Ignore digits common to all full file paths. This leaves behind the files' 
    # geological time label.
    geological_times = []
    filepath_digits = np.array(filepath_digits).T
    for digit_array in filepath_digits:
        if not all(digit == digit_array[0] for digit in digit_array):
            geological_times.append(digit_array)

    # If files have geological time labels, allocate indices to the current filename order, 
    # and sort files from recent to deep geological time.
    if geological_times:
        sorted_geological_times = sorted(
            enumerate(geological_times[0]), 
            key=lambda x: x[1]
        )
        sorted_geological_time_indices = [geo_time[0] for geo_time in sorted_geological_times]
        filenames_sorted = [fnames[index] for index in sorted_geological_time_indices]
    else:
        # If given filenames do not have a time label, return them as is.
        filenames_sorted = fnames
    return filenames_sorted


def collection_sorter(fnames, string_identifier):
    """If multiple file collections or plate reconstruction models are downloaded from
    a single zip folder, only return the needed model."""
    studyname = re.findall(r'[A-Za-z]+|\d+', string_identifier)[0]
    newfnames = []
    for files in fnames:
        if studyname not in files:
            continue
        newfnames.append(files)
    return newfnames


def ignore_macOSX(fnames):
    """For Mac users: filters out duplicate filenames extracted from the __MACOSX folder."""
    for fname in fnames:
        if fname.find("__MACOSX") != -1:
            fnames.remove(fname)
    return fnames


class DataServer(object):
    """Uses Pooch to download geological feature data from plate reconstruction models and other studies
    that are stored on web servers (e.g. EarthByte's webDAV server). Downloaded files are kept in
    a 'gplately' cache folder. 

    Currently,DataServer supports the following plate reconstruction models:
    +-----------------------+-----------------------------------------------------------+-------------------+
    | Plate reconstruction  | Paper citation                                            | String identifier |
    | model                 |                                                           |                   |
    +-----------------------+-----------------------------------------------------------+-------------------+
    | Muller et al. 2019    | Müller, R. D., Zahirovic, S., Williams, S. E.,            | "Muller2019"      |
    |                       | Cannon, J., Seton, M., Bower, D. J., Tetley, M. G.,       |                   |
    |                       | Heine, C., Le Breton, E., Liu, S., Russell, S. H. J.,     |                   |
    |                       | Yang, T., Leonard, J., and Gurnis, M. (2019),             |                   |
    |                       | A global plate model including lithospheric deformation   |                   |
    |                       | along major rifts and orogens since the Triassic.         |                   |
    |                       | Tectonics, vol. 38, https://doi.org/10.1029/2018TC005462. |                   |
    +-----------------------+-----------------------------------------------------------+-------------------+
    | Merdith et al. 2021   | Andrew S. Merdith, Simon E. Williams, Alan S. Collins,    | "Merdith2021"     |
    |                       | Michael G. Tetley, Jacob A. Mulder, Morgan L. Blades,     |                   |
    |                       | Alexander Young, Sheree E. Armistead, John Cannon,        |                   |
    |                       | Sabin Zahirovic, R. Dietmar Müller, (2021).               |                   |
    |                       | Extending full-plate tectonic models into deep time:      |                   |
    |                       | Linking the Neoproterozoic and the Phanerozoic,           |                   |
    |                       | Earth-Science Reviews, Volume 214, 2021, 103477,          |                   |
    |                       | ISSN 0012-8252,                                           |                   |
    |                       | https://doi.org/10.1016/j.earscirev.2020.103477.          |                   |
    +-----------------------+-----------------------------------------------------------+-------------------+
    |                       |                                                           |                   |
    +-----------------------+-----------------------------------------------------------+-------------------+


    Methods
    -------
    get_plate_reconstruction_files
        Downloads the `rotation_model`, `topology_features`, and `static_polygons` needed to create an
        instance of the gplately.reconstruct.PlateReconstruction object.
    get_topology_geometries
        Downloads the `coastlines`, `continents` and `COBs` needed to create an instance of the
        gplately.plot.PlotTopologies object.
    get_netcdf_rasters
        Downloads netCDF (.nc) and .grd rasters, as well as .tif images

    Examples
    --------
    Calling the object: 
        # string identifier to access the Muller et al. 2019 model
        gDownload = gplately.download.DataServer("Muller2019")

    """
    def __init__(self, file_collection=None):
        if file_collection is None:
            raise ValueError(
                "Please supply a file collection to fetch."
            )
        self.file_collection = file_collection

    def get_plate_reconstruction_files(self):
        """Downloads and constructs a rotation model, a pygplates.FeatureCollection and a set of 
        static polygons needed to call the gplately.PlateReconstruction object.

        Notes
        -----
        This method accesses the plate reconstruction model requested in the gplately.DataServer
        object. For example, if the object was called as:

            gDownload = gplately.download.DataServer("Muller2019")

        the method will download:
            - a rotation file
            - GPML topology features
            - static polygons
        from the Muller et al. (2019) plate reconstruction model.
        """
        database = {

            "Muller2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics.zip"], 
            "Seton2012" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Rotations/Seton_etal_ESR2012_2012.1.rot",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Plate_polygons/Seton_etal_ESR2012_PP_2012.1.gpml",
                          None], 
            #"Merdith2021" : ["https://zenodo.org/record/4320873/files/SM2-Merdith_et_al_1_Ga_reconstruction.zip?download=1"],
            "Merdith2021" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/AltPlateReconstructions.zip"],
            "Matthews2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Matthews_etal_2016_Global_Plate_Model_GPC.zip"], 
            "Merdith2017" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2017_GR.zip"], 

        }

        # Set to true if we find the given collection in our database
        found_collection = False
        for collection, url in database.items():

            # Only continue if the user's chosen collection exists in our database
            if self.file_collection.lower() == collection.lower():
                found_collection = True
                rotation_filenames = []
                topology_features = pygplates.FeatureCollection()
                static_polygons = []

                if len(url) == 1:
                    fnames = collection_sorter(
                        fetch_from_web(url[0]), self.file_collection
                    )
                    rotation_filenames = collect_file_extension(
                        fnames, [".rot"]
                    )
                    rotation_model = pygplates.RotationModel(rotation_filenames)

                    topology_filenames = collect_file_extension(
                        str_in_filename(fnames, strings_to_ignore=["__MACOSX"]),
                        [".gpml", ".gpmlz"]
                    )
                    for file in topology_filenames:
                        topology_features.add(pygplates.FeatureCollection(file))

                    static_polygons = check_gpml_or_shp(
                        str_in_filename(
                            fnames, 
                            strings_to_include=["Static", "StaticPolygon", "Static_Polygon"],
                            strings_to_ignore=["__MACOSX"]
                        )
                    )
                else:
                    for file in url[0]:
                        rotation_filenames.append(
                            collect_file_extension(
                                fetch_from_web(file),
                                [".rot"])
                        )
                        rotation_model = pygplates.RotationModel(rotation_filenames)

                    for file in url[1]:
                        topology_filenames.append(
                            collect_file_extension(
                                fetch_from_web(file),
                                [".gpml"])
                        )
                        for file in topology_filenames:
                            topology_features.add(
                                pygplates.FeatureCollection(file)
                            )

                    for file in url[2]:
                        static_polygons.append(
                            check_gpml_or_shp(
                                str_in_filename(
                                    fetch_from_web(url[0]), 
                                    strings_to_include=["Static", "StaticPolygon", "Static_Polygon"],
                                    strings_to_ignore=["__MACOSX"]
                                )
                            )   
                        )
                break

        if not rotation_filenames:
            print("No .rot files in %s. No rotation model created." %self.file_collection)
        if not topology_filenames:
            print("No topology features in %s. No FeatureCollection created." %self.file_collection)
        if not static_polygons:
            print("No static polygons in %s." %self.file_collection)

        return rotation_model, topology_features, static_polygons


    def get_topology_geometries(self):
        """Downloads coastline, continent and continent-ocean boundary geometries needed to call
        the gplately.PlotTopologies object.

        Notes
        -----
        This method accesses the plate reconstruction model requested in the gplately.DataServer
        object. For example, if the object was called as:

            gDownload = gplately.download.DataServer("Muller2019")

        the method will download:
            - Coastlines: present-day coastlines cookie-cut using static polygons
            - Continents: cookie-cutting polygons for non-oceanic regions (continents, 
                          intra-oceanic arcs, etc.)
            - COBs: COB line segments
        from the Muller et al. (2019) plate reconstruction model.
        """
        database = {

            "Muller2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics.zip"], 
            "Seton2012" : [["https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Coastlines/Seton_etal_ESR2012_Coastline_2012.1.gpml",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Coastlines/Seton_etal_ESR2012_Coastline_2012.1_polyline.dbf",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Coastlines/Seton_etal_ESR2012_Coastline_2012.1_polyline.kml",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Coastlines/Seton_etal_ESR2012_Coastline_2012.1_polyline.prj",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Coastlines/Seton_etal_ESR2012_Coastline_2012.1_polyline.shp",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Coastlines/Seton_etal_ESR2012_Coastline_2012.1_polyline.shx"],
                           None,
                           ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Continent-ocean_boundaries/Seton_etal_ESR2012_ContinentOceanBoundaries_2012.1.dbf",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Continent-ocean_boundaries/Seton_etal_ESR2012_ContinentOceanBoundaries_2012.1.gpml",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Continent-ocean_boundaries/Seton_etal_ESR2012_ContinentOceanBoundaries_2012.1.kml",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Continent-ocean_boundaries/Seton_etal_ESR2012_ContinentOceanBoundaries_2012.1.prj",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Continent-ocean_boundaries/Seton_etal_ESR2012_ContinentOceanBoundaries_2012.1.shp",
                           "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Continent-ocean_boundaries/Seton_etal_ESR2012_ContinentOceanBoundaries_2012.1.shx"]], 
            "Merdith2021" : ["https://zenodo.org/record/4320873/files/SM2-Merdith_et_al_1_Ga_reconstruction.zip?download=1"],
            "Matthews2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Matthews_etal_2016_Global_Plate_Model_GPC.zip"], 
            "Merdith2017" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2017_GR.zip"],              
        }

        # Set to true if we find the given collection in our database
        found_collection = False
        for collection, url in database.items():

            # Only continue if the user's chosen collection exists in our database
            if self.file_collection.lower() == collection.lower():
                found_collection = True
                coastlines = []
                continents = []
                COBs = []

                if len(url) == 1:
                    fnames = fetch_from_web(url[0])
                    coastlines = check_gpml_or_shp(
                        str_in_filename(
                            fnames, 
                            strings_to_include=["coastline"], 
                            strings_to_ignore=["__MACOSX"])
                    )
                    continents = check_gpml_or_shp(
                        str_in_filename(
                            fnames, 
                            strings_to_include=["continent"], 
                            strings_to_ignore=["__MACOSX"])
                    )
                    COBs = check_gpml_or_shp(
                        str_in_filename(
                            fnames, 
                            strings_to_include=["cob"],
                            strings_to_ignore=["__MACOSX"])
                    )
                    files = coastlines, continents, COBs

                else:
                    for file in url[0]:
                        coastlines.append(str_in_filename(
                            fetch_from_web(file), 
                            strings_to_include=["coastline"])
                        )
                        coastlines = check_gpml_or_shp(coastlines)

                    for file in url[1]:
                        continents.append(str_in_filename(
                            fetch_from_web(file), 
                            strings_to_include=["continent"])
                        )
                        continents = check_gpml_or_shp(continents)

                    for file in url[2]:
                        COBs.append(str_in_filename(
                            fetch_from_web(file), 
                            strings_to_include=["cob"])
                        )
                        COBs = check_gpml_or_shp(COBs)

                    files = coastlines, continents, COBs
                break

        if not coastlines:
            print("No coastlines in %s." %self.file_collection)
        if not continents:
            print("No continents in %s." %self.file_collection)
        if not COBs:
            print("No continent-ocean boundaries in %s." %self.file_collection)
            
        return files


    def get_rasters(self, raster_id_string=None, time=None):
        """Downloads rasters and images from the web into the "gplately" cache.

        Currently supports the following rasters and images:
        +--------------+------------------------+---------------------------------------+-------------------+
        | Raster/image | Description            | Citation                              | String identifier |
        | name         |                        |                                       |                   |
        +--------------+------------------------+---------------------------------------+-------------------+
        | Muller et    | Seafloor age grid      | Müller, R. D., Zahirovic, S.,         | "Muller2019_nc"   |
        | al. 2019     | netCDF (.nc) rasters,  | Williams, S. E., Cannon, J.,          |                   |
        |              | as well as JPEG (.jpg) | Seton, M., Bower, D. J.,              | "Muller2019_jpg"  |
        |              | and PNG (.png) image   | Tetley, M. G., Heine, C.,             |                   |
        |              | equivalents.           | Le Breton, E., Liu, S.,               | "Muller2019_png"  |
        |              |                        | Russell, S. H. J., Yang, T.,          |                   |
        |              |                        | Leonard, J., and Gurnis, M. (2019),   |                   |
        |              |                        | A global plate model including        |                   |
        |              |                        | lithospheric deformation along major  |                   |
        |              |                        | rifts and orogens since the Triassic. |                   |
        |              |                        | Tectonics, vol. 38,                   |                   |
        |              |                        | https://doi.org/10.1029/2018TC005462. |                   |
        +--------------+------------------------+---------------------------------------+-------------------+
        | ETOPO1       | A 1-arc minute global  | doi:10.7289/V5C8276M                  | "ETOPO1_grd"      |
        |              | relief model combining |                                       |                   |
        |              | land topography and    |                                       | "ETOPO1_tif"      |
        |              | ocean bathymetry.      |                                       |                   |
        |              | Available in           |                                       |                   |
        |              | netCDF (.grd) and TIFF |                                       |                   |
        |              | (.tif)                 |                                       |                   |
        +--------------+------------------------+---------------------------------------+-------------------+
        |              |                        |                                       |                   |
        +--------------+------------------------+---------------------------------------+-------------------+

        Parameters
        ----------
        raster_id_string : str, default=None
            A string to identify which raster to download.
        time : int, default=None
            For time-reconstructed rasters, specify which time to extract.

        Returns
        -------
        raster_filenames : list of str
            A list containing the full path to the cached raster(s).

        Raises
        ------
        ValueError
            if a raster_id_string is not supplied.

        Notes
        -----
        By default, if rasters exist for each timestep (typically every 1 Ma), get_netcdf_rasters will return
        the full list of all rasters. For example, Muller et al. 2019 has 251 rasters for 0-250Ma. If the
        `time` parameter is passed, only the raster for that timestep is returned. This will be useful in time
        for loops.

        Examples
        --------
        1) Downloading the Muller et al. 2019 netCDF raster for 0-10 Ma:

                gdownload = gplately.download.DataServer("Muller2019")
                for time in arange(0,10):
                    time_raster = gdownload.get_raster("Muller2019_nc", time=time)

        2) Plotting a TIF global relief model on a GeoAxis instance

                etopo = gdownload.get_rasters("ETOPO1_tif")
                etopo_tif_img = image.imread(etopo)
                fig = plt.figure(figsize=(18,14), dpi=600)
                ax = fig.add_subplot(111, projection=ccrs.Mollweide(central_longitude = -150))
                ax.imshow(etopo_tif_img, extent=[-180,180,-90,90], transform=ccrs.PlateCarree())

        """
        if raster_id_string is None:
            raise ValueError(
                "Please specify which raster data to fetch."
            )

        database = {

            "Muller2019_nc" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_netCDF.zip"],
            "Muller2019_jpeg" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_jpgs.zip"],
            "Muller2019_png" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_pngs.zip"],
            "ETOPO1_grd" : ["https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"],
            "ETOPO1_tif" : ["https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/image/color_etopo1_ice_low.tif.gz"]
        }

        archive_formats = tuple([".gz", ".xz", ".bz2"])
        # Set to true if we find the given collection in database
        found_collection = False
        raster_filenames = []

        for collection, zip_url in database.items():
            if raster_id_string.lower() == collection.lower():
                found_collection = True
                raster_filenames = fetch_from_web(zip_url[0])
                if collection.endswith("nc"):
                    raster_filenames = order_filenames_by_time(
                        collect_file_extension(
                            raster_filenames, [".nc", ".grd"])
                    )
                elif collection.endswith("jpeg"):
                    raster_filenames = order_filenames_by_time(
                        collect_file_extension(
                            raster_filenames, [".jpg"])
                    )
                elif collection.endswith("png"):
                    raster_filenames = order_filenames_by_time(
                        collect_file_extension(
                            raster_filenames, [".png"])
                    )
                if time is not None:
                    raster_filenames = order_filenames_by_time(raster_filenames)[time]
                break

        if found_collection is False:
            raise ValueError("%s not in collection database." % (raster_id_string))
        return raster_filenames


    def get_feature_data(self, feature_data_id_string=None):
        """Downloads geological feature data from the web into the "gplately" cache.

        Currently supports the following feature data:
        +-------------------+---------+------------------------------------+-------------------------+
        | Feature data type | Formats | Paper citation                     | String identifier       |
        +-------------------+---------+------------------------------------+-------------------------+
        | Large Igneous     | .gpmlz  | Johansson, L., Zahirovic, S.,      | "LIP_VolcanicProvinces" |
        | Province products |         | and Müller, R. D., In Prep,        |                         |
        | from Johansson    |         | The interplay between the          |                         |
        | et al. (2018)     |         | eruption and weathering of         |                         |
        |                   |         | Large Igneous Provinces and        |                         |
        |                   |         | the deep-time carbon cycle:        |                         |
        |                   |         | Geophysical Research Letters.      |                         |
        +-------------------+---------+------------------------------------+-------------------------+
        | Large Igneous     | .gpmlz  | Whittaker, J. M., Afonso, J. C.,   | "LIP_VolcanicProvinces" |
        | Provinces         | .shp    | Masterton, S., Müller, R. D.,      |                         |
        | interpreted to    |         | Wessel, P., Williams, S. E.,       |                         |
        | be plume products |         | & Seton, M. (2015).                |                         |
        | from Whittaker    |         |  Long-term interaction between     |                         |
        | et al. (2015).    |         | mid-ocean ridges and mantle        |                         |
        |                   |         | plumes. Nature Geoscience, 8(6),   |                         |
        |                   |         | 479-483. doi:10.1038/ngeo2437.     |                         |
        +-------------------+---------+------------------------------------+-------------------------+
        | Seafloor tectonic | .gpml   | Matthews, K.J., M¸ller, R.D.,      | "SeafloorFabric"        |
        | fabric (fracture  |         | Wessel, P. and Whittaker, J.M.,    |                         |
        | zones, discordant |         | 2011. The tectonic fabric of the   |                         |
        | zones, V-shaped   |         | ocean basins. Journal of           |                         |
        | structures,       |         | Geophysical Research, 116(B12):    |                         |
        | unclassified      |         | B12109, DOI: 10.1029/2011JB008413. |                         |
        | V-anomalies,      |         |                                    |                         |
        | propagating ridge |         |                                    |                         |
        | lineations and    |         |                                    |                         |
        | extinct ridges)   |         |                                    |                         |
        | from Matthews     |         |                                    |                         |
        | et al. (2011)     |         |                                    |                         |
        +-------------------+---------+------------------------------------+-------------------------+
        | Present day       | .gpmlz  | Whittaker, J., Afonso, J.,         | "Hotspots"              |
        | surface hotspot/  |         | Masterton, S., Müller, R., Wessel, |                         |
        | plume locations   |         | P., Williams, S., and Seton, M.,   |                         | 
        | from Whittaker et |         | 2015, Long-term interaction between|                         |
        | al. (2013)        |         | mid-ocean ridges and mantle plumes:|                         |
        |                   |         | Nature Geoscience, v. 8, no. 6,    |                         |
        |                   |         | p. 479-483, doi:10.1038/ngeo2437.  |                         |
        +-------------------+---------+------------------------------------+-------------------------+
        
        Parameters
        ----------
        feature_data_id_string : str, default=None
            A string to identify which feature data to download to the cache. See table above.

        Returns
        -------
        feature_data_filenames : list of str
            A list containing the full path to the requested feature data. This is ready to be turned
            into a pygplates.FeatureCollection.

        Raises
        ------
        ValueError
            If a feature_data_id_string is not provided.
        """
        if feature_data_id_string is None:
            raise ValueError(
                "Please specify which feature data to fetch."
            )

        database = {

            "SeafloorFabric" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/SeafloorFabric.zip"],
            "LIP_VolcanicProvinces" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/IgneousProvinces.zip"],
            "Hotspots" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/Hotspots.zip"]
        }

        found_collection = False
        for collection, zip_url in database.items():
            if feature_data_id_string.lower() == collection.lower():
                found_collection = True
                feature_data_filenames = collect_file_extension(
                    fetch_from_web(zip_url[0]), [".gpml", ".gpmlz"]
                )
                break
        return feature_data_filenames
    