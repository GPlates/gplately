import re as _re

# ============ LINKS FOR DOWNLOADING PLATE RECONSTRUCTION FILES ON-THE-FLY ==============
""" This auxiliary script contains links to download plate model data from web servers
such as EarthByte's webDAV server. These links are accessed by GPlately's DataServer object
and downloaded into the "gplately" folder in your machine's cache via Pooch.
"""

def _find_needed_collection(
    collection_identifier,
    data_dictionary,
    times=None
    ):
    """Search through link databases and time arrays (or single integer
    time) for requested download links."""

    all_urls = []
    for collection, url in data_dictionary.items():
        if collection_identifier.lower() == collection.lower():
            if times is not None:
                if isinstance(times, int):
                    all_urls.append(url[0].format(str(times)))
                    return all_urls
                elif isinstance(times, list):
                    for t in times:
                        url_current_time = url[0].format(str(t))
                        all_urls.append(url_current_time)
                    return all_urls
            else:
                return url


def _studyname(study_name):
    name = _re.findall(r'[A-Za-z]+|\d+', study_name)[0]
    return name

class DataCollection(object):
    """GPlately's collection of plate model data is a dictionary where
    the plate model's identifier string is the key, while values are 
    lists containing any relevant file download links."""

    def __init__(self, file_collection):
        """Uses a string to identify the needed plate model, taken from
        <gplately.data.DataServer>."""
        self.file_collection = file_collection


    def netcdf4_age_grids(self, times):

        age_grid_links = {

            "Muller2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-{}.nc"],
            "Muller2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Agegrids/Muller_etal_2016_AREPS_Agegrids_v1.17/Muller_etal_2016_AREPS_v1.17_netCDF/Muller_etal_2016_AREPS_v1.17_AgeGrid-{}.nc"],
        }

        links_to_download = _find_needed_collection(
            self.file_collection, 
            age_grid_links,
            times)

        return links_to_download


    def plate_reconstruction_files(self):

        database = {

            "Cao2020" : ["https://zenodo.org/record/3854549/files/1000Myr_synthetic_tectonic_reconstructions.zip"],
            "Muller2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics.zip"], 
            "Muller2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17.zip"],
            "Mather2021" : ["https://zenodo.org/record/5769002/files/plate_model.zip"],
            "Seton2012" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR.zip"],
            "Merdith2021" : ["https://zenodo.org/record/4485738/files/SM2_4485738_V2.zip"],
            "Matthews2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Matthews_etal_2016_Global_Plate_Model_GPC.zip"], 
            "Merdith2017" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2017_GR.zip"], 
            "Li2008" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Li_etal_2008_RodiniaModel.zip"],
            "Pehrsson2015" : ["https://www.geolsoc.org.uk/~/media/Files/GSL/shared/Sup_pubs/2015/18822_7.zip"],
            "TorsvikCocks2017" : ["http://www.earthdynamics.org/earthhistory/bookdata/CEED6.zip"],
            "Young2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Young_etal_2018_GeoscienceFrontiers/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel.zip"], 
            "Scotese2008" : ["https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip"],      
            "Golonka2007" : ["https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip"],
            "Clennett2020_M2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_M2019.zip"],
            "Clennett2020_S2013" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_S2013.zip"],

        }

        return database


    def rotation_strings_to_ignore(self):

        strings = [
            "OLD",
            "__MACOSX",
            "DO_NOT",
        ]    
        return strings


    def dynamic_polygon_strings_to_include(self):

        strings = [
            "plate_boundaries",
            "PlateBoundaries",
            "Transform",
            "Divergence",
            "Convergence",
            "Topologies",
            "Topology",
            #"ContinentOceanBoundaries",
            "Seton_etal_ESR2012_Coastline_2012",
            "Deforming_Mesh",
            "boundaries",
            "Clennett_etal_2020_Plates", # For Clennett 2020 (M2019)

        ]
        return strings 


    def dynamic_polygon_strings_to_ignore(self):

        strings = [
            "OLD",
            "__MACOSX",
            "DO_NOT",
        ]
        return strings


    def static_polygon_strings_to_include(self):

        strings = [
            "Static",
            "StaticPolygon",
            "Static_Polygon",
            "RodiniaBlocks_WithPlateIDColumnAndIDs",
            "PlatePolygons.shp",
            "CEED6_TERRANES.shp",
            "CEED6_MICROCONTINENTS.shp",
            "CEED6_LAND.gpml",
            "Scotese_2008_PresentDay_ContinentalPolygons", # Scotese 2008
            "Golonka_2007_PresentDay_ContinentalPolygons.shp", # Golonka 2007
        ]
        return strings


    def static_polygon_strings_to_ignore(self):

        strings = [

            "DO_NOT",
            "OLD",
            "__MACOSX"
        ]
        return strings


    def topology_geometries(self):

        database = {

            "Cao2020" : ["https://zenodo.org/record/3854549/files/1000Myr_synthetic_tectonic_reconstructions.zip"],
            "Muller2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics.zip"], 
            "Muller2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17.zip"],
            "Mather2021" : ["https://zenodo.org/record/5769002/files/plate_model.zip"],
            "Seton2012" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR.zip"],
            "Merdith2021" : ["https://zenodo.org/record/4485738/files/SM2_4485738_V2.zip"],
            "Matthews2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Matthews_etal_2016_Global_Plate_Model_GPC.zip"], 
            "Merdith2017" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2017_GR.zip"],  
            "Pehrsson2015" : ["https://www.geolsoc.org.uk/~/media/Files/GSL/shared/Sup_pubs/2015/18822_7.zip"],
            "TorsvikCocks2017" : ["http://www.earthdynamics.org/earthhistory/bookdata/CEED6.zip"],
            "Young2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Young_etal_2018_GeoscienceFrontiers/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel.zip"],
            "Scotese2008" : ["https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip"],
            "Golonka2007" : ["https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip"],
            "Clennett2020_M2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_M2019.zip"],
            "Clennett2020_S2013" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_S2013.zip"],

        }
        return database


    def coastline_strings_to_include(self):

        strings = [

            "coastline",
            "CEED6_LAND.gpml" # for TorsvikCocks2017
        ]
        return strings


    def coastline_strings_to_ignore(self):

        strings = [

            "DO_NOT",
            "OLD",
            "__MACOSX"
        ]
        return strings


    def continent_strings_to_include(self):

        strings = [

            "continent",
            "COBfile_1000_0_Toy_introversion",
            "continental",
            "Scotese_2008_PresentDay_ContinentalPolygons.shp", # Scotese 2008
            "Terrane"
        ]
        return strings


    def continent_strings_to_ignore(self):

        strings = [

            "DO_NOT",
            "OLD",
            "__MACOSX",
            "Continent-ocean_boundaries"
        ]
        return strings


    def COB_strings_to_include(self):

        strings = [

            "cob",
            "ContinentOceanBoundaries",
        ]
        return strings


    def COB_strings_to_ignore(self):

        strings = [

            "DO_NOT",
            "OLD",
            "__MACOSX"
        ]
        return strings


    def rasters(self):

        database = {

            "ETOPO1_grd" : ["https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"],
            "ETOPO1_tif" : ["https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/image/color_etopo1_ice_low.tif.gz"],
        }
        return database


    def feature_data(self):
        """Assorted feature data from EarthByte's GPlates 2.3 sample data repository."""

        database = {

            "SeafloorFabric" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/SeafloorFabric.zip"],
            "LIP_VolcanicProvinces" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/IgneousProvinces.zip"],
            "Hotspots" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/Hotspots.zip"]
        }
        return database





















import pooch as _pooch
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
from pooch import Decompress as _Decompress
from matplotlib import image as _image
from .data import DataCollection
import gplately as _gplately
import pygplates as _pygplates
import re as _re
import numpy as _np

def _fetch_from_web(url):
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


def _collect_file_extension(fnames, file_extension):
    """Searches cached directory for filenames with a specified extension(s)."""
    sorted_fnames = []
    file_extension=tuple(file_extension)
    for file in fnames:
        if file.endswith(file_extension):
            sorted_fnames.append(file)
    return sorted_fnames


def _str_in_folders(fnames, strings_to_include=None, strings_to_ignore=None):
    """Collect and ignore file paths with strings to include and/or ignore
    from their parent directories."""
    sorted_fnames = []
    for i, fname in enumerate(fnames):
        parent_directory = '/'.join(fname.split("/")[:-1])
        if strings_to_ignore is not None:
            check = [s for s in strings_to_ignore if s in parent_directory]
            if check:
                continue
        if strings_to_include is not None:
            if any(x.lower() in fname.lower() for x in strings_to_include):
                sorted_fnames.append(fname)
        else:
            sorted_fnames.append(fname)
    #if sorted_fnames:
        #return sorted_fnames
    #else:
    return sorted_fnames


def _str_in_folder(fnames, strings_to_include=None, strings_to_ignore=None):
    fnames_to_ignore = []
    fnames_to_include = []
    sorted_fnames = []
    for i, fname in enumerate(fnames):
        parent_directory = '/'.join(fname.split("/")[:-1])
        if strings_to_ignore is not None:
            for s in strings_to_ignore:
                if s in parent_directory:
                    fnames_to_ignore.append(fname)
            sorted_fnames = list(set(fnames) - set(fnames_to_ignore))

    if strings_to_include is not None:
        for fname in sorted_fnames:
            parent_directory = '/'.join(fname.split("/")[:-1])
            for s in strings_to_include:
                if s in parent_directory:
                    fnames_to_include.append(fname)
        sorted_fnames = list(set(sorted_fnames).intersection(set(fnames_to_include)))
    return sorted_fnames


def _str_in_filename(fnames, strings_to_include=None, strings_to_ignore=None):
    sorted_fnames = []
    if strings_to_ignore is not None:
        for f in fnames:
            f = f.split("/")[-1]
            check = [s for s in strings_to_ignore if s.lower() in f.lower()]
    if strings_to_include is not None:
        for s in strings_to_include:
            for f in fnames:
                fname = f.split("/")[-1]
                if s.lower() in fname.lower():
                    sorted_fnames.append(f)
    return sorted_fnames


def _check_gpml_or_shps(fnames):
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


def _check_gpml_or_shp(fnames):
    """For topology features, returns GPML by default. Searches for ESRI Shapefiles 
    instead if GPML files not found."""
    sorted_fnames = []
    for file in fnames:
        if file.endswith(".gpml"):
            sorted_fnames.append(file)
        elif file.endswith(".shp"):
            sorted_fnames.append(file)
    return sorted_fnames


def _remove_hash(fname):
    """Removes hashes (32 character file IDs) from cached filenames."""
    split_paths = fname.split("-")
    cache_path = split_paths[0][:-32]
    new_path = cache_path + "-".join(split_paths[1:])
    return new_path


def _order_filenames_by_time(fnames):
    """Orders filenames in a list from present day to deeper geological time if they
    are labelled by time."""
    # Collect all digits in each filename.
    filepath_digits=[]
    for i, file in enumerate(fnames):
        digits = []
        for element in _re.split('([0-9]+)', _remove_hash(file)):
            if element.isdigit():
                digits.append(int(str(element)))
        filepath_digits.append(digits)

    # Ignore digits common to all full file paths. This leaves behind the files' 
    # geological time label.
    geological_times = []
    filepath_digits = _np.array(filepath_digits).T
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


def _collection_sorter(fnames, string_identifier):
    """If multiple file collections or plate reconstruction models are downloaded from
    a single zip folder, only return the needed model. 

    The plate models that need separating are listed."""

    needs_sorting = [
        "merdith2021",
        "scotese2008",
        "golonka2007",
        "clennett2020"
    ]
    if string_identifier.lower() in needs_sorting:
        studyname = _re.findall(r'[A-Za-z]+|\d+', string_identifier)[0]
        newfnames = []
        for files in fnames:
            if studyname not in files:
                continue
            newfnames.append(files)
        return newfnames
    else:
        return fnames


def _match_filetype_to_extension(filetype):
    extensions = []
    if filetype == "netCDF":
        extensions.append(".nc")
    elif filetype == "jpeg":
        extensions.append(".jpg")
    elif filetype == "png":
        extensions.append(".png")
    elif filetype == "TIFF":
        extensions.append(".tif")
    return extensions


class DataServer(object):
    """Uses Pooch to download geological feature data from plate reconstruction models 
    and other studies that are stored on web servers (e.g. EarthByte's webDAV server). 
    Downloaded files are kept in a 'gplately' cache folder. 

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
        Downloads the `rotation_model`, `topology_features`, and `static_polygons` 
        needed to create an instance of the gplately.reconstruct.PlateReconstruction object.
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
        self.data_collection = DataCollection(self.file_collection)


    def get_plate_reconstruction_files(self):
        """Downloads and constructs a rotation model, a pygplates.FeatureCollection 
        and a set of static polygons needed to call the gplately.PlateReconstruction 
        object.

        Notes
        -----
        This method accesses the plate reconstruction model requested in the 
        gplately.DataServer object. For example, if the object was called as:

            gDownload = gplately.download.DataServer("Muller2019")

        the method will download:
            - a rotation file
            - GPML topology features
            - static polygons
        from the Muller et al. (2019) plate reconstruction model.
        """

        rotation_filenames = []
        rotation_model = []
        topology_filenames = []
        topology_features = _pygplates.FeatureCollection()
        static_polygons = []

        # Locate all plate reconstruction files from GPlately's DataCollection
        database = DataCollection.plate_reconstruction_files(self)

        # Set to true if we find the given collection in our database
        found_collection = False
        for collection, url in database.items():

            # Only continue if the user's chosen collection exists in our database
            if self.file_collection.lower() == collection.lower():
                found_collection = True
                if len(url) == 1:
                    fnames = _collection_sorter(
                        _fetch_from_web(url[0]), self.file_collection
                    )
                    rotation_filenames = _str_in_folder(
                        _collect_file_extension(fnames, [".rot"]),
                        strings_to_ignore=DataCollection.rotation_strings_to_ignore(self)
                    )

                    #print(rotation_filenames)
                    rotation_model = _pygplates.RotationModel(rotation_filenames)

                    topology_filenames = _collect_file_extension(
                        _str_in_folder(
                            _str_in_filename(fnames, 
                                strings_to_include=DataCollection.dynamic_polygon_strings_to_include(self)
                            ), 
                            strings_to_ignore=DataCollection.dynamic_polygon_strings_to_ignore(self)
                        ),
                        [".gpml", ".gpmlz"]
                    )
                    #print(topology_filenames)
                    for file in topology_filenames:
                        topology_features.add(_pygplates.FeatureCollection(file))

                    static_polygons = _check_gpml_or_shp(
                        _str_in_folder(
                            _str_in_filename(fnames, 
                                strings_to_include=DataCollection.static_polygon_strings_to_include(self)
                            ),
                            strings_to_ignore=DataCollection.static_polygon_strings_to_ignore(self)
                        )
                    )
                    #print(static_polygons)
                else:
                    for file in url[0]:
                        rotation_filenames.append(
                            _collect_file_extension(
                                _fetch_from_web(file), [".rot"])
                        )
                        rotation_model = _pygplates.RotationModel(rotation_filenames)

                    for file in url[1]:
                        topology_filenames.append(
                            _collect_file_extension(
                                _fetch_from_web(file), [".gpml"])
                        )
                        for file in topology_filenames:
                            topology_features.add(
                                _pygplates.FeatureCollection(file)
                            )

                    for file in url[2]:
                        static_polygons.append(
                            _check_gpml_or_shp(
                                _str_in_folder(
                                    _str_in_filename(_fetch_from_web(url[0]), 
                                        strings_to_include=DataCollection.static_polygon_strings_to_include(self)
                                    ),    
                                        strings_to_ignore=DataCollection.static_polygon_strings_to_ignore(self)
                                )
                            )   
                        )
                break

        if found_collection is False:
            raise ValueError("{} is not in GPlately's DataServer.".format(self.file_collection))

        if not rotation_filenames:
            print("No .rot files in {}. No rotation model created.".format(self.file_collection))
        if not topology_filenames:
            print("No topology features in {}. No FeatureCollection created - unable to plot trenches, ridges and transforms.".format(self.file_collection))
        if not static_polygons:
            print("No static polygons in {}.".format(self.file_collection))

        return rotation_model, topology_features, static_polygons


    def get_topology_geometries(self):
        """Downloads coastline, continent and continent-ocean boundary geometries from the 
        requested plate model. These are needed to call the <gplately.plot.PlotTopologies> 
        object.

        Returns
        -------
        continents, coastlines, COBs : instance of <pygplates.FeatureCollection>

        Notes
        -----
        This method accesses the plate reconstruction model requested when calling the 
        <gplately.data.DataServer> object. For example, if the object was called as:

            gDownload = gplately.download.DataServer("Muller2019")

        the method will download:
            - Coastlines: present-day coastlines cookie-cut using static polygons
            - Continents: cookie-cutting polygons for non-oceanic regions (continents, 
                          intra-oceanic arcs, etc.)
            - COBs: COB line segments
        from the Muller et al. (2019) plate reconstruction model. They are returned as 
        individual pyGPlates Feature Collections. 
        """

        # Locate all topology geometries from GPlately's DataCollection
        database = DataCollection.topology_geometries(self)

        coastlines = []
        continents = []
        COBs = []
        
        # Find the requested plate model data collection
        found_collection = False
        for collection, url in database.items():

            if self.file_collection.lower() == collection.lower():
                found_collection = True

                if len(url) == 1:
                    fnames = _collection_sorter(
                        _fetch_from_web(url[0]), self.file_collection
                    )
                    coastlines = _check_gpml_or_shp(
                        _str_in_folder(
                            _str_in_filename(
                                fnames,
                                strings_to_include=DataCollection.coastline_strings_to_include(self)
                            ), 
                            strings_to_ignore=DataCollection.coastline_strings_to_ignore(self)
                        )
                    )
                    continents = _check_gpml_or_shp(
                        _str_in_folder(
                            _str_in_filename(
                                fnames, 
                                strings_to_include=DataCollection.continent_strings_to_include(self)
                            ), 
                            strings_to_ignore=DataCollection.continent_strings_to_ignore(self)
                        )
                    )
                    COBs = _check_gpml_or_shp(
                        _str_in_folder(
                            _str_in_filename(
                                fnames,
                                strings_to_include=DataCollection.COB_strings_to_include(self)
                            ), 
                            strings_to_ignore=DataCollection.COB_strings_to_ignore(self)
                        )
                    )
                else:
                    for file in url[0]:
                        coastlines.append(_str_in_filename(
                            _fetch_from_web(file), 
                            strings_to_include=["coastline"])
                        )
                        coastlines = _check_gpml_or_shp(coastlines)

                    for file in url[1]:
                        continents.append(_str_in_filename(
                            _fetch_from_web(file), 
                            strings_to_include=["continent"])
                        )
                        continents = _check_gpml_or_shp(continents)

                    for file in url[2]:
                        COBs.append(_str_in_filename(
                            _fetch_from_web(file), 
                            strings_to_include=["cob"])
                        )
                        COBs = _check_gpml_or_shp(COBs)
                break

        if found_collection is False:
            raise ValueError("{} is not in GPlately's DataServer.".format(self.file_collection))

        if not coastlines:
            print("No coastlines in {}.".format(self.file_collection))
        else:
            #print(coastlines)
            coastlines_featurecollection = _pygplates.FeatureCollection()
            for coastline in coastlines:
                coastlines_featurecollection.add(_pygplates.FeatureCollection(coastline))
        
        if not continents:
            print("No continents in {}.".format(self.file_collection))
        else:
            #print(continents)
            continents_featurecollection = _pygplates.FeatureCollection()
            for continent in continents:
                continents_featurecollection.add(_pygplates.FeatureCollection(continent))
        
        if not COBs:
            print("No continent-ocean boundaries in {}.".format(self.file_collection))
        else:
            #print(COBs)
            COBs_featurecollection = _pygplates.FeatureCollection()
            for COB in COBs:
                COBs_featurecollection.add(_pygplates.FeatureCollection(COB))
        
        geometries = coastlines, continents, COBs
        return geometries


    def get_age_grids(self, times=None):
        """Downloads age grids from plate reconstruction files on GPlately's DataServer 
        into the "gplately" cache.

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
        |              | 0-250 Ma               | Leonard, J., and Gurnis, M. (2019),   |                   |
        |              |                        | A global plate model including        |                   |
        |              |                        | lithospheric deformation along major  |                   |
        |              |                        | rifts and orogens since the Triassic. |                   |
        |              |                        | Tectonics, vol. 38,                   |                   |
        |              |                        | https://doi.org/10.1029/2018TC005462. |                   |
        +--------------+------------------------+---------------------------------------+-------------------+
        | Muller et    | Seafloor age grid      |
        | al. 2016     | netCDF (.nc) rasters,  |
        |              | as well as JPEG (.jpg) |
        |              | ang PNG (.png) image   |
        |              | equivalents.   
        |              | 
        |              | 0-240 Ma       
        +--------------+------------------------+---------------------------------------+-------------------+
        
        Parameters
        ----------
        times : int, or list of int, default=None
            Request an age grid from one or more reconstruction times, e.g. from 0-5 Ma
            requires times=np.arange(0,5).

        Returns
        -------
        raster_array : ndarray
            A masked array containing the read netCDF4 grid, ready for plotting.

        Raises
        -----
        ValueError
            If `times` (a list of reconstruction times to extract the age grids from) is 
            not passed.

        Notes
        -----
        Once the requested age grid(s) are downloaded to the cache, the grid(s) are read 
        by GPlately's Raster object and returned as a masked array
        objects. 

        """
        if times is None:
            raise ValueError("Please supply a list of times.")

        age_grids = []
        age_grid_links = DataCollection.netcdf4_age_grids(self, times)
        for link in age_grid_links:
            age_grid_file = _fetch_from_web(link)
            age_grid = _gplately.grids.read_netcdf_grid(age_grid_file)
            age_grids.append(age_grid)

        if not age_grids:
            raise ValueError("{} netCDF4 age grids not found.".format())

        if len(age_grids) == 1:
            return age_grids[0]
        else: 
            return age_grids


    def get_raster(self, raster_id_string=None, filetype=None):
        """Downloads assorted rasters and images from the web into the "gplately" cache.

        Currently supports the following rasters and images:
        +--------------+------------------------+---------------------------------------+-------------------+
        | Raster/image | Description            | Citation                              | String identifier |
        | name         |                        |                                       |                   |
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
        raster_id_string : str, defaultNone
            A string to identify which raster to download.
        filetype : str, default None
            A string to request an age grid of a particular filetype. Currently supports
                - netCDF
                - JPEG
                - PNG

        Returns
        -------
        raster_filenames : list of str
            A list containing the full path to the cached raster(s).

        Raises
        ------
        ValueError
            if a raster_id_string is not supplied.
            if a filetype is not supplied.

        Notes
        -----
        Rasters obtained by this method are (so far) only reconstructed to present-day. 
        """
        if raster_id_string is None:
            raise ValueError(
                "Please specify which raster to download."
            )
        if filetype is None:
            raise ValueError(
                "Please specify which raster filetype to download (i.e. 'tif')."
            )
        filetype = "."+filetype

        archive_formats = tuple([".gz", ".xz", ".bz2"])
        # Set to true if we find the given collection in database
        found_collection = False
        raster_filenames = []
        database = DataCollection.rasters(self)

        for collection, zip_url in database.items():
            # Isolate the raster name and the file type
            raster_name = collection.split("_")[0]
            raster_type = "."+collection.split("_")[-1]
            if (raster_id_string.lower() == raster_name.lower()
                and filetype.lower() == raster_type.lower()
                ):
                raster_filenames = _fetch_from_web(zip_url[0])
                found_collection = True
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

        database = DataCollection.feature_data(self)

        found_collection = False
        for collection, zip_url in database.items():
            if feature_data_id_string.lower() == collection.lower():
                found_collection = True
                feature_data_filenames = _collect_file_extension(
                    _fetch_from_web(zip_url[0]), [".gpml", ".gpmlz"]
                )
                break
        return feature_data_filenames
    