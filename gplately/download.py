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
    """Uses Pooch to download plate reconstruction feature data from plate models 
    and other studies that are stored on web servers (e.g. EarthByte's webDAV server). 
    
    All requested files are downloaded once into a 'gplately' cache folder. When 
    workflows requiring these files are rerun, they will be accessed from the cache. 

    Currently, DataServer supports the following plate reconstruction models:

        - Muller et al. 2019 : string identifier "Muller2019"
            Citation: Müller, R. D., Zahirovic, S., Williams, S. E., Cannon, J., Seton, M., 
            Bower, D. J., Tetley, M. G., Heine, C., Le Breton, E., Liu, S., Russell, S. H. J., 
            Yang, T., Leonard, J., and Gurnis, M., 2019, A global plate model including 
            lithospheric deformation along major rifts and orogens since the Triassic: 
            Tectonics, v. 38, no. Fifty Years of Plate Tectonics: Then, Now, and Beyond.
                
        - Muller et al. 2016 : string identifier "Muller2016"
            Citation: Müller R.D., Seton, M., Zahirovic, S., Williams, S.E., Matthews, K.J.,
             Wright, N.M., Shephard, G.E., Maloney, K.T., Barnett-Moore, N., Hosseinpour, M., 
             Bower, D.J., Cannon, J., InPress. Ocean basin evolution and global-scale plate 
             reorganization events since Pangea breakup, Annual Review of Earth and Planetary 
             Sciences, Vol 44, 107-138. DOI: 10.1146/annurev-earth-060115-012211.

        - Merdith et al. 2021 : string identifier "Merdith2021"
            Citation: Merdith et al. (in review), 'A continuous, kinematic full-plate motion model
             from 1 Ga to present'

        - Cao et al. 2020 : string identifier "Cao2020"
            Citation: Toy Billion-year reconstructions from Cao et al (2020). 
            Coupled Evolution of Plate Tectonics and Basal Mantle Structure Tectonics, 
            doi: 10.1029/2020GC009244

        - Mather et al. 2021 : string identifier "Mather2021"
            Citation: Mather, B., Müller, R.D.,; Alfonso, C.P., Seton, M., 2021, Kimberlite eruption 
            driven by slab flux and subduction angle. DOI: 10.5281/zenodo.5769002

        - Seton et al. 2012 : string identifier "Seton2012"
            Citation: M. Seton, R.D. Müller, S. Zahirovic, C. Gaina, T.H. Torsvik, G. Shephard, A. Talsma, 
            M. Gurnis, M. Turner, S. Maus, M. Chandler, Global continental and ocean basin reconstructions 
            since 200 Ma, Earth-Science Reviews, Volume 113, Issues 3-4, July 2012, Pages 212-270, 
            ISSN 0012-8252, 10.1016/j.earscirev.2012.03.002.

        - Matthews et al. 2016 : string identifier "Matthews2016"
            Citation: Matthews, K.J., Maloney, K.T., Zahirovic, S., Williams, S.E., Seton, M.,
            and Müller, R.D. (2016). Global plate boundary evolution and kinematics since the 
            late Paleozoic, Global and Planetary Change, 146, 226-250. 
            DOI: 10.1016/j.gloplacha.2016.10.002

        - Merdith et al. 2017 : string identifier "Merdith2017"
            Citation: Merdith, A., Collins, A., Williams, S., Pisarevskiy, S., Foden, J., Archibald, D. 
            and Blades, M. et al. 2016. A full-plate global reconstruction of the Neoproterozoic. 
            Gondwana Research. 50: pp. 84-134. DOI: 10.1016/j.gr.2017.04.001

        - Li et al. 2008 : string identifier "Li2008"
            Citation: Rodinia reconstruction from Li et al (2008), doi: 10.1016/j.precamres.2007.04.021

        - Pehrsson et al. 2015 : string identifier "Pehrsson2015"
            Citation:

        - Torsvik and Cocks et al. 2017 : string identifier "TorsvikCocks2017"
            Citation:

        - Young et al. 2019 : string identifier "Young2019"
            Citation:

        - Scotese et al. 2008 : string identifier "Scotese2008"
            Citation: 

        - Golonka et al. 2007 : string identifier "Golonka2007"
            Citation:

        - Clennett et al. 2020 (based on Muller et al. 2019) : string identifier "Clennett2020_M2019"
            Citation:

        - Clennett et al. 2020 (rigid topological model based on Shephard et al, 2013) : string identifier "Clennett2020_S2013"
            Citation:


    When calling the object, just supply one of these string identifiers. The object
    will look through the requested plate model when calling the following methods:

    Methods
    -------
    get_plate_reconstruction_files
        Downloads the `rotation_model`, `topology_features`, and `static_polygons` 
        needed to create an instance of the <gplately.reconstruction.PlateReconstruction> 
        object. Returns them as <pygplates.RotationModel> and <pygplates.FeatureCollection> 
        instances respectively.
    get_topology_geometries
        Downloads the `coastlines`, `continents` and `COBs` needed to create an instance of the
        <gplately.plot.PlotTopologies> object. These geometries are returned as 
        <pygplates.FeatureCollection> objects.
    get_age_grids
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
        """Downloads and constructs a rotation model, a set of dynamic polygons and
        and a set of static polygons needed to call the gplately.PlateReconstruction 
        object.

        Returns
        -------
        rotation_model : instance of <pygplates.RotationModel>
            A rotation model to query equivalent and/or relative topological plate rotations
            from a time in the past relative to another time in the past or to present day.
        topology_features : instance of <pygplates.FeatureCollection>
            Point, polyline and/or polygon feature data in motion through geological time.
        static_polygons : instance of <pygplates.FeatureCollection>
            Present-day polygons whose shapes do not change through geological time. They are
            used to cookie-cut dynamic polygons into identifiable topological plates (assigned 
            an ID) according to their present-day locations.

        Notes
        -----
        This method accesses the plate reconstruction model requested in the 
        gplately.DataServer object. For example, if the object was called as:

            gDownload = gplately.download.DataServer("Muller2019")

        the method will download:
            - a rotation file
                returned as a <pygplates.RotationModel>
            - GPML topology features
                returned as a <pygplates.FeatureCollection>
            - static polygons
                returned as a <pygplates.FeatureCollection>

        from the Muller et al. (2019) plate reconstruction model. If the requested plate 
        model does not have a certain files, say `static_polygons`, a message will be printed 
        to alert the user. 

        Once the reconstruction objects are returned, they can be passed into:

            model = gplately.reconstruction.PlateReconstruction(
                rotation_model
                topology_features,
                coastlines,
                COBs)

        the <gplately.reconstruction.PlateReconstruction> object reconstructs features
        to certain geological times (see documentation for more info).
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
            Shapely polygons and polylines resolved from .shp and/or .gpml topology files,
            ready for reconstruction to a particular geological time and plotting.

        Notes
        -----
        This function searches for the plate model requested when calling the 
        <gplately.data.DataServer> object. For example, if the object was called as:

            gDownload = gplately.download.DataServer("Muller2019")

        the method will attempt to download:
            - Coastlines: present-day coastlines cookie-cut using static polygons
            - Continents: cookie-cutting polygons for non-oceanic regions (continents, 
                          intra-oceanic arcs, etc.)
            - COBs: COB line segments
        from the Muller et al. (2019) plate reconstruction model. If found, these files are 
        returned as individual pyGPlates Feature Collections. 

        If the requested plate model does not have a certain geometry, say `continents`, a
        message will be printed to alert the user. 

        Once the continents, coastlines and COBs Feature Collections are returned, they can 
        be passed into:

            gPlot = gplately.plot.PlotTopologies(
                <gplately.reconstruction.PlateReconstruction>,
                time,
                continents,
                coastlines,
                COBs)

        to reconstruct features to a certain geological time. The <gplately.plot.PlotTopologies> 
        object provides simple methods to plot these geometries along with trenches, ridges and 
        transforms (see documentation for more info).
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

        Currently, DataServer supports the following rasters and images:

        - Muller et al. 2019
            Time range: 0-250 Ma
            Seafloor age grid rasters in netCDF format.

        - Muller et al. 2016
            Time range: 0-240 Ma
            Seafloor age grid rasters in netCDF format. 

        Note that this function will download the age grid(s) from the plate model 
        requested via the string identifer passed into the DataServer object. For example, 
        if the object was called as:

            gDownload = gplately.download.DataServer("Muller2019")

        the method will attempt to download age grid(s) from the Muller et al. (2019) plate
        reconstruction model from the geological time(s) requested in the `times` parameter. 
        If found, these age grids are returned as masked arrays. 
        
        Parameters
        ----------
        times : int, or list of int, default=None
            Request an age grid from one or more reconstruction times, e.g. from 0-5 Ma
            requires times=np.arange(0,5). If a single integer is passed, a single raster
            masked array is returned. If a list of integers is passed, a list of raster
            masked arrays is returned.

        Returns
        -------
        raster_array : ndarray
            A masked array containing the read netCDF4 grid, ready for plotting or for
            passing into the <gplately.grid.Raster> object for raster manipulation.

        Raises
        -----
        ValueError
            If `times` (a list of reconstruction times to extract the age grids from) is 
            not passed.

        Notes
        -----
        Once requested age grid(s) are downloaded to the gplately cache once, they are not 
        re-downloaded if the same workflow (or even a different one!) requires them. Rather,
        DataServer fetches them from the cache.
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
        """Downloads assorted rasters and images from the web that are not associated with
        a plate reconstruction model in DataServer into the "gplately" cache.

        Currently supports the following rasters and images:

        - ETOPO1 
            Filetypes available : TIF, netCDF (GRD)
            string identifiers : "ETOPO1_grd", "ETOPO1_tif"
            A 1-arc minute global relief model combining lang topography and ocean bathymetry.
            Available in netCDF (in .grd) and TIFF (.tif) format. 


        Parameters
        ----------
        raster_id_string : str, default=None
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

        - Large igneous provinces from Johansson et al. (2018)
            Formats: .gpmlz
            String identifier: "LIP_VolcanicProvinces"
            Citation: Johansson, L., Zahirovic, S., and Müller, R. D., In Prep, The 
            interplay between the eruption and weathering of Large Igneous Provinces and 
            the deep-time carbon cycle: Geophysical Research Letters.

        - Large igneous province products interpreted as plume products from Whittaker 
        et al. (2015).
            Formats: .gpmla, .shp
            String identifier: "LIP_VolcanicProvinces"
            Citation: Whittaker, J. M., Afonso, J. C., Masterton, S., Müller, R. D., 
            Wessel, P., Williams, S. E., & Seton, M. (2015). Long-term interaction between 
            mid-ocean ridges and mantle plumes. Nature Geoscience, 8(6), 479-483. 
            doi:10.1038/ngeo2437.

        - Seafloor tectonic fabric (fracture zones, discordant zones, V-shaped structures, 
        unclassified V-anomalies, propagating ridge lineations and extinct ridges) from 
        Matthews et al. (2011)
            Formats: .gpml
            String identifier: "SeafloorFabric"
            Citation: Matthews, K.J., Müller, R.D., Wessel, P. and Whittaker, J.M., 2011. The 
            tectonic fabric of the ocean basins. Journal of Geophysical Research, 116(B12): 
            B12109, DOI: 10.1029/2011JB008413. 

        - Present day surface hotspot/plume locations from Whittaker et al, (2013)
            Formats: .gpmlz
            String identifier: "Hotspots"
            Citation: Whittaker, J., Afonso, J., Masterton, S., Müller, R., Wessel, P., 
            Williams, S., and Seton, M., 2015, Long-term interaction between mid-ocean ridges and 
            mantle plumes: Nature Geoscience, v. 8, no. 6, p. 479-483, doi:10.1038/ngeo2437.

        
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
    