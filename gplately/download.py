""" Functions for downloading assorted plate reconstruction data to use with GPlately's
main objects. Files are stored in the user's cache and can be reused after being
downloaded once. 

These data have been created and used in plate reconstruction models and studies, and
are available from public web servers (like EarthByte's webDAV server, or the GPlates 
2.3 sample dataset library).

"""
import pooch as _pooch
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import utils as _utils
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
from pooch import Decompress as _Decompress
from matplotlib import image as _image
from gplately.data import DataCollection
import gplately as _gplately
import pygplates as _pygplates
import re as _re
import os as _os
import numpy as _np
import urllib.request as _request
import hashlib as _hashlib
import pathlib as _pathlib
import shutil as _shutil
import requests as _requests

def _test_internet_connection(url):
    """Test whether a connection to the required web server
    can be made given a `url`.
    
    Returns `False` the `url` is incorrect, and/or if there
    is no internet connection."""
    try:
        _request.urlopen(url) 
        return True
    except:
        return False
    

def _determine_processor(url):
    """Set instructions for how to process/unpack a file depending on
    its filetype. The unpacked file paths will have an .unzip, or 
    .decomp, or no file extension in their processed form."""
    archive_formats = tuple([".gz", ".xz", ".bz2"])
    if url.endswith(".zip"):
        processor=_Unzip()
        ext = ".unzip"
    elif url.endswith(archive_formats):
        processor=_Decompress()
        ext = ".decomp"
    else:
        processor = None
        ext = ""
    return processor, ext


def _path_of_cached_file(url):
    """Determine the full path to the cache where the file in `url` 
    will be downloaded to."""
    cached_filename = _pooch.utils.unique_file_name(url)
    path = _pooch.utils.cache_location(
        _os_cache('gplately'), 
        env=None, 
        version=None
    )
    _pooch.utils.make_local_storage(path)
    full_path = path.resolve() / cached_filename
    return full_path


def _extract_processed_files(processed_directory):
    """Return a list of all full filenames from a given directory 
    in the GPlately cache.
    """
    if _os.path.isdir(processed_directory):
        fnames = []
        for root, dirs, files in _os.walk(processed_directory):
            for file in files:
                fnames.append(_os.path.join(root,file))
        return(fnames)
    elif _os.path.isfile(str(processed_directory)):
        return(processed_directory)

    
def _get_url_etag(url):
    """Obtain the E-Tag of a web server URL. 
    
    The E-Tag identifies a resource under a URL. If the resource
    is modified, a new E-Tag is generated. DataServer uses the
    E-Tag to determine whether local copies of plate model files
    available from a web server need to be updated to match the
    version on the web server.
    """
    # Determine E-Tag of the URL
    etag = str(_requests.head(url).headers.get("ETag"))

    # Determine the filename of an E-Tag txt file
    md5 = _hashlib.md5(url.encode()).hexdigest()
    fname = _pooch.utils.parse_url(url)["path"].split("/")[-1]
    fname = fname[-(255 - len(md5) - 1) :]
    unique_name = f"{md5}-{fname}".split(".")[0]+"-ETAG.txt"
    cachepath = str(
        _pathlib.Path(
            _os.path.expanduser(str(_os_cache('gplately'))))
    )
    text_path = cachepath+"/"+unique_name 
    return(etag, text_path)
    
    
def _save_url_etag_to_txt(etag, text_path):
    """Write an E-Tag to a text file.
    """       
    # Write E-Tag to a text file on the GPlately cache. 
    text_file = open(text_path, "w")
    text_file.write(etag)
    text_file.close()
    
    
def _first_time_download_from_web(url):
    """
    # Provided a web connection to a server can be established,
    download the files from the URL into the GPlately cache.
    """
    if _test_internet_connection(url):
        fnames = _retrieve(
                url=url,
                known_hash=None,  
                downloader=_HTTPDownloader(progressbar=True),
                path=_os_cache('gplately'),
                processor=_determine_processor(url)[0]
        )
        # Get the URL's E-Tag for the first time
        etag, textfilename = _get_url_etag(url)
        _save_url_etag_to_txt(etag, textfilename)
        return(fnames, etag, textfilename)
    
    
def download_from_web(url, download_changes=True):
    """Download a file from a `url` into the `gplately` cache.
    
    Notes
    -----
    After the file belonging to the given `url` is downloaded 
    to the `gplately` cache once, subsequent runs of 
    `download_from_web` with this `url` will not redownload 
    the file as long as:
    
    * The file has not been updated on the web server,
    * The file has not been removed from the `gplately` cache.
    
    Instead, the file will be re-accessed from the `gplately` 
    cache it was downloaded to. 
    
    However, if the file has been updated on the web server,
    `download_from_web` overwrites the cached file with the 
    updated version. The following messages will be displayed 
    to the user:
    
        "Checking whether the requested files need to be updated..."
        "Yes - updating requested files..."
        "Requested files downloaded to the GPlately cache folder!"
    
    If ever a connection to the web server (and the file(s)) in 
    `url` is unsuccessful, this is likely because:

    * An internet connection could not be established; or
    * The `url` passed to `download_from_web` is incorrect

    In either case, `download_from_web` attempts to find a version 
    of the requested file(s) in `url` already stored in the 
    `gplately` cache (assuming it has been downloaded from the same 
    `url` once before). This version may not match the one on the web
    server. If a copy of the file(s) cannot be found in the `gplately`
    cache, a `ConnectionError` is raised.
    
    Parameters
    ----------
    url : str
        The full URL used to download a file from a public web server
        like webDAV.
    download_changes : bool, default=True
        Permit the re-downloading/update of the file from `url` if 
        it has been updated on the web server since the last download.
        
    Returns
    -------
    fnames : list of str
        A list of strings representing the full paths to all cached data
        downloaded from the given `url`.

    Raises
    ------
    ConnectionError
        If a connection to the web server and file(s) in the given `url` is
        unsuccessful (because there is no internet access, and/or the `url`
        is incorrect) and no version of the requested file(s) have been
        cached before. In this case, nothing is returned. 
    """
    # Identify the final filename from the given url within the GPlately cache
    full_path = _path_of_cached_file(url)
    
    # If the files are not yet on the cache,
    if not _os.path.isfile(str(full_path)):
        
        # ...and if a connection to the web server can be established,
        # download files from the URL and create a textfile for this URL's E-Tag
        if _test_internet_connection(url):
            fnames, etag, textfilename = _first_time_download_from_web(url)
            print("Requested files downloaded to the GPlately cache folder!")
            return(fnames)
    
        # ... if a connection to the web server cannot be established
        else:
            raise ConnectionError(
                    "A connection to {} could not be made. Please check your internet connection and/or ensure the URL is correct. No file from the given URL has been cached to {} yet - nothing has been returned.".format(url, full_path.parent))

    # If the files have been downloaded before...
    else:
        #... and if a connection to the web server can be made...
        if _test_internet_connection(url):
            
            # If the newest version of the files in `url` must be cached 
            # at all times, perform E-Tag comparisons:
            if download_changes:

                # Get the path to the file's E-Tag textfile
                local_etag_txtfile = str(full_path).split(".")[0]+"-ETAG.txt"

                # If an e-tag text file does not exist, erase the cached files
                # and download the latest version from the web server. This, in turn,
                # creates an e-tag textfile for this version.
                if not _os.path.isfile(local_etag_txtfile):
                    _os.remove(str(full_path))
                    fnames, etag, local_etag_txtfile = _first_time_download_from_web(url)
                    return(fnames)

                # If the e-tag textfile exists for the local files,
                else: 
                    print("Checking whether the requested files need to be updated...")

                    # Determine the local file's URL e-tag from the textfile
                    with open(local_etag_txtfile) as f:
                        local_etag = str(f.readlines()[0])

                    # Get the e-tag of the web server URL at current time
                    remote_etag, remote_etag_textfile = _get_url_etag(url)

                    # If the local and remote e-tags are unequal, the web-server URL 
                    # contains an updated version of the cached files.
                    if str(remote_etag) != str(local_etag):
                        print("Yes - updating requested files...")
                        
                        # Update the e-tag textfile with this newly-identified URL e-tag
                        _save_url_etag_to_txt(remote_etag, remote_etag_textfile)

                        # Re-download the file, and process it if need-be.
                        with _pooch.utils.temporary_file(path=str(full_path.parent)) as tmp:
                            downloader = _HTTPDownloader(progressbar=True)
                            downloader(url, tmp, _pooch) 
                            _shutil.move(tmp, str(full_path))
                            processor=_determine_processor(url)[0]
                            processor(str(full_path), "update", None)
                            # determine_processor gives the file its processed filename
                            processed_file = str(full_path)+_determine_processor(url)[1]

                        print("Requested files downloaded to the GPlately cache folder!")
                        return(_extract_processed_files(processed_file))

                    # If the e-tags are equal, the local and remote files are the same.
                    # Just return the file(s) as-is.
                    else:
                        print("Requested files are up-to-date!")
                        return(_extract_processed_files(
                            str(full_path)+_determine_processor(url)[1]))
            
            # If file versioning doesn't matter, just keep returning the cached files.
            else:
                fnames, etag, local_etag_txtfile = _first_time_download_from_web(url)
                return(fnames)
                
        # If a connection to the web server could not be made, and the files exist in
        # the GPlately cache, just return the files as-is.
        else:
            print("No connection to {} established. The requested file(s) (potentially older versions) exist in the GPlately cache ({}) and have been returned.".format(url, full_path.parent))
            return(_extract_processed_files(
                    str(full_path)+_determine_processor(url)[1]))



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
    if strings_to_include is not None:
        for f in fnames:
            f_splitted = f.split("/")[-1]
            check = [s for s in strings_to_include if s.lower() in f_splitted.lower()]
            if check:
                sorted_fnames.append(f)
    else:
        sorted_fnames = fnames
    
    if strings_to_ignore is not None:
        more_sorted = []
        for f in sorted_fnames:
            f_splitted = f.split("/")[-1]
            check = [s for s in strings_to_ignore if s.lower() in f_splitted.lower()]
            if not check:
                more_sorted.append(f)
        return(more_sorted)
    else:
        return(sorted_fnames)


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
        "clennett2020",
        "johansson2018",
        "whittaker2015"
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
    """Uses [Pooch](https://www.fatiando.org/pooch/latest/) to download plate reconstruction 
    feature data from plate models and other studies that are stored on web servers 
    (e.g. EarthByte's [webDAV server](https://www.earthbyte.org/webdav/ftp/Data_Collections/)). 
    
    If the `DataServer` object and its methods are called for the first time, i.e. by:

        # string identifier to access the Muller et al. 2019 model
        gDownload = gplately.download.DataServer("Muller2019")

    all requested files are downloaded into the user's 'gplately' cache folder only _once_. If the same
    object and method(s) are re-run, the files will be re-accessed from the cache provided they have not been 
    moved or deleted. 

    Currently, `DataServer` supports a number of plate reconstruction models. To call the object,
    supply a `file_collection` string from one of the following models:

    * __[Müller et al. 2019](https://www.earthbyte.org/muller-et-al-2019-deforming-plate-reconstruction-and-seafloor-age-grids-tectonics/):__ 

        file_collection = `Muller2019`
    
        Information
        -----------
        * Downloadable files: a `rotation_model`, `topology_features`, `static_polygons`, `coastlines`, `continents`, `COBs`, and
        seafloor `age_grids` from 0 to 250 Ma. 
        * Maximum reconstruction time: 250 Ma

        Citations
        ---------
        Müller, R. D., Zahirovic, S., Williams, S. E., Cannon, J., Seton, M., 
        Bower, D. J., Tetley, M. G., Heine, C., Le Breton, E., Liu, S., Russell, S. H. J., 
        Yang, T., Leonard, J., and Gurnis, M. (2019), A global plate model including 
        lithospheric deformation along major rifts and orogens since the Triassic. 
        Tectonics, vol. 38, https://doi.org/10.1029/2018TC005462.
            

    * __Müller et al. 2016__:

        file_collection = `Muller2016`

        Information
        -----------
        * Downloadable files: a `rotation_model`, `topology_features`, `static_polygons`, `coastlines`, and
        seafloor `age_grids` from 0-230 Ma. 
        * Maximum reconstruction time: 230 Ma

        Citations
        ---------
        * Müller R.D., Seton, M., Zahirovic, S., Williams, S.E., Matthews, K.J.,
        Wright, N.M., Shephard, G.E., Maloney, K.T., Barnett-Moore, N., Hosseinpour, M., 
        Bower, D.J., Cannon, J., InPress. Ocean basin evolution and global-scale plate 
        reorganization events since Pangea breakup, Annual Review of Earth and Planetary 
        Sciences, Vol 44, 107-138. DOI: 10.1146/annurev-earth-060115-012211.


    * __[Merdith et al. 2021](https://zenodo.org/record/4485738#.Yhrm8hNBzA0)__: 

        file_collection = `Merdith2021`

        Information
        -----------
        * Downloadable files: a `rotation_model`, `topology_features`, `static_polygons`, `coastlines`
        and `continents`.
        * Maximum reconstruction time: 1 Ga (however, `PlotTopologies` correctly visualises up to 410 Ma) 

        Citations: 
        ----------
        * Merdith et al. (in review), 'A continuous, kinematic full-plate motion model
        from 1 Ga to present'. 
        * Andrew Merdith. (2020). Plate model for 'Extending Full-Plate Tectonic Models 
        into Deep Time: Linking the Neoproterozoic and the Phanerozoic ' (1.1b) [Data set]. 
        Zenodo. https://doi.org/10.5281/zenodo.4485738


    * __Cao et al. 2020__: 

        file_collection = `Cao2020`

        Information
        -----------
        * Downloadable files: `rotation_model`, `topology_features`, `static_polygons`, `coastlines`
        and `continents`.
        * Maximum reconstruction time: 1 Ga

        Citations
        ---------
        * Toy Billion-year reconstructions from Cao et al (2020). 
        Coupled Evolution of Plate Tectonics and Basal Mantle Structure Tectonics, 
        doi: 10.1029/2020GC009244


    - __Clennett et al. 2020__ : 

        file_collection = `Clennett2020`
        
        Information
        -----------
        * Downloadable files: `rotation_model`, `topology_features`, `static_polygons`, `coastlines`
        and `continents`
        * Maximum reconstruction time: 170 Ma

        Citations
        ---------
        * Mather, B., Müller, R.D.,; Alfonso, C.P., Seton, M., 2021, Kimberlite eruption 
        driven by slab flux and subduction angle. DOI: 10.5281/zenodo.5769002


    - __Seton et al. 2012__:

        file_collection = `Seton2012`

        Information
        -----------
        * Downloadable files: `rotation_model`, `topology_features`, `coastlines`,
        `COBs`, and paleo-age grids (0-200 Ma)
        * Maximum reconstruction time: 200 Ma

        Citations
        ---------
        * M. Seton, R.D. Müller, S. Zahirovic, C. Gaina, T.H. Torsvik, G. Shephard, A. Talsma, 
        M. Gurnis, M. Turner, S. Maus, M. Chandler, Global continental and ocean basin reconstructions 
        since 200 Ma, Earth-Science Reviews, Volume 113, Issues 3-4, July 2012, Pages 212-270, 
        ISSN 0012-8252, 10.1016/j.earscirev.2012.03.002.


    - __Matthews et al. 2016__: 

        file_collection = `Matthews2016`

        Information
        -----------
        * Downloadable files: `rotation_model`, `topology_features`, `static_polygons`, `coastlines`,
        and `continents`
        * Maximum reconstruction time(s): 410-250 Ma, 250-0 Ma

        Citations
        ---------
        * Matthews, K.J., Maloney, K.T., Zahirovic, S., Williams, S.E., Seton, M.,
        and Müller, R.D. (2016). Global plate boundary evolution and kinematics since the 
        late Paleozoic, Global and Planetary Change, 146, 226-250. 
        DOI: 10.1016/j.gloplacha.2016.10.002


    - __Merdith et al. 2017__: 

        file_collection = `Merdith2017`

        Information
        -----------
        * Downloadable files: `rotation_files` and `topology_features`
        * Maximum reconstruction time: 410 Ma

        Citations
        ---------
        * Merdith, A., Collins, A., Williams, S., Pisarevskiy, S., Foden, J., Archibald, D. 
        and Blades, M. et al. 2016. A full-plate global reconstruction of the Neoproterozoic. 
        Gondwana Research. 50: pp. 84-134. DOI: 10.1016/j.gr.2017.04.001


    - __Li et al. 2008__: 

        file_collection = `Li2008`

        Information
        -----------
        * Downloadable files: `rotation_model` and `static_polygons`
        * Maximum reconstruction time: 410 Ma

        Citations
        ---------
        * Rodinia reconstruction from Li et al (2008), Assembly, configuration, and break-up 
        history of Rodinia: A synthesis. Precambrian Research. 160. 179-210. 
        DOI: 10.1016/j.precamres.2007.04.021.


    - __Pehrsson et al. 2015__: 

        file_collection = `Pehrsson2015`

        Information
        -----------
        * Downloadable files: `rotation_model` and `static_polygons`
        * Maximum reconstruction time: N/A

        Citations
        ---------
        * Pehrsson, S.J., Eglington, B.M., Evans, D.A.D., Huston, D., and Reddy, S.M., (2015),
        Metallogeny and its link to orogenic style during the Nuna supercontinent cycle. Geological 
        Society, London, Special Publications, 424, 83-94. DOI: https://doi.org/10.1144/SP424.5


    - __Torsvik and Cocks et al. 2017__: 

        file_collection = `TorsvikCocks2017`

        Information
        -----------
        * Downloadable files: `rotation_model`, and `coastlines`
        * Maximum reconstruction time: 410 Ma

        Citations
        ---------
        * Torsvik, T., & Cocks, L. (2016). Earth History and Palaeogeography. Cambridge: 
        Cambridge University Press. doi:10.1017/9781316225523


    - __Young et al. 2019__: 

        file_collection = `Young2019`

        Information
        -----------
        * Downloadable files: `rotation_model`, `topology_features`, `static_polygons`, `coastlines`
        and `continents`.
        * Maximum reconstruction time: 410-250 Ma, 250-0 Ma
        
        Citations
        ---------
        * Young, A., Flament, N., Maloney, K., Williams, S., Matthews, K., Zahirovic, S.,
        Müller, R.D., (2019), Global kinematics of tectonic plates and subduction zones since the late 
        Paleozoic Era, Geoscience Frontiers, Volume 10, Issue 3, pp. 989-1013, ISSN 1674-9871,
        DOI: https://doi.org/10.1016/j.gsf.2018.05.011.


    - __Scotese et al. 2008__: 

        file_collection = `Scotese2008`

        Information
        -----------
        * Downloadable files: `rotation_model`, `static_polygons`, and `continents`
        * Maximum reconstruction time: 
        
        Citations
        ---------
        * Scotese, C.R. 2008. The PALEOMAP Project PaleoAtlas for ArcGIS, Volume 2, Cretaceous 
        paleogeographic and plate tectonic reconstructions. PALEOMAP Project, Arlington, Texas.


    - __Golonka et al. 2007__: 

        file_collection = `Golonka2007`

        Information
        -----------
        * Downloadable files: `rotation_model`, `static_polygons`, and `continents`
        * Maximum reconstruction time: 410 Ma
        
        Citations
        ---------
        * Golonka, J. 2007. Late Triassic and Early Jurassic palaeogeography of the world. 
        Palaeogeography, Palaeoclimatology, Palaeoecology 244(1–4), 297–307.


    - __Clennett et al. 2020 (based on Müller et al. 2019)__: 

        file_collection = `Clennett2020_M2019`

        Information
        -----------
        * Downloadable files: `rotation_model`, `topology_features`, `continents` and `coastlines`
        * Maximum reconstruction time: 250 Ma
        
        Citations
        ---------
        * Clennett, E.J., Sigloch, K., Mihalynuk, M.G., Seton, M., Henderson, M.A., Hosseini, K.,
        Mohammadzaheri, A., Johnston, S.T., Müller, R.D., (2020), A Quantitative Tomotectonic Plate 
        Reconstruction of Western North America and the Eastern Pacific Basin. Geochemistry, Geophysics, 
        Geosystems, 21, e2020GC009117. DOI: https://doi.org/10.1029/2020GC009117


    - __Clennett et al. 2020 (rigid topological model based on Shephard et al, 2013)__: 

        file_collection = `Clennett2020_S2013`

        Information
        -----------
        * Downloadable files: `rotation_model`, `topology_features`, `continents` and `coastlines`
        * Maximum reconstruction time: 
        
        Citations
        ---------
        * Clennett, E.J., Sigloch, K., Mihalynuk, M.G., Seton, M., Henderson, M.A., Hosseini, K.,
        Mohammadzaheri, A., Johnston, S.T., Müller, R.D., (2020), A Quantitative Tomotectonic Plate 
        Reconstruction of Western North America and the Eastern Pacific Basin. Geochemistry, Geophysics, 
        Geosystems, 21, e2020GC009117. DOI: https://doi.org/10.1029/2020GC009117

    """
    def __init__(self, file_collection):

        self.file_collection = file_collection
        self.data_collection = DataCollection(self.file_collection)


    def get_plate_reconstruction_files(self):
        """Downloads and constructs a `rotation model`, a set of `topology_features` and
        and a set of `static_polygons` needed to call the `PlateReconstruction` object.

        Returns
        -------
        rotation_model : instance of <pygplates.RotationModel>
            A rotation model to query equivalent and/or relative topological plate rotations
            from a time in the past relative to another time in the past or to present day.
        topology_features : instance of <pygplates.FeatureCollection>
            Point, polyline and/or polygon feature data that are reconstructable through 
            geological time.
        static_polygons : instance of <pygplates.FeatureCollection>
            Present-day polygons whose shapes do not change through geological time. They are
            used to cookie-cut dynamic polygons into identifiable topological plates (assigned 
            an ID) according to their present-day locations.

        Notes
        -----
        This method accesses the plate reconstruction model ascribed to the `file_collection` string passed 
        into the `DataServer` object. For example, if the object was called with `"Muller2019"`:

            gDownload = gplately.download.DataServer("Muller2019")
            rotation_model, topology_features, static_polygons = gDownload.get_plate_reconstruction_files()

        the method will download a `rotation_model`, `topology_features` and `static_polygons` from the 
        Müller et al. (2019) plate reconstruction model. Once the reconstruction objects are returned, 
        they can be passed into:

            model = gplately.reconstruction.PlateReconstruction(rotation_model, topology_features, static_polygons)

        * Note: If the requested plate model does not have a certain file(s), a message will be printed 
        to alert the user. For example, using `get_plate_reconstruction_files()`
        for the Torsvik and Cocks (2017) plate reconstruction model yields the printed message:

                No topology features in TorsvikCocks2017. No FeatureCollection created - unable to 
                plot trenches, ridges and transforms.
                No continent-ocean boundaries in TorsvikCocks2017.

        """

        rotation_filenames = []
        rotation_model = []
        topology_filenames = []
        topology_features = _pygplates.FeatureCollection()
        static_polygons= _pygplates.FeatureCollection()
        static_polygon_filenames = []

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
                        download_from_web(url[0]), self.file_collection
                    )
                    rotation_filenames = _collect_file_extension(
                        _str_in_folder(
                            _str_in_filename(fnames,
                                strings_to_ignore=DataCollection.rotation_strings_to_ignore(self)
                            ),
                        strings_to_ignore=DataCollection.rotation_strings_to_ignore(self)
                        ),
                        [".rot"]
                    )
                    #print(rotation_filenames)
                    rotation_model = _pygplates.RotationModel(rotation_filenames)

                    topology_filenames = _collect_file_extension(
                        _str_in_folder(
                            _str_in_filename(fnames, 
                                strings_to_include=DataCollection.dynamic_polygon_strings_to_include(self),
                                strings_to_ignore=DataCollection.dynamic_polygon_strings_to_ignore(self)
                            ), 
                            strings_to_ignore=DataCollection.dynamic_polygon_strings_to_ignore(self)
                        ),
                        [".gpml", ".gpmlz"]
                    )
                    #print(topology_filenames)
                    for file in topology_filenames:
                        topology_features.add(_pygplates.FeatureCollection(file))

                    static_polygon_filenames = _check_gpml_or_shp(
                        _str_in_folder(
                            _str_in_filename(fnames, 
                                strings_to_include=DataCollection.static_polygon_strings_to_include(self),
                                strings_to_ignore=DataCollection.static_polygon_strings_to_ignore(self)
                            ),
                            strings_to_ignore=DataCollection.static_polygon_strings_to_ignore(self)
                        )
                    )
                    #print(static_polygon_filenames)
                    for stat in static_polygon_filenames:
                        static_polygons.add(_pygplates.FeatureCollection(stat))

                else:
                    for file in url[0]:
                        rotation_filenames.append(
                            _collect_file_extension(
                                download_from_web(file), [".rot"])
                        )
                        rotation_model = _pygplates.RotationModel(rotation_filenames)

                    for file in url[1]:
                        topology_filenames.append(
                            _collect_file_extension(
                                download_from_web(file), [".gpml"])
                        )
                        for file in topology_filenames:
                            topology_features.add(
                                _pygplates.FeatureCollection(file)
                            )

                    for file in url[2]:
                        static_polygon_filenames.append(
                            _check_gpml_or_shp(
                                _str_in_folder(
                                    _str_in_filename(download_from_web(url[0]), 
                                        strings_to_include=DataCollection.static_polygon_strings_to_include(self)
                                    ),    
                                        strings_to_ignore=DataCollection.static_polygon_strings_to_ignore(self)
                                )
                            )   
                        )
                        for stat in static_polygon_filenames:
                            static_polygons.add(_pygplates.FeatureCollection(stat))
                break

        if found_collection is False:
            raise ValueError("{} is not in GPlately's DataServer.".format(self.file_collection))

        if not rotation_filenames:
            print("No .rot files in {}. No rotation model created.".format(self.file_collection))
            rotation_model = []
        if not topology_filenames:
            print("No topology features in {}. No FeatureCollection created - unable to plot trenches, ridges and transforms.".format(self.file_collection))
            topology_features = []
        if not static_polygons:
            print("No static polygons in {}.".format(self.file_collection))
            static_polygons = []

        return rotation_model, topology_features, static_polygons


    def get_topology_geometries(self):
        """Uses Pooch to download coastline, continent and COB (continent-ocean boundary)
        Shapely geometries from the requested plate model. These are needed to call the `PlotTopologies`
        object and visualise topological plates through time.

        Returns
        -------
        coastlines : instance of <pygplates.FeatureCollection>
            Present-day global coastline Shapely polylines cookie-cut using static polygons. Ready for
            reconstruction to a particular geological time and for plotting.

        continents : instance of <pygplates.FeatureCollection>
            Cookie-cutting Shapely polygons for non-oceanic regions (continents, inta-oceanic arcs, etc.)
            ready for reconstruction to a particular geological time and for plotting.

        COBs : instance of <pygplates.FeatureCollection>
            Shapely polylines resolved from .shp and/or .gpml topology files that represent the 
            locations of the boundaries between oceanic and continental crust.
            Ready for reconstruction to a particular geological time and for plotting.

        Notes
        -----
        This method accesses the plate reconstruction model ascribed to the `file_collection` 
        string passed into the `DataServer` object. For example, if the object was called with
        `"Muller2019"`:

            gDownload = gplately.download.DataServer("Muller2019")
            coastlines, continents, COBs = gDownload.get_topology_geometries()

        the method will attempt to download `coastlines`, `continents` and `COBs` from the Müller
        et al. (2019) plate reconstruction model. If found, these files are returned as individual 
        pyGPlates Feature Collections. They can be passed into:

            gPlot = gplately.plot.PlotTopologies(gplately.reconstruction.PlateReconstruction, time, continents, coastlines, COBs)

        to reconstruct features to a certain geological time. The `PlotTopologies`
        object provides simple methods to plot these geometries along with trenches, ridges and 
        transforms (see documentation for more info). Note that the `PlateReconstruction` object 
        is a parameter.

        * Note: If the requested plate model does not have a certain geometry, a
        message will be printed to alert the user. For example, if `get_topology_geometries()` 
        is used with the `"Matthews2016"` plate model, the workflow will print the following 
        message: 

                No continent-ocean boundaries in Matthews2016.
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
                    # Some plate models do not have reconstructable geometries i.e. Li et al. 2008
                    if url[0] is None:
                        break
                    else:
                        fnames = _collection_sorter(
                            download_from_web(url[0]), self.file_collection
                        )
                        coastlines = _check_gpml_or_shp(
                            _str_in_folder(
                                _str_in_filename(
                                    fnames,
                                    strings_to_include=DataCollection.coastline_strings_to_include(self),
                                    strings_to_ignore=DataCollection.coastline_strings_to_ignore(self)
                                ), 
                                strings_to_ignore=DataCollection.coastline_strings_to_ignore(self)
                            )
                        )
                        continents = _check_gpml_or_shp(
                            _str_in_folder(
                                _str_in_filename(
                                    fnames, 
                                    strings_to_include=DataCollection.continent_strings_to_include(self),
                                    strings_to_ignore=DataCollection.continent_strings_to_ignore(self)
                                ), 
                                strings_to_ignore=DataCollection.continent_strings_to_ignore(self)
                            )
                        )
                        COBs = _check_gpml_or_shp(
                            _str_in_folder(
                                _str_in_filename(
                                    fnames,
                                    strings_to_include=DataCollection.COB_strings_to_include(self),
                                    strings_to_ignore=DataCollection.COB_strings_to_ignore(self)
                                ), 
                                strings_to_ignore=DataCollection.COB_strings_to_ignore(self)
                            )
                        )
                else:
                    for file in url[0]:
                        if url[0] is not None:
                            coastlines.append(_str_in_filename(
                                download_from_web(file), 
                                strings_to_include=["coastline"])
                            )
                            coastlines = _check_gpml_or_shp(coastlines)
                        else:
                            coastlines = []

                    for file in url[1]:
                        if url[1] is not None:
                            continents.append(_str_in_filename(
                                download_from_web(file), 
                                strings_to_include=["continent"])
                            )
                            continents = _check_gpml_or_shp(continents)
                        else:
                            continents = []

                    for file in url[2]:
                        if url[2] is not None:
                            COBs.append(_str_in_filename(
                                download_from_web(file), 
                                strings_to_include=["cob"])
                            )
                            COBs = _check_gpml_or_shp(COBs)
                        else:
                            COBs = []
                break

        if found_collection is False:
            raise ValueError("{} is not in GPlately's DataServer.".format(self.file_collection))

        if not coastlines:
            print("No coastlines in {}.".format(self.file_collection))
            coastlines_featurecollection = []
        else:
            #print(coastlines)
            coastlines_featurecollection = _pygplates.FeatureCollection()
            for coastline in coastlines:
                coastlines_featurecollection.add(_pygplates.FeatureCollection(coastline))
        
        if not continents:
            print("No continents in {}.".format(self.file_collection))
            continents_featurecollection = []
        else:
            #print(continents)
            continents_featurecollection = _pygplates.FeatureCollection()
            for continent in continents:
                continents_featurecollection.add(_pygplates.FeatureCollection(continent))
        
        if not COBs:
            print("No continent-ocean boundaries in {}.".format(self.file_collection))
            COBs_featurecollection = []
        else:
            #print(COBs)
            COBs_featurecollection = _pygplates.FeatureCollection()
            for COB in COBs:
                COBs_featurecollection.add(_pygplates.FeatureCollection(COB))
        
        geometries = coastlines_featurecollection, continents_featurecollection, COBs_featurecollection
        return geometries


    def get_age_grid(self, time):
        """Downloads seafloor and paleo-age grids from the plate reconstruction model (`file_collection`)
        passed into the `DataServer` object. Stores grids in the "gplately" cache.

        Currently, `DataServer` supports the following age grids:

        * __Muller et al. 2019__

            * `file_collection` = `Muller2019`
            * Time range: 0-250 Ma
            * Seafloor age grid rasters in netCDF format.

        * __Muller et al. 2016__
            
            * `file_collection` = `Muller2016`
            * Time range: 0-240 Ma
            * Seafloor age grid rasters in netCDF format. 

        * __Seton et al. 2012__

            * `file_collection` = `Seton2012`
            * Time range: 0-200 Ma
            * Paleo-age grid rasters in netCDF format.

        
        Parameters
        ----------
        time : int, or list of int, default=None
            Request an age grid from one (an integer) or multiple reconstruction times (a
            list of integers).

        Returns
        -------
        raster_array : MaskedArray
            A masked array containing the netCDF4 age grid ready for plotting or for
            passing into GPlately's `Raster` object for raster manipulation.

        Raises
        -----
        ValueError
            If `time` (a single integer, or a list of integers representing reconstruction
            times to extract the age grids from) is not passed.

        Notes
        -----
        The first time that `get_age_grid` is called for a specific time(s), the age grid(s) 
        will be downloaded into the GPlately cache once. Upon successive calls of `get_age_grid`
        for the same reconstruction time(s), the age grids will not be re-downloaded; rather, 
        they are re-accessed from the same cache provided the age grid(s) have not been moved or deleted. 

        Examples
        --------
        if the `DataServer` object was called with the `Muller2019` `file_collection` string:

            gDownload = gplately.download.DataServer("Muller2019")

        `get_age_grid` will download seafloor age grids from the Müller et al. (2019) plate 
        reconstruction model for the geological time(s) requested in the `time` parameter. 
        If found, these age grids are returned as masked arrays. 

        For example, to download  Müller et al. (2019) seafloor age grids for 0Ma, 1Ma and
        100 Ma:

            age_grids = gDownload.get_age_grid([0, 1, 100])
            
        """
        age_grids = []
        age_grid_links = DataCollection.netcdf4_age_grids(self, time)
        for link in age_grid_links:
            age_grid_file = download_from_web(link)
            age_grid = _gplately.grids.read_netcdf_grid(age_grid_file)
            age_grids.append(age_grid)

        if not age_grids:
            raise ValueError("{} netCDF4 age grids not found.".format(self.file_collection))

        if len(age_grids) == 1:
            return age_grids[0]
        else: 
            return age_grids


    def get_spreading_rate_grid(self, time):
        """Downloads seafloor spreading rate grids from the plate reconstruction 
        model (`file_collection`) passed into the `DataServer` object. Stores 
        grids in the "gplately" cache.

        Currently, `DataServer` supports spreading rate grids from the following plate
        models:

        * __Clennett et al. 2020__

            * `file_collection` = `Clennett2020`
            * Time range: 0-250 Ma
            * Seafloor spreading rate grids in netCDF format.

        
        Parameters
        ----------
        time : int, or list of int, default=None
            Request a spreading grid from one (an integer) or multiple reconstruction 
            times (a list of integers).

        Returns
        -------
        raster_array : MaskedArray
            A masked array containing the netCDF4 spreading rate grid ready for 
            plotting or for passing into GPlately's `Raster` object.

        Raises
        -----
        ValueError
            If `time` (a single integer, or a list of integers representing reconstruction
            times to extract the spreading rate grids from) is not passed.

        Notes
        -----
        The first time that `get_spreading_rate_grid` is called for a specific time(s), 
        the spreading rate grid(s) will be downloaded into the GPlately cache once. 
        Upon successive calls of `get_spreading_rate_grid` for the same reconstruction 
        time(s), the grids will not be re-downloaded; rather, they are re-accessed from 
        the same cache location provided they have not been moved or deleted. 

        Examples
        --------
        if the `DataServer` object was called with the `Clennett2020` `file_collection` string:

            gDownload = gplately.download.DataServer("Clennett2020")

        `get_spreading_rate_grid` will download seafloor spreading rate grids from the 
        Clennett et al. (2020) plate reconstruction model for the geological time(s) 
        requested in the `time` parameter. When found, these spreading rate grids are 
        returned as masked arrays. 

        For example, to download Clennett et al. (2020) seafloor spreading rate grids for 
        0Ma, 1Ma and 100 Ma as MaskedArray objects:

            spreading_rate_grids = gDownload.get_spreading_rate_grid([0, 1, 100])
            
        """
        spreading_rate_grids = []
        spreading_rate_grid_links = DataCollection.netcdf4_spreading_rate_grids(self, time)
        for link in spreading_rate_grid_links:
            spreading_rate_grid_file = download_from_web(link)
            spreading_rate_grid = _gplately.grids.read_netcdf_grid(spreading_rate_grid_file)
            spreading_rate_grids.append(spreading_rate_grid)

        if not spreading_rate_grids:
            raise ValueError("{} netCDF4 seafloor spreading rate grids not found.".format(self.file_collection))

        if len(spreading_rate_grids) == 1:
            return spreading_rate_grids[0]
        else: 
            return spreading_rate_grids


    def get_raster(self, raster_id_string=None):
        """Downloads assorted raster data that are not associated with the plate 
        reconstruction models supported by GPlately's `DataServer`. Stores rasters in the 
        "gplately" cache.

        Currently, `DataServer` supports the following rasters and images:

        * __[ETOPO1](https://www.ngdc.noaa.gov/mgg/global/)__: 
            * Filetypes available : TIF, netCDF (GRD)
            * `raster_id_string` = `"ETOPO1_grd"`, `"ETOPO1_tif"` (depending on the requested format)
            * A 1-arc minute global relief model combining lang topography and ocean bathymetry.
            * Citation: doi:10.7289/V5C8276M


        Parameters
        ----------
        raster_id_string : str, default=None
            A string to identify which raster to download.

        Returns
        -------
        raster_filenames : ndarray
            An ndarray of the cached raster. This can be plotted using `matplotlib.pyplot.imshow` on
            a `cartopy.mpl.GeoAxis` GeoAxesSubplot (see example below).

        Raises
        ------
        ValueError
            * if a `raster_id_string` is not supplied.

        Notes
        -----
        Rasters obtained by this method are (so far) only reconstructed to present-day. 

        Examples
        --------
        To download ETOPO1 and plot it on a Mollweide projection:

            import gplately
            import numpy as np
            import matplotlib.pyplot as plt
            import cartopy.crs as ccrs

            gdownload = gplately.DataServer("Muller2019")
            etopo1 = gdownload.get_raster("ETOPO1_tif")
            fig = plt.figure(figsize=(18,14), dpi=300)
            ax = fig.add_subplot(111, projection=ccrs.Mollweide(central_longitude = -150))
            ax2.imshow(etopo1, extent=[-180,180,-90,90], transform=ccrs.PlateCarree()) 

        """
        from matplotlib import image
        if raster_id_string is None:
            raise ValueError(
                "Please specify which raster to download."
            )
        #filetype = "."+"_".split(raster_id_string)[-1]

        archive_formats = tuple([".gz", ".xz", ".bz2"])
        # Set to true if we find the given collection in database
        found_collection = False
        raster_filenames = []
        database = DataCollection.rasters(self)

        for collection, zip_url in database.items():
            # Isolate the raster name and the file type
            #raster_name = collection.split("_")[0]
            #raster_type = "."+collection.split("_")[-1]
            if (raster_id_string.lower() == collection.lower()):
                raster_filenames = download_from_web(zip_url[0])
                found_collection = True
                break

        if found_collection is False:
            raise ValueError("{} not in collection database.".format(raster_id_string))
        else:
            raster_matrix = image.imread(raster_filenames)
        return raster_matrix


    def get_feature_data(self, feature_data_id_string=None):
        """Downloads assorted geological feature data from web servers (i.e. 
        [GPlates 2.3 sample data](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/))
        into the "gplately" cache.

        Currently, `DataServer` supports the following feature data:

        * __Large igneous provinces from Johansson et al. (2018)__

            Information
            -----------
            * Formats: .gpmlz
            * `feature_data_id_string` = `Johansson2018`

            Citations
            ---------
            * Johansson, L., Zahirovic, S., and Müller, R. D., In Prep, The 
            interplay between the eruption and weathering of Large Igneous Provinces and 
            the deep-time carbon cycle: Geophysical Research Letters.


        - __Large igneous province products interpreted as plume products from Whittaker 
        et al. (2015)__.

            Information
            -----------
            * Formats: .gpmlz, .shp
            * `feature_data_id_string` = `Whittaker2015`
            
            Citations
            ---------
            * Whittaker, J. M., Afonso, J. C., Masterton, S., Müller, R. D., 
            Wessel, P., Williams, S. E., & Seton, M. (2015). Long-term interaction between 
            mid-ocean ridges and mantle plumes. Nature Geoscience, 8(6), 479-483. 
            doi:10.1038/ngeo2437.


        - __Seafloor tectonic fabric (fracture zones, discordant zones, V-shaped structures, 
        unclassified V-anomalies, propagating ridge lineations and extinct ridges) from 
        Matthews et al. (2011)__

            Information
            -----------
            * Formats: .gpml
            * `feature_data_id_string` = `SeafloorFabric`

            Citations
            ---------
            * Matthews, K.J., Müller, R.D., Wessel, P. and Whittaker, J.M., 2011. The 
            tectonic fabric of the ocean basins. Journal of Geophysical Research, 116(B12): 
            B12109, DOI: 10.1029/2011JB008413. 


        - __Present day surface hotspot/plume locations from Whittaker et al. (2013)__

            Information
            -----------
            * Formats: .gpmlz
            * `feature_data_id_string` = `Hotspots`

            Citation
            --------
            * Whittaker, J., Afonso, J., Masterton, S., Müller, R., Wessel, P., 
            Williams, S., and Seton, M., 2015, Long-term interaction between mid-ocean ridges and 
            mantle plumes: Nature Geoscience, v. 8, no. 6, p. 479-483, doi:10.1038/ngeo2437.

        
        Parameters
        ----------
        feature_data_id_string : str, default=None
            A string to identify which feature data to download to the cache (see list of supported
            feature data above).

        Returns
        -------
        feature_data_filenames : instance of <pygplates.FeatureCollection>, or list of instance <pygplates.FeatureCollection>
            If a single set of feature data is downloaded, a single pyGPlates `FeatureCollection` 
            object is returned. Otherwise, a list containing multiple pyGPlates `FeatureCollection` 
            objects is returned (like for `SeafloorFabric`). In the latter case, feature reconstruction 
            and plotting may have to be done iteratively.

        Raises
        ------
        ValueError
            If a `feature_data_id_string` is not provided.

        Examples
        --------
        For examples of plotting data downloaded with `get_feature_data`, see GPlately's sample 
        notebook 05 - Working With Feature Geometries [here](https://github.com/GPlates/gplately/blob/master/Notebooks/05-WorkingWithFeatureGeometries.ipynb).
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
                feature_data_filenames = _collection_sorter(
                    _collect_file_extension(
                    download_from_web(zip_url[0]), [".gpml", ".gpmlz"]
                    ),
                    collection
                )

                break

        feat_data = _pygplates.FeatureCollection()
        if len(feature_data_filenames) == 1:
                feat_data.add(_pygplates.FeatureCollection(feature_data_filenames[0]))
                return feat_data
        else:    
            feat_data=[]
            for file in feature_data_filenames:
                feat_data.append(_pygplates.FeatureCollection(file))
            return feat_data
    