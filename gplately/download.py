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

"""
Functions for downloading assorted plate reconstruction data to use with GPlately's
main objects. Files are stored in the user's cache and can be reused after being
downloaded once.

These data have been created and used in plate reconstruction models and studies, and
are available from public web servers (like EarthByte's webDAV server, or the GPlates
2.3 sample dataset library).

"""
import hashlib as _hashlib
import os as _os
import pathlib as _pathlib
import re as _re
import shutil as _shutil
import urllib.request as _request
from typing import Union

import numpy as _np
import numpy as np
import pooch as _pooch
import pygplates
import requests as _requests
from matplotlib import image
from plate_model_manager import PlateModelManager, PresentDayRasterManager
from pooch import Decompress as _Decompress
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve

from .data import _feature_data, _rasters
from .grids import Raster


class DownloadWarning(RuntimeWarning):
    pass


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
        processor = _Unzip()
        ext = ".unzip"
    elif url.endswith(archive_formats):
        processor = _Decompress()
        ext = ".decomp"
    else:
        processor = None
        ext = ""
    return processor, ext


def path_of_cached_file(url, model_name=None):
    """Determine the absolute path where the file(s) from `url`
    will be downloaded to.

    Parameters
    ----------
    url : str
        The full download URL for the file passed as a string.

    model_name : str, default None
        An optional substring that ideally describes/labels the file.
        This will help name the file(s) when it is downloaded to the
        gplately cache (this directory can be
        found using `gplately.download.path_to_cache()`).
    """

    cached_filename = _pooch.utils.unique_file_name(url)
    cached_filename = _remove_hash(cached_filename)
    path = path_to_cache()

    processor_to_use, processor_extension = _determine_processor(url)

    fn = _parse_url_for_filenames(url)
    assert isinstance(fn, str)
    # If the requested files need processing (i.e. zip, gz folders)
    if processor_extension:
        # Are they from plate models? These typically are the .zip folders for plate models
        if model_name:
            cached_filename = str(path) + "/" + model_name + processor_extension + "/"
            unprocessed_path = str(path) + "/" + model_name
            # cached_filename = cached_filename = str(path) + '/' + model_name

        # If not from plate models but need processing, i.e. ETOPO1
        else:

            cached_filename = (
                str(path) + "/" + "gplately_" + fn + processor_extension + "/"
            )
            unprocessed_path = str(path) + "/" + "gplately_" + fn
            # cached_filename = "gplately_"+_parse_url_for_filenames(url)

    # If the requested files do not need processing, like standalone .nc files:
    else:
        if model_name:
            cached_filename = str(path) + "/" + model_name + "_" + fn
            unprocessed_path = None
        else:
            cached_filename = str(path) + "/" + "gplately_" + fn
            unprocessed_path = None

    _pooch.utils.make_local_storage(path)
    full_path = path.resolve() / cached_filename

    return full_path, unprocessed_path


def _extract_processed_files(processed_directory):
    """Return a list of all full filenames from a given directory
    in the GPlately cache.
    """
    if _os.path.isdir(processed_directory):
        fnames = []
        for root, dirs, files in _os.walk(processed_directory):
            for file in files:
                fnames.append(_os.path.join(root, file))
        return fnames
    elif _os.path.isfile(str(processed_directory)):
        return processed_directory


def path_to_cache():
    """Determine the absolute path to the system gplately cache."""
    path = _pooch.utils.cache_location(_os_cache("gplately"), env=None, version=None)
    return path


def clear_cache():
    """Clear the `gplately` cache directory.

    The absolute path of this directory can be found by running
    [gplately.download.path_to_cache()](file:///Users/laurenilano/gplately/api/gplately/download.html#gplately.download.path_to_cache).

    Caution - when called, all files in /path/to/caches/gplately will
    be deleted. This action cannot be undone.
    """
    cache_path = path_to_cache()
    _shutil.rmtree(str(cache_path))
    _pooch.utils.make_local_storage(str(cache_path))
    return


def _parse_url_for_filenames(url, return_hash=False):
    # Determine the filename of an E-Tag txt file
    md5 = _hashlib.md5(url.encode()).hexdigest()
    fname = _os.path.basename(_pooch.utils.parse_url(url)["path"])
    fname = fname[-(255 - len(md5) - 1) :]
    if return_hash:
        return str(fname), str(md5)
    else:
        return str(fname)


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
    parsed_fname, filehash = _parse_url_for_filenames(url, return_hash=True)

    unique_name = filehash + "-ETAG.txt"

    cachepath = str(_pathlib.Path(_os.path.expanduser(str(_os_cache("gplately")))))

    text_path = _os.path.join(cachepath, unique_name)
    return (etag, text_path)


def _save_url_etag_to_txt(etag, text_path):
    """Write an E-Tag to a text file."""
    # Write E-Tag to a text file on the GPlately cache.
    text_file = open(text_path, "w")
    text_file.write(etag)
    text_file.close()


def _match_url_to_extension(url):
    url = str(url)
    if url.endswith(".nc"):
        return ".nc"
    elif url.endswith(".jpg"):
        return ".jpg"
    elif url.endswith(".png"):
        return ".png"
    elif url.endswith(".tif"):
        return ".tif"


def _first_time_download_from_web(url, model_name=None, verbose=True):
    """
    # Provided a web connection to a server can be established,
    download the files from the URL into the GPlately cache.
    """
    logger = _pooch.get_logger()
    log_level = logger.level
    if _test_internet_connection(url):

        if not verbose:
            logger.setLevel("WARNING")

        # The filename pooch saves the requested file is derived from
        # one of four permutations:
        # 1. File is from a plate model and needs processing (i.e. zip --> unzip)
        # 2. File is from a plate model and does not need processing (i.e. .nc age grids)
        # 3. File is not from a plate model but needs processing (i.e. ETOPO, .grd.gz --> .decomp)
        # 4. File is not from a plate model and does not need processing
        processor_to_use, processor_extension = _determine_processor(url)

        fn = _parse_url_for_filenames(url)
        assert isinstance(fn, str)

        # If the requested files need processing (i.e. zip, gz folders)
        if processor_extension:
            # Are they from plate models? These typically are the .zip folders for plate models
            if model_name:
                # Download the files with a naming structure like:
                # /path/to/cache/gplately/model_name+processor_extension
                used_fname = model_name
                fnames = _retrieve(
                    url=url,
                    known_hash=None,
                    downloader=_HTTPDownloader(progressbar=verbose),
                    fname=used_fname,
                    path=_os_cache("gplately"),
                    processor=processor_to_use,
                )
            # If not from plate models but need processing, i.e. ETOPO1
            else:
                # Download the files with a naming structure like:
                # /path/to/cache/gplately/file_name-as_inteded_in_url+processor_extension
                used_fname = "gplately_" + fn
                fnames = _retrieve(
                    url=url,
                    known_hash=None,
                    downloader=_HTTPDownloader(progressbar=verbose),
                    fname=used_fname,
                    path=_os_cache("gplately"),
                    processor=processor_to_use,
                )
        # If the requested files do not need processing, like standalone .nc files:
        else:
            # Are they from plate models? These typically are age or spreading rate grids
            if model_name:
                # Download the files with a naming structure like:
                # /path/to/cache/gplately/file_name-as_inteded_in_url+processor_extension
                used_fname = model_name + "_" + fn
                fnames = _retrieve(
                    url=url,
                    known_hash=None,
                    downloader=_HTTPDownloader(progressbar=verbose),
                    fname=used_fname,
                    path=_os_cache("gplately"),
                    processor=processor_to_use,
                )
            # If not from plate models and do not need processing,
            else:
                used_fname = "gplately_" + fn
                fnames = _retrieve(
                    url=url,
                    known_hash=None,
                    downloader=_HTTPDownloader(progressbar=verbose),
                    fname=used_fname,
                    path=_os_cache("gplately"),
                    processor=processor_to_use,
                )

        if not verbose:
            logger.setLevel(log_level)

        # Get the URL's E-Tag for the first time
        etag, textfilename = _get_url_etag(url)
        _save_url_etag_to_txt(etag, textfilename)
        return (fnames, etag, textfilename, used_fname)


def download_from_web(url, verbose=True, download_changes=True, model_name=None):
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
    verbose : bool
        Choose whether to print user alerts regarding file availability,
        data server/internet connection status etc.
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

    #   NOTE: We need a way to verify the existence of requested file(s) in the gplately
    #   cache to determine whether a file needs to be installed, updated, or re-accessed
    #   from the cache. Every time a file is installed for the first time,
    #   DataServer creates a directory called `full_path`. Its existence verifies the
    #   existence
    #
    #   The nature of `full_path` is dependent on the file-type:
    #
    #   .zip files will be downloaded and expanded in an inside folder:
    #
    #   /path/to/cache/gplately/fname.zip.unzip/
    #
    #   Thus, for zips, `full_path` is a directory that ends in ".zip.unzip":
    #
    #   For example: /Users/laurenilano/Library/Caches/gplately/Muller2019.zip.unzip/

    #   Other types of files that need processing, like .gz --> .decomp, aren't
    #   expanded in an internal folder. This is also the case for files that do not
    #   need processing, e.g. ".nc" files. In these cases, `full_path` is the exact
    #   directory that the cached file is saved to.
    #
    #
    #   For example: /Users/laurenilano/Library/Caches/gplately/Muller_etal_2019_Tectonics_v2.0_AgeGrid-100.nc
    #
    #   `full_path` is an empty directory for non-zips, and is the parent directory of
    #   unzipped contents in ".zip" URLs.
    #
    #   Why do we need `full_path`?
    #   We search the top-level gplately cache directory for the `full_path` directory as it is
    #   installed with the requested files. Its existence verifies the
    #   existence of the requested file(s), and thus to decide whether to install the
    #   files or re-access existing cached versions. This also helps with E-Tag versioning
    #   in instances where the download URL remains the same but its contents may have changed
    #   since the file(s) were last cached.

    full_path, unprocessed_path = path_of_cached_file(url, model_name)

    # If the file required processing (zips make a directory to unzip in, and .gz for example
    # makes a file just saved to the top-level directory), and the directory or file is not
    # yet on the cache,
    if _determine_processor(url)[1] and not (
        _os.path.isdir(str(full_path)) or _os.path.isfile(str(full_path))
    ):

        # ...and if a connection to the web server can be established,
        # download files from the URL and create a textfile for this URL's E-Tag
        if _test_internet_connection(url):
            fnames, etag, textfilename, used_fname = _first_time_download_from_web(
                url, model_name=model_name, verbose=verbose
            )  # type: ignore
            if verbose:
                print("Requested files downloaded to the GPlately cache folder!")
            return fnames

        # ... if a connection to the web server cannot be established
        else:
            raise ConnectionError(
                "A connection to {} could not be made. Please check your internet connection and/or ensure the URL is correct. No file from the given URL has been cached to {} yet - nothing has been returned.".format(
                    url, full_path.parent
                )
            )

    # If the file does not require processing, it did not open up a directory, so check isfile,
    # and if the file is not yet on the cache,
    elif not _determine_processor(url)[1] and not _os.path.isfile(str(full_path)):
        # ...and if a connection to the web server can be established,
        # download files from the URL and create a textfile for this URL's E-Tag
        if _test_internet_connection(url):
            fnames, etag, textfilename, used_fname = _first_time_download_from_web(
                url, model_name=model_name, verbose=verbose
            )  # type: ignore
            if verbose:
                print("Requested files downloaded to the GPlately cache folder!")
            return fnames

        # ... if a connection to the web server cannot be established
        else:
            raise ConnectionError(
                "A connection to {} could not be made. Please check your internet connection and/or ensure the URL is correct. No file from the given URL has been cached to {} yet - nothing has been returned.".format(
                    url, full_path.parent
                )
            )

    # If the files have been downloaded before...
    else:
        # ... and if a connection to the web server can be made...
        if _test_internet_connection(url):

            _, local_etag_txtfile = _get_url_etag(url)

            # If the newest version of the files in `url` must be cached
            # at all times, perform E-Tag comparisons:
            if download_changes:

                # Walk through the top-level cache directory to find an E-Tag textfile unique to the URL
                etag_exists = False

                cache_path = str(path_to_cache())
                if _os.path.isfile(local_etag_txtfile):
                    etag_exists = True

                # If an e-tag text file does not exist, erase the cached files
                # and download the latest version from the web server. This, in turn,
                # creates an e-tag textfile for this version.
                if not etag_exists:
                    if _os.path.isdir(full_path):
                        _shutil.rmtree(str(full_path))
                    elif _os.path.isfile(full_path):
                        _os.remove(full_path)

                    if unprocessed_path:
                        if _os.path.isdir(unprocessed_path):
                            _shutil.rmtree(str(unprocessed_path))
                        elif _os.path.isfile(unprocessed_path):
                            _os.remove(unprocessed_path)

                    fnames, etag, local_etag_txtfile, used_fname = (
                        _first_time_download_from_web(
                            url, model_name=model_name, verbose=verbose
                        )
                    )  # type: ignore
                    return fnames

                # If the e-tag textfile exists for the local files,
                else:
                    if verbose:
                        print(
                            "Checking whether the requested files need to be updated..."
                        )

                    # Determine the local file's URL e-tag from the textfile
                    with open(local_etag_txtfile) as f:
                        local_etag = str(f.readlines()[0])

                    # Get the e-tag of the web server URL at current time
                    remote_etag, remote_etag_textfile = _get_url_etag(url)

                    # If the local and remote e-tags are unequal, the web-server URL
                    # contains an updated version of the cached files.
                    if str(remote_etag) != str(local_etag):
                        if verbose:
                            print("Yes - updating requested files...")

                        # Update the e-tag textfile with this newly-identified URL e-tag
                        _save_url_etag_to_txt(remote_etag, local_etag_txtfile)

                        # Delete existing version of the files...
                        # If it didn't need processing, i.e. 'unzipping', just delete as-is
                        if _os.path.isdir(full_path):
                            _shutil.rmtree(str(full_path))
                        elif _os.path.isfile(full_path):
                            _os.remove(full_path)

                        # If it's the kind of file that needs processing, delete the
                        # unprocessed version so we can re-download it
                        if unprocessed_path:
                            if _os.path.isdir(unprocessed_path):
                                _shutil.rmtree(str(unprocessed_path))
                            elif _os.path.isfile(unprocessed_path):
                                _os.remove(unprocessed_path)

                        # Treat as if downloading the file(s) from the URL for the first time
                        fnames, etag, local_etag_txtfile, used_fname = (
                            _first_time_download_from_web(
                                url, model_name=model_name, verbose=verbose
                            )
                        )  # type: ignore

                        if verbose:
                            print(
                                "Updated requested files downloaded to the GPlately cache folder!"
                            )
                        return fnames

                    # If the e-tags are equal, the local and remote files are the same.
                    # Just return the file(s) as-is.
                    else:
                        if verbose:
                            print("Requested files are up-to-date!")

                        # If files were processed once, return the processed files.
                        if _determine_processor(url):
                            if str(full_path).endswith(_determine_processor(url)[1]):
                                return _extract_processed_files((str(full_path)))
                            else:
                                return _extract_processed_files(
                                    str(full_path) + _determine_processor(url)[1]
                                )
                        # If not, return as-is.
                        else:
                            return _extract_processed_files(
                                str(full_path) + _match_url_to_extension(url)  # type: ignore
                            )

            # If file versioning doesn't matter, just keep returning the cached files.
            else:
                fnames, etag, local_etag_txtfile = _first_time_download_from_web(  # type: ignore
                    url, model_name
                )  # type: ignore
                return fnames

        # If a connection to the web server could not be made, and the files exist in
        # the GPlately cache, just return the files as-is.
        else:
            print(
                "No connection to {} established. The requested file(s) (potentially older versions) exist in the GPlately cache ({}) and have been returned.".format(
                    url, full_path.parent
                )
            )
            # print(str(full_path)+_determine_processor(url)[1])
            return _extract_processed_files(str(full_path))
            # This created zip.unzip.unzip, so i deleted it but not sure if this will affect other files.
            # return(_extract_processed_files(str(full_path)+_determine_processor(url)[1]))


def _collect_file_extension(fnames, file_extension):
    """Searches cached directory for filenames with a specified extension(s)."""
    sorted_fnames = []
    file_extension = tuple(file_extension)
    for file in fnames:
        if file.endswith(file_extension):
            sorted_fnames.append(file)
    return sorted_fnames


def _str_in_folder(fnames, strings_to_include=None, strings_to_ignore=None):
    fnames_to_ignore = []
    fnames_to_include = []
    sorted_fnames = []
    for i, fname in enumerate(fnames):
        parent_directory = _os.path.dirname(fname)
        if strings_to_ignore is not None:
            for s in strings_to_ignore:
                if s in parent_directory:
                    fnames_to_ignore.append(fname)
            sorted_fnames = list(set(fnames) - set(fnames_to_ignore))

    if strings_to_include is not None:
        for fname in sorted_fnames:
            parent_directory = _os.path.dirname(fname)
            for s in strings_to_include:
                if s in parent_directory:
                    fnames_to_include.append(fname)
        sorted_fnames = list(set(sorted_fnames).intersection(set(fnames_to_include)))
    return sorted_fnames


def _str_in_filename(
    fnames,
    strings_to_include=None,
    strings_to_ignore=None,
    file_collection=None,
    file_collection_sensitive=False,
):
    out = []

    def filter_func(fname):
        basename = _os.path.basename(fname)
        keep = False
        if strings_to_include is None:
            keep = True
        else:
            # If a file collection was passed to the string to include, there is at least one file specific to
            # this model that must be included. Such a file should be presented in the respective
            # strings_to_include list in data.py with format:

            #    "file_collection string_to_include"

            # That is, a whitespace must be placed between the file collection and the string to include.
            # The file collection must be identical to the string allocated to the key.

            # For example, strings_to_include = ["Muller2022 1000_0_rotfile_Merdith_et_al_optimised.rot"]

            # In this example, "Muller2022" and the strings_to_include list from data.py are passed to this function
            # when sorting through rotation files.
            # The list is looped through - if the current string has "Muller2022" (case insensitive) in it,
            # we will only pass through the filename following "Muller2022", i.e. the optmised plate model.
            # All other rotation files bundled in the webDAV zip (including the published Merdith et al. 2021 rot files)
            # are excluded from the filter.

            # If no strings in the list include the passed file collection, we have one of two options, depending on whether
            # file_collection_sensitive is True or False.

            # If it is set to True, that means that we should only treat strings_to_include as True if and only if the
            #  passed file collection was found in the strings_to_include list. Otherwise, we have to treat strings_to_include
            # as if it was NoneType, and therefore place no filter for the files we accept through (that is, accept all files).

            # If it is set to False, that means that we should treat strings_to_include as True always, irrespective of
            # whether the passed file collection was found in the strings_to_include list. An example is the static polygon
            # filetype - this relies on strings_to_include being True no matter what.

            # For example, Merdith2021, Muller2019 would have file_collection_sensitive = False because these
            # models currently don't have any files that MUST be excluded for their own instance, but MUST
            # be included for other model instances.

            # Conversely, Muller2022 would have file_collection_sensitive = True because it requires all published Merdith2021
            # rot models to be ignored (in favour of the optimised model). However, we do not want to ignore Merdith2021 rot
            # models when we are using DataServer to collect Merdith2021 files.

            if file_collection is not None:

                # If the file collection is in the provided list of strings to include...
                strings_with_file_collection = [
                    s
                    for s in strings_to_include
                    if file_collection.lower() in s.lower()
                ]
                if strings_with_file_collection:

                    # Include the string, and break out.
                    for s in strings_with_file_collection:
                        if s.split(" ")[-1].lower() in basename.lower():
                            keep = True
                            break

                # If there is a file collection passed, but none of the strings to include include the file collection,
                else:
                    # If we no longer require strings_to_include, treat as if strings_to_include is False, and just pass
                    # all files through.
                    if file_collection_sensitive is True:
                        keep = True

                    # If we still need strings_to_include, treat as if strings_to_include is True, and pass only required
                    # files through.
                    else:
                        for s in strings_to_include:
                            if s.lower() in basename.lower():
                                keep = True
                                break

            # If a file collection is not passed, but strings_to_include exists, only pass through those requested.
            else:
                for s in strings_to_include:
                    if s.lower() in basename.lower():
                        keep = True
                        break

        if strings_to_ignore is not None:
            for s in strings_to_ignore:
                if s.lower() in basename.lower():
                    keep = False
                    break
        return keep

    return list(filter(filter_func, fnames))


def _check_gpml_or_shp(fnames):
    """For topology features, returns GPML by default. Searches for ESRI Shapefiles
    instead if GPML files not found."""
    sorted_fnames = []
    for file in fnames:
        if file.endswith(".gpml") or file.endswith(".gpmlz"):
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
    filepath_digits = []
    for i, file in enumerate(fnames):
        digits = []
        for element in _re.split("([0-9]+)", _remove_hash(file)):
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
            enumerate(geological_times[0]), key=lambda x: x[1]
        )
        sorted_geological_time_indices = [
            geo_time[0] for geo_time in sorted_geological_times
        ]
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
        "whittaker2015",
    ]
    if string_identifier.lower() in needs_sorting:
        studyname = _re.findall(r"[A-Za-z]+|\d+", string_identifier)[0]
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


def get_raster(raster_id_string=None, verbose=True):
    """Downloads assorted raster data that are not associated with the plate
    reconstruction models supported by GPlately's `DataServer`. Stores rasters in the
    "gplately" cache.

    Currently, gplately supports the following rasters and images:

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
    a gplately.Raster object
        A gplately.Raster object containing the raster data. The gridded data can be extracted
        into a numpy ndarray or MaskedArray by appending `.data` to the variable assigned to `get_raster()`.

        For example:

            graster = gplately.download.get_raster(raster_id_string, verbose)

            graster_data = graster.data

        where `graster_data` is a numpy ndarray. This array can be visualised using
        `matplotlib.pyplot.imshow` on a `cartopy.mpl.GeoAxis` GeoAxesSubplot
        (see example below).

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

        etopo1 = gplately.download.get_raster("ETOPO1_tif")
        etopo1_data = etopo1.data

        fig = plt.figure(figsize=(18,14), dpi=300)
        ax = fig.add_subplot(111, projection=ccrs.Mollweide(central_longitude = -150))
        etopo1.imshow(ax)

    """
    from matplotlib import image

    if raster_id_string is None:
        raise ValueError("Please specify which raster to download.")
    # filetype = "."+"_".split(raster_id_string)[-1]

    archive_formats = tuple([".gz", ".xz", ".bz2"])
    grid_extensions = tuple([".grd", ".nc"])

    # Set to true if we find the given collection in database
    found_collection = False
    raster_filenames = []
    database = _rasters()

    for collection, zip_url in database.items():
        # Isolate the raster name and the file type
        # raster_name = collection.split("_")[0]
        # raster_type = "."+collection.split("_")[-1]
        if raster_id_string.lower() == collection.lower():
            raster_filenames = download_from_web(zip_url[0], verbose)
            found_collection = True
            break

    if found_collection is False:
        raise ValueError("{} not in collection database.".format(raster_id_string))
    else:
        # If the downloaded raster is a grid, process it with the gplately.Raster object
        if any(
            grid_extension in raster_filenames for grid_extension in grid_extensions
        ):
            raster = Raster(data=raster_filenames)

        # Otherwise, the raster is an image; use imread to process
        else:
            raster_matrix = image.imread(raster_filenames)  # type: ignore
            raster = Raster(data=raster_matrix)

        if raster_id_string.lower() == "etopo1_tif":
            raster.lats = raster.lats[::-1]
        if raster_id_string.lower() == "etopo1_grd":
            raster._data = raster._data.astype(float)  # type: ignore

    return raster


def get_feature_data(feature_data_id_string=None, verbose=True):
    """Downloads assorted geological feature data from web servers (i.e.
    [GPlates 2.3 sample data](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/))
    into the "gplately" cache.

    Currently, gplately supports the following feature data:

    --------------

    | **Feature data string identifier** | **Description**                                                                                                                                                                              |
    |------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | Johansson2018                      | Large igneous provinces from  Johansson et al. (2018)                                                                                                                                        |
    | Whittaker2015                      | Large igneous province products  interpreted as plume products  from Whittaker et al. (2015).                                                                                                |
    | SeafloorFabric                     | Seafloor tectonic fabric  (fracture zones, discordant zones,  V-shaped structures, unclassified  V-anomalies, propagating ridge  lineations and extinct ridges)  from Matthews et al. (2011) |
    | Hotspots                           | Present day surface hotspot/plume  locations from Whittaker et al. (2013)                                                                                                                    |

    ---------------

    Detailed descriptions can be found below:

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
        raise ValueError("Please specify which feature data to fetch.")

    database = _feature_data()

    found_collection = False
    feature_data_filenames = []
    for collection, zip_url in database.items():
        if feature_data_id_string.lower() == collection.lower():
            found_collection = True
            feature_data_filenames = _collection_sorter(
                _collect_file_extension(
                    download_from_web(zip_url[0], verbose), [".gpml", ".gpmlz"]
                ),
                collection,
            )

            break

    if found_collection is False:
        raise ValueError(
            "{} are not in GPlately's DataServer.".format(feature_data_id_string)
        )

    feat_data = pygplates.FeatureCollection()
    if len(feature_data_filenames) == 1:
        feat_data.add(pygplates.FeatureCollection(feature_data_filenames[0]))
        return feat_data
    else:
        feat_data = []
        for file in feature_data_filenames:
            feat_data.append(pygplates.FeatureCollection(file))
        return feat_data


class DataServer(object):
    """
    Download the plate reconstruction models from the `EarthByte server <https://repo.gplates.org/webdav/pmm/>`__.

    The :class:`DataServer` object downloads the model files to the ``GPlately cache folder``.
    If the same model is requested again, a new :class:`DataServer` instance will retrieve the files from the cache --
    provided they haven't been moved or deleted.

    .. seealso::

        - `This table <https://gplates.github.io/gplately/sphinx-latest/html/use_cases.html#id1>`__ provides a list of available plate reconstruction models.
        - Visit this `EarthByte web page <https://www.earthbyte.org/category/resources/data-models/global-regional-plate-motion-models/>`__ for more information about these plate models.
        - Call :meth:`gplately.auxiliary.get_data_server_cache_path` to see the path to the ``GPlately cache folder``.
    """

    def __init__(self, file_collection, data_dir=None, verbose=True):
        """Constructor. Create a :class:`DataServer` object.

        Example
        -------
        .. code-block:: python

            # create a DataServer object for the Cao2024 model (https://zenodo.org/records/11536686)
            data_server = gplately.download.DataServer("Cao2024")

        Parameters
        ----------
        file_collection: str
            The model name of interest.

        verbose: bool, default=True
            Toggle print messages regarding server/internet connection status, file availability, etc.
        """

        if not data_dir:
            _data_dir = path_to_cache()
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

    @property
    def cache_path(self):
        """The location of DataServer cache on your computer."""
        return path_to_cache()

    @property
    def rotation_model(self):
        """A pygplates.RotationModel object for the plate reconstruction model."""
        if self._rotation_model is None and self.pmm:
            self._rotation_model = pygplates.RotationModel(
                self.pmm.get_rotation_model()
            )
            self._rotation_model.reconstruction_identifier = self.file_collection
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
            Müller et al. (2019) plate reconstruction model and create a :class:`gplately.PlateReconstruction` object.

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

            The example code below will attempt to download ``coastlines``, ``continents`` and ``COBs`` from the Müller
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

            :attr:`DataServer.cache_path`

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
        To download  Müller et al. (2019) seafloor age grids for 0Ma, 1Ma and 100 Ma:

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
            except:
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

            :attr:`DataServer.cache_path`

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
            except:
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
                data_dir=str(path_to_cache())
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

    def get_feature_data(self, feature_data_id_string=None):
        """Downloads assorted geological feature data from web servers (i.e.
        `GPlates 2.3 sample data <https://www.earthbyte.org/gplates-2-3-software-and-data-sets/>`__) into the "gplately" cache.

        Currently, ``DataServer`` supports the following feature data:

        * Large igneous provinces from Johansson et al. (2018)

            - Information

                * Formats: .gpmlz
                * `feature_data_id_string` = `Johansson2018`

            - Citations

                Johansson, L., Zahirovic, S., and Müller, R. D., In Prep, The
                interplay between the eruption and weathering of Large Igneous Provinces and
                the deep-time carbon cycle: Geophysical Research Letters.


        - Large igneous province products interpreted as plume products from Whittaker et al. (2015).

            - Information

                * Formats: .gpmlz, .shp
                * `feature_data_id_string` = `Whittaker2015`

            - Citations

                Whittaker, J. M., Afonso, J. C., Masterton, S., Müller, R. D.,
                Wessel, P., Williams, S. E., & Seton, M. (2015). Long-term interaction between
                mid-ocean ridges and mantle plumes. Nature Geoscience, 8(6), 479-483. doi:10.1038/ngeo2437.


        - Seafloor tectonic fabric (fracture zones, discordant zones, V-shaped structures, unclassified V-anomalies, propagating ridge lineations and extinct ridges) from Matthews et al. (2011)

            - Information

                * Formats: .gpml
                * `feature_data_id_string` = `SeafloorFabric`

            - Citations

                Matthews, K.J., Müller, R.D., Wessel, P. and Whittaker, J.M., 2011. The
                tectonic fabric of the ocean basins. Journal of Geophysical Research, 116(B12):
                B12109, DOI: 10.1029/2011JB008413.


        - Present day surface hotspot/plume locations from Whittaker et al. (2013)

            - Information

                * Formats: .gpmlz
                * `feature_data_id_string` = `Hotspots`

            - Citation

                Whittaker, J., Afonso, J., Masterton, S., Müller, R., Wessel, P.,
                Williams, S., and Seton, M., 2015, Long-term interaction between mid-ocean ridges and
                mantle plumes: Nature Geoscience, v. 8, no. 6, p. 479-483, doi:10.1038/ngeo2437.


        Parameters
        ----------
        feature_data_id_string : str, default=None
            A string to identify which feature data to download to the cache (see list of supported feature data above).

        Returns
        -------
        feature_data_filenames : instance of <pygplates.FeatureCollection>, or list of instance <pygplates.FeatureCollection>
            If a single set of feature data is downloaded, a single pyGPlates ``FeatureCollection``
            object is returned. Otherwise, a list containing multiple pyGPlates ``FeatureCollection``
            objects is returned (like for ``SeafloorFabric``). In the latter case, feature reconstruction
            and plotting may have to be done iteratively.

        Raises
        ------
        ValueError
            If a ``feature_data_id_string`` is not provided.

        Examples
        --------
        For examples of plotting data downloaded with ``get_feature_data``, see GPlately's sample notebook 05 - Working With Feature Geometries
        `here <https://github.com/GPlates/gplately/blob/master/Notebooks/05-WorkingWithFeatureGeometries.ipynb>`__.
        """
        if feature_data_id_string is None:
            raise ValueError("Please specify which feature data to fetch.")

        database = _feature_data()

        found_collection = False
        feature_data_filenames = []
        for collection, zip_url in database.items():
            if feature_data_id_string.lower() == collection.lower():
                found_collection = True
                feature_data_filenames = _collection_sorter(
                    _collect_file_extension(
                        download_from_web(zip_url[0], self.verbose), [".gpml", ".gpmlz"]
                    ),
                    collection,
                )

                break

        if found_collection is False:
            raise ValueError(
                "{} are not in GPlately's DataServer.".format(feature_data_id_string)
            )

        feat_data = pygplates.FeatureCollection()
        if len(feature_data_filenames) == 1:
            feat_data.add(pygplates.FeatureCollection(feature_data_filenames[0]))
            return feat_data
        else:
            feat_data = []
            for file in feature_data_filenames:
                feat_data.append(pygplates.FeatureCollection(file))
            return feat_data
