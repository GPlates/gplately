""" A module that utilises Pooch functionalities for accessing web-based data on-the-fly. Pooch is a data manager that
downloads files only when necessary and stores them in your machine's local cache.

Methods in this module require two parameters: a URL to a downloadable zip file containing the needed data; and the
file extension format needed. They return the paths (as a list of strings) to files in the zip folder with the specified
extension. These paths can be extracted and used in regular GPlately functions.

To see the full paths of the filenames extracted by functions in this module, print() the lists returned by the fetch_X
modules. To see just the extracted filenames, pass "show_obtained_files" as TRUE when calling the methods in this module.

Methods
-------
fetch_filenames
    This can be used to find any file contained in a web-based zip folder. Useful for finding ".rot" files. Note: if using
    this to locate ".shp" files, it will not automatically differentiate between different geometries and static polygons.
    Use the functions detailed below to find specific ".shp" files.

fetch_coastlines
    Finds and returns files in a web-based zip folder with the substring "coastline" and the user-given extension.
fetch_continents
    Finds and returns files in a web-based zip folder with the substring "continent" and the user-given extension.
fetch_COBs
    Finds and returns files in a web-based zip folder with the substring "COB" and the user-given extension.
fetch_static_polygons
    Finds and returns files in a web-based zip folder with the substring "static" and the user-given extension.

Sources
-------
Methods in this module take inspiration from and utilise:
Pooch : https://github.com/fatiando/pooch, DOI: https://doi.org/10.21105/joss.01943)
GPlatesReconstruction Model : S. Williams (2021) https://github.com/siwill22/GPlatesReconstructionModel/blob/dc29746d073050987fc91258b7786f1074df8efa/gprm/datasets/
Reconstructions.py#L174
"""
import os

from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip


class DataServer(object):

    def fetch_filenames(self, zip_url=None, file_ext=None, show_obtained_files=False):
        """ Obtains files from a zip web URL with the specified file extension format.

        Uses Pooch to download a zip file on-the-fly. Stores all zip file contents in the machine's local cache.
        Returns a list of paths to all files in the zip folder with the specified file extension format.

        Parameters
        ----------
        zip_url : str, default=None
            A single string pointing to a zip file on the web that contains topology features to extract.
        file_ext : str, default=None
            The desired file extension format to extract
        show_obtained_files : bool, default=False
            If set to true, prints the filename(s) of the extracted file(s) without the cache directory.

        Returns
        -------
        feature_filenames : list
            A list containing paths to all files in the zip folder with the specified file extension format.

        Raises
        ------
        ValueError
            If a URL and/or file extension are not supplied.
            If there are no files in the zip file with the specified file extension.
        """

        if zip_url is None:
            raise ValueError("Please supply a URL to a zipfile containing feature topologies.")
        if file_ext is None:
            raise ValueError("Please supply a file extension to fetch.")

        # Use pooch to download topologies into local cache
        fnames = _retrieve(
            url=zip_url,
            known_hash=None,
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gplately'),
            processor=_Unzip(),
        )
        dirname = os.path.split(fnames[0])[0]
        feature_filenames = []
        for subdir, dirs, files in os.walk(dirname):
            for file in files:
                if file.endswith(file_ext):
                    feature_filenames.append(subdir+"/"+file)
                    if show_obtained_files:
                        print(file)

        if len(feature_filenames) == 0:
            raise ValueError("No files with extension %s found." % (file_ext))
        return feature_filenames

    def fetch_coastlines(self, zip_url=None, file_ext=None, show_obtained_files=False):
        """Obtains coastline geometries from a zip web URL with the specified file extension format.

        Uses Pooch to download a zip file on-the-fly. Stores all zip file contents in the machine's local cache.
        Returns a list of paths to all coastline geometry files in the zip folder with the specified file extension format.

        Parameters
        ----------
        zip_url : str, default=None
            A single string pointing to a zip file on the web that contains coastline features to extract.
        file_ext : str, default=None
            The desired file extension format to extract
        show_obtained_files : bool, default=False
            If set to true, prints the filename(s) of the extracted file(s) without the cache directory.

        Returns
        -------
        feature_filenames : list
            A list containing paths to all coastline files in the zip folder with the specified file extension format.

        Raises
        ------
        ValueError
            If a URL and/or file extension are not supplied.
            If there are no coastlines in the zip file with the specified file extension.
        """

        if zip_url is None:
            raise ValueError("Please supply a URL to a zipfile containing feature topologies.")
        if file_ext is None:
            raise ValueError("Please supply a file extension to fetch.")

        # Use pooch to download topologies into local cache
        fnames = _retrieve(
            url=zip_url,
            known_hash=None,
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gplately'),
            processor=_Unzip(),
        )
        dirname = os.path.split(fnames[0])[0]
        feature_filenames = []
        for subdir, dirs, files in os.walk(dirname):
            for file in files:
                if file.endswith(file_ext):
                    if file.lower().find("coastline") != -1:
                        feature_filenames.append(subdir+"/"+file)
                        if show_obtained_files:
                            print(file)

        if len(feature_filenames) == 0:
            raise ValueError("No files with extension %s found." % (file_ext))

        return feature_filenames

    def fetch_continents(self, zip_url=None, file_ext=None, show_obtained_files=False):
        """Obtains continental geometries from a zip web URL with the specified file extension format.

        Uses Pooch to download a zip file on-the-fly. Stores all zip file contents in the machine's local cache.
        Returns a list of paths to all continental geometries in the zip folder with the specified file extension format.

        Parameters
        ----------
        zip_url : str, default=None
            A single string pointing to a zip file on the web that contains continental features to extract.
        file_ext : str, default=None
            The desired file extension format to extract
        show_obtained_files : bool, default=False
            If set to true, prints the filename(s) of the extracted file(s) without the cache directory.

        Returns
        -------
        feature_filenames : list
            A list containing paths to all continental geometries in the zip folder with the specified file extension format.

        Raises
        ------
        ValueError
            If a URL and/or file extension are not supplied.
            If there are no files in the zip file with the specified file extension.
        """

        if zip_url is None:
            raise ValueError("Please supply a URL to a zipfile containing feature topologies.")
        if file_ext is None:
            raise ValueError("Please supply a file extension to fetch.")

        # Use pooch to download topologies into local cache
        fnames = _retrieve(
            url=zip_url,
            known_hash=None,
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gplately'),
            processor=_Unzip(),
        )
        dirname = os.path.split(fnames[0])[0]
        feature_filenames = []
        for subdir, dirs, files in os.walk(dirname):
            for file in files:
                if file.endswith(file_ext):
                    if file.lower().find("continent") != -1:
                        feature_filenames.append(subdir+"/"+file)
                        if show_obtained_files:
                            print(file)

        if len(feature_filenames) == 0:
            raise ValueError("No files with extension %s found." % (file_ext))
        return feature_filenames

    def fetch_COBs(self, zip_url=None, file_ext=None, show_obtained_files=False):
        """Obtains continent-ocean boundary (COB) geometries from a zip web URL with the specified file extension format.

        Uses Pooch to download a zip file on-the-fly. Stores all zip file contents in the machine's local cache.
        Returns a list of paths to all COB geometries in the zip folder with the specified file extension format.

        Parameters
        ----------
        zip_url : str, default=None
            A single string pointing to a zip file on the web that contains COB features to extract.
        file_ext : str, default=None
            The desired file extension format to extract
        show_obtained_files : bool, default=False
            If set to true, prints the filename(s) of the extracted file(s) without the cache directory.

        Returns
        -------
        feature_filenames : list
            A list containing paths to all COB feature files in the zip folder with the specified file extension format.

        Raises
        ------
        ValueError
            If a URL and/or file extension are not supplied.
            If there are no files in the zip file with the specified file extension.
        """

        if zip_url is None:
            raise ValueError("Please supply a URL to a zipfile containing feature topologies.")
        if file_ext is None:
            raise ValueError("Please supply a file extension to fetch.")

        # Use pooch to download topologies into local cache
        fnames = _retrieve(
            url=zip_url,
            known_hash=None,
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gplately'),
            processor=_Unzip(),
        )
        dirname = os.path.split(fnames[0])[0]
        feature_filenames = []
        for subdir, dirs, files in os.walk(dirname):
            for file in files:
                if file.endswith(file_ext):
                    if file.lower().find("cob") != -1:
                        feature_filenames.append(subdir+"/"+file)
                        if show_obtained_files:
                            print(file)

        if len(feature_filenames) == 0:
            raise ValueError("No files with extension %s found." % (file_ext))
        return feature_filenames

    def fetch_static_polygons(self, zip_url=None, file_ext=None, show_obtained_files=False):
        """Obtains static polygon files from a zip web URL with the specified file extension format.

        Uses Pooch to download a zip file on-the-fly. Stores all zip file contents in the machine's local cache.
        Returns a list of paths to all static polygon files in the zip folder with the specified file extension format.

        Parameters
        ----------
        zip_url : str, default=None
            A single string pointing to a zip file on the web that contains static polygons to extract.
        file_ext : str, default=None
            The desired file extension format to extract
        show_obtained_files : bool, default=False
            If set to true, prints the filename(s) of the extracted file(s) without the cache directory.

        Returns
        -------
        feature_filenames : list
            A list containing paths to all static polygon files in the zip folder with the specified file extension format.

        Raises
        ------
        ValueError
            If a URL and/or file extension are not supplied.
            If there are no files in the zip file with the specified file extension.
        """

        if zip_url is None:
            raise ValueError("Please supply a URL to a zipfile containing feature topologies.")
        if file_ext is None:
            raise ValueError("Please supply a file extension to fetch.")

        # Use pooch to download topologies into local cache
        fnames = _retrieve(
            url=zip_url,
            known_hash=None,
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gplately'),
            processor=_Unzip(),
        )
        dirname = os.path.split(fnames[0])[0]
        feature_filenames = []
        for subdir, dirs, files in os.walk(dirname):
            for file in files:
                if file.endswith(file_ext):
                    if file.lower().find("static") != -1:
                        feature_filenames.append(subdir+"/"+file)
                        if show_obtained_files:
                            print(file)

        if len(feature_filenames) == 0:
            raise ValueError("No files with extension %s found." % (file_ext))
        return feature_filenames
