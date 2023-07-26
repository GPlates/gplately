import asyncio
import concurrent.futures
import functools
import glob
import json
import os
import shutil
import time
from datetime import datetime, timedelta
from pathlib import Path

import platformdirs
import requests

from . import network_requests, pygplates


class PlateModel:
    """Class to manage a plate model"""

    def __init__(self, model_name, data_dir=None, force_fresh=False):
        """The valid model names can be found at https://www.earthbyte.org/webdav/ftp/gplately/models.json"""
        self.model_name = model_name
        self.meta_filename = "metadata.json"
        self.expiry_format = "%Y/%m/%d, %H:%M:%S"
        self.model = None
        if not data_dir:
            self.data_dir = platformdirs.user_data_dir("gplately")
        else:
            self.data_dir = data_dir
        models_file = f"{self.data_dir}/models.json"
        models = None

        # async and concurrent things
        self.executor = concurrent.futures.ThreadPoolExecutor(max_workers=15)
        self.loop = asyncio.new_event_loop()
        self.run = functools.partial(self.loop.run_in_executor, self.executor)
        asyncio.set_event_loop(self.loop)

        # force refresh models.json
        if os.path.isfile(models_file) and force_fresh:
            os.remove(models_file)

        # check the local models cfg first. if not too old, use it
        if os.path.isfile(models_file):
            if (time.time() - os.path.getmtime(models_file)) < 6 * 60 * 60:  # 6 hours
                with open(models_file) as f:
                    models = json.load(f)

        # and then try the network
        if not models:
            try:
                r = requests.get(
                    "https://www.earthbyte.org/webdav/ftp/gplately/models.json"
                )
                Path(self.data_dir).mkdir(parents=True, exist_ok=True)
                models = r.json()
                with open(models_file, "w+") as of:
                    of.write(r.text)
            except requests.exceptions.ConnectionError:
                print("No network connection!")
                if os.path.isfile(models_file):
                    with open(models_file) as f:
                        models = json.load(f)
        if models:
            if model_name in models:
                self.model = models[model_name]
                # print(models[model_name])
            else:
                raise Exception(f"Fatal: invalid model name: {model_name}")
        else:
            raise Exception("Fatal: failed to load models.")

    def __del__(self):
        self.loop.close()

    def get_data_dir(self):
        return self.data_dir

    def get_avail_layers(self):
        """get all available layers in this plate model"""
        if not self.model:
            raise Exception("Fatal: No model configuration found!")
        return list(self.model["Layers"].keys())

    def get_rotation_model(self):
        """return a pygplates.RotationModel object"""
        rotation_folder = self.download_layer_files("Rotations")
        rotation_files = glob.glob(f"{rotation_folder}/*.rot")
        rotation_files.extend(glob.glob(f"{rotation_folder}/*.grot"))
        # print(rotation_files)
        return pygplates.RotationModel(rotation_files)

    def get_coastlines(self):
        """return coastlines feature collection"""
        return self.get_layer("Coastlines")

    def get_static_polygons(self):
        """return StaticPolygons feature collection"""
        return self.get_layer("StaticPolygons")

    def get_continental_polygons(self):
        """return ContinentalPolygons feature collection"""
        return self.get_layer("ContinentalPolygons")

    def get_topologies(self):
        """return Topologies feature collection"""
        return self.get_layer("Topologies")

    def get_COBs(self):
        """return COBs feature collection"""
        return self.get_layer("COBs")

    def get_layer(self, layer_name):
        """get a layer by name

        :param layer_name: layer name
        :returns: pygplates.FeatureCollection
        """
        file_extensions = [
            "gpml",
            "gpmlz",
            "gpml.gz",
            "dat",
            "pla",
            "shp",
            "geojson",
            "json",
            ".gpkg",
            "gmt",
            "vgp",
        ]
        layer_folder = self.download_layer_files(layer_name)
        files = []
        for ext in file_extensions:
            files.extend(glob.glob(f"{layer_folder}/*.{ext}"))

        fc = pygplates.FeatureCollection()
        for f in files:
            if os.path.basename(f) != self.meta_filename:
                fc.add(pygplates.FeatureCollection(f))
        return fc

    def download_layer_files(self, layer_name, dst_path=None, force=False):
        """given the layer name, download the layer files.
        The layer files are in a .zip file. download and unzip it.

        :param layer_name: such as "Rotations","Coastlines", "StaticPolygons", "ContinentalPolygons", "Topologies", etc
        :param force: delete the local files and download again

        :returns: the folder path which contains the layer files

        """
        print(f"downloading {layer_name}")
        download_flag = False
        meta_etag = None

        # find layer file url. two parts. one is the rotation, the other is all other geometry layers
        if layer_name in self.model:
            layer_file_url = self.model[layer_name]
        elif "Layers" in self.model and layer_name in self.model["Layers"]:
            layer_file_url = self.model["Layers"][layer_name]
        else:
            raise Exception(f"Fatal: No {layer_name} files in configuration file!")

        now = datetime.now()

        # make sure the model folder is a folder indeed; otherwise rename it.
        if dst_path:
            model_folder = f"{dst_path}/{self.model_name}"
        else:
            model_folder = f"{self.data_dir}/{self.model_name}"
        if os.path.isfile(model_folder):
            new_name = model_folder + now.strftime("-%Y-%m-%d-%H-%M-%S")
            os.rename(model_folder, new_name)
            print(
                f"{model_folder} is a file. This should not happen. Rename the file name to {new_name}"
            )

        layer_folder = f"{model_folder}/{layer_name}"

        if os.path.isdir(layer_folder) and force:
            shutil.rmtree(layer_folder)

        # first check if the layer folder exists
        if os.path.isdir(layer_folder):
            metadata_file = f"{layer_folder}/{self.meta_filename}"
            download_flag, meta_etag = self._check_redownload_need(
                metadata_file, layer_file_url
            )
        else:
            # if the layer folder does not exist
            download_flag = True

        if not download_flag:
            print("The local files are still good. Will not download again.")
        else:
            new_etag = network_requests.fetch_file(
                layer_file_url,
                model_folder,
                etag=meta_etag,
                auto_unzip=True,
            )

            if new_etag != meta_etag:
                # save metadata
                metadata = {
                    "layer_name": layer_name,
                    "url": layer_file_url,
                    "expiry": (now + timedelta(hours=12)).strftime(self.expiry_format),
                    "etag": new_etag,
                }
                with open(f"{layer_folder}/{self.meta_filename}", "w+") as f:
                    json.dump(metadata, f)

        return layer_folder

    def download_all_layers(self, dst_path=None, force=False):
        """download all layers. Call download_layer_files() on every layer"""

        async def f():
            tasks = []
            if "Rotations" in self.model:
                tasks.append(
                    self.run(self.download_layer_files, "Rotations", dst_path, force)
                )
            if "Layers" in self.model:
                for layer in self.model["Layers"]:
                    tasks.append(
                        self.run(self.download_layer_files, layer, dst_path, force)
                    )

            # print(tasks)
            await asyncio.wait(tasks)

        self.loop.run_until_complete(f())

    def _check_redownload_need(self, metadata_file, url):
        """check the metadata file and decide if redownload is necessary

        :param metadata_file: metadata file path
        :param url: url for the target file

        :returns download_flag, etag: a flag indicates if redownload is neccesarry and old etag if needed.
        """
        download_flag = False
        meta_etag = None
        if os.path.isfile(metadata_file):
            with open(metadata_file, "r") as f:
                meta = json.load(f)
                if "url" in meta:
                    meta_url = meta["url"]
                    if meta_url != url:
                        # if the data url has changed, re-download
                        download_flag = True
                else:
                    download_flag = True

                # if the url is the same, now check the expiry date
                if not download_flag:
                    if "expiry" in meta:
                        try:
                            meta_expiry = meta["expiry"]
                            expiry_date = datetime.strptime(
                                meta_expiry, self.expiry_format
                            )
                            now = datetime.now()
                            if now > expiry_date:
                                download_flag = True  # expired
                        except ValueError:
                            download_flag = True  # invalid expiry date
                    else:
                        download_flag = True  # no expiry date in metafile

                    if download_flag and "etag" in meta:
                        meta_etag = meta["etag"]
        else:
            download_flag = True  # if metadata_file does not exist

        return download_flag, meta_etag

    def get_avail_time_dependent_raster_names(self):
        """return the names of all time dependent rasters which have been configurated in this model."""
        if not "TimeDepRasters" in self.model:
            return []
        else:
            return [name for name in self.model["TimeDepRasters"]]

    def download_time_dependent_rasters(self, raster_name, dst_path=None, times=None):
        """download time dependent rasters, such agegrids

        :param raster_name: raster name, such as AgeGrids. see the models.json
        :param dst_path: where to save the files
        :param times: if not given, download from begin to end with 1My interval
        """
        if (
            "TimeDepRasters" in self.model
            and raster_name in self.model["TimeDepRasters"]
        ):

            async def f():
                nonlocal times
                nonlocal dst_path
                tasks = []
                if not dst_path:
                    dst_path = f"{self.get_data_dir()}/{self.model_name}/{raster_name}"
                if not times:
                    times = range(self.model["SmallTime"], self.model["BigTime"])
                for time in times:
                    tasks.append(
                        self.run(
                            self.download_raster,
                            self.model["TimeDepRasters"][raster_name].format(time),
                            dst_path,
                        )
                    )

                # print(tasks)
                await asyncio.wait(tasks)

            self.loop.run_until_complete(f())

        else:
            raise Exception(
                f"Unable to find {raster_name} configuration in this model {self.model_name}."
            )

    def download_raster(self, url, dst_path):
        """download a single raster file from "url" and save the file in "dst_path"
        a metadata file will also be created for the raster file in folder f"{dst_path}/metadata"

        """
        print(f"downloading {url}")
        filename = url.split("/")[-1]
        metadata_folder = f"{dst_path}/metadata"
        metadata_file = f"{metadata_folder}/{filename}.json"
        download_flag, etag = self._check_redownload_need(metadata_file, url)
        # only redownload when necessary
        if download_flag:
            new_etag = network_requests.fetch_file(
                url,
                dst_path,
                etag=etag,
                auto_unzip=True,
            )
            if etag != new_etag:
                # save metadata file
                metadata = {
                    "url": url,
                    "expiry": (datetime.now() + timedelta(hours=12)).strftime(
                        self.expiry_format
                    ),
                    "etag": new_etag,
                }
                Path(metadata_folder).mkdir(parents=True, exist_ok=True)
                with open(metadata_file, "w+") as f:
                    json.dump(metadata, f)
