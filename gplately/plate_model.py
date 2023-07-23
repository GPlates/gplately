import glob
import json
import os
import shutil, time
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
            self.data_dir = platformdirs.user_cache_dir("gplately")
        else:
            self.data_dir = data_dir
        models_file = f"{self.data_dir}/models.json"
        models = None

        if force_fresh:
            os.remove(models_file)

        # check the local models cfg first
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

    def get_layer(self, layer_name):
        """get a layer by name"""
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

    def download_layer_files(self, layer_name, force=False):
        """given the layer name, download the layer files

        :param layer_name: such as "Rotations","Coastlines", "StaticPolygons", "ContinentalPolygons", "Topologies", etc
        :param force: delete the local files and download again

        :returns: the folder path which contains the layer files

        """

        download_flag = False
        download_with_etag = False
        meta_etag = None
        if layer_name in self.model:
            layer_file_url = self.model[layer_name]
        elif "Layers" in self.model and layer_name in self.model["Layers"]:
            layer_file_url = self.model["Layers"][layer_name]
        else:
            raise Exception(f"Fatal: No {layer_name} files in configuration file!")

        now = datetime.now()

        # make sure the model folder is a folder indeed; otherwise rename it.
        model_folder = f"{self.data_dir}/{self.model_name}"
        if os.path.isfile(model_folder):
            new_name = model_folder + now.strftime("-%Y-%m-%d-%H-%M-%S")
            os.rename(model_folder, new_name)
            print(
                f"{model_folder} is a file. This should not happen. Rename the file name to {new_name}"
            )

        # first check if the "Rotations" folder exists
        layer_folder = f"{model_folder}/{layer_name}"

        if force:
            shutil.rmtree(layer_folder)

        if os.path.isdir(layer_folder):
            metadata_file = f"{layer_folder}/{self.meta_filename}"
            if os.path.isfile(metadata_file):
                with open(metadata_file, "r") as f:
                    meta = json.load(f)
                    if "url" in meta:
                        meta_url = meta["url"]
                        if meta_url != layer_file_url:
                            # if the data url has changed, re-download
                            download_flag = True
                    if "expiry" in meta:
                        meta_expiry = meta["expiry"]
                        expiry_date = datetime.strptime(meta_expiry, self.expiry_format)
                        now = datetime.now()
                        if now > expiry_date:
                            download_with_etag = True
                        # now.strftime("%m/%d/%Y, %H:%M:%S")
                    if "etag" in meta:
                        meta_etag = meta["etag"]
            else:
                # if the metadata.json does not exist
                download_flag = True

        else:
            # if the "Rotations" folder does not exist
            download_flag = True

        if not download_flag and not download_with_etag:
            print("The local files are still good. Will not download again.")
        else:
            new_etag = None
            if download_flag:
                new_etag = network_requests.fetch_file(
                    layer_file_url,
                    model_folder,
                    auto_unzip=True,
                )
            elif download_with_etag and meta_etag:
                new_etag = network_requests.fetch_file(
                    layer_file_url,
                    model_folder,
                    etag=meta_etag,
                    auto_unzip=True,
                )

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
