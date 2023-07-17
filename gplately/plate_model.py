import platformdirs
import pygplates
import requests
from pathlib import Path
import os, json, glob

from datetime import datetime

from . import network_requests


class PlateModel:
    """Class to manage a plate model"""

    def __init__(self, model_name, data_dir=None):
        self.model_name = model_name
        self.model = None
        if not data_dir:
            self.data_dir = platformdirs.user_cache_dir("gplately")
        else:
            self.data_dir = data_dir
        models_file = f"{self.data_dir}/models.json"
        models = None
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
                print(models[model_name])
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
        if "Rotations" not in self.model:
            raise Exception("Fatal: No rotation files in configuration file!")
        download_flag = False
        download_with_etag = False
        meta_etag = None
        rotation_file_url = self.model["Rotations"]
        now = datetime.now()

        # make sure the model folder is a folder indeed
        model_folder = f"{self.data_dir}/{self.model_name}"
        if os.path.isfile(model_folder):
            new_name = model_folder + now.strftime("-%Y-%m-%d-%H-%M-%S")
            os.rename(model_folder, new_name)
            print(
                f"{model_folder} is a file. This should not happen. Rename the file name to {new_name}"
            )

        # first check if the "Rotations" folder exists
        rotations_folder = f"{model_folder}/Rotations"
        if os.path.isdir(rotations_folder):
            metadata_file = f"{rotations_folder}/metadata.json"
            if os.path.isfile(metadata_file):
                with open(metadata_file, "r") as f:
                    meta = json.load(f)
                    if "url" in meta:
                        meta_url = meta["url"]
                        if meta_url != rotation_file_url:
                            # if the data url has changed, re-download
                            download_flag = True
                    if "expiry" in meta:
                        meta_expiry = meta["expiry"]
                        expiry_date = datetime.strptime(
                            meta_expiry, "%Y/%m/%d, %H:%M:%S"
                        )
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

        if download_flag:
            network_requests.fetch_file(
                rotation_file_url,
                model_folder,
                auto_unzip=True,
            )
        elif download_with_etag and meta_etag:
            network_requests.fetch_file(
                rotation_file_url,
                model_folder,
                etag=meta_etag,
                auto_unzip=True,
            )

        rotation_files = glob.glob(f"{rotations_folder}/*.rot")
        return pygplates.RotationModel(rotation_files)
