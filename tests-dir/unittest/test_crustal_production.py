#!/usr/bin/env python3


from common import *

import gplately

print(gplately.__file__)
from plate_model_manager.utils import download

from gplately.utils import crustal_production

data_dir = "./crustal_production_agegrid"
if __name__ == "__main__":
    downloader = download.FileDownloader(
        "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc",
        f"{data_dir}/.metadata.json",
        f"{data_dir}",
    )
    # only re-download when necessary
    if downloader.check_if_file_need_update():
        downloader.download_file_and_update_metadata()
    else:
        print(f"The local age grid is still good. No need to download again!")

    for i in range(3, 10):
        crustal_production_rate = crustal_production.compute_crustal_production(
            f"{data_dir}/Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc", i
        )
        print(f"{i}: {crustal_production_rate}")
