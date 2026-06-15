import logging
import os
import sys
from pathlib import Path
from urllib.request import urlretrieve

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false

import xarray as xr
import matplotlib.pyplot as plt
from plate_model_manager import PresentDayRasterManager

OUTPUT_DIR = "output"

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

MODEL_REPO_DIR = "plate-model-repo"

logger = logging.getLogger("gplately-unittest-logger")
fhandler = logging.FileHandler(filename="gplately-unittest.log", mode="a")
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
logger.setLevel(logging.INFO)


def get_logger():
    return logger


def save_fig(filename):
    output_file = f"{OUTPUT_DIR}/{Path(filename).stem}.png"
    plt.gcf().savefig(output_file, dpi=120, bbox_inches="tight")  # transparent=True)
    print(f"Done! The {output_file} has been saved.")
    plt.close(plt.gcf())


def get_test_local_code_flag():
    if "GPLATELY_TEST_LOCAL_FLAG" in os.environ and (
        os.environ["GPLATELY_TEST_LOCAL_FLAG"].lower() == "false"
        or os.environ["GPLATELY_TEST_LOCAL_FLAG"].lower() == "0"
    ):
        return False
    return True


if get_test_local_code_flag():
    sys.path.insert(0, f"{os.path.dirname(os.path.realpath(__file__))}/../..")


def _get_topo_raster():
    data_dir = Path(__file__).resolve().parent / "unittest-data"
    topo_file = PresentDayRasterManager(data_dir=data_dir).get_raster("topography")

    try:
        raster_xr = xr.open_dataarray(topo_file)
    except ValueError:
        # Fallback for NetCDF files with multiple variables.
        dataset = xr.open_dataset(topo_file)
        first_var = next(iter(dataset.data_vars))
        raster_xr = dataset[first_var]
    return raster_xr
