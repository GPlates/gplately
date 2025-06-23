#!/usr/bin/env python3
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from common import MODEL_REPO_DIR, save_fig
from matplotlib import image
from plate_model_manager import PresentDayRasterManager

import gplately

print(gplately.__file__)

# gplately.turn_on_debug_logging()


def main(show=True, anchor_pid=0):
    model = gplately.auxiliary.get_plate_reconstruction(
        "Muller2019", model_repo_dir="plate-model-repo"
    )

    etopo = gplately.Raster(
        data=image.imread(PresentDayRasterManager().get_raster("ETOPO1_tif")),
        plate_reconstruction=model,
    )
    etopo.lats = etopo.lats[::-1]

    etopo_downscaled = etopo.resample(0.5, 0.5)
    assert etopo_downscaled

    fig = plt.figure(figsize=(10, 5), dpi=96)
    # plot 0Ma
    ax_1 = fig.add_subplot(221, projection=ccrs.PlateCarree())
    etopo_downscaled.imshow(ax_1, interpolation="none")
    ax_1.set_title(f"0 Ma")

    # reconstruct to 50 Ma
    r_50 = etopo_downscaled.reconstruct(
        50, fill_value="white", anchor_plate_id=anchor_pid
    )

    # plot 50 Ma
    ax_2 = fig.add_subplot(222, projection=ccrs.PlateCarree())
    r_50.imshow(ax_2, interpolation="none")  # type: ignore
    ax_2.set_title(f"50 Ma - APID:{anchor_pid}")

    # reconstruct to 100 Ma
    r_100 = etopo_downscaled.reconstruct(
        100, fill_value="white", anchor_plate_id=anchor_pid
    )

    # plot 100 Ma
    ax_3 = fig.add_subplot(223, projection=ccrs.PlateCarree())
    r_100.imshow(ax_3, interpolation="none")  # type: ignore
    ax_3.set_title(f"100 Ma - APID:{anchor_pid}")

    # reconstruct to 200 Ma
    r_200 = etopo_downscaled.reconstruct(
        200, fill_value="white", anchor_plate_id=anchor_pid
    )

    # plot 200 Ma
    ax_4 = fig.add_subplot(224, projection=ccrs.PlateCarree())
    r_200.imshow(ax_4, interpolation="none")  # type: ignore
    ax_4.set_title(f"200 Ma - APID:{anchor_pid}")

    if show:
        plt.show()
    else:
        save_fig(f"{__file__[:-3]}-{anchor_pid}.py")


if __name__ == "__main__":
    if len(sys.argv) == 2:
        if sys.argv[1] == "save":
            main(show=False)
        else:
            main(show=True, anchor_pid=int(sys.argv[1]))
    elif len(sys.argv) == 3:
        main(show=False, anchor_pid=int(sys.argv[1]))
    else:
        main(show=True)
