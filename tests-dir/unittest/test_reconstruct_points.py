#!/usr/bin/env python3
import os
import sys

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from common import MODEL_REPO_DIR, save_fig

import gplately
from gplately.auxiliary import get_gplot

print(gplately.__file__)


def main(show=True):
    gplot = get_gplot("Muller2019", time=0)
    model = gplot.plate_reconstruction

    pt_lons = np.array([140.0, 47, 13, 78])
    pt_lats = np.array([-30.0, 22, 42, 23])
    gpts = gplately.Points(model, pt_lons, pt_lats)

    fig = plt.figure(figsize=(16, 8))

    ax_1 = fig.add_subplot(121, projection=ccrs.Mollweide(0))
    ax_1.set_global()  # type: ignore
    gplot.plot_coastlines(ax_1, color="0.8")
    ax_1.plot(pt_lons, pt_lats, "o", transform=ccrs.PlateCarree())
    ax_1.set_title(f"{int(gplot.time)} Ma")  # type: ignore

    ax_2 = fig.add_subplot(122, projection=ccrs.Mollweide(0))
    ax_2.set_global()  # type: ignore
    r_time = 100
    gplot.time = r_time
    gplot.plot_coastlines(ax_2, color="0.8")
    rlons, rlats = gpts.reconstruct(r_time, return_array=True)  # type: ignore
    ax_2.plot(rlons, rlats, "o", transform=ccrs.PlateCarree())
    ax_2.set_title(f"{r_time} Ma")

    fig.tight_layout()
    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
