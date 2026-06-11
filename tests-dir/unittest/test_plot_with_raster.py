#!/usr/bin/env python3
import os
import sys

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false
import xarray as xr

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

from gplately.auxiliary import get_gplot

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from common import MODEL_REPO_DIR, save_fig, _get_topo_raster

import gplately

print(gplately.__file__)


def main(show=True):
    gplot = get_gplot("Muller2025", model_repo_dir=MODEL_REPO_DIR, time=55)

    agegrid = gplately.Raster(
        data=gplot.plate_reconstruction.plate_model.get_raster(
            "AgeGrids", time=gplot.time
        )
    )

    fig = plt.figure(figsize=(12, 6), dpi=96)
    ax1 = fig.add_subplot(111, projection=ccrs.Mollweide(180))

    gplot.plot_continents(ax1, facecolor="0.8")
    gplot.plot_coastlines(ax1, color="0.5")
    gplot.plot_ridges(ax1, color="red")
    gplot.plot_trenches(ax1, color="k")
    gplot.plot_subduction_teeth(ax1, color="k")

    im = gplot.plot_grid(
        ax1,
        agegrid.data,
        cmap="YlGnBu",
        shading=_get_topo_raster(),
        # shading=True,
        vmin=0,
        vmax=200,
    )

    gplot.plot_plate_motion_vectors(
        ax1, spacingX=10, spacingY=10, normalise=True, zorder=10, alpha=0.5
    )

    fig.colorbar(im, orientation="horizontal", shrink=0.4, pad=0.05, label="Age (Ma)")
    plt.title(f"{gplot.time} Ma")

    if show:
        plt.show()
    else:
        save_fig(__file__)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
