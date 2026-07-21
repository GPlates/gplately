#!/usr/bin/env python3
# import matplotlib

# matplotlib.use("QtAgg")

import os
import sys

from gplately.mapping.cartopy_plot import _create_a_basic_cartopy_ax
from gplately.mapping.pygmt_plot import PygmtPlotEngine
from plate_model_manager import PresentDayRasterManager

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
# pyright: reportMissingImports=false

import matplotlib.pyplot as plt  # pyright: ignore[reportMissingModuleSource]
from common import MODEL_REPO_DIR, save_fig

import gplately
import pygmt  # pyright: ignore[reportMissingImports]

from gplately.auxiliary import get_pygmt_basemap_figure, get_gplot

print(gplately.__file__)

topo_file = PresentDayRasterManager(data_dir="./unittest-data").get_raster("topography")


def test_plot_netcdf_cartopy(show):
    gplot = get_gplot(
        "muller2025",
        model_repo_dir=MODEL_REPO_DIR,
        time=0,
    )
    ax = _create_a_basic_cartopy_ax(
        figsize=(12, 12),
    )

    gplot.plot_grid_from_netCDF(
        ax,
        filename=topo_file,
        # var_name="spreading_rate", https://repo.gplates.org/webdav/PlateModel_Age_SR_Grids/Alfonso_etal_2024_modClennettMuller/archive/2025-09-15/Agegrids/
        shading=topo_file,
        cmap="terrain",
    )
    ax.set_title(" NetCDF Plot with Cartopy")
    if show:
        # LOOK HERE! 👀👀 👇👇
        # If the figure did not show up, you need to set your matplotlib plotting backend properly.
        # On Windows, you may install PyQt and do
        # import matplotlib
        # matplotlib.use('QtAgg')

        # if you are interested in finding what backends available on your computer and what is your current backend, do the following
        # import matplotlib.rcsetup as rcsetup
        # print(rcsetup.all_backends) # get all available backends
        # import matplotlib
        # matplotlib.get_backend() # your current backend
        #
        plt.show()
    else:
        save_fig(__file__)


def test_plot_netcdf_pygmt(show):
    gplot = get_gplot(
        "muller2025",
        model_repo_dir=MODEL_REPO_DIR,
        time=0,
        plot_engine=PygmtPlotEngine(),
    )

    with pygmt.config(MAP_TITLE_OFFSET="-7p"):
        fig = get_pygmt_basemap_figure(
            region="d",
            frame=["af", "+tNetCDF Plot with PyGMT"],
        )

        gplot.plot_grid_from_netCDF(
            fig,
            filename=topo_file,
            # var_name="seafloor_age", #https://repo.gplates.org/webdav/PlateModel_Age_SR_Grids/Alfonso_etal_2024_modClennettMuller/archive/2025-09-15/Agegrids/
            shading=topo_file,
            cmap="geo",
        )

    if show:
        fig.show(crop="+m0.4c")
    else:
        output_file = "./output/test_plot_netcdf_pygmt.pdf"
        fig.savefig(output_file, crop="+m0.4c")  # type: ignore
        print(f"Done! The {output_file} has been saved.")


def main(show=True, use_pygmt=False):
    if use_pygmt:
        test_plot_netcdf_pygmt(show)
    else:
        test_plot_netcdf_cartopy(show)


if __name__ == "__main__":
    show = True
    use_pygmt = False
    if "save" in sys.argv:
        show = False
    if "pygmt" in sys.argv:
        use_pygmt = True

    main(show=show, use_pygmt=use_pygmt)
