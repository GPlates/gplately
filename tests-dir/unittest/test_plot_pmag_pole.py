#!/usr/bin/env python3
# import matplotlib

# matplotlib.use("QtAgg")

import os
import sys

from gplately.mapping.cartopy_plot import _create_a_basic_cartopy_ax
from gplately.mapping.pygmt_plot import PygmtPlotEngine

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
# pyright: reportMissingImports=false

import cartopy.crs as ccrs  # pyright: ignore[reportMissingImports]
import matplotlib.pyplot as plt  # pyright: ignore[reportMissingModuleSource]
from common import MODEL_REPO_DIR, save_fig

import gplately
import pygmt  # pyright: ignore[reportMissingImports]
import cartopy.feature as cfeature

from gplately.auxiliary import get_pygmt_basemap_figure, get_gplot

print(gplately.__file__)


def test_plot_pole_cartopy(show):
    gplot = get_gplot(
        "muller2025",
        model_repo_dir=MODEL_REPO_DIR,
        time=0,
    )
    ax = _create_a_basic_cartopy_ax(
        figsize=(12, 12),
        projection=ccrs.Orthographic(central_longitude=-80, central_latitude=-20),
    )
    ax.add_feature(cfeature.OCEAN, facecolor="lightblue")
    ax.add_feature(cfeature.LAND, facecolor="whitesmoke")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor="gray")
    ax.add_feature(cfeature.BORDERS, linewidth=0.4, edgecolor="silver")
    gplot.plot_pole(
        ax,
        lon=-80,
        lat=-20,
        a95=15,
        color="red",
    )
    ax.set_title("Paleomagnetic Pole Plot with Cartopy")
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


def test_plot_pole_pygmt(show):
    gplot = get_gplot(
        "muller2025",
        model_repo_dir=MODEL_REPO_DIR,
        time=0,
        plot_engine=PygmtPlotEngine(),
    )

    with pygmt.config(MAP_TITLE_OFFSET="-7p"):
        fig = get_pygmt_basemap_figure(
            projection="G-80/-20/15c",
            region="d",
            frame=["af", "+tPaleomagnetic Pole Plot with PyGMT"],
        )

        fig.coast(
            land="gray90",
            water="lightblue",
            shorelines="1/0.3p,gray50",
            borders="1/0.4p,gray70",
        )
        gplot.plot_pole(
            fig,
            lon=-80,
            lat=-20,
            a95=15,
            color="green",
        )

    if show:
        fig.show(crop="+m0.4c")
    else:
        output_file = "./output/test_plot_pole_pygmt.pdf"
        fig.savefig(output_file, crop="+m0.4c")
        print(f"Done! The {output_file} has been saved.")


def main(show=True, use_pygmt=False):
    if use_pygmt:
        test_plot_pole_pygmt(show)
    else:
        test_plot_pole_cartopy(show)


if __name__ == "__main__":
    show = True
    use_pygmt = False
    if "save" in sys.argv:
        show = False
    if "pygmt" in sys.argv:
        use_pygmt = True

    main(show=show, use_pygmt=use_pygmt)
