#!/usr/bin/env python3
# import matplotlib

# matplotlib.use("QtAgg")

import os
import sys

from gplately.mapping.cartopy_plot import _create_a_basic_cartopy_ax
from gplately.auxiliary import get_gplot, get_pygmt_basemap_figure
from gplately.mapping.pygmt_plot import PygmtPlotEngine
import pygmt

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import matplotlib.pyplot as plt
from common import MODEL_REPO_DIR, save_fig

import gplately

print(gplately.__file__)

MODEL_NAME = "Muller2025"

age = 55


def plot_with_pygmt(show=True):
    gplot = get_gplot(
        MODEL_NAME,
        model_repo_dir=MODEL_REPO_DIR,
        time=age,
        plot_engine=PygmtPlotEngine(),
    )
    fig = get_pygmt_basemap_figure(
        projection="N180/10c", region="d", frame=["xafg30", "yafg30"]
    )

    gplot.plot_all_topological_sections(
        fig,
        plot_subduction_teeth=True,
        other_kwargs={
            "color": "lightgrey",
            "linewidth": 0.5,
            "gmtlabel": "Miscellaneous",
        },
        ridge_kwargs={"color": "red", "linewidth": 0.7, "gmtlabel": "Ridges"},
        transform_kwargs={"color": "green", "linewidth": 0.7, "gmtlabel": "Transforms"},
        trench_kwargs={
            "color": "blue",
            "linewidth": 0.7,
            "gmtlabel": "Subduction Zones",
        },
    )
    fig.text(
        text=f"{age} Ma",
        position="TC",
        no_clip=True,
        font="12p,Helvetica,black",
        offset="j0/-0.5c",
    )
    with pygmt.config(FONT_ANNOT_PRIMARY="5p"):
        fig.legend(position="jBL+o-1.1c/-.55c+w2.0c", box="+gwhite+p0.25p")
    fig.show(crop="+m0.4c")
    fig.savefig("./output/test-pygmt-plot-topologies.pdf", crop="+m0.4c")


def plot_with_cartopy(show=True):
    gplot = get_gplot(
        MODEL_NAME,
        model_repo_dir=MODEL_REPO_DIR,
        time=age,
    )
    ax = _create_a_basic_cartopy_ax()

    gplot.plot_all_topological_sections(
        ax,
        plot_subduction_teeth=True,
        other_kwargs={"color": "lightgrey", "linewidth": 0.8},
        ridge_kwargs={"color": "red", "linewidth": 1.0},
        transform_kwargs={"color": "green", "linewidth": 1.0},
        trench_kwargs={"color": "blue", "linewidth": 1.0},
    )

    ax.set_title(f"{age} Ma")

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


def main(show=True, use_pygmt=False):
    if use_pygmt:
        plot_with_pygmt(show=show)
    else:
        plot_with_cartopy(show=show)


if __name__ == "__main__":
    show = True
    use_pygmt = False
    if "save" in sys.argv:
        show = False
    if "pygmt" in sys.argv:
        use_pygmt = True

    main(show=show, use_pygmt=use_pygmt)
