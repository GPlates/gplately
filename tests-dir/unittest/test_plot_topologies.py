#!/usr/bin/env python3
# import matplotlib

# matplotlib.use("QtAgg")

import os
import sys

from gplately.mapping.cartopy_plot import _create_a_basic_cartopy_ax
from gplately.auxiliary import get_gplot

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import matplotlib.pyplot as plt
from common import MODEL_REPO_DIR, save_fig

import gplately

print(gplately.__file__)

MODEL_NAME = "Muller2025"


def main(show=True):
    age = 55
    gplot = get_gplot(MODEL_NAME, model_repo_dir=MODEL_REPO_DIR, time=55)
    ax = _create_a_basic_cartopy_ax()

    other_kwargs = {"color": "lightgrey", "linewidth": 0.8}
    ridge_kwargs = {"color": "red", "linewidth": 1.0}
    transform_kwargs = {"color": "green", "linewidth": 1.0}
    trench_kwargs = {"color": "blue", "linewidth": 1.0}

    gplot.plot_all_topological_sections(
        ax,
        plot_subduction_teeth=True,
        other_kwargs=other_kwargs,
        ridge_kwargs=ridge_kwargs,
        transform_kwargs=transform_kwargs,
        trench_kwargs=trench_kwargs,
    )

    # gplot.plot_all_topological_sections(ax)
    # gplot.plot_misc_transforms(ax=ax, color="red", linewidth=0.5)
    # gplot.misc_transforms
    # gplot.get_misc_transforms()
    # gplot.get_transforms()
    # gplot.plot_all_topological_sections(ax, color="0.5", linewidth=0.5)
    # gplot.plot_ridges(ax, color="red")
    # gplot.plot_transforms(ax, color="goldenrod")
    # gplot.plot_subduction_teeth(ax, color="blue")

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


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
