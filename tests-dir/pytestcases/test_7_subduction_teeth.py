"""Tests for the `gplately.plot.SubductionTeeth` class."""
import os
from itertools import product

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pytest
from gplately import PlotTopologies

from conftest import (
    projections,
    reconstruction_times,
    gplately_plot_topologies_object as gplot,
)


@pytest.mark.parametrize(
    "time,projection",
    product(reconstruction_times, projections),
)
def test_plot_topologies_subduction_teeth(
    time: float,
    gplot: PlotTopologies,
    projection: ccrs.Projection,
):
    """Test the `PlotTopologies.plot_subduction_teeth` method."""
    gplot.time = time
    fig, ax = plt.subplots(subplot_kw={"projection": projection})
    ax.set_global()
    teeth = gplot.plot_subduction_teeth(
        ax,
        size=5.5,
        # Test a few keyword args
        color="red",
        zorder=-1,
        alpha=0.3,
    )
    assert len(teeth.left_projected) > 0
    assert len(teeth.right_projected) > 0
    assert teeth.figure is fig
    assert teeth.ax is ax
    if os.environ.get("SAVEFIG", False):
        outdir = os.path.join(os.path.dirname(__file__), "test_outputs")
        os.makedirs(outdir, exist_ok=True)
        filename = os.path.join(
            outdir,
            (
                "test_plot_topologies_subduction_teeth_"
                + f"{type(projection).__name__}_"
                + f"{time:0.0f}Ma.png"
            ),
        )
        fig.savefig(
            filename,
            dpi=200,
            bbox_inches="tight",
        )
    plt.close(fig)
