#!/usr/bin/env python3

# This test script generates a sample plot using the PygmtPlotEngine.

import os

import pygmt

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

from gplately.auxiliary import get_gplot, get_pygmt_basemap_figure
from gplately.mapping.pygmt_plot import PygmtPlotEngine

if __name__ == "__main__":
    model_name = "muller2019"
    reconstruction_name = 55

    gplot = get_gplot(
        model_name,
        "plate-model-repo",
        time=reconstruction_name,
        plot_engine=PygmtPlotEngine(),
    )
    fig = get_pygmt_basemap_figure(projection="N180/10c", region="d")

    gplot.plot_grid(
        fig, "AgeGrids", cmap="create-age-grids-video/agegrid.cpt", nan_transparent=True
    )

    # fig.coast(shorelines=True)

    gplot.plot_coastlines(
        fig,
        edgecolor="none",
        facecolor="gray",
        linewidth=0.1,
        central_meridian=180,
        gmtlabel="Coastlines",
    )
    gplot.plot_topological_plate_boundaries(
        fig,
        edgecolor="black",
        linewidth=0.25,
        central_meridian=180,
        gmtlabel="Plate Boundaries",
    )
    gplot.plot_ridges(fig, pen="0.5p,black", gmtlabel="Ridges")
    gplot.plot_transforms(fig, pen="0.5p,green", gmtlabel="Transforms")
    gplot.plot_subduction_teeth(fig, color="blue", gmtlabel="Subduction Zones")

    try:
        gplot.plot_plate_motion_vectors(fig)
    except NotImplementedError as e:
        print(e)

    fig.text(
        text=f"{reconstruction_name}Ma",
        position="TC",
        no_clip=True,
        font="12p,Helvetica,black",
        offset="j0/-0.5c",
    )
    with pygmt.config(FONT_ANNOT_PRIMARY=4):
        fig.legend(position="jBL+o-1.0/0", box="+gwhite+p0.25p")

    # fig.show(width=1200)
    output_file = "./output/test-pygmt-plot.pdf"
    fig.savefig(output_file)
    print(f"The figure has been saved to: {output_file}.")
