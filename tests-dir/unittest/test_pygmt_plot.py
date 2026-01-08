#!/usr/bin/env python3
import os

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

from gplately import Raster
from gplately.auxiliary import get_gplot, get_pygmt_basemap_figure
from gplately.mapping.pygmt_plot import PygmtPlotEngine

if __name__ == "__main__":
    gplot = get_gplot(
        "muller2019", "plate-model-repo", time=55, plot_engine=PygmtPlotEngine()
    )
    fig = get_pygmt_basemap_figure(projection="N180/10c", region="d")

    age_grid_raster = Raster(
        data=gplot.plate_reconstruction.plate_model.get_raster("AgeGrids", 55),
        plate_reconstruction=gplot.plate_reconstruction,
        extent=(-180, 180, -90, 90),
    )
    gplot.plot_grid(fig, age_grid_raster, nan_transparent=True)

    # fig.coast(shorelines=True)

    gplot.plot_topological_plate_boundaries(
        fig,
        edgecolor="black",
        linewidth=0.25,
        central_meridian=180,
        gmtlabel="plate boundaries",
    )
    gplot.plot_coastlines(
        fig, edgecolor="none", facecolor="gray", linewidth=0.1, central_meridian=180
    )
    gplot.plot_ridges(fig, pen="0.5p,red", gmtlabel="ridges")
    gplot.plot_transforms(fig, pen="0.5p,red", gmtlabel="transforms")
    gplot.plot_subduction_teeth(fig, color="blue", gmtlabel="subduction zones")

    try:
        gplot.plot_plate_motion_vectors(fig)
    except NotImplementedError as e:
        print(e)

    fig.text(
        text="55Ma (Merdith2021)",
        position="TC",
        no_clip=True,
        font="12p,Helvetica,black",
        offset="j0/-0.5c",
    )
    fig.legend(position="jBL+o-2.7/0", box="+gwhite+p0.5p")

    # fig.show(width=1200)
    fig.savefig("./output/test-pygmt-plot.pdf")
