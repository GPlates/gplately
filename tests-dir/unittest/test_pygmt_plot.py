#!/usr/bin/env python3

from gplately.auxiliary import get_gplot
from gplately.mapping.pygmt_plot import PygmtPlotEngine, get_pygmt_basemap_figure

if __name__ == "__main__":
    gplot = get_gplot(
        "merdith2021", "plate-model-repo", time=55, plot_engine=PygmtPlotEngine()
    )
    fig = get_pygmt_basemap_figure(projection="N180/10c", region="d")
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

    fig.text(
        text="55Ma (Merdith2021)",
        position="TC",
        no_clip=True,
        font="12p,Helvetica,black",
        offset="j0/-0.5c",
    )
    fig.legend(position="jBL+o-2.7/0", box="+gwhite+p0.5p")

    # fig.show(width=1200)
    fig.savefig("test-pygmt-plot.pdf")
