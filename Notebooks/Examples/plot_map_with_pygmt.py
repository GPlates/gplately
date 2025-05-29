#!/usr/bin/env python3

from gplately.auxiliary import get_gplot, get_pygmt_basemap_figure
from gplately.mapping.pygmt_plot import PygmtPlotEngine

# for now, the pygmt integration is still pretty basic.
# please create GitHub issues and let us know how we can enhance the pygmt integration.
# we are grateful for your constructive feedbacks. Thank you.
# https://github.com/GPlates/gplately/issues

if __name__ == "__main__":
    # tell PlotTopologies object to use the PygmtPlotEngine
    gplot = get_gplot(
        "merdith2021", "plate-model-repo", time=55, plot_engine=PygmtPlotEngine()
    )
    # you need to know how to specify projection and region in pygmt way
    fig = get_pygmt_basemap_figure(projection="N180/10c", region="d")

    # now you can plot some features with the PlotTopologies object
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
        # plotting grid has not been implemented in PygmtPlotEngine yet.
        gplot.plot_grid(fig, None)
    except NotImplementedError as e:
        # print(e)
        pass

    try:
        # plotting velocities has not been implemented in PygmtPlotEngine yet.
        gplot.plot_plate_motion_vectors(fig)
    except NotImplementedError as e:
        # print(e)
        pass

    # use pygmt directly to plot title and legend
    fig.text(
        text="55Ma (Merdith2021)",
        position="TC",
        no_clip=True,
        font="12p,Helvetica,black",
        offset="j0/-0.5c",
    )
    fig.legend(position="jBL+o-2.7/0", box="+gwhite+p0.5p")

    # fig.show(width=1200)
    out_f = "test-pygmt-plot.pdf"
    fig.savefig(out_f)
    print(f"the file {out_f} has been saved.")
