#!/usr/bin/env python3

# This test script generates a sample plot using the PygmtPlotEngine.

import os
from pathlib import Path

# pyright: reportMissingImports=false

import pygmt
import xarray as xr

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

from gplately.auxiliary import get_gplot, get_pygmt_basemap_figure
from gplately.mapping.pygmt_plot import PygmtPlotEngine
from gplately import Raster
from plate_model_manager import PresentDayRasterManager

if __name__ == "__main__":
    model_name = "muller2025"
    reconstruction_time = 55

    gplot = get_gplot(
        model_name,
        "plate-model-repo",
        time=reconstruction_time,
        plot_engine=PygmtPlotEngine(),
    )
    fig = get_pygmt_basemap_figure(
        projection="N180/10c", region="d", frame=["xafg30", "yafg30"]
    )

    # reconstruct the topography grid for the specified reconstruction time,
    # and use it to create an illumination grid for shading
    topo_file = PresentDayRasterManager(data_dir="./unittest-data").get_raster(
        "topography"
    )
    topo_grid = Raster(
        data=topo_file, plate_reconstruction=gplot.plate_reconstruction
    ).reconstruct(time=reconstruction_time)

    gplot.plot_grid(
        fig,
        "AgeGrids",
        cmap="create-age-grids-video/agegrid.cpt",
        nan_transparent=True,
        # shading=True,
        # shading="+a315+ne0.6",
        shading=topo_grid.to_data_array(),
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

    gplot.plot_plate_motion_vectors(fig, normalise=True, color="red")

    fig.text(
        text=f"{reconstruction_time} Ma",
        position="TC",
        no_clip=True,
        font="12p,Helvetica,black",
        offset="j0/-0.5c",
    )
    with pygmt.config(FONT_ANNOT_PRIMARY=4):
        fig.legend(position="jBL+o-1.0/0", box="+gwhite+p0.25p")

    fig.show(width=1200, crop="+m0.4c")
    output_file = "./output/test-pygmt-plot.pdf"
    fig.savefig(output_file, crop="+m0.4c")
    print(f"The figure has been saved to: {output_file}.")
