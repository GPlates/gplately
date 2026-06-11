#!/usr/bin/env python3

# This test script generates a sample plot using the PygmtPlotEngine.

import os
from pathlib import Path
from urllib.request import urlretrieve

import pygmt  # pyright: ignore[reportMissingImports]
import xarray as xr  # pyright: ignore[reportMissingImports]

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

from gplately.auxiliary import get_gplot, get_pygmt_basemap_figure
from gplately.mapping.pygmt_plot import PygmtPlotEngine


def _get_topo_raster():
    data_dir = Path(__file__).resolve().parent / "test-pygmt-plot-data"
    data_dir.mkdir(parents=True, exist_ok=True)
    topo_file = data_dir / "topo15-3601x1801.nc"
    topo_url = "https://repo.gplates.org/webdav/mchin/data/topo15-3601x1801.nc"

    if not topo_file.exists():
        print(f"Downloading {topo_file.name} ...")
        urlretrieve(topo_url, topo_file)

    try:
        raster_xr = xr.open_dataarray(topo_file)
    except ValueError:
        # Fallback for NetCDF files with multiple variables.
        dataset = xr.open_dataset(topo_file)
        first_var = next(iter(dataset.data_vars))
        raster_xr = dataset[first_var]
    return raster_xr


if __name__ == "__main__":
    model_name = "muller2019"
    reconstruction_time = 55

    gplot = get_gplot(
        model_name,
        "plate-model-repo",
        time=reconstruction_time,
        plot_engine=PygmtPlotEngine(),
    )
    fig = get_pygmt_basemap_figure(projection="N180/10c", region="d")

    # Warning: the topography raster is present-day,
    # so it may not be appropriate to plot it with a past reconstruction time.
    # This is just for testing the plotting of raster data with pygmt.
    illumination = pygmt.grdgradient(grid=_get_topo_raster(), radiance=[315, 45])
    gplot.plot_grid(
        fig,
        "AgeGrids",
        cmap="create-age-grids-video/agegrid.cpt",
        nan_transparent=True,
        # shading=True,
        # shading="+a315+ne0.6",
        shading=illumination,
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
