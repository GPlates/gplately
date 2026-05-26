#!/usr/bin/env python3
# import matplotlib

# matplotlib.use("QtAgg")

import os
import sys

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"
import cartopy.crs as ccrs  # pyright: ignore[reportMissingImports]
import matplotlib.pyplot as plt  # pyright: ignore[reportMissingModuleSource]
from common import save_fig

import gplately
import pygplates  # pyright: ignore[reportMissingModuleSource]

print(gplately.__file__)

point_geom = pygplates.PointOnSphere(-10, -20)  # type: ignore
line_geom = pygplates.PolylineOnSphere(  # type: ignore
    [pygplates.PointOnSphere(0, 0), pygplates.PointOnSphere(10, 30), pygplates.PointOnSphere(20, 40)]  # type: ignore
)
polygon_geom = pygplates.PolygonOnSphere(
    [
        pygplates.PointOnSphere(5, 5),  # type: ignore
        pygplates.PointOnSphere(10, 25),  # type: ignore
        pygplates.PointOnSphere(35, 30),  # type: ignore
        pygplates.PointOnSphere(45, -60),  # type: ignore
    ]
)
multipoint_geom = pygplates.MultiPointOnSphere(  # type: ignore
    [pygplates.PointOnSphere(-20, -20), pygplates.PointOnSphere(40, 0), pygplates.PointOnSphere(50, 30)]  # type: ignore
)
feature1 = pygplates.Feature()  # type: ignore
feature1.set_geometry(point_geom)
feature2 = pygplates.Feature()  # type: ignore
feature2.set_geometry(line_geom)
feature3 = pygplates.Feature()  # type: ignore
feature3.set_geometry(polygon_geom)
feature4 = pygplates.Feature()  # type: ignore
feature4.set_geometry(multipoint_geom)

test_features = [feature1, feature2, feature3, feature4]


def test_plot_pygplates_features_cartopy(show):
    fig = plt.figure(figsize=(10, 5), dpi=96)
    ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=0))

    ax.set_global()  # type: ignore
    ax.gridlines(  # type: ignore
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.8,
        color="gray",
        alpha=0.6,
        linestyle="--",
    )

    cartopy_plot_engine = gplately.CartopyPlotEngine()
    cartopy_plot_engine.plot_pygplates_features(ax, test_features)

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


def test_plot_pygplates_features_pygmt(show):
    import pygmt  # pyright: ignore[reportMissingImports]

    from gplately.auxiliary import get_pygmt_basemap_figure

    with pygmt.config(MAP_TITLE_OFFSET="-7p"):
        fig = get_pygmt_basemap_figure(
            projection="N0/10c",
            region="d",
            frame=["xafg30", "yafg30", "+tPyGMT Plot of Pygplates Features"],
        )

    pygmt_plot_engine = gplately.PygmtPlotEngine()
    pygmt_plot_engine.plot_pygplates_features(
        fig,
        test_features,
        edgecolor="blue",
        facecolor="none",
        gmtlabel="Lines and Polygons",
        pointlabel="Points",
    )
    fig.legend(position="jBL+o-1.1c/-.55c", box="+gwhite+p0.25p")

    if show:
        fig.show(crop="+m0.4c")
    else:
        output_file = "./output/test_plot_pygplates_features_pygmt.pdf"
        fig.savefig(output_file, crop="+m0.4c")
        print(f"Done! The {output_file} has been saved.")


def main(show=True):
    test_plot_pygplates_features_cartopy(show)
    test_plot_pygplates_features_pygmt(show)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "save":
        main(show=False)
    else:
        main(show=True)
