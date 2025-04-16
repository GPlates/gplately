Examples
========

.. contents::
   :local:
   :depth: 2

Use PMM with pyGPlates 
----------------------

.. code-block:: python
    :linenos:
    :emphasize-lines: 14,15,24

    import pygplates

    from plate_model_manager import PlateModelManager

    pm_manager = PlateModelManager()
    model = pm_manager.get_model("Muller2019")

    # create a point feature at (0,0)
    point_feature = pygplates.Feature()
    point_feature.set_geometry(pygplates.PointOnSphere(0, 0))

    # assign plate ID
    point_feature_with_PID = pygplates.partition_into_plates(
        model.get_static_polygons(),  # 👈👀 LOOK HERE
        model.get_rotation_model(),  # 👈👀 LOOK HERE
        [point_feature],
    )

    # Reconstruct the point features.
    reconstructed_feature_geometries = []
    time = 140
    pygplates.reconstruct(
        point_feature_with_PID,
        model.get_rotation_model(),  # 👈👀 LOOK HERE
        reconstructed_feature_geometries,
        time,
    )
    print(reconstructed_feature_geometries[0].get_reconstructed_geometry().to_lat_lon())

.. seealso::
    `use PMM with pyGPlates example notebook`_

.. _use PMM with pyGPlates example notebook: https://github.com/GPlates/pygplates-tutorials/blob/master/notebooks/working-with-plate-model-manager.ipynb


Use PMM with GPlately 
---------------------

.. code-block:: python
    :linenos:
    :emphasize-lines: 16,17,21,22,26,28

    from gplately import (
    PlateModelManager,
    PlateReconstruction,
    PlotTopologies,
    PresentDayRasterManager,
    Raster,
    )

    model = PlateModelManager().get_model(
        "Muller2019",  # model name
        data_dir="plate-model-repo",  # the folder to save the model files
    )

    recon_model = PlateReconstruction(
        model.get_rotation_model(),
        topology_features=model.get_layer("Topologies"),
        static_polygons=model.get_layer("StaticPolygons"),
    )
    gplot = PlotTopologies(
        recon_model,
        coastlines=model.get_layer("Coastlines"),
        COBs=model.get_layer("COBs"),
        time=55,
    )
    # get present-day topography raster
    raster = Raster(PresentDayRasterManager().get_raster("topography"))
    # get paleo-agegrid raster at 100Ma from Muller2019 model
    agegrid = Raster(model.get_raster("AgeGrids", time=100))

.. seealso::
    `use PMM with GPlately example`_

.. _use PMM with GPlately example: https://github.com/GPlates/gplately/blob/master/Notebooks/Examples/working_with_plate_model_manager.py