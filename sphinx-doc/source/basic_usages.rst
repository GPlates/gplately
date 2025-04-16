Basic Usages
============

.. contents::
   :local:
   :depth: 2

Download .rot file(s) from a model
----------------------------------

The Python code below downloads the ``rotation`` files from the ``Cao2024`` model into a local folder ``plate-models-data-dir`` 
and displays the file paths on screen.

.. code-block:: python
    :linenos:
    :emphasize-lines: 4,6

    from plate_model_manager import PlateModelManager

    mgr = PlateModelManager()
    cao2024_model = mgr.get_model("Cao2024", data_dir="plate-models-data-dir")
    # download the rotation files and print their local paths
    print(cao2024_model.get_rotation_model())

.. seealso::
    See the :ref:`list-all-models` section for how to get a list of available models.

Download a layer from a model
-----------------------------

The Python code below downloads the ``Coastlines`` layer from the ``Cao2024`` model 
and displays the file paths on screen.

.. code-block:: python
    :linenos:
    :emphasize-lines: 4,6
   
    from plate_model_manager import PlateModelManager

    mgr = PlateModelManager()
    cao2024_model = mgr.get_model("Cao2024", data_dir="plate-models-data-dir")
    # download Coastlines from model Cao2024 and display the local path
    print(cao2024_model.get_layer("Coastlines"))

.. seealso::
    See the :ref:`list-all-layers` section for how to get a list of available layers.

.. _list-all-layers:

List all layer names in a model
-------------------------------

The Python code below displays a list of available layers in the model Cao2024 on the screen.

.. code-block:: python
    :linenos:
    :emphasize-lines: 4,6
   
    from plate_model_manager import PlateModelManager

    mgr = PlateModelManager()
    cao2024_model = mgr.get_model("Cao2024", data_dir="plate-models-data-dir")
    # display a list of available layers in model Cao2024
    print(cao2024_model.get_avail_layers())

.. _list-all-models:

List all available model names
------------------------------

The Python code below displays a list of available plate models.

.. code-block:: python
    :linenos:
    :emphasize-lines: 3

    from plate_model_manager import PlateModelManager

    print(PlateModelManager().get_available_model_names())