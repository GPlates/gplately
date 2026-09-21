Command Line Interface (CLI)
============================

GPlately comes with a collection of useful command-line tools, each designed as a subcommand of GPlately. 

.. note::

   - Run ``gplately -h`` to see a list of available subcommands.
   - Run ``gplately one-of-the-subcommands -h`` to see the usage of the given subcommand. For example, ``gplately list -h``.

.. contents::
   :local:
   :depth: 2

list
----

Show a list of available plate reconstruction models. 
If a model name is given, show the details about this model.

👉 list all available plate reconstruction models

.. code:: console

   $ gplately list
   
👉 show details about model ``merdith2021``

.. code:: console

   $ gplately list -m merdith2021

If you are using GPlately Docker image, run these commands.

.. code:: console

   $ docker run gplates/gplately gplately list
   $ docker run gplates/gplately gplately list -m merdith2021

combine
-------

Combine multiple feature collections into one. 

👉 combine three feature collections into one file ``output_file.gpmlz``

.. code:: console

   $ gplately combine input_file_1.shp input_file_2.gpmlz input_file_3.gpml output_file.gpmlz   

filter
------

Filter feature collection by various criteria and save the result to the `output_file`.  

👉 only keep features whose name contains "Africa" or "North America" 

.. code:: console

   $ gplately filter input_file output_file -n Africa "North America"
  
👉 only keep features whose plate ID is one of "701 714 715 101"

.. code:: console

   $ gplately filter input_file output_file -p 701 714 715 101

👉 only keep features whose birth age is older than 500 Myr

.. code:: console

   $ gplately filter input_file output_file --min-birth-age 500

👉 only keep features whose birth age is younger than 500 Myr

.. code:: console

   $ gplately filter input_file output_file --max-birth-age 500
   
👉 only keep features whose name contains "Africa" or "North America" and plate ID is one of "701 714 715 101" and birth age is older than 500 Myr

.. code:: console

   $ gplately filter input_file output_file -n Africa "North America" -p 701 714 715 101 --min-birth-age 500
   
👉 only keep gpml:Basin features

.. code:: console

   $ gplately filter input_file output_file -t gpml:Basin
   
👉 only keep gpml:Basin and gpml:IslandArc features

.. code:: console

   $ gplately filter input_file output_file -t "gpml:IslandArc|gpml:Basin"
   
.. note::

   If you are using Docker, prefix ``docker run gplates/gplately`` to the command, such as ``docker run gplates/gplately gplately filter input_file output_file -t gpml:Basin``.

.. seealso::

   Check out `this shell script <https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_feature_filter.sh>`__ for more "gplately filter" examples. 

reset-feature-type (rft)
-------------------------

Reset the feature type for the selected features. Also available as ``reset_feature_type`` (the
old underscore-separated name still works, kept for backward compatibility).

👉 change all gpml:ClosedContinentalBoundary to gpml:UnclassifiedFeature

.. code:: console

   $ gplately reset-feature-type -s gpml:ClosedContinentalBoundary -t gpml:UnclassifiedFeature input_file output_file

👉 change all gpml:ContinentalFragment and gpml:Coastline to gpml:UnclassifiedFeature

.. code:: console

   $ gplately reset-feature-type -s "gpml:ContinentalFragment|gpml:Coastline" -t gpml:UnclassifiedFeature input_file output_file

👉 change all feature types to gpml:UnclassifiedFeature

.. code:: console

   $ gplately reset-feature-type -s ".*" -t gpml:UnclassifiedFeature input_file output_file

.. note::

  If you are using Docker, prefix ``docker run gplates/gplately`` to the command, such as ``docker run gplates/gplately gplately reset-feature-type -s ".*" -t gpml:UnclassifiedFeature input_file output_file``.

.. seealso::

  Check out `this shell script <https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_reset_feature_type.sh>`__ for more "gplately reset-feature-type" examples.

agegrid (ag)
------------

Generate seafloor age grids for a plate reconstruction model. 

👉 create age grids from 10 Ma to 0 Ma with 1 Myr increment for the `merdith2021` reconstruction model

.. code:: console

   $ gplately ag output -m merdith2021 -e 0 -s 10
   
👉 create age grids from 10 Ma to 0 Ma with 1 Myr increment using the specified reconstruction files

.. code:: console

   $ gplately ag rotations.rot topologies.gpmlz output -c continental_polygons.gpmlz -e 0 -s 10
   
.. note::

   To use a configuration file instead of specifying the command line arguments,
   create a TOML file (e.g. ``gplately-cli-config.toml``) and run the command below.

   .. code:: console

      $ gplately ag output-dir --config gplately-cli-config.toml

   Some example configuration files are available in the ``tests-dir/unittest/`` folder of the GPlately repository.
   And you can also run ``gplately get-cli-config-example`` to get an example configuration file.

paleobathymetry (pb)
---------------------

Run the *simple_paleobathymetry* workflow end to end: seafloor age → basement depth (Step 1),
distance to the nearest passive continental margin (Step 2), predicted sediment thickness (Step 3),
and isostatically-compensated paleobathymetry (Step 4). Optionally also merges in pyBacktrack's
present-day paleobathymetry (Step 5) to cover submerged continental crust and crust that has since
subducted.

👉 generate paleobathymetry grids from 10 Ma to 0 Ma for the `muller2025` reconstruction model

.. code:: console

   $ gplately pb output -m muller2025 --proximity-features cobs.gpml -e 0 -s 10

Distances route around continents by default, as in the original workflow, taking the
continent geometries from the plate model's Coastlines layer (or from
``--continent-obstacles``). Pass ``--no-route-around-continents`` for straight-line
great-circle distances, which are much faster but are not what the workflow computes.

👉 straight-line distances, and merge in pyBacktrack's present-day paleobathymetry

.. code:: console

   $ gplately pb output -m muller2025 --proximity-features cobs.gpml -e 0 -s 10 --no-route-around-continents --pybacktrack

.. note::

   The Plate Model Manager does not deliver passive-margin continent-ocean-boundary (COB) line
   segments for most models, so ``--proximity-features`` is usually required (unless you generate
   dynamically-contoured passive margins with ``gplately generate-passive-margins`` and pass those
   instead).

.. seealso::

   ``gplately generate-distance-grids`` / ``gplately generate-sediment-grids`` run Steps 2 and 3
   individually; ``gplately generate-passive-margins`` can generate a dynamically-contoured
   proximity target instead of a static COB file.

generate-distance-grids (gdg)
------------------------------

For each ocean point in a seafloor-age grid, reconstruct it backward through time and compute its
lifetime-mean distance to the nearest passive continental margin (Step 2 of the paleobathymetry
workflow, see ``gplately paleobathymetry``).

👉 generate distance-to-margin grids from 10 Ma to 0 Ma

.. code:: console

   $ gplately gdg output -m muller2025 --proximity-features cobs.gpml -e 0 -s 10

👉 straight-line distances instead of routing around continents (routing is the default)

.. code:: console

   $ gplately gdg output -m muller2025 --proximity-features cobs.gpml -e 0 -s 10 --no-route-around-continents

generate-sediment-grids (gsg)
------------------------------

Combine a seafloor-age grid with the distance-to-passive-margin grids from
``gplately generate-distance-grids`` into predicted compacted sediment-thickness grids
(Dutkiewicz et al., 2017) (Step 3 of the paleobathymetry workflow).

👉 predict sediment thickness from previously-generated distance grids

.. code:: console

   $ gplately gsg output -m muller2025 -e 0 -s 10 --distance-grids-dir distances/

generate-passive-margins (gpm)
---------------------------------

Dynamically contour continents through time and split each contour into passive-margin segments,
as an alternative to a static continent-ocean-boundary (COB) line-segment file for
``gplately generate-distance-grids``/``gplately paleobathymetry``.

👉 generate passive margins from 100 Ma to 0 Ma for the `muller2025` reconstruction model

.. code:: console

   $ gplately gpm output -m muller2025 -e 0 -s 100

fix-crossovers (fc)
---------------------

Fixes crossovers in rotation file(s). Also available as ``fix_crossovers``.

👉 fix crossovers in two rotation files with a threshold of 0.01 degrees and ignore plate IDs 201 and 701

.. code:: console

   $ gplately fix-crossovers -d -c 0.01 -i 201 701 -- input_rotations1.rot input_rotations2.rot


remove-rotations (rr)
------------------------

Remove one or more plate IDs from a rotation model (consisting of one or more rotation files).
Also available as ``remove_rotations``.

👉 remove plate IDs 70, 4, 3 and 1 from a rotation file

.. code:: console

   $ gplately remove-rotations -p 70 4 3 1 -o removed_ref_frames_ -- rotations.rot


cleanup-topologies (ct)
--------------------------

Remove any regular features not referenced by topological features. Also available as
``cleanup_topologies``.

👉 remove all features which are not referenced by any topological feature from topologies.gpml

.. code:: console

   $ gplately cleanup-topologies -o cleanup-topologies- -- topologies.gpml


convert-xy-to-gplates (cxg)
------------------------------

Converts geometry in one or more input ASCII files (such as '.xy' files) to output files suitable
for loading into GPlates. Also available as ``convert_xy_to_gplates``.

👉 convert two .xy files into a shapefile

.. code:: console

   $ gplately convert-xy-to-gplates -e shp -- input1.xy input2.xy


diagnose-rotations (dr)
--------------------------

Diagnose one or more rotation files to check for inconsistencies. Also available as
``diagnose_rotations``.

👉 check two rotation files and print the diagnostic results on screen

.. code:: console

   $ gplately diagnose-rotations input_rotations1.rot input_rotations2.rot


resolve-topologies (rt)
--------------------------

Resolve topological plate polygons (and deforming networks) and save (to separate files) the resolved topologies,
and their boundary sections as subduction zones, mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges).
Also available as ``resolve_topologies``.

👉 resolve topologies at 10 Ma

.. code:: console

   $ gplately resolve-topologies -r rotations1.rot rotations2.rot -m topologies1.gpml topologies2.gpml -t 10

rotation-tools (rots)
------------------------

Calculate stage rotations between consecutive finite rotations in plate pairs. Also available as
``rotation_tools``.

👉 calculate stage rotations for moving plate 701 relative to the fixed plate 0

.. code:: console

   $ gplately rotation-tools -p 701 0 -o stage_ -- rotations.rot


separate-ridge-transform-segments (srts)
-------------------------------------------

Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments.
Also available as ``separate_ridge_transform_segments``.

👉 pick out the ridge and transform features from the file spreading_features.gpml

.. code:: console

   $ gplately separate-ridge-transform-segments -r rotations.rot -d 45 -s _ridges -t _transforms -- spreading_features.gpml

.. note::

   There's also PlateReconstruction.crustal_production_destruction_rate() and PlateReconstruction.divergent_convergent_plate_boundaries()
   for users interested in convergent/divergent plate boundaries and crustal production/destruction - it's more accurate
   because it avoids the issue of whether features are labelled as mid-ocean ridges (or missing left/right plate IDs, etc)
   and the default ridge-transform separation angle of 70 degrees is just an empirically observed compromise to
   suit most cases (but not all cases).

subduction-convergence (sc)
------------------------------

Find the convergence rates along trenches (subduction zones) over time. Also available as
``subduction_convergence``.

👉 calculate the convergence rates along subduction zones from 200 Ma to 0 Ma

.. code:: console

   $ gplately subduction-convergence -r rotations.rot -m topologies.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- convergence


gpmdb
-----

Retrieve the paleomagnetic data from the `GPMDB website <http://www.gpmdb.net>`__, create GPlates-compatible VGP features and save them in a .gpmlz file. 

👉 download the paleomagnetic data and generate GPlates-compatible VGP features using the `zahirovic2022` reconstruction model

.. code:: console

   $ gplately gpmdb -m zahirovic2022 -o vgp.gpmlz
   

rotate-grid (rtg)
---------------------

Rotate a grid (or all grids in a folder) between plate-model reference frames. Also available as
``rotate_grid``.
The source and target rotation models can each be a named plate model or local rotation files (the two options are mutually exclusive).
If the reconstruction time is not given explicitly, it is deduced from the filename (e.g. ``paleobathymetry_103Ma.nc`` → 103 Ma).

👉 rotate ``input.nc`` at 100 Ma from the Alfonso2024 mantle frame (anchor 0) to the Alfonso2024 pmag frame (anchor 701701)

.. code:: console

   $ gplately rotate-grid input.nc output.nc --from-model Alfonso2024 --to-model Alfonso2024 --from-anchor 0 --to-anchor 701701 --time 100

👉 rotate all .nc files in a directory; reconstruction times are deduced from filenames such as ``paleobathymetry_103Ma.nc``

.. code:: console

   $ gplately rotate-grid input_dir output_dir --from-model Alfonso2024 --to-model Alfonso2024 --from-anchor 0 --to-anchor 701701

👉 rotate using local rotation files instead of a named model

.. code:: console

   $ gplately rotate-grid input.nc output.nc --from-rotation-files from.rot --to-rotation-files to.rot --time 100


get-cli-config-example (gcce)
---------------------------------

Print an example CLI configuration in TOML format to stdout (for use with --config); also
available as ``get_cli_config_example``. Redirect it to save, e.g.

.. code:: console

   $ gplately get-cli-config-example > my-gplately-cli-config.toml