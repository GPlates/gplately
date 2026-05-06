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

reset_feature_type
------------------

Reset the feature type for the selected features. 

👉 change all gpml:ClosedContinentalBoundary to gpml:UnclassifiedFeature

.. code:: console

   $ gplately reset_feature_type -s gpml:ClosedContinentalBoundary -t gpml:UnclassifiedFeature input_file output_file
   
👉 change all gpml:ContinentalFragment and gpml:Coastline to gpml:UnclassifiedFeature

.. code:: console

   $ gplately reset_feature_type -s "gpml:ContinentalFragment|gpml:Coastline" -t gpml:UnclassifiedFeature input_file output_file

👉 change all feature types to gpml:UnclassifiedFeature    

.. code:: console

   $ gplately reset_feature_type -s ".*" -t gpml:UnclassifiedFeature input_file output_file
    
.. note::

  If you are using Docker, prefix ``docker run gplates/gplately`` to the command, such as ``docker run gplates/gplately gplately reset_feature_type -s ".*" -t gpml:UnclassifiedFeature input_file output_file``.

.. seealso::

  Check out `this shell script <https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_reset_feature_type.sh>`__ for more "gplately reset_feature_type" examples. 

agegrid (ag)
------------

Generate age grids for a plate reconstruction model. 

👉 create age grids from 10 Ma to 0 Ma with 1 Myr increment for the `merdith2021` reconstruction model

.. code:: console

   $ gplately ag output -m merdith2021 -e 0 -s 10
   
👉 create age grids from 10 Ma to 0 Ma with 1 Myr increment using the specified reconstruction files

.. code:: console

   $ gplately ag rotations.rot topologies.gpmlz output -c continental_polygons.gpmlz -e 0 -s 10
   

fix_crossovers
--------------

Fixes crossovers in rotation file(s).

👉 fix crossovers in two rotation files with a threshold of 0.01 degrees and ignore plate IDs 201 and 701

.. code:: console

   $ gplately fix_crossovers -d -c 0.01 -i 201 701 -- input_rotations1.rot input_rotations2.rot
   

remove_rotations
----------------

Remove one or more plate IDs from a rotation model (consisting of one or more rotation files). 

👉 remove plate IDs 70, 4, 3 and 1 from a rotation file

.. code:: console

   $ gplately remove_rotations -p 70 4 3 1 -o removed_ref_frames_ -- rotations.rot
  

cleanup_topologies
------------------

Remove any regular features not referenced by topological features. 

👉 remove all features which are not referenced by any topological feature from topologies.gpml

.. code:: console

   $ gplately cleanup_topologies -o cleanup_topologies_ -- topologies.gpml
   

convert_xy_to_gplates
---------------------

Converts geometry in one or more input ASCII files (such as '.xy' files) to output files suitable for loading into GPlates. 

👉 convert two .xy files into a shapefile

.. code:: console

   $ gplately convert_xy_to_gplates -e shp -- input1.xy input2.xy
  

diagnose_rotations
------------------

Diagnose one or more rotation files to check for inconsistencies. 

👉 check two rotation files and print the diagnostic results on screen

.. code:: console

   $ gplately diagnose_rotations input_rotations1.rot input_rotations2.rot
   

resolve_topologies
------------------

Resolve topological plate polygons (and deforming networks) and save (to separate files) the resolved topologies, 
and their boundary sections as subduction zones, mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges). 

👉 resolve topologies at 10 Ma

.. code:: console

   $ gplately resolve_topologies -r rotations1.rot rotations2.rot -m topologies1.gpml topologies2.gpml -t 10

rotation_tools
--------------

Calculate stage rotations between consecutive finite rotations in plate pairs. 

👉 calculate stage rotations for moving plate 701 relative to the fixed plate 0

.. code:: console

   $ gplately rotation_tools -p 701 0 -o stage_ -- rotations.rot
   

separate_ridge_transform_segments
---------------------------------

Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments. 

👉 pick out the ridge and transform features from the file spreading_features.gpml

.. code:: console

   $ gplately separate_ridge_transform_segments -r rotations.rot -d 45 -s _ridges -t _transforms -- spreading_features.gpml
    

subduction_convergence
----------------------

Find the convergence rates along trenches (subduction zones) over time. 

👉 calculate the convergence rates along subduction zones from 200 Ma to 0 Ma

.. code:: console

   $ gplately subduction_convergence -r rotations.rot -m topologies.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- convergence
   

gpmdb
-----

Retrieve the paleomagnetic data from the `GPMDB website <http://www.gpmdb.net>`__, create GPlates-compatible VGP features and save them in a .gpmlz file. 

👉 download the paleomagnetic data and generate GPlates-compatible VGP features using the `zahirovic2022` reconstruction model

.. code:: console

   $ gplately gpmdb -m zahirovic2022 -o vgp.gpmlz
   

rotate_grid (rg2)
-----------------

Rotate a grid (or all grids in a folder) between plate-model reference frames.
The source and target rotation models can each be a named plate model or local rotation files (the two options are mutually exclusive).
If the reconstruction time is not given explicitly, it is deduced from the filename (e.g. ``paleobathymetry_103Ma.nc`` → 103 Ma).

👉 rotate ``input.nc`` at 100 Ma from the Alfonso2024 mantle frame (anchor 0) to the Alfonso2024 pmag frame (anchor 701701)

.. code:: console

   $ gplately rotate_grid input.nc output.nc --from-model Alfonso2024 --to-model Alfonso2024 --from-anchor 0 --to-anchor 701701 --time 100

👉 rotate all .nc files in a directory; reconstruction times are deduced from filenames such as ``paleobathymetry_103Ma.nc``

.. code:: console

   $ gplately rotate_grid input_dir output_dir --from-model Alfonso2024 --to-model Alfonso2024 --from-anchor 0 --to-anchor 701701

👉 rotate using local rotation files instead of a named model

.. code:: console

   $ gplately rotate_grid input.nc output.nc --from-rotation-files from.rot --to-rotation-files to.rot --time 100

