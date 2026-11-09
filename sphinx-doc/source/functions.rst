Functions
=========

This page lists the assorted functions of the GPlately Python package.


.. contents::
   :local:
   :depth: 3

Paleobathymetry
-----------------
.. autosummary::
   :nosignatures:
   :toctree: generated

   gplately.simple_paleobathymetry

Step 1 -- seafloor age to basement depth:

.. autosummary::
   :nosignatures:
   :toctree: generated

   gplately.age_to_basement_depth
   gplately.AGE_DEPTH_MODELS

Step 2 -- distance to the nearest passive continental margin:

.. autosummary::
   :nosignatures:
   :toctree: generated

   gplately.generate_distance_grids
   gplately.generate_input_points_grid
   gplately.generate_passive_margins
   gplately.passive_margin_polylines

Steps 3 and 4 -- sediment thickness, and its isostatic correction:

.. autosummary::
   :nosignatures:
   :toctree: generated

   gplately.generate_sediment_thickness_grids
   gplately.dutkiewicz_2017_sediment_thickness
   gplately.DUTKIEWICZ_2017_SEDIMENT_THICKNESS
   gplately.sediment_isostatic_correction
   gplately.paleobathymetry

Step 5 -- merging in pyBacktrack's paleobathymetry:

.. autosummary::
   :nosignatures:
   :toctree: generated

   gplately.merge_pybacktrack_paleobathymetry


Reconstruction
--------------

.. autosummary::
   :nosignatures:
   :toctree: generated

   gplately.reconstruct_grid
   gplately.reconstruct_points
   gplately.reconstruct_points_with_model_files
   gplately.reverse_reconstruct_points
   gplately.reverse_reconstruct_points_with_model_files

Resolve Topologies
------------------

.. autosummary::
   :nosignatures:
   :toctree: generated

   gplately.resolve_topologies
   gplately.resolve_topological_snapshot
   gplately.resolve_topologies_into_features
   gplately.resolve_topological_snapshot_into_features

Colour Maps
-----------
.. autosummary::
   :nosignatures:
   :toctree: generated

   gplately.plot.get_age_grid_cmap
   gplately.plot.get_spreading_rate_cmap
   gplately.plot.get_topo_cmap
   

Miscellaneous
-------------

.. autosummary::
   :nosignatures:
   :toctree: generated

   gplately.auxiliary
   gplately.ridge_spreading_rate
   gplately.subduction_convergence
   gplately.load_feature_collection
   