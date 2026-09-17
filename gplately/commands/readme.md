## command-line interface (CLI)

GPlately comes with a suite of useful command line tools. These tools are designed as GPlately subcommands. Run `gplately -h` to show the list of tools in a terminal window.

📌 Subcommand names use hyphens (e.g. `reset-feature-type`), and most also have a short alias (e.g. `rft`) -- run `gplately <subcommand> -h` to see a subcommand's short alias. The old underscore-separated names (e.g. `reset_feature_type`) still work too, kept for backward compatibility.

- [__list__](#-list) -- show all available reconstruction models
- [__combine__](#-combine) -- combine feature collection files
- [__filter__](#-filter) -- filter feature collection by various criteria 
- [__reset-feature-type (rft)__](#-reset-feature-type-rft) -- change feature type
- [__agegrid (ag)__](#-agegrid-ag) -- generate age grids
- [__fix-crossovers (fc)__](#-fix-crossovers-fc) -- fix crossovers
- [__remove-rotations (rr)__](#-remove-rotations-rr) -- remove rotations by plate ID
- [__cleanup-topologies (ct)__](#-cleanup-topologies-ct) -- remove unreferenced features
- [__convert-xy-to-gplates (cxg)__](#-convert-xy-to-gplates-cxg) -- convert .xy files to a GPlates-compatible file
- [__diagnose-rotations (dr)__](#-diagnose-rotations-dr) -- check rotation files for inconsistencies
- [__resolve-topologies (rt)__](#-resolve-topologies-rt) -- resolve topologies at given times
- [__rotation-tools (rots)__](#-rotation-tools-rots) -- calculate stage rotations  
- [__separate-ridge-transform-segments (srts)__](#-separate-ridge-transform-segments-srts) -- pick out ridge and transform features
- [__subduction-convergence (sc)__](#-subduction-convergence-sc) -- calculate the convergence rates along subduction zones
- [__gpmdb__](#-gpmdb) -- download the paleomagnetic data and create GPlates-compatible VGP features
- [__rotate-grid (rtg)__](#-rotate-grid-rtg) -- rotate a grid between plate-model reference frames

### 🟢 **list**

  Show a list of available plate reconstruction models. Run `gplately list -h` to see the details of this subcommand.

  Examples:

  - `gplately list`
    (list all available plate reconstruction models)

  - `gplately list -m merdith2021`
    (show details about model merdith2021)

  If you are using GPlately Docker image

  - `docker run gplates/gplately gplately list`
  - `docker run gplates/gplately gplately list -m merdith2021`

### 🟢 **combine**

  Combine multiple feature collections into one. Run `gplately combine -h` to see the details of this subcommand.

  Example:

  - `gplately combine input_file_1.shp input_file_2.gpmlz input_file_3.gpml output_file.gpmlz`
    (combine three feature collection files and save to the "output_file.gpmlz")

### 🟢 **filter**

  Filter feature collection by various criteria. Run `gplately filter -h` to see the details of this subcommand.

  Examples: 

  - `gplately filter input_file output_file -n Africa "North America"`
    (get features whose name contains "Africa" or "North America")

  - `gplately filter input_file output_file -p 701 714 715 101`
    (get features whose plate ID is one of 701 714 715 101)
    
  - `gplately filter input_file output_file --min-birth-age 500`
    (get features whose birth age is older than 500Myr)
    
  - `gplately filter input_file output_file --max-birth-age 500`
    (get features whose birth age is younger than 500Myr)
    
  - `gplately filter input_file output_file -n Africa "North America" -p 701 714 715 101 --min-birth-age 500`
    (get features whose name contains "Africa" or "North America" and plate ID is one of 701 714 715 101 and birth age is older than 500Myr)
    
  - `gplately filter input_file output_file -t gpml:Basin`
    (get all gpml:Basin features)
    
  - `gplately filter input_file output_file -t "gpml:IslandArc|gpml:Basin"`
    (get all gpml:Basin and gpml:IslandArc features)

  If you are using Docker, prefix `docker run gplates/gplately ` to the command, such as `docker run gplates/gplately gplately filter input_file output_file -t gpml:Basin`.

  Check out [this shell script](https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_feature_filter.sh) for more `gplately filter` examples. 

### 🟢 **reset-feature-type (rft)**

  Reset the feature type for the selected features. Also available as `reset_feature_type`. Run `gplately reset-feature-type -h` to see the details of this subcommand.

  Examples: 

  - `gplately reset-feature-type -s gpml:ClosedContinentalBoundary -t gpml:UnclassifiedFeature input_file output_file`
    (change all gpml:ClosedContinentalBoundary to gpml:UnclassifiedFeature)
        
  - `gplately reset-feature-type -s "gpml:ContinentalFragment|gpml:Coastline" -t gpml:UnclassifiedFeature input_file output_file`
    (change all gpml:ContinentalFragment and gpml:Coastline to gpml:UnclassifiedFeature)
        
  - `gplately reset-feature-type -s ".*" -t gpml:UnclassifiedFeature input_file output_file` 
    (change all feature types to gpml:UnclassifiedFeature)     

  If you are using Docker, prefix `docker run gplates/gplately ` to the command, such as `docker run gplates/gplately gplately reset-feature-type -s ".*" -t gpml:UnclassifiedFeature input_file output_file`.

  Check out [this shell script](https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_reset_feature_type.sh) for more `gplately reset-feature-type` examples. 

### 🟢 **agegrid (ag)**

  Generate age grids for a plate reconstruction model. Run `gplately agegrid -h` to see the details of this subcommand.

  Examples:

  - `gplately ag output -m merdith2021 -e 0 -s 10`
    (create age grids from 10Ma to 0Ma with 1Myr increment for the merdith2021 reconstruction mode)

  - `gplately ag rotations.rot topologies.gpmlz output -c continental_polygons.gpmlz -e 0 -s 10`
    (create age grids from 10Ma to 0Ma with 1Myr increment using the specified reconstruction files)

### 🟢 **fix-crossovers (fc)**

  Fixes crossovers in rotation file(s). Also available as `fix_crossovers`. Run `gplately fix-crossovers -h` to see the details of this subcommand.

  Example:

  - `gplately fix-crossovers -d -c 0.01 -i 201 701 -- input_rotations1.rot input_rotations2.rot`
    (fix crossovers in two rotation files with a threshold 0.01 degree and ignore plate ID 201 and 701)

### 🟢 **remove-rotations (rr)**

  Remove one or more plate IDs from a rotation model (consisting of one or more rotation files). Also available as `remove_rotations`. Run `gplately remove-rotations -h` to see the details of this subcommand.

  Example:

  - `gplately remove-rotations -p 70 4 3 1 -o removed_ref_frames_ -- rotations.rot`
    (remove plate IDs 70,4,3 and 1 from a rotation file)

### 🟢 **cleanup-topologies (ct)**

  Remove any regular features not referenced by topological features. Also available as `cleanup_topologies`. Run `gplately cleanup-topologies -h` to see the details of this subcommand.

  Example:

  - `gplately cleanup-topologies -o cleanup-topologies- -- topologies.gpml`
    (remove all features which are not referenced by any topological feature from topologies.gpml)

### 🟢 **convert-xy-to-gplates (cxg)**

  Converts geometry in one or more input ascii files (such as '.xy' files) to output files suitable for loading into GPlates. Also available as `convert_xy_to_gplates`. Run `gplately convert-xy-to-gplates -h` to see the details of this subcommand.

  Example:

  - `gplately convert-xy-to-gplates -e shp -- input1.xy input2.xy`
    (convert two .xy file into a shapefile)

### 🟢 **diagnose-rotations (dr)**

  Diagnose one or more rotation files to check for inconsistencies. Also available as `diagnose_rotations`. Run `gplately diagnose-rotations -h` to see the details of this subcommand.

  Example:

  - `gplately diagnose-rotations input_rotations1.rot input_rotations2.rot`
    (check two rotation files and print the diagnostic results on screen)

### 🟢 **resolve-topologies (rt)**

  Resolve topological plate polygons (and deforming networks) and saves (to separate files) the resolved topologies, and their boundary sections as subduction zones, mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges). Also available as `resolve_topologies`. Run `gplately resolve-topologies -h` to see the details of this subcommand.

  Example:

  - `gplately resolve-topologies -r rotations1.rot rotations2.rot -m topologies1.gpml topologies2.gpml -t 10`
    (resolve topologies at 10Ma)


### 🟢 **rotation-tools (rots)**

  Calculate stage rotations between consecutive finite rotations in plate pairs. Also available as `rotation_tools`. Run `gplately rotation-tools -h` to see the details of this subcommand.

  Example:

  - `gplately rotation-tools -p 701 0 -o stage_ -- rotations.rot`
    (calculate stage rotations for moving plate 701 relative to the fixed plate 0)

### 🟢 **separate-ridge-transform-segments (srts)**

  Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments. Also available as `separate_ridge_transform_segments`. Run `gplately separate-ridge-transform-segments -h` to see the details of this subcommand.

  Example:

  - `gplately separate-ridge-transform-segments -r rotations.rot -d 45 -s _ridges -t _transforms -- spreading_features.gpml`
    (pick out ridge and transform features from the file spreading_features.gpml)

### 🟢 **subduction-convergence (sc)**

  Find the convergence rates along trenches (subduction zones) over time. Also available as `subduction_convergence`. Run `gplately subduction-convergence -h` to see the details of this subcommand.

  Example:

  - `gplately subduction-convergence -r rotations.rot -m topologies.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- convergence`
    (calculate the convergence rates along subduction zones from 200Ma to 0Ma)

### 🟢 **gpmdb**

  Retrieve the paleomagnetic data from https://www.gpmdb.net, create GPlates-compatible VGP features and save the VGP features in a .gpmlz file. Run `gplately gpmdb -h` to see the details of this subcommand.

  Example:

  - `gplately gpmdb -m zahirovic2022 -o vgp.gpmlz`
    (download the paleomagnetic data and generate GPlates-compatible VGP features using the zahirovic2022 reconstruction model)


### 🟢 **rotate-grid (rtg)**

  Rotate a grid (or all grids in a folder) between plate-model reference frames. Also available as `rotate_grid`. Run `gplately rotate-grid -h` to see the details of this subcommand.

  Examples:

  - `gplately rotate-grid input.nc output.nc --from-model Alfonso2024 --to-model Alfonso2024 --from-anchor 0 --to-anchor 701701 --time 100`
    (rotate input.nc at 100 Ma from the Alfonso2024 mantle frame to the Alfonso2024 pmag frame)

  - `gplately rotate-grid input_dir output_dir --from-model Alfonso2024 --to-model Alfonso2024 --from-anchor 0 --to-anchor 701701`
    (rotate all .nc files in input_dir; reconstruction times are deduced from filenames such as paleobathymetry_103Ma.nc)

  - `gplately rotate-grid input.nc output.nc --from-rotation-files from.rot --to-rotation-files to.rot --time 100`
    (rotate using local rotation files instead of a named model)
