## command-line interface (CLI)

GPlately comes with a suite of useful command line tools. These tools are designed as GPlately subcommands. Run `gplately -h` to show the list of tools in a terminal window.

- [__list__](#list) -- show all available reconstruction models
- [__combine__](#combine) -- combine feature collection files
- [__filter__](#filter) -- filter feature collection by various criteria 
- [__reset_feature_type__](#reset_feature_type) -- change feature type
- [__agegrid (ag)__](#agegrid-(ag)) -- generate age grids
- [__fix_crossovers__](#fix_crossovers) -- fix crossovers
- [__remove_rotations__](#remove_rotations) -- remove rotations by plate ID
- [__cleanup_topologies__](#cleanup_topologies) -- remove unreferenced features
- [__convert_xy_to_gplates__](#convert_xy_to_gplates) -- convert .xy files to a GPlates-compatible file
- [__diagnose_rotations__](#diagnose_rotations) -- check rotation files for inconsistencies
- [__resolve_topologies__](#resolve_topologies) -- resolve topologies at given times
- [__rotation_tools__](#rotation_tools) -- calculate stage rotations  
- [__separate_ridge_transform_segments__](#separate_ridge_transform_segments) -- pick out ridge and transform features
- [__subduction_convergence__](#subduction_convergence) -- calculate the convergence rates along subduction zones
- [__gpmdb__](#gpmdb) -- download the paleomagnetic data and create GPlates-compatible VGP features 

### **list**

  Show a list of available plate reconstruction models. Run `gplately list -h` to see the details of this subcommand.

  Examples:

  - `gplately list`
    (list all available plate reconstruction models)

  - `gplately list -m merdith2021`
    (show details about model merdith2021)

  If you are using GPlately Docker image

  - `docker run gplates/gplately gplately list`
  - `docker run gplates/gplately gplately list -m merdith2021`

### **combine**

  Combine multiple feature collections into one. Run `gplately combine -h` to see the details of this subcommand.

  Example:

  - `gplately combine input_file_1.shp input_file_2.gpmlz input_file_3.gpml output_file.gpmlz`
    (combine three feature collection files and save to the "output_file.gpmlz")

### **filter**

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

### **reset_feature_type**

  Reset the feature type for the selected features. Run `gplately reset_feature_type -h` to see the details of this subcommand.

  Examples: 

  - `gplately reset_feature_type -s gpml:ClosedContinentalBoundary -t gpml:UnclassifiedFeature input_file output_file`
    (change all gpml:ClosedContinentalBoundary to gpml:UnclassifiedFeature)
        
  - `gplately reset_feature_type -s "gpml:ContinentalFragment|gpml:Coastline" -t gpml:UnclassifiedFeature input_file output_file`
    (change all gpml:ContinentalFragment and gpml:Coastline to gpml:UnclassifiedFeature)
        
  - `gplately reset_feature_type -s ".*" -t gpml:UnclassifiedFeature input_file output_file` 
    (change all feature types to gpml:UnclassifiedFeature)     

  If you are using Docker, prefix `docker run gplates/gplately ` to the command, such as `docker run gplates/gplately gplately reset_feature_type -s ".*" -t gpml:UnclassifiedFeature input_file output_file`.

  Check out [this shell script](https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_reset_feature_type.sh) for more `gplately reset_feature_type` examples. 

### **agegrid (ag)**

  Generate age grids for a plate reconstruction model. Run `gplately agegrid -h` to see the details of this subcommand.

  Examples:

  - `gplately ag output -m merdith2021 -e 0 -s 10`
    (create age grids from 10Ma to 0Ma with 1Myr increment for the merdith2021 reconstruction mode)

  - `gplately ag rotations.rot topologies.gpmlz output -c continental_polygons.gpmlz -e 0 -s 10`
    (create age grids from 10Ma to 0Ma with 1Myr increment using the specified reconstruction files)

### **fix_crossovers**

  Fixes crossovers in rotation file(s). Run `gplately fix_crossovers -h` to see the details of this subcommand.

  Example:

  - `gplately fix_crossovers -d -c 0.01 -i 201 701 -- input_rotations1.rot input_rotations2.rot`
    (fix crossovers in two rotation files with a threshold 0.01 degree and ignore plate ID 201 and 701)

### **remove_rotations**

  Remove one or more plate IDs from a rotation model (consisting of one or more rotation files). Run `gplately remove_rotations -h` to see the details of this subcommand.

  Example:

  - `gplately remove_rotations -p 70 4 3 1 -o removed_ref_frames_ -- rotations.rot`
    (remove plate IDs 70,4,3 and 1 from a rotation file)

### **cleanup_topologies**

  Remove any regular features not referenced by topological features. Run `gplately cleanup_topologies -h` to see the details of this subcommand.

  Example:

  - `gplately cleanup_topologies -o cleanup_topologies_ -- topologies.gpml`
    (remove all features which are not referenced by any topological feature from topologies.gpml)

### **convert_xy_to_gplates**

  Converts geometry in one or more input ascii files (such as '.xy' files) to output files suitable for loading into GPlates. Run `gplately convert_xy_to_gplates -h` to see the details of this subcommand.

  Example:

  - `gplately convert_xy_to_gplates -e shp -- input1.xy input2.xy`
    (convert two .xy file into a shapefile)

### **diagnose_rotations**

  Diagnose one or more rotation files to check for inconsistencies. Run `gplately diagnose_rotations -h` to see the details of this subcommand.

  Example:

  - `gplately diagnose_rotations input_rotations1.rot input_rotations2.rot`
    (check two rotation files and print the diagnostic results on screen)

### **resolve_topologies**

  Resolve topological plate polygons (and deforming networks) and saves (to separate files) the resolved topologies, and their boundary sections as subduction zones, mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges). Run `gplately resolve_topologies -h` to see the details of this subcommand.

  Example:

  - `gplately resolve_topologies -r rotations1.rot rotations2.rot -m topologies1.gpml topologies2.gpml -t 10`
    (resolve topologies at 10Ma)


### **rotation_tools**

  Calculate stage rotations between consecutive finite rotations in plate pairs. Run `gplately rotation_tools -h` to see the details of this subcommand.

  Example:

  - `gplately rotation_tools -p 701 0 -o stage_ -- rotations.rot`
    (calculate stage rotations for moving plate 701 relative to the fixed plate 0)

### **separate_ridge_transform_segments**

  Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments. Run `gplately separate_ridge_transform_segments -h` to see the details of this subcommand.

  Example:

  - `gplately separate_ridge_transform_segments -r rotations.rot -d 45 -s _ridges -t _transforms -- spreading_features.gpml`
    (pick out ridge and transform features from the file spreading_features.gpml)

### **subduction_convergence**

  Find the convergence rates along trenches (subduction zones) over time. Run `gplately subduction_convergence -h` to see the details of this subcommand.

  Example:

  - `gplately subduction_convergence -r rotations.rot -m topologies.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- convergence`
    (calculate the convergence rates along subduction zones from 200Ma to 0Ma)

### **gpmdb**

  Retrieve the paleomagnetic data from https://www.gpmdb.net, create GPlates-compatible VGP features and save the VGP features in a .gpmlz file. Run `gplately gpmdb -h` to see the details of this subcommand.

  Example:

  - `gplately gpmdb -m zahirovic2022 -o vgp.gpmlz`
    (download the paleomagnetic data and generate GPlates-compatible VGP features using the zahirovic2022 reconstruction model)
