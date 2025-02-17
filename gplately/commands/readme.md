## command-line interface (CLI)

GPlately comes with a suite of useful command line tools. These tools are designed as GPlately subcommands. Run `gplately -h` to show the list of tools.

### **list**

  Display a list of available plate models from GPlates server. These model names can then be used by the Plate Model Manager to download model files over the Internet. Run `gplately list -h` for details.

  Examples:

    - `gplately list`
      (list all available plate models from GPlates server)

    - `gplately list -m merdith2021`
      (show details about model merdith2021)

    If you are using GPlately Docker image

    - `docker run gplates/gplately gplately list`
    - `docker run gplates/gplately gplately list -m merdith2021`

### **combine**

  Combine multiple feature collections into one. Run `gplately combine -h` for details.

### **filter**

  Filter feature collection by various criteria. Run `gplately filter -h` for details.

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

    See https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_feature_filter.sh for more examples. 

### **reset_feature_type**

  Reset the feature type for the selected features. Run `gplately reset_feature_type -h` for details.

  Examples: 

    - `gplately reset_feature_type -s gpml:ClosedContinentalBoundary -t gpml:UnclassifiedFeature input_file output_file`
        (change all gpml:ClosedContinentalBoundary to gpml:UnclassifiedFeature)
        
    - `gplately reset_feature_type -s "gpml:ContinentalFragment|gpml:Coastline" -t gpml:UnclassifiedFeature input_file output_file`
        (change all gpml:ContinentalFragment and gpml:Coastline to gpml:UnclassifiedFeature)
        
    - `gplately reset_feature_type -s ".*" -t gpml:UnclassifiedFeature input_file output_file` 
        (change all feature types to gpml:UnclassifiedFeature)     

    If you are using Docker, prefix `docker run gplates/gplately ` to the command, such as `docker run gplates/gplately gplately reset_feature_type -s ".*" -t gpml:UnclassifiedFeature input_file output_file`.

  See https://github.com/GPlates/gplately/blob/master/tests-dir/unittest/test_reset_feature_type.sh for more examples. 

### **agegrid (ag)**

  Create age grids for a plate model. Run `gplately agegrid -h` for details.

### **fix_crossovers**

  Loads one or more input rotation files, fixes any crossovers and saves the rotations to output rotation files. Run `gplately fix_crossovers -h` for details.

### **remove_rotations**

  Remove one or more plate IDs from a rotation model (consisting of one or more rotation files). Run `gplately remove_rotations -h` for details.

### **cleanup_topologies**

  Remove any regular features not referenced by topological features. Run `gplately cleanup_topologies -h` for details.

### **convert_xy_to_gplates**

  Converts geometry in one or more input ascii files (such as '.xy' files) to output files suitable for loading into GPlates. Run `gplately convert_xy_to_gplates -h` for details.

### **diagnose_rotations**

  Diagnose one or more rotation files to check for inconsistencies. Run `gplately diagnose_rotations -h` for details.

### **resolve_topologies**

  Resolve topological plate polygons (and deforming networks) and saves (to separate files) the resolved topologies, and their boundary sections as subduction zones, mid-ocean ridges (ridge/transform) and others (not subduction zones or mid-ocean ridges). Run `gplately resolve_topologies -h` for details.

### **rotation_tools**

  Calculate stage rotations between consecutive finite rotations in plate pairs. Run `gplately rotation_tools -h` for details.

### **separate_ridge_transform_segments**

  Split the geometries of isochrons and mid-ocean ridges into ridge and transform segments. Run `gplately separate_ridge_transform_segments -h` for details.

### **subduction_convergence**

  Find the convergence rates along trenches (subduction zones) over time. Run `gplately subduction_convergence -h` for details.

### **gpmdb**

  Retrieve paleomagnetic data from https://www.gpmdb.net, create GPlates-compatible VGP features and save the VGP features in a .gpmlz file. Run `gplately gpmdb -h` for details.
