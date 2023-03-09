Readme

Muller et al. (2019), Journal: Tectonics
A Global Plate Model Including Lithospheric Deformation Along Major Rifts and Orogens Since the Triassic

Last updated: 20 December 2019 (Licensing information updated on 16 January 2020)
v 1.0 - Submitted set of files
v 1.1 - Resubmitted set of files, with fixes to Mediterranean 
v 1.2 - Final submitted files with tidied absolute reference frame 
v 1.3 - Fixes to Arctic interpolation issue, Andes subduction polarity, plate boundaries at SAM-AFR-NAM boundary, and other minor issues (internal r1448)
v 1.4 - Fixing topology break along Andes in the last 15 Myr
v 2.0 - Major update of model, with new absolute reference frame, synchronisation with Young et al. relative plate motions, and improvements to interpolation of stretching factors (see release notes for details). SVN r1630. 

The latest version of this supplementary dataset can be downloaded from:
https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/

This directory contains a number of files and folders:
* GPML files in root folder contain the geometries (lines, topological lines, points, polygons) that participate in the rigid and deforming topologies
* ROT files in root folder contain the Euler rotations, which must be in a merged rotation layer in GPlates to function properly
* DeformingMeshPoints: Folder containing low- (7 level mesh) and high- (9 level mesh) resolution of deforming points. 
* ProjectFiles
	- Muller_etal_2019_7_Point_Density_Project_File_v2.0.gproj: The default global project file to use for quick visualisations
	- Muller_etal_2019_9_Point_Density_Project_File_v2.0.gproj: High resolution deforming points used to make the figures (Note: These require a lot of memory to load, and can take about 30 minutes to open in GPlates.)
	- Muller_etal_2019_v2.0.gproj: Plate motion model with deforming networks but no crustal deformation points (Note: quickest option to load)
* StaticGeometries
	- AgeGridInput: Isochrons, ridges, and COBs for generating (paleo-) seafloor age-grids 
	- Coastlines: Present-day coastlines cookie-cut using the Static Polygons
	- COBLineSegments: COB line segments not used in the age-gridding
	- ContinentalPolygons: Cookie-cutting polygons just for the non-oceanic regions (continents,  intra-oceanic arcs, etc.)
	- ContinentalStencils: 
	- StaticPolygons: Cookie-cutting global polygons for reconstruction raster and vector data

To load the files in GPlates do the following:
1.  Open GPlates
2.  Pull down the GPlates File menu and select the operation Open Project
3.  Point to the 7_Point_Density_Project_File.gproj file in the ProjectFiles folder
4.  Click Open 

Alternatively, drag and drop the project file onto the globe. 
 
** Release notes for model v2.0 ** 

Relative and absolute plate motions
* Full synchronisation of east Asia reconstructions of Young et al. (2018) from 250 to  130 Ma, including the closure of the Mongol-Okhotsk Ocean
* Improved kinematics of the back-arc basin history along western North America (250 to 180 Ma) so that the seafloor spreading isochrons are consistent with the implied plate motions
* Pacific absolute plate motion fixes following Torsvik et al. (2019)
* Absolute optimised frame, with simplified hierarchy at 5 Myr intervals with optimisation weighting parameters of
    - 0-80 Ma:  TR=1, NR=1, HS=1;
    - 80-170 Ma: TR=1, NR=0.5;
    - 170-250 Ma: TR=1, NR=0.2
	Where TR=trench migration optimisation, NR=net rotation optimisation and HS=hotspot track fitting, in an optimisation framework as described in Tetley et al. (2019), with numbers referring to relative weights of the optimisation parameters. These parameters were iteratively determined to provide the most optimal absolute plate motion model, given the relative plate motion revisions listed above.

Interpolation of deformation grids
* Using nearneighbor rather than spherical interpolation (removing areas that have had no deformation)
* Applying rigid blocks within deformed regions  (e.g. Tarim craton, Khorat platform, etc.), as well as oceanic age-grid, as masks to remove anomalous areas of extrapolation
* Previous versions had a tiny amount of negative values in stretching factors, but these grids are now clipped to values above zero

References

Müller, R. D., Zahirovic, S., Williams, S. E., Cannon, J., Seton, M., Bower, D. J., Tetley, M. G., Heine, C., Le Breton, E., Liu, S., Russell, S. H. J., Yang, T., Leonard, J., and Gurnis, M., 2019, A global plate model including lithospheric deformation along major rifts and orogens since the Triassic: Tectonics, v. 38, no. Fifty Years of Plate Tectonics: Then, Now, and Beyond.

Tetley, M.G., Williams, S.E., Gurnis, M., Flament, N. and Müller, R.D., 2019. Constraining absolute plate motions since the Triassic. Journal of Geophysical Research: Solid Earth, 124, 7231-7258.

Torsvik, T. H., Steinberger, B., Shephard, G. E., Doubrovine, P. V., Gaina, C., Domeier, M., Conrad, C. P., and Sager, W. W., 2019, Pacific‐Panthalassic reconstructions: Overview, errata and the way forward: Geochemistry, Geophysics, Geosystems, v. 20, no. 7, p. 3659-3689.

Young, A., Flament, N., Maloney, K., Williams, S., Matthews, K., Zahirovic, S., and Müller, R. D., 2018, Global kinematics of tectonic plates and subduction zones since the late Paleozoic Era: Geoscience Frontiers.


Play around with the GPlates buttons to make an animation, select features, draw features, etc.  
For more information, read the GPlates manual which can be downloaded from www.gplates.org or http://www.earthbyte.org/Resources/earthbyte_gplates.html

This plate reconstruction model is released under the Creative Commons "Attribution-ShareAlike 4.0 International" license. This allows free use for everyone (including researchers, members of the public, commercial use, etc.) with attribution. Modification of the plate reconstruction model and associated files in this folder is allowed, but any modifications must retain the same license (i.e., will need to be made free to use under the same terms). More information on the license can be found here: https://creativecommons.org/licenses/by-sa/4.0/legalcode

Any questions, please email: dietmar.muller@sydney.edu.au, sabin.zahirovic@sydney.edu.au or maria.seton@sydney.edu.au


