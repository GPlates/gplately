#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
import re as _re

import numpy as _np
from numpy import size as _size

# ============ LINKS FOR DOWNLOADING PLATE RECONSTRUCTION FILES ON-THE-FLY ==============
""" This auxiliary script contains links to download plate model data from web servers
such as EarthByte's webDAV server. These links are accessed by GPlately's DataServer object
and downloaded into the "gplately" folder in your machine's cache via Pooch.
"""


def _find_needed_collection(collection_identifier, data_dictionary, time=None):
    """Search through link databases and time arrays (or single integer
    time) for requested download links."""

    all_urls = []
    for collection, url in data_dictionary.items():
        if collection_identifier.lower() == collection.lower():

            if time is not None:

                valid_times_dict = DataCollection(
                    collection_identifier
                ).plate_model_valid_reconstruction_times()

                for collection, valid_times in list(valid_times_dict.items()):
                    if collection_identifier.lower() == collection.lower():
                        min_time = valid_times[0]
                        max_time = valid_times[1]

                if isinstance(time, list) or isinstance(time, _np.ndarray):

                    for t in time:
                        if t < min_time or t > max_time:
                            all_urls.append(None)
                        else:
                            url_current_time = url[0].format(str(int(t)))
                            all_urls.append(url_current_time)
                    return all_urls

                else:
                    if time < min_time or time > max_time:
                        all_urls.append(None)
                    else:
                        all_urls.append(url[0].format(str(int(time))))
                        return all_urls
            else:
                return url


def _studyname(study_name):
    """Locate the citation surname for a particular identifier string
    (i.e. "Muller" for "Muller2019")."""
    name = _re.findall(r"[A-Za-z]+|\d+", study_name)[0]
    return name


def _rasters():

    database = {
        "ETOPO1_grd": [
            "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"
        ],
        "ETOPO1_tif": [
            "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/image/color_etopo1_ice_low.tif.gz"
        ],
    }
    return database


def _feature_data():
    """Assorted feature data from EarthByte's GPlates 2.3 sample data repository."""

    database = {
        "SeafloorFabric": [
            "https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/SeafloorFabric.zip"
        ],
        "Johansson2018": [
            "https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/IgneousProvinces.zip"
        ],
        "Whittaker2015": [
            "https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/IgneousProvinces.zip"
        ],
        "Hotspots": [
            "https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/Hotspots.zip"
        ],
    }
    return database


class DataCollection(object):
    """GPlately's collection of plate model data is a dictionary where
    the plate model's identifier string is the key, and values are
    lists containing any relevant file download links."""

    def __init__(self, file_collection):
        """Uses a string to identify the needed plate model, taken from
        <gplately.data.DataServer>."""
        # Allow strings with capitalisation anywhere.
        database = [
            model_name.lower() for model_name in self.plate_reconstruction_files()
        ]
        if file_collection.lower() not in database:
            raise ValueError(
                "Enter a valid plate model identifier, e.g. Muller2019, Seton2012, etc."
            )

        self.file_collection = file_collection.capitalize()

    def netcdf4_age_grids(self, time):

        age_grid_links = {
            "Muller2019": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-{}.nc"
            ],
            "Muller2016": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Agegrids/Muller_etal_2016_AREPS_Agegrids_v1.17/Muller_etal_2016_AREPS_v1.17_netCDF/Muller_etal_2016_AREPS_v1.17_AgeGrid-{}.nc"
            ],
            "Seton2012": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR/Seton_etal_2012_ESR_Agegrids/netCDF_0-200Ma/agegrid_{}.nc"
            ],
            "Clennett2020": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennet_AgeGrids_0.1d_masked/seafloor_age_mask_{}.0Ma.nc"
            ],
        }

        links_to_download = _find_needed_collection(
            self.file_collection, age_grid_links, time
        )

        return links_to_download

    def netcdf4_spreading_rate_grids(self, time):

        spread_grid_links = {
            "Clennett2020": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_SpreadRate_Grids/rategrid_final_mask_{}.nc"
            ]
        }

        links_to_download = _find_needed_collection(
            self.file_collection, spread_grid_links, time
        )

        return links_to_download

    def plate_reconstruction_files(self):

        database = {
            "Cao2020": [
                "https://zenodo.org/record/3854549/files/1000Myr_synthetic_tectonic_reconstructions.zip"
            ],
            "Muller2019": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics_Updated.zip"
            ],
            "Muller2016": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17.zip"
            ],
            "Clennett2020": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Global_Model_WD_Internal_Release_2019_v2_Clennett_NE_Pacific.zip"
            ],
            "Seton2012": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR.zip"
            ],
            # "Merdith2021" : ["https://zenodo.org/record/4485738/files/SM2_4485738_V2.zip"],
            "Merdith2021": [
                "https://earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2021_ESR/SM2-Merdith_et_al_1_Ga_reconstruction_v1.1.zip"
            ],
            "Matthews2016": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Matthews_etal_2016_Global_Plate_Model_GPC.zip"
            ],
            "Merdith2017": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2017_GR.zip"
            ],
            "Li2008": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Li_etal_2008_RodiniaModel.zip"
            ],
            "Pehrsson2015": [
                "https://www.geolsoc.org.uk/~/media/Files/GSL/shared/Sup_pubs/2015/18822_7.zip"
            ],
            "TorsvikCocks2017": [
                "http://www.earthdynamics.org/earthhistory/bookdata/CEED6.zip"
            ],
            "Young2019": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Young_etal_2018_GeoscienceFrontiers/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel.zip"
            ],
            "Scotese2008": [
                "https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip"
            ],
            "Golonka2007": [
                "https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip"
            ],
            "Clennett2020_M2019": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_M2019.zip"
            ],
            "Clennett2020_S2013": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_S2013.zip"
            ],
            "Muller2008": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller2008/Global_Model_Rigid_Internal_Release_2010.zip"
            ],
            "Scotese2016": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Scotese2016/PALEOMAP_GlobalPlateModel.zip"
            ],
            "Shephard2013": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Shephard_etal_2013_ESR.zip"
            ],
            "Muller2022": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2022_SE/Muller_etal_2022_SE_1Ga_Opt_PlateMotionModel_v1.1.zip"
            ],
            "Cao2023": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Cao_etal_2023/1.8Ga_model_submit.zip"
            ],
            "Cao2023_Opt": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Cao_etal_2023_Opt/Cao1800Opt.zip"
            ],
        }

        return database

    def plate_model_valid_reconstruction_times(self):

        database = {
            "Cao2020": [0, 1000],
            "Muller2019": [0, 250],
            "Muller2016": [0, 230],
            "Clennett2020": [0, 170],
            "Seton2012": [0, 200],
            "Merdith2021": [0, 1000],
            "Matthews2016": [0, 410],
            "Merdith2017": [0, 410],
            "Li2008": [0, 410],
            # "Pehrsson2015" : [25], (First implement continuous rotation)
            "TorsvikCocks2017": [0, 410],
            "Young2019": [0, 410],
            "Scotese2008": [0, 410],
            "Golonka2007": [0, 410],
            "Clennett2020_M2019": [0, 170],
            "Clennett2020_S2013": [0, 170],
            "Scotese2016": [0, 410],
            "Shephard2013": [0, 200],
            "Muller2008": [0, 141],  # GPlates static polygons reconstruct to this time
            "Muller2022": [0, 1000],
            "Cao2023": [0, 1800],
            "Cao2023_Opt": [0, 1800],
        }
        return database

    def rotation_strings_to_include(self):

        strings = [
            "Muller2022 1000_0_rotfile_Merdith_et_al_optimised.rot",  # For Muller et al. 2022
        ]
        return strings

    def rotation_strings_to_ignore(self):

        strings = ["OLD", "__MACOSX", "DO_NOT", "Blocks_crossing_Poles"]
        return strings

    def dynamic_polygon_strings_to_include(self):

        strings = [
            "plate_boundaries",
            "PlateBoundaries",
            "Transform",
            "Divergence",
            "Convergence",
            "Topologies",
            "Topology",
            "_PP_",  # for Seton 2012
            # "ContinentOceanBoundaries",
            # "Seton_etal_ESR2012_Coastline_2012",
            "Deforming_Mesh",
            "Deforming",
            "Flat_Slabs",
            "Feature_Geometries",
            "boundaries",
            "Clennett_etal_2020_Plates",  # For Clennett 2020 (M2019)
            "Clennett_2020_Plates",  # For topologies in Clennett et al 2020 (Pacific)
            "Clennett_2020_Terranes",  # For topologies in Clennett et al 2020 (Pacific)
            "Angayucham",
            "Farallon",
            "Guerrero",
            "Insular",
            "Intermontane",
            "Kula",
            "North_America",
            "South_America",
            "Western_Jurassic",
            "Clennett_2020_Isochrons",
            "Clennett_2020_Coastlines",
            "Clennett_2020_NAm_boundaries",
            "Shephard_etal_ESR2013_Global_EarthByte_2013",  # For Shephard et al. 2013
            "1800-1000Ma-plate-boundary_new_valid_time_and_subduction_polarity.gpml",  # for Cao2023
        ]
        return strings

    def dynamic_polygon_strings_to_ignore(self):

        strings = [
            "OLD",
            "__MACOSX",
            "DO_NOT",
            "9_Point",  # Muller et al 2019
            "9_Point_Density",  # Clennett et al 2020
            "Density",  # Clennett et al 2020
            "Inactive_Meshes_and_Topologies",  # Clennett et al 2020
            "ContinentOceanBoundaries",  # Seton 2012
            "Seton_etal_ESR2012_Coastline_2012",  # Seton 2012
            "PALEOMAP_PoliticalBoundaries",  # Scotese 2016
            "SimplifiedFiles",  # Muller et al. 2019 (updated)
            "1000-410_poles",  # Merdith
        ]
        return strings

    def static_polygon_strings_to_include(self):

        strings = [
            "StaticPolygon",
            "StaticPolygons",
            "Static_Polygon",
            "StaticPlatePolygons_",
            "RodiniaBlocks_WithPlateIDColumnAndIDs",
            # "PlatePolygons.shp",
            "CEED6_TERRANES.shp",
            "CEED6_MICROCONTINENTS.shp",
            "CEED6_LAND.gpml",
            "Scotese_2008_PresentDay_ContinentalPolygons",  # Scotese 2008
            "Golonka_2007_PresentDay_ContinentalPolygons.shp",  # Golonka 2007
            "PALEOMAP_PlatePolygons.gpml",  # For Scotese 2016
        ]
        return strings

    def static_polygon_strings_to_ignore(self):

        strings = [
            "DO_NOT",
            "OLD",
            "__MACOSX",
            "Global_Model_WD_Internal_Release_2019_v2_Clennett_NE_Pacific/StaticGeometries/StaticPolygons/Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons.shp",  # Clennett 2020
        ]
        return strings

    def topology_geometries(self):

        database = {
            "Cao2020": [
                "https://zenodo.org/record/3854549/files/1000Myr_synthetic_tectonic_reconstructions.zip"
            ],
            "Muller2019": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics_Updated.zip"
            ],
            "Muller2016": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17.zip"
            ],
            "Clennett2020": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Global_Model_WD_Internal_Release_2019_v2_Clennett_NE_Pacific.zip"
            ],
            "Seton2012": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR.zip"
            ],
            # "Merdith2021" : ["https://zenodo.org/record/4485738/files/SM2_4485738_V2.zip"],
            "Merdith2021": [
                "https://earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2021_ESR/SM2-Merdith_et_al_1_Ga_reconstruction_v1.1.zip"
            ],
            "Matthews2016": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Matthews_etal_2016_Global_Plate_Model_GPC.zip"
            ],
            "Merdith2017": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2017_GR.zip"
            ],
            "Li2008": [None],
            "Pehrsson2015": [
                "https://www.geolsoc.org.uk/~/media/Files/GSL/shared/Sup_pubs/2015/18822_7.zip"
            ],
            "TorsvikCocks2017": [
                "http://www.earthdynamics.org/earthhistory/bookdata/CEED6.zip"
            ],
            "Young2019": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Young_etal_2018_GeoscienceFrontiers/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel.zip"
            ],
            "Scotese2008": [
                "https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip"
            ],
            "Golonka2007": [
                "https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip"
            ],
            "Clennett2020_M2019": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_M2019.zip"
            ],
            "Clennett2020_S2013": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_S2013.zip"
            ],
            "Muller2008": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller2008/Global_Model_Rigid_Internal_Release_2010.zip"
            ],
            "Scotese2016": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Scotese2016/PALEOMAP_GlobalPlateModel.zip"
            ],
            "Shephard2013": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Shephard_etal_2013_ESR.zip"
            ],
            "Muller2022": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2022_SE/Muller_etal_2022_SE_1Ga_Opt_PlateMotionModel_v1.1.zip"
            ],
            "Cao2023": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Cao_etal_2023/1.8Ga_model_submit.zip"
            ],
            "Cao2023_Opt": [
                "https://www.earthbyte.org/webdav/ftp/Data_Collections/Cao_etal_2023_Opt/Cao1800Opt.zip"
            ],
        }
        return database

    def coastline_strings_to_include(self):

        strings = [
            "coastline",
            "CEED6_LAND.gpml",  # for TorsvikCocks2017
            "PALEOMAP_PoliticalBoundaries",  # For Scotese 2016
            "coast",  # for Cao2023
        ]
        return strings

    def coastline_strings_to_ignore(self):

        strings = [
            "DO_NOT",
            "OLD",
            "__MACOSX",
            "Clennett_2020_Coastlines",  # Clennett et al. 2020
            "COB_polygons_and_coastlines_combined_1000_0_Merdith_etal",  # Muller et al. 2022
        ]
        return strings

    def continent_strings_to_include(self):

        strings = [
            "continent",
            "COBfile_1000_0_Toy_introversion",
            "continental",
            "Scotese_2008_PresentDay_ContinentalPolygons.shp",  # Scotese 2008
            # "Terrane",
        ]
        return strings

    def continent_strings_to_ignore(self):

        strings = [
            "DO_NOT",
            "OLD",
            "__MACOSX",
            "Continent-ocean_boundaries",
            "COB",
        ]
        return strings

    def COB_strings_to_include(self):

        strings = [
            "cob",
            "ContinentOceanBoundaries",
            "COBLineSegments",
        ]
        return strings

    def COB_strings_to_ignore(self):

        strings = [
            "DO_NOT",
            "OLD",
            "__MACOSX",
        ]
        return strings
