
# ============ LINKS FOR DOWNLOADING PLATE RECONSTRUCTION FILES ON-THE-FLY ==============
""" This auxiliary script contains links to download plate model data from web servers
such as EarthByte's webDAV server. These links are accessed by GPlately's DataServer object
and downloaded into the "gplately" folder in your machine's cache via Pooch.
"""

def _find_needed_collection(
    collection_identifier,
    data_dictionary,
    times=None
    ):
    """Search through link databases and time arrays for requested 
    download links."""

    all_urls = []
    for collection, url in data_dictionary.items():
        if collection_identifier.lower() == collection.lower():
            if times is not None:
                for t in times:
                    url_current_time = url[0].format(str(t))
                    all_urls.append(url_current_time)
                return all_urls
            else:
                return url


class DataCollection(object):
    """GPlately's collection of plate model data is a dictionary where
    the plate model's identifier string is the key, while values are 
    lists containing any relevant file download links."""

    def __init__(self, file_collection):
        """Uses a string to identify the needed plate model, taken from
        <gplately.data.DataServer>."""
        self.file_collection = file_collection


    def netcdf4_age_grids(self, times):

        age_grid_links = {

            "Muller2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Agegrids/Muller_etal_2016_AREPS_Agegrids_v1.17/Muller_etal_2016_AREPS_v1.17_netCDF/Muller_etal_2016_AREPS_v1.17_AgeGrid-{}.nc"]
        }

        links_to_download = _find_needed_collection(
            self.file_collection, 
            age_grid_links,
            times)

        return links_to_download


    def jpeg_agegrids(self, times):

        age_grid_links = {

            "Muller2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_jpgs.zip"],
            "Muller2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Agegrids/Muller_etal_2016_AREPS_Agegrids_v1.17/Muller_etal_2016_AREPS_v1.17_jpgs/Muller_etal_2016_AREPS_v1.17_AgeGrid-{}.jpg"],
        }

        links_to_download = _find_needed_collection(
            self.file_collection, 
            age_grid_links,
            times)

        return links_to_download


    def png_agegrids(self, times):

        age_grid_links = {

            "Muller2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Agegrids/Muller_etal_2016_AREPS_Agegrids_v1.17/Muller_etal_2016_AREPS_v1.17_pngs/Muller_etal_2016_AREPS_v1.17_AgeGrid-{}.png"],
        }

        links_to_download = _find_needed_collection(
            self.file_collection, 
            age_grid_links,
            times)

        return links_to_download


    def plate_reconstruction_files(self):

        database = {

            "Cao2020" : ["https://zenodo.org/record/3854549/files/1000Myr_synthetic_tectonic_reconstructions.zip"],
            "Muller2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics.zip"], 
            "Muller2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17.zip"],
            "Mather2021" : ["https://zenodo.org/record/5769002/files/plate_model.zip"],
            "Seton2012" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR.zip"],
            "Merdith2021" : ["https://zenodo.org/record/4485738/files/SM2_4485738_V2.zip"],
            "Matthews2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Matthews_etal_2016_Global_Plate_Model_GPC.zip"], 
            "Merdith2017" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2017_GR.zip"], 
            "Li2008" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Li_etal_2008_RodiniaModel.zip"],
            "Pehrsson2015" : ["https://www.geolsoc.org.uk/~/media/Files/GSL/shared/Sup_pubs/2015/18822_7.zip"]

        }

        return database


    def rotation_strings_to_ignore(self):

        strings = [
            "OLD",
            "__MACOSX",
            "DO_NOT",
        ]    
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
        ]
        return strings 


    def dynamic_polygon_strings_to_ignore(self):

        strings = [
            "OLD",
            "__MACOSX",
            "DO_NOT",
        ]
        return strings


    def static_polygon_strings_to_include(self):

        strings = [
            "Static",
            "StaticPolygon",
            "Static_Polygon",
            "RodiniaBlocks_WithPlateIDColumnAndIDs",
            "PlatePolygons.shp",
        ]
        return strings


    def static_polygon_strings_to_ignore(self):

        strings = [

            "DO_NOT",
            "OLD",
            "__MACOSX"
        ]
        return strings


    def topology_geometries(self):

        database = {

            "Cao2020" : ["https://zenodo.org/record/3854549/files/1000Myr_synthetic_tectonic_reconstructions.zip"],
            "Muller2019" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics.zip"], 
            "Muller2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17.zip"],
            "Mather2021" : ["https://zenodo.org/record/5769002/files/plate_model.zip"],
            "Seton2012" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR.zip"],
            "Merdith2021" : ["https://zenodo.org/record/4485738/files/SM2_4485738_V2.zip"],
            "Matthews2016" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Matthews_etal_2016_Global_Plate_Model_GPC.zip"], 
            "Merdith2017" : ["https://www.earthbyte.org/webdav/ftp/Data_Collections/Merdith_etal_2017_GR.zip"],  
            "Pehrsson2015" : ["https://www.geolsoc.org.uk/~/media/Files/GSL/shared/Sup_pubs/2015/18822_7.zip"]            
        }
        return database


    def coastline_strings_to_include(self):

        strings = [

            "coastline"
        ]
        return strings


    def coastline_strings_to_ignore(self):

        strings = [

            "DO_NOT",
            "OLD",
            "__MACOSX"
        ]
        return strings


    def continent_strings_to_include(self):

        strings = [

            "continent",
            "COBfile_1000_0_Toy_introversion",
            "continental",
        ]
        return strings


    def continent_strings_to_ignore(self):

        strings = [

            "DO_NOT",
            "OLD",
            "__MACOSX"
        ]
        return strings


    def COB_strings_to_include(self):

        strings = [

            "cob",
        ]
        return strings


    def COB_strings_to_ignore(self):

        strings = [

            "DO_NOT",
            "OLD",
            "__MACOSX"
        ]
        return strings


    def rasters(self):

        database = {

            "ETOPO1_grd" : ["https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"],
            "ETOPO1_tif" : ["https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/image/color_etopo1_ice_low.tif.gz"],
        }
        return database


    def feature_data(self):
        """Assorted feature data from EarthByte's GPlates 2.3 sample data repository."""

        database = {

            "SeafloorFabric" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/SeafloorFabric.zip"],
            "LIP_VolcanicProvinces" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/IgneousProvinces.zip"],
            "Hotspots" : ["https://www.earthbyte.org/webdav/ftp/earthbyte/GPlates/GPlates2.3_GeoData/Individual/Hotspots.zip"]
        }
        return database
