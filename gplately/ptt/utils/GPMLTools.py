# Copyright (C) 2014  Michael Tetley
# EarthByte Group, University of Sydney
# Contact email: michael.tetley@sydney.edu.au
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.



#### GPML Tools ####

from __future__ import print_function

import pygplates as pgp
import datetime
import time
import os
import numpy as np


# Filter GPML by selected criteria and output new GPML file of filtered data
def filterGPML(**kwargs):

    # Start the clock
    start = time.time()

    filterProperties = ["inputFile", "outputFile", "filterSequence", "rPlateID", "cPlateID", "ageAppearWindow", "ageDisappearWindow",
                        "ageExistsWindow", "boundingBox", "featureType", "geometryType", "featureID", "featureName", "feature_truncate_age", "inverse", "cascade"]

    # Process supplied arguments and assign values to variables

    # Inverse is set to False by default
    inverse = False

    # Cascade is set to True  by default
    cascade = True

    for parameter, value in kwargs.items():

        if parameter in filterProperties:
            if parameter == filterProperties[0]:
                inputFile = value
            elif parameter == filterProperties[1]:
                outputFile = value
            elif parameter == filterProperties[2]:
                filterSequence = value
            elif parameter == filterProperties[3]:
                rPlateID = value
            elif parameter == filterProperties[4]:
                cPlateID = value
            elif parameter == filterProperties[5]:
                ageAppearWindow = value
            elif parameter == filterProperties[6]:
                ageDisappearWindow = value
            elif parameter == filterProperties[7]:
                ageExistsWindow = value


                if ageExistsWindow[1] > ageExistsWindow[0]:
                    print(" ")
                    print("ERROR - Age exists window end age older than begin age: " + str(ageExistsWindow[1]))
                    

            elif parameter == filterProperties[8]:
                boundingBox = value

                if pgp.LatLonPoint.is_valid_longitude(boundingBox[0]) is False:
                    print(" ")
                    print("ERROR - Bounding box longitude is not valid: " + str(boundingBox[0]))
                    
                if pgp.LatLonPoint.is_valid_longitude(boundingBox[1]) is False:
                    print(" ")
                    print("ERROR - Bounding box longitude is not valid: " + str(boundingBox[1]))
                    
                if pgp.LatLonPoint.is_valid_latitude(boundingBox[2]) is False:
                    print(" ")
                    print("ERROR - Bounding box latitude is not valid: " + str(boundingBox[2]))
                    
                if pgp.LatLonPoint.is_valid_latitude(boundingBox[3]) is False:
                    print(" ")
                    print("ERROR - Bounding box latitude is not valid: " + str(boundingBox[3]))
                    

            elif parameter == filterProperties[9]:
                featureTemp = value
                featureType = []

                if "ALL" in featureTemp:
                    featureType = ["Isochron", "MidOceanRidge", "PassiveContinentalBoundary"]
                if "ISO" in featureTemp:
                    featureType.append("Isochron")
                if "MOR" in featureTemp:
                    featureType.append("MidOceanRidge")
                if "PCB" in featureTemp:
                    featureType.append("PassiveContinentalBoundary")

            elif parameter == filterProperties[10]:
                geometryType = value

                if "ALL" in geometryType:
                    geometryType = ["PolylineOnSphere", "PolygonOnSphere", "PointOnSphere", "MultiPointOnSphere"]

            elif parameter == filterProperties[11]:
                featureID = value
            elif parameter == filterProperties[12]:
                featureName = value
            elif parameter == filterProperties[13]:
                feature_truncate_age = value
            elif parameter == filterProperties[14]:
                inverse = value
            elif parameter == filterProperties[15]:
                cascade = value


        else:
            print(" ")
            print("ERROR - Filter criteria not found: " + str(parameter))
            print(" ")
            


    date = datetime.date.today()

    #if inputFile != "none":
        #output = pgp.FeatureCollection()

    featureCollection = pgp.FeatureCollectionFileFormatRegistry()

    print(" ")
    print("--------------------------------------------")
    print(" ### GPMLTools - filterGPML ###")

    # Check for existing output directory and create it if not found
    if not os.path.exists("output"):
        os.makedirs("output")
        print(" ")
        print("Housekeeping:")
        print("    No output folder found. Folder 'output' created.")

    # Check for existing output file with same name and remove if found
    if os.path.isfile("output/output.gpml"):
        os.remove("output/output.gpml")
        print(" ")
        print("Housekeeping:")
        print("    Previous 'output.gpml' found in destination folder. File removed for new filter sequence.")


    try:
        feature = featureCollection.read(inputFile)
        feature1 = featureCollection.read(inputFile)
        print(" ")
        print("Data handling:")
        print("    Successfully loaded data file:  '" + str(inputFile) + "'")
        print("       - File contains " + str(len(feature)) + " features.")

    except pgp.OpenFileForReadingError:
        print(" ")
        print("    ERROR - File read error in: '" + inputFile + "'. Is this a valid GPlates file?")
        
    except pgp.FileFormatNotSupportedError:
        print(" ")
        print("    ERROR - File format not supported: '" + inputFile + "'. Please check the file name and try again")
        



    #Filter data

    print(" ")
    print("Filter sequence:")

    previousFilter = 0

    f1_result = pgp.FeatureCollection()
    f2_result = pgp.FeatureCollection()
    f3_result = pgp.FeatureCollection()
    f4_result = pgp.FeatureCollection()
    f5_result = pgp.FeatureCollection()
    f6_result = pgp.FeatureCollection()
    f7_result = pgp.FeatureCollection()
    f8_result = pgp.FeatureCollection()
    f9_result = pgp.FeatureCollection()
    f10_result = pgp.FeatureCollection()
    f11a_result = pgp.FeatureCollection()
    f11b_result = pgp.FeatureCollection()


    for filter_ in filterSequence:

        if previousFilter == 0:
            data_ = feature
        elif previousFilter == 1:
            data_ = f1_result
        elif previousFilter == 2:
            data_ = f2_result
        elif previousFilter == 3:
            data_ = f3_result
        elif previousFilter == 4:
            data_ = f4_result
        elif previousFilter == 5:
            data_ = f5_result
        elif previousFilter == 6:
            data_ = f6_result
        elif previousFilter == 7:
            data_ = f7_result
        elif previousFilter == 8:
            data_ = f8_result
        elif previousFilter == 9:
            data_ = f9_result
        elif previousFilter == 10:
            data_ = f10_result
        elif previousFilter == 11:
            data_ = f11_result



        # Filter by reconstruction plate ID
        if filter_ == 1:

            for feature in data_:

                if cascade == False:

                    for property in feature:

                        filter_property = property.get_name()

                        if filter_property.get_name() == "reconstructionPlateId":
                            selected_filter_property = property.get_value()

                            for property in feature:

                                filter_property = property.get_name()

                                if filter_property.get_name() == "conjugatePlateId":
                                    second_filter_property = property.get_value()

                                    if inverse == False:
                         
                                        # Isolate criteria match and process
                                        if str(selected_filter_property) == str(rPlateID[0]) and str(second_filter_property) == str(cPlateID[0]):

                                            # Append filtered data to associated Feature Collection
                                            f1_result.add(feature)

                                    elif inverse == True:

                                        if int(str(selected_filter_property)) not in rPlateID or int(str(second_filter_property)) not in cPlateID:

                                            # Append filtered data to associated Feature Collection
                                            f1_result.add(feature)



                elif cascade == True:

                    for property in feature:

                        filter_property = property.get_name()

                        if filter_property.get_name() == "reconstructionPlateId":
                            selected_filter_property = property.get_value()

                            if inverse == False:

                                for plateID in rPlateID:

                                    if str(selected_filter_property) == str(plateID):

                                        # Append filtered data to associated Feature Collection
                                        f1_result.add(feature)

                            elif inverse == True:

                                if int(str(selected_filter_property)) not in rPlateID:
        							
        							# Append filtered data to associated Feature Collection
                                    f1_result.add(feature)

                        

            if cascade == False:
                print("Oooooh, you found the secret command...")
                print(" ")
                print("    1. Filtering data by reconstruction plate ID: " + str(rPlateID) + " and conjugate plate ID: " + str(cPlateID))
                print("       - Found " + str(len(f1_result)) + " feature(s).")

                cascade = True

            else:

                print(" ")
                print("    1. Filtering data by reconstruction plate ID(s): " + str(rPlateID))
                print("       - Found " + str(len(f1_result)) + " feature(s).")

            print(" ")

            previousFilter = 1


        # Filter by conjugate plate ID
        if filter_ == 2:

            for feature in data_:
                for property in feature:

                    filter_property = property.get_name()

                    if filter_property.get_name() == "conjugatePlateId":
                        selected_filter_property = property.get_value()

                        # Isolate criteria match and process
                        if inverse == False:

                            for plateID in cPlateID:

                                if str(selected_filter_property) == str(plateID):

                                    # Append filtered data to associated Feature Collection
                                    f2_result.add(feature)

                        elif inverse == True:

                            if int(str(selected_filter_property)) not in cPlateID:
                                
                                # Append filtered data to associated Feature Collection
                                f2_result.add(feature)

            print("    2. Filtering data by conjugate plate ID(s): " + str(cPlateID))
            print("       - Found " + str(len(f2_result)) + " feature(s).")
            print(" ")

            previousFilter = 2



        # Filter by age of appearance
        if filter_ == 3:

            if ageAppearWindow[0] == "DP":
                ageAppearWindow[0] = float("inf")

            for feature in data_:

                begin_time, end_time = feature.get_valid_time()

                if begin_time <= ageAppearWindow[0] and begin_time >= ageAppearWindow[1]:
                    f3_result.add(feature)

            print("    3. Filtering data by age of appearance window: " + str(ageAppearWindow[0]) + " - " + str(ageAppearWindow[1]) + " Ma")
            print("       - Found " + str(len(f3_result)) + " feature(s).")
            print(" ")

            previousFilter = 3



        # Filter by age of disappearance
        if filter_ == 4:

            if ageDisappearWindow[1] == "DF":
                    ageDisappearWindow[1] = float("-inf")

            for feature in data_:

                begin_time, end_time = feature.get_valid_time()

                if end_time <= ageDisappearWindow[0] and end_time >= ageDisappearWindow[1]:
                    f4_result.add(feature)

            print("    4. Filtering data by age of disappearance window: " + str(ageDisappearWindow[0]) + " - " + str(ageDisappearWindow[1]) + " Ma")
            print("       - Found " + str(len(f4_result)) + " feature(s).")
            print(" ")

            previousFilter = 4



        # Filter by geographic selection / polygon
        if filter_ == 5:

            for feature in data_:
                for property in feature:

                    filter_property = property.get_name()

                    if filter_property.get_name() == "centerLineOf":
                        selected_filter_property = property.get_value()

                        points = selected_filter_property.get_value().get_base_curve().get_polyline().get_points_view()

                        for point in points:
                            point_latlong = pgp.convert_point_on_sphere_to_lat_lon_point(point)

                            if point_latlong.get_longitude() >= float(boundingBox[0]) and point_latlong.get_longitude() <= float(boundingBox[1])\
                                    and point_latlong.get_latitude() >= float(boundingBox[2]) and point_latlong.get_latitude() <= float(boundingBox[3]):

                                # If point is found within bounding box, add feature and break loop (search next feature)
                                f5_result.add(feature)
                                break

            print("    5. Filtering data by geographic bounding box: " + str(boundingBox[0]) + "/" + str(boundingBox[1]) + "/" + str(boundingBox[2]) + "/" + str(boundingBox[3]))
            print("       - Found " + str(len(f5_result)) + " feature(s).")
            print(" ")

            previousFilter = 5



        # Filter by age exists window
        if filter_ == 6:

            for feature in data_:

                begin_time, end_time = feature.get_valid_time()

                if begin_time >= ageExistsWindow[0] and end_time <= ageExistsWindow[1]:
                    f6_result.add(feature)
                elif begin_time >= ageExistsWindow[0] and end_time <= ageExistsWindow[0] and end_time >= ageExistsWindow[1]:
                    f6_result.add(feature)
                elif begin_time <= ageExistsWindow[0] and end_time >= ageExistsWindow[1]:
                    f6_result.add(feature)
                elif begin_time <= ageExistsWindow[0] and end_time >= ageExistsWindow[1] and end_time <= ageExistsWindow[1]:
                    f6_result.add(feature)

            print("    6. Filtering data by age of existence window: " + str(ageExistsWindow[0]) + " - " + str(ageExistsWindow[1]) + " Ma")
            print("       - Found " + str(len(f6_result)) + " feature(s).")
            print(" ")

            previousFilter = 6



        # Filter by feature type
        if filter_ == 7:

            iso_count = 0
            mor_count = 0
            pcb_count = 0

            for feature in data_:

                if "Isochron" in featureType:
                    if str(feature.get_feature_type()) == "gpml:Isochron":
                        f7_result.add(feature)
                        iso_count += 1
                if "MidOceanRidge" in featureType:
                    if str(feature.get_feature_type()) == "gpml:MidOceanRidge":
                        f7_result.add(feature)
                        mor_count += 1
                if "PassiveContinentalBoundary" in featureType:
                    if str(feature.get_feature_type()) == "gpml:PassiveContinentalBoundary":
                        f7_result.add(feature)
                        pcb_count += 1


            print("    7. Filtering data by feature type(s): " + str(featureType))

            if "Isochron" in featureType:
                print("       - Found " + str(iso_count) + " Isochron(s).")
            if "MidOceanRidge" in featureType:
                print("       - Found " + str(mor_count) + " MidOceanRidge(s).")
            if "PassiveContinentalBoundary" in featureType:
                print("       - Found " + str(pcb_count) + " PassiveContinentalBoundary(s).")

            print(" ")

            previousFilter = 7



        # Filter by geometry type
        if filter_ == 8:

            polylineCount = 0
            polygonCount = 0
            pointCount = 0
            multiPointCount = 0

            for feature in data_:

                geometries = feature.get_geometry()

                for geometry in geometryType:

                    if str(geometry) in str(geometries):
                        f8_result.add(feature)

                        if str(geometry) == "PolylineOnSphere":
                            polylineCount += 1
                        if str(geometry) == "PolygonOnSphere":
                            polygonCount += 1
                        if str(geometry) == "PointOnSphere":
                            pointCount += 1
                        if str(geometry) == "MultiPointOnSphere":
                            multiPointCount += 1

            print("    8. Filtering data by feature geometries present: " + str(geometryType))

            if "PolylineOnSphere" in geometryType:
                print("       - Found " + str(polylineCount) + " PolylineOnSphere(s).")
            if "PolygonOnSphere" in geometryType:
                print("       - Found " + str(polygonCount) + " PolygonOnSphere(s).")
            if "PointOnSphere" in geometryType:
                print("       - Found " + str(pointCount) + " PointOnSphere(s).")
            if "MultiPointOnSphere" in geometryType:
                print("       - Found " + str(multiPointCount) + " MultiPointOnSphere(s).")

            print(" ")

            previousFilter = 8



        # Filter by feature ID
        if filter_ == 9:

            for feature in data_:

                for id in featureID:
                    if str(feature.get_feature_id()).lower() == str(id).lower():
                        f9_result.add(feature)

            print("    9. Filtering data by feature ID: " + str(featureID))
            print("       - Found " + str(len(f9_result)) + " feature(s).")
            print(" ")

            previousFilter = 9



        # Filter by feature name (case insensitive)
        if filter_ == 10:

            for feature in data_:

                feature_name = feature.get_name()

                for names in featureName:

                    if names.lower() in feature_name.lower():
                        f10_result.add(feature)

            print("    10. Filtering data by feature name: " + str(featureName))
            print("       - Found " + str(len(f10_result)) + " feature(s).")
            print(" ")

            previousFilter = 10



        # Truncate file by age boundary
        if filter_ == 11:

            for feature in data_:

                begin_time, end_time = feature.get_valid_time()

                
                if begin_time > feature_truncate_age and end_time > feature_truncate_age:

                    f11a_result.add(feature)


                elif begin_time > feature_truncate_age and end_time <= feature_truncate_age:

                    # Special case if SubductionZone - need to incorporate start age
                    if str(feature.get_feature_type()) == "gpml:SubductionZone":

                        feature.set_valid_time(begin_time, feature_truncate_age + 0.1)
                        f11a_result.add(feature)
                        
                    else:

                        feature.set_valid_time(begin_time, feature_truncate_age + 0.1)
                        f11a_result.add(feature)


            # Limitation of pyGPlates: need to loop through copy of feature
            for feature in feature1:

                begin_time, end_time = feature.get_valid_time()

        
                if begin_time > feature_truncate_age and end_time <= feature_truncate_age:

                    # Special case if SubductionZone - need to incorporate start age
                    if str(feature.get_feature_type()) == "gpml:SubductionZone":

                        #print(feature.get_name())
                        sz_age_array = []

                        for property in feature:

                            # Create array or properties and add 'gpml:subductionZoneAge' if there isn't one already
                            sz_age_array.append(str(property.get_name()))

                        
                        if "gpml:subductionZoneAge" not in sz_age_array:
                                
                            #print(feature.get_name())
                            feature.add(pgp.PropertyName.create_gpml('subductionZoneAge'), pgp.XsDouble(begin_time))


                        feature.set_valid_time(feature_truncate_age, end_time)

                        f11b_result.add(feature)


                    else:

                        feature.set_valid_time(feature_truncate_age, end_time)
                        f11b_result.add(feature)


                elif begin_time == feature_truncate_age:

                    f11b_result.add(feature)


                elif begin_time <= feature_truncate_age and end_time < feature_truncate_age:

                    f11b_result.add(feature)

                
            print("    11. File truncated by age boundary: " + str(feature_truncate_age) + " Ma")
            print("       - Created " + str(len(f11a_result)) + " feature(s) older than truncation boundary.")
            print("       - Created " + str(len(f11b_result)) + " feature(s) younger than truncation boundary.")
            print(" ")

            previousFilter = 11




    # output new feature collection from filtered data to file

    if previousFilter == 11:

        iso_output1 = eval("f" + str(previousFilter) + "a_result")
        iso_output2 = eval("f" + str(previousFilter) + "b_result")

        if len(iso_output1) != 0:
            
            outputFeatureCollection1 = pgp.FeatureCollectionFileFormatRegistry()
            outputFeatureCollection1.write(iso_output1, "output/>" + str(feature_truncate_age) + "Ma_" + str(outputFile))

            print(" ")
            print("Output file 1:")
            print("    ../output/" + ">" + str(feature_truncate_age) + "Ma_"+ str(outputFile))
            print(" ")
            print("Process took " + str(round(time.time() - start, 2)) + " seconds.")
            print(" ")

        if len(iso_output2) != 0:
            
            outputFeatureCollection2 = pgp.FeatureCollectionFileFormatRegistry()
            outputFeatureCollection2.write(iso_output2, "output/<" + str(feature_truncate_age) + "Ma_" + str(outputFile))

            print("Output file 2:")
            print("    ../output/" + "<" + str(feature_truncate_age) + "Ma_"+ str(outputFile))
            print(" ")
            print("Process took " + str(round(time.time() - start, 2)) + " seconds.")
            print("--------------------------------------------")


    else:

        iso_output = eval("f" + str(previousFilter) + "_result")

        outputFeatureCollection = pgp.FeatureCollectionFileFormatRegistry()
        outputFeatureCollection.write(iso_output, outputFile)

        print("Output file:")
        print(str(outputFile))
        print(" ")
        print("Process took " + str(round(time.time() - start, 2)) + " seconds.")
        print("--------------------------------------------")
