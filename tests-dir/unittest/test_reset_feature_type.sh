#!/bin/bash
export DISABLE_GPLATELY_DEV_WARNING=true

TEST_DATA_DIR="test-reset-feature-type-data"
mkdir -p $TEST_DATA_DIR


# download the test file
IN_FILE="$TEST_DATA_DIR/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"
FILE_URL="https://repo.gplates.org/webdav/mchin/data/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"

if ! test -f $IN_FILE; then
    if command -v curl; then
        curl -o $IN_FILE "$FILE_URL"
    else
        if command -v wget; then
            wget "$FILE_URL"
        fi
    fi
fi

if ! test -f $IN_FILE; then
    echo "download test file from" "$FILE_URL" "to folder $TEST_DATA_DIR first and then try again!"
fi

# case 1
echo "Reset all gpml:ClosedContinentalBoundary features to gpml:UnclassifiedFeature."
gplately reset_feature_type -s gpml:ClosedContinentalBoundary -t gpml:UnclassifiedFeature $IN_FILE $TEST_DATA_DIR/ClosedContinentalBoundary-to-UnclassifiedFeature.gpml

# case 2
echo "Reset all gpml:ContinentalFragment and gpml:Coastline features to gpml:UnclassifiedFeature."
gplately reset_feature_type -s "gpml:ContinentalFragment|gpml:Coastline" -t gpml:UnclassifiedFeature $IN_FILE $TEST_DATA_DIR/ContinentalFragment-Coastline-to-UnclassifiedFeature.gpml

# case  3
echo "Reset all features to gpml:UnclassifiedFeature and use the same feature IDs in the new file."
gplately reset_feature_type -s ".*" -t gpml:UnclassifiedFeature --keep-feature-id $IN_FILE $TEST_DATA_DIR/all-to-UnclassifiedFeature.gpml

# case 4
echo "It is by design if you see an error message below. Do not use --verify-information-model if you need not to guarantee GPGIM integrity."
gplately reset_feature_type -s ".*" -t gpml:Transform --verify-information-model $IN_FILE $TEST_DATA_DIR/should-not-exist.gpml

# case 5
echo "It is by design if you see error messages below. Users need to make sure using the correct feature type."
gplately reset_feature_type -s ".*" -t gpml:BadFeatureType  $IN_FILE $TEST_DATA_DIR/should-not-exist.gpml

echo "Testing reset_feature_type is done. Open the output files in GPlates and check the results."
