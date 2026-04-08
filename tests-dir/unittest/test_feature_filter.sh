#!/bin/bash
set -euo pipefail  # Exit on error, undefined vars, pipe failures

export DISABLE_GPLATELY_DEV_WARNING=true

TEST_DATA_DIR="./output/test-feature-filter-data"
mkdir -p "$TEST_DATA_DIR"

# Trap to show which command failed
trap 'echo "ERROR: Command failed on line $LINENO"; exit 1' ERR

# download the test file
IN_FILE="$TEST_DATA_DIR/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"
FILE_URL="https://repo.gplates.org/webdav/mchin/data/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"

if ! test -f "$IN_FILE"; then
    echo "Downloading test file..."
    if command -v curl >/dev/null 2>&1; then
        curl -o "$IN_FILE" "$FILE_URL" || { echo "ERROR: curl failed to download test file"; exit 1; }
    elif command -v wget >/dev/null 2>&1; then
        wget -O "$IN_FILE" "$FILE_URL" || { echo "ERROR: wget failed to download test file"; exit 1; }
    else
        echo "ERROR: Neither curl nor wget found. Download test file manually from: $FILE_URL to $IN_FILE"
        exit 1
    fi
fi

if ! test -f "$IN_FILE"; then
    echo "ERROR: Test file not found at $IN_FILE"
    echo "Download from: $FILE_URL"
    exit 1
fi

echo "Test file ready: $IN_FILE"
echo "Starting tests..."

# get features whose name contains "Africa" or "North America"
echo "Test 1: Filter by name (Africa, North America)..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/africa_north_america.gpmlz" -n Africa "North America"

# get features whose plate ID is one of 701 714 715 101
echo "Test 2: Filter by plate ID (701, 714, 715, 101)..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/pid_701_714_715_101.gpmlz" -p 701 714 715 101

# get features whose birth age is older than 500Myr
echo "Test 3: Filter by min birth age (500 Myr)..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/birth_age_older_500.gpmlz" --min-birth-age 500

# get features whose birth age is younger than 500Myr
echo "Test 4: Filter by max birth age (500 Myr)..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/birth_age_younger_500.gpmlz" --max-birth-age 500

# get features whose name conains "Africa" or "North America" and plate ID is one of 701 714 715 101
echo "Test 5: Filter by name AND plate ID..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/africa_north_america_pid_701_714_715_101.gpmlz" -n Africa "North America" -p 701 714 715 101

# get features whose name conains "Africa" or "North America" and plate ID is one of 701 714 715 101 and birth age is older than 500Myr
echo "Test 6: Filter by name AND plate ID AND min birth age..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/africa_north_america_pid_701_714_715_101_birth_age_500.gpmlz" -n Africa "North America" -p 701 714 715 101 --min-birth-age 500

# only Africa is saved because the name is --case-sensitive
echo "Test 7: Case-sensitive name filter..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/africa.gpmlz" -n Africa "North america" --case-sensitive

# only Africa is saved because the name is --exact-match
echo "Test 8: Exact name match filter..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/africa_1.gpmlz" -n Africa "North America" --exact-match

# get features whose name does not conain "Africa" or "North America"
echo "Test 9: Exclude by name (Africa, North America)..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/no_africa_north_america.gpmlz" --exclude-names Africa "North America"

# get features whose name does not conain "Africa" or "North America" and plate ID is not in 401 702
echo "Test 10: Exclude by name AND plate ID..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/no_africa_north_america_no_401_702.gpmlz" --exclude-names Africa "North America" --exclude-pids 401 702

# get all gpml:Basin features
echo "Test 11: Filter by feature type (gpml:Basin)..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/basins.gpmlz" -t gpml:Basin 

# get all gpml:Basin + gpml:IslandArc features
echo "Test 12: Filter by multiple feature types..."
gplately filter "$IN_FILE" "$TEST_DATA_DIR/basin_islandarc.gpmlz" -t "gpml:IslandArc|gpml:Basin"

echo ""
echo "✓ All tests passed successfully!"
