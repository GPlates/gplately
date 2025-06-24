export DISABLE_GPLATELY_DEV_WARNING=true

TEST_DATA_DIR="./output/test-feature-filter-data"
mkdir -p $TEST_DATA_DIR


# download the test file
IN_FILE="$TEST_DATA_DIR/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"
FILE_URL="https://repo.gplates.org/webdav/mchin/data/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz"

if ! test -f $IN_FILE; then
    if command -v curl; then
        curl -o $IN_FILE "$FILE_URL"
    else
        if command -v wget; then
            wget -O $IN_FILE "$FILE_URL"
        fi
    fi
fi

if ! test -f $IN_FILE; then
    echo "download test file from" "$FILE_URL" "to folder $TEST_DATA_DIR first and then try again!"
fi


# get features whose name contains "Africa" or "North America"

gplately filter "$IN_FILE" $TEST_DATA_DIR/africa_north_america.gpmlz -n Africa "North America"

# get features whose plate ID is one of 701 714 715 101

gplately filter "$IN_FILE" $TEST_DATA_DIR/pid_701_714_715_101.gpmlz -p 701 714 715 101

# get features whose birth age is older than 500Myr

gplately filter "$IN_FILE" $TEST_DATA_DIR/birth_age_older_500.gpmlz --min-birth-age 500

# get features whose birth age is younger than 500Myr

gplately filter "$IN_FILE" $TEST_DATA_DIR/birth_age_younger_500.gpmlz --max-birth-age 500

# get features whose name conains "Africa" or "North America" and plate ID is one of 701 714 715 101

gplately filter "$IN_FILE" $TEST_DATA_DIR/africa_north_america_pid_701_714_715_101.gpmlz -n Africa "North America" -p 701 714 715 101

# get features whose name conains "Africa" or "North America" and plate ID is one of 701 714 715 101 and birth age is older than 500Myr

gplately filter "$IN_FILE" $TEST_DATA_DIR/africa_north_america_pid_701_714_715_101_birth_age_500.gpmlz -n Africa "North America" -p 701 714 715 101 --min-birth-age 500

# only Africa is saved because the name is --case-sensitive

gplately filter "$IN_FILE" $TEST_DATA_DIR/africa.gpmlz -n Africa "North america" --case-sensitive

# only Africa is saved because the name is --exact-match

gplately filter "$IN_FILE" $TEST_DATA_DIR/africa_1.gpmlz -n Africa "North America" --exact-match

# get features whose name does not conain "Africa" or "North America"

gplately filter "$IN_FILE" $TEST_DATA_DIR/no_africa_north_america.gpmlz --exclude-names Africa "North America"

# get features whose name does not conain "Africa" or "North America" and plate ID is not in 401 702

gplately filter "$IN_FILE" $TEST_DATA_DIR/no_africa_north_america_no_401_702.gpmlz --exclude-names Africa "North America" --exclude-pids 401 702

# get all gpml:Basin features

gplately filter "$IN_FILE" $TEST_DATA_DIR/basins.gpmlz -t gpml:Basin 

# get all gpml:Basin + gpml:IslandArc features

gplately filter "$IN_FILE" $TEST_DATA_DIR/basin_islandarc.gpmlz -t "gpml:IslandArc|gpml:Basin"
