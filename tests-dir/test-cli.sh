#!/usr/bin/env bash

set -euo pipefail

gplately list

gplately get_cli_config_example 

gplately gpmdb -m Zahirovic2022 -o test-gpmdb.gpmlz

pmm download muller2025 plate-model-repo

gplately paleobathymetry test-paleobathymetry \
    -m muller2025 -f plate-model-repo \
    -e 0 -s 1 --time-step 1 -r 10 --max-reconstruction-time 5

gplately gdg test-distance-grids \
    -m muller2025 -f plate-model-repo \
    -e 0 -s 1 --time-step 1 -r 10 --max-reconstruction-time 5

gplately gsg test-sediment-grids \
    -m muller2025 -f plate-model-repo \
    -e 0 -s 1 --time-step 1 -r 10 \
    --distance-grids-dir test-distance-grids

gplately combine \
    plate-model-repo/muller2025/Coastlines/shapes_coasts.gpmlz \
    plate-model-repo/muller2025/ContinentalPolygons/shapes_continents.gpmlz \
    test-gplates-combine.gpmlz

gplately reset_feature_type -s gpml:ClosedContinentalBoundary -t gpml:UnclassifiedFeature \
    test-gplates-combine.gpmlz \
    test-gplates-reset.gpmlz

./unittest/test_feature_filter.sh
./unittest/test_reset_feature_type.sh
./unittest/test-seafloor-gridding.sh

rm test-gpmdb.gpmlz test-gplates-combine.gpmlz test-gplates-reset.gpmlz
rm -r test-paleobathymetry test-distance-grids test-sediment-grids
rm ptt.log gplately.log
echo "All CLI tests passed!"