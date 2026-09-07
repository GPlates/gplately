#!/usr/bin/env bash

set -euo pipefail

gplately ag seafloor-gridding-test-output/test-1-gplately-cli-config --config gplately-cli-config.toml

gplately ag seafloor-gridding-test-output/test-2-isochron-seafloor-gridding-1 --config isochron-seafloor-gridding-1.toml

gplately ag seafloor-gridding-test-output/test-3-isochron-seafloor-gridding-2 --config isochron-seafloor-gridding-2.toml

gplately ag seafloor-gridding-test-output/test-4-topology-seafloor-gridding-1 --config topology-seafloor-gridding-1.toml

gplately ag seafloor-gridding-test-output/test-5-topology-seafloor-gridding-2 --config topology-seafloor-gridding-2.toml