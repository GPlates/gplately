import sys

sys.path.insert(0, "../")
import time
import gplately
from platformdirs import *

print(user_cache_dir("gplately", ""))

st = time.time()
spt = time.process_time()


gdownload = gplately.download.DataServer("Muller2019")
(
    rotation_model,
    topology_features,
    static_polygons,
) = gdownload.get_plate_reconstruction_files()
# coastlines, continents, COBs = gdownload.get_topology_geometries()


et = time.time()
ept = time.process_time()


print(f"time: {et - st}")
print(f"process time: {ept - spt}")
