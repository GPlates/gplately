import sys
import time

sys.path.insert(0, "../")
from gplately import cache, network

st = time.time()
spt = time.process_time()

print(cache.get_user_cache_dir())

print(
    network.get_etag(
        "https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics_Updated.zip"
    )
)

print(
    cache.get(
        "https://www.earthbyte.org/webdav/ftp/Data_Collections/Gibbons_etal_2013_JGR.zip"
    )
)

print(
    cache.get(
        "https://www.earthbyte.org/webdav/ftp/Data_Collections/Zahirovic_etal_2016_ESR_AgeGrid/netCDF/EarthByte_Zahirovic_etal_2016_ESR_r888_AgeGrid-0.nc"
    )
)

et = time.time()
ept = time.process_time()


print(f"time: {et - st}")
print(f"process time: {ept - spt}")
