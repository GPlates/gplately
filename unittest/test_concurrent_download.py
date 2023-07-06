import sys
import time

sys.path.insert(0, "../")
from gplately import network, network_requests

test_url = "https://www.earthbyte.org/webdav/ftp/Data_Collections/Zahirovic_etal_2016_ESR_AgeGrid/netCDF/EarthByte_Zahirovic_etal_2016_ESR_r888_AgeGrid-0.nc"
# test_url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Zahirovic_etal_2016_ESR_AgeGrid/jpegs.zip"

auto_unzip = False


def test():
    st = time.time()
    spt = time.process_time()

    for i in range(20):
        print(i)
        network.fetch_file(
            f"https://www.earthbyte.org/webdav/ftp/Data_Collections/Zahirovic_etal_2016_ESR_AgeGrid/jpegs/EarthByte_Zahirovic_etal_2016_ESR_r888_AgeGrid-{i}.jpg",
            "./download-directly/",
            auto_unzip=auto_unzip,
        )

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")


def test_concurrent():
    st = time.time()
    spt = time.process_time()
    urls = []
    paths = []
    for i in range(20):
        urls.append(
            f"https://www.earthbyte.org/webdav/ftp/Data_Collections/Zahirovic_etal_2016_ESR_AgeGrid/jpegs/EarthByte_Zahirovic_etal_2016_ESR_r888_AgeGrid-{i}.jpg",
        )
        paths.append("./download-concurrent/")

    network.fetch_files(
        urls,
        paths,
        auto_unzip=auto_unzip,
    )

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")


test()
test_concurrent()
