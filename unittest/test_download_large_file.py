import sys
import time

sys.path.insert(0, "../")
from gplately import network, network_requests

test_url = "http://212.183.159.230/100MB.zip"
# test_url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Zahirovic_etal_2016_ESR_AgeGrid/jpegs.zip"

auto_unzip = False


def test():
    st = time.time()
    spt = time.process_time()

    network.fetch_file(test_url, "./download-directly", auto_unzip=auto_unzip)

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")


def test_requests():
    st = time.time()
    spt = time.process_time()

    network_requests.fetch_large_file(
        test_url, "./download-with-requests", auto_unzip=auto_unzip
    )

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")


def test_aiohttp():
    st = time.time()
    spt = time.process_time()

    network.fetch_large_file(test_url, "./download-with-aiohttp", auto_unzip=auto_unzip)

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")


test()
test_requests()
test_aiohttp()
