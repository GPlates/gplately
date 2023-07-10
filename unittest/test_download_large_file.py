import sys
import time

sys.path.insert(0, "../")
from gplately import network_aiohttp, network_requests

# concurrent download from www.earthbyte.org does not work because of Cloud provider's network traffic control
# use "http://212.183.159.230/100MB.zip" to test. you can see the performance improvement
# the performance improvement depends on the server/network's configuration

test_url = "http://212.183.159.230/100MB.zip"
# test_url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Zahirovic_etal_2016_ESR_AgeGrid/jpegs.zip"

auto_unzip = False


def test_download_directly():
    """Download the large file directly with requests"""
    st = time.time()
    spt = time.process_time()

    print("Start download directly ...")

    network_requests.fetch_file(test_url, "./download-directly", auto_unzip=auto_unzip)

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")

    print("End download directly ...")


def test_requests():
    """requests + ThreadPoolExecutor + asyncio"""
    st = time.time()
    spt = time.process_time()

    print("Start download with requests+executor ...")

    network_requests.fetch_large_file(
        test_url, "./download-with-requests-executor", auto_unzip=auto_unzip
    )

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")

    print("End download with requests+executor ...")


def test_aiohttp():
    """Download with aiohttp+asyncio"""
    st = time.time()
    spt = time.process_time()

    print("Start download with aiohttp ...")

    network_aiohttp.fetch_large_file(
        test_url, "./download-with-aiohttp", auto_unzip=auto_unzip
    )

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")

    print("End download with aiohttp ...")


test_download_directly()
test_requests()
test_aiohttp()
