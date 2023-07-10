import sys
import time

sys.path.insert(0, "../")
from gplately import network_aiohttp, network_requests

test_urls = [
    f"https://www.earthbyte.org/webdav/ftp/Data_Collections/Zahirovic_etal_2016_ESR_AgeGrid/jpegs/EarthByte_Zahirovic_etal_2016_ESR_r888_AgeGrid-{i}.jpg"
    for i in range(20)
]

auto_unzip = False


def test_with_for_loop():
    """requests + "for loop" """
    st = time.time()
    spt = time.process_time()

    print("Start test_with_for_loop ... ")

    for url in test_urls:
        network_requests.fetch_file(
            url,
            "./download-with-for-loop/",
            auto_unzip=auto_unzip,
        )

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")

    print("End test_with_for_loop ... ")


def test_concurrent_aiohttp():
    """asyncio + aiohttp"""
    st = time.time()
    spt = time.process_time()
    paths = "./download-concurrently-with-aiohttp/"

    print("Start test_concurrent_aiohttp ... ")

    network_aiohttp.fetch_files(
        test_urls,
        paths,
        auto_unzip=auto_unzip,
    )

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")

    print("End test_concurrent_aiohttp ... ")


def test_concurrent_executor():
    """requests + ThreadPoolExecutor + asyncio"""
    st = time.time()
    spt = time.process_time()

    paths = "./download-concurrently-with-executor/"

    print("Start test_concurrent_executor ... ")

    network_requests.fetch_files(
        test_urls,
        paths,
        auto_unzip=auto_unzip,
    )

    et = time.time()
    ept = time.process_time()

    print(f"time: {et - st}")
    print(f"process time: {ept - spt}")

    print("End test_concurrent_executor ... ")


test_with_for_loop()
test_concurrent_aiohttp()
test_concurrent_executor()
