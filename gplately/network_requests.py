import asyncio
import concurrent.futures
import functools
import io
import os
import zipfile
from pathlib import Path

import requests

from typing import List, Type


def get_etag(url):
    """return the etag of the given url. could be none"""
    return requests.head(url).headers.get("ETag")


def fetch_file(
    url: str,
    filepath: str,
    filename: str = None,
    etag: str = None,
    auto_unzip: bool = True,
):
    """download a file from "url" and save to "filepath" """

    if etag:
        headers = {"If-None-Match": etag}
    else:
        headers = {}

    if os.path.isfile(filepath):
        raise Exception(
            f"The 'filepath' is in fact a file. The 'filepath' should be a folder path(non-exist is fine). {filepath}"
        )
    Path(filepath).mkdir(parents=True, exist_ok=True)

    r = requests.get(url, allow_redirects=True, headers=headers)
    # print(r.headers)

    if r.status_code == 304:
        print(url)
        print(
            "The file has not been changed since it was downloaded last time.Do nothing and return."
        )
    elif r.status_code == 200:
        if auto_unzip and url.endswith(".zip"):
            # unzip zip file
            z = zipfile.ZipFile(io.BytesIO(r.content))
            z.extractall(filepath)
        else:
            _save_file(url, filepath, filename, r.content)
    else:
        raise Exception(f"HTTP request failed with code {r.status_code}.")
    new_etag = r.headers.get("ETag").replace(
        "-gzip", ""
    )  # remove the content-encoding awareness thing

    return new_etag


def _fetch_range(url, index: int, chunk_size: int, data: list):
    """"""
    headers = {
        "Range": f"bytes={index*chunk_size}-{(index+1)*chunk_size-1}",
        "Accept-Encoding": "identity",
    }

    r = requests.get(url, headers=headers)
    data[index].write(r.content)


async def _fetch_large_file(
    run, url, file_size: int, data: list, chunk_size=2 * 1000 * 1000
):
    """"""
    num_chunks = file_size // chunk_size + 1
    data_array = [io.BytesIO() for i in range(num_chunks)]
    tasks = [
        run(
            _fetch_range,
            url,
            i,
            chunk_size,
            data_array,
        )
        for i in range(num_chunks)
    ]

    await asyncio.wait(tasks)

    for i in range(num_chunks):
        data_array[i].seek(0)
        data[0].write(data_array[i].read())


def fetch_large_file(
    url: str,
    filepath: str,
    filename: str = None,
    etag: str = None,
    auto_unzip: bool = True,
):
    """use multi-thread to fetch a large file.
        Be careful when use this function. You cannot get partial content if the content is gzip encoded.
        So the file might be larger than the one download directly.
        It is useful when downloading large .zip file.

    :param url: the file url
    :param file_size: the file size. if none, get via header request inside this function

    :returns:
    """

    # check file size and etag
    headers = {"Accept-Encoding": "identity"}
    r = requests.head(url, headers=headers)
    # print(r.headers)

    file_size = r.headers.get("Content-Length")
    if not file_size:
        raise Exception(
            "Unable to find the size of the file. Call fetch_file() instead."
        )
    else:
        file_size = int(file_size)

    new_etag = r.headers.get("ETag")
    if new_etag:
        new_etag = new_etag.replace(
            "-gzip", ""
        )  # remove the content-encoding awareness thing
        if new_etag == etag:
            print(url)
            print(
                "The file has not been changed since it was downloaded last time. Do nothing and return."
            )
            return new_etag

    # create folder to keep the file
    if os.path.isfile(filepath):
        raise Exception(
            f"The 'filepath' is in fact a file. The 'filepath' should be a folder path(non-exist is fine). {filepath}"
        )
    Path(filepath).mkdir(parents=True, exist_ok=True)

    # set up concurrent functions
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=15)
    loop = asyncio.new_event_loop()
    run = functools.partial(loop.run_in_executor, executor)

    asyncio.set_event_loop(loop)

    data = [io.BytesIO()]

    try:
        loop.run_until_complete(_fetch_large_file(run, url, file_size, data))
    finally:
        loop.close()

    data[0].seek(0)
    # save the file
    if auto_unzip and url.endswith(".zip"):
        # unzip zip file
        zipfile.ZipFile(data[0]).extractall(filepath)
    else:
        _save_file(url, filepath, filename, data[0].read())

    return new_etag


def _save_file(url, filepath, filename, data):
    """"""
    Path(filepath).mkdir(parents=True, exist_ok=True)
    if not filename:
        filename = url.split("/")[-1]  # use the filename in the url
    if os.path.isfile(f"{filepath}/{filename}"):
        print(f"Warning: overwriting {filename}")
    with open(f"{filepath}/{filename}", "wb+") as of:
        of.write(data)
