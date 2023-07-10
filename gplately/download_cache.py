import hashlib, uuid, os
import platformdirs

from . import network_aiohttp


def get_user_cache_dir():
    """return the cache dir which is different across operating systems"""
    return platformdirs.user_cache_dir("gplately", "")


def get(url: str, auto_unzip: bool = True, large_file: bool = False):
    """get the file at url. return a folder path for zip file. return a file path for other files."""

    url_id = uuid.UUID(hex=hashlib.md5(url.encode("UTF-8")).hexdigest())
    cache_path = f"{get_user_cache_dir()}/{url_id}"
    meta_file = f"{cache_path}.meta"

    # get the etag from meta file
    current_etag = None
    if os.path.isfile(meta_file):
        with open(meta_file, "r") as f:
            for line in f:
                if line.startswith("etag="):
                    current_etag = line[5:-1]

    if large_file:
        etag = network_aiohttp.fetch_large_file(
            url, cache_path, etag=current_etag, auto_unzip=auto_unzip
        )
    else:
        etag = network_aiohttp.fetch_file(
            url, cache_path, etag=current_etag, auto_unzip=auto_unzip
        )

    with open(f"{cache_path}.meta", "w+") as f:
        f.write(f"url={url}\n")
        if etag:
            f.write(f"etag={etag}\n")
    return cache_path


def get_all(urls, auto_unzip: bool = True):
    """get the file at url. return a folder path for zip file. return a file path for other files."""

    filepaths = []
    etags = []
    for url in urls:
        url_id = uuid.UUID(hex=hashlib.md5(url.encode("UTF-8")).hexdigest())
        cache_path = f"{get_user_cache_dir()}/{url_id}"
        meta_file = f"{cache_path}.meta"

        # get the etag from meta file
        current_etag = None
        if os.path.isfile(meta_file):
            with open(meta_file, "r") as f:
                for line in f:
                    if line.startswith("etag="):
                        current_etag = line[5:-1]
        filepaths.append(cache_path)
        etags.append(current_etag)

    new_etags = network_aiohttp.fetch_files(
        urls, filepaths, etags=etags, auto_unzip=auto_unzip
    )
    print(new_etags)

    for idx, filepath in enumerate(filepaths):
        with open(f"{filepath}.meta", "w+") as f:
            f.write(f"url={urls[idx]}\n")
            if len(new_etags) > idx:
                etag = new_etags[idx]
                if etag:
                    f.write(f"etag={etag}\n")
    return filepaths
