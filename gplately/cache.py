import hashlib, uuid, os
import platformdirs

from . import network


def get_user_cache_dir():
    """return the cache dir which is different across operating systems"""
    return platformdirs.user_cache_dir("gplately", "")


def get(url: str):
    """get the file at url. return a folder path for zip file. return a file path for other files."""

    url_id = uuid.UUID(hex=hashlib.md5(url.encode("UTF-8")).hexdigest())
    cache_path = f"{get_user_cache_dir()}/{url_id}"
    meta_file = f"{cache_path}.meta"

    current_etag = None
    if os.path.isfile(meta_file):
        with open(meta_file, "r") as f:
            for line in f:
                if line.startswith("etag="):
                    current_etag = line[5:-1]

    if url.endswith(".zip"):
        etag = network.fetch_file(url, cache_path, etag=current_etag)
    else:
        etag = network.fetch_file(url, cache_path, etag=current_etag)

    with open(f"{cache_path}.meta", "w+") as f:
        f.write(f"url={url}\n")
        if etag:
            f.write(f"etag={etag}\n")
    return cache_path
