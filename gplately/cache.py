import hashlib, uuid
import platformdirs

from . import network


def get_user_cache_dir():
    """return the cache dir which is different across operating systems"""
    return platformdirs.user_cache_dir("gplately", "")


def get(url: str):
    url_id = uuid.UUID(hex=hashlib.md5(url.encode("UTF-8")).hexdigest())
    cache_path = f"{get_user_cache_dir()}/{url_id}"
    etag = network.download(url, cache_path)
    with open(f"{cache_path}.meta", "w+") as f:
        f.write(f"url={url}\n")
        if etag:
            f.write(f"etag={etag}\n")
    return cache_path
