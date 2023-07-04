import requests


def get_etag(url):
    """return the etag of the given url. could be none"""
    return requests.head(url).headers.get("ETag")


def download(url, filename):
    """download a file from "url" and save to "filename"

    Parameters
    ----------
    url : str
        url
    filename : str
        local file path

    Returns
    ----------
    ETag

    """

    r = requests.get(url, allow_redirects=True)
    open(filename, "wb").write(r.content)
    return r.headers.get("ETag")
