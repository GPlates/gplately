import requests, zipfile, io, os
from pathlib import Path


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

    r = requests.get(url, allow_redirects=True, headers=headers)
    # print(r.headers)

    if os.path.isfile(filepath):
        raise Exception(
            f"The 'filepath' is in fact a file. The 'filepath' should be a folder path(non-exist is fine). {filepath}"
        )
    Path(filepath).mkdir(parents=True, exist_ok=True)

    if r.status_code == 304:
        print("The file has not been changed since it was downloaded last time.")
    elif r.status_code == 200:
        if auto_unzip and url.endswith(".zip"):
            # unzip zip file
            z = zipfile.ZipFile(io.BytesIO(r.content))
            Path(filepath).mkdir(parents=True, exist_ok=True)
            z.extractall(filepath)
        else:
            if not filename:
                filename = url.split("/")[-1]  # use the filename in the url
            if os.path.isfile(f"{filepath}/{filename}"):
                print(f"Warning: overwriting {filename}")
            with open(f"{filepath}/{filename}", "wb+") as of:
                of.write(r.content)
    else:
        raise Exception(f"HTTP request failed with code {r.status_code}.")

    new_etag = r.headers.get("ETag").replace(
        "-gzip", ""
    )  # remove the content-encoding awareness thing
    return new_etag
