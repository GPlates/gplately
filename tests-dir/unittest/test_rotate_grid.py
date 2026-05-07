#!/usr/bin/env python3
"""
Test script for the ``gplately rotate_grid`` CLI command.

Steps
-----
1. Download paleobathymetry grids (100 – 110 Ma) from the Zahirovic 2022
   *Optimised Mantle Frame* WebDAV folder into ``Zahirovic_mantle/``.
2. Download the matching grids from the *Palaeomagnetic Frame* WebDAV folder
   into ``Zahirovic_pmag/`` (used as reference).
3. Run ``gplately rotate_grid`` to rotate the mantle-frame grids into the
   palaeomagnetic reference frame (anchor plate 701701), writing results to
   ``Zahirovic2022_PMAG_rotated/``.
4. Compare the rotated grids against the reference pmag grids.  The values
   will not be identical (different rotation pipeline / numerical precision)
   but they should be strongly correlated and the mean absolute difference
   should be small relative to the data range.

Usage
-----
    cd tests-dir/unittest
    python test_rotate_grid.py
"""

import os
import re
import subprocess
import sys
from pathlib import Path

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

import numpy as np
import requests
import xarray as xr

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

MANTLE_WEBDAV = (
    "https://repo.gplates.org/webdav/PlateModel_Age_SR_Grids/"
    "Zahirovic_etal_2022_GDJ/02_AgegridsUsingTopologies/"
    "OptimisedMantleFrame/PaleobathymetryGrids/GDH1/"
)

PMAG_WEBDAV = (
    "https://repo.gplates.org/webdav/PlateModel_Age_SR_Grids/"
    "Zahirovic_etal_2022_GDJ/02_AgegridsUsingTopologies/"
    "PaleomagneticFrame/PaleobathymetryGrids/GDH1/"
)

MANTLE_DIR = Path("Zahirovic_mantle")
PMAG_DIR = Path("Zahirovic_pmag")
ROTATED_DIR = Path("Zahirovic2022_PMAG_rotated")

TIME_MIN = 100  # Ma
TIME_MAX = 110  # Ma

# Acceptance thresholds for the automated comparison
MIN_CORRELATION = 0.95  # Pearson r must be at least this
MAX_MEAN_ABS_DIFF_FRACTION = 0.10  # mean |diff| / value-range ≤ 10 %

DOWNLOAD_CHUNK_SIZE = 64 * 1024  # 64 KB per streaming chunk

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_TIME_RE = re.compile(r"(\d+(?:\.\d+)?)Ma", re.IGNORECASE)


def _time_from_filename(fname):
    """Return the reconstruction time (Ma) encoded in *fname*, or None."""
    m = _TIME_RE.search(str(fname))
    return float(m.group(1)) if m else None


def _list_nc_files(base_url):
    """Return a list of .nc filenames found in an HTTP-browseable directory."""
    resp = requests.get(base_url, timeout=60)
    resp.raise_for_status()
    # Extract href values pointing at .nc files from the HTML directory listing
    return re.findall(r'href="([^"?#]+\.nc)"', resp.text)


# ---------------------------------------------------------------------------
# Step 1 & 2: Download grids from WebDAV
# ---------------------------------------------------------------------------


def download_grids(base_url, dest_dir):
    """Download all .nc files for TIME_MIN – TIME_MAX Ma into *dest_dir*."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    all_files = _list_nc_files(base_url)
    target = []
    for f in all_files:
        t = _time_from_filename(f)
        if t is not None and TIME_MIN <= t <= TIME_MAX:
            target.append(f)
    if not target:
        raise RuntimeError(
            f"No .nc files found for {TIME_MIN}–{TIME_MAX} Ma in {base_url}"
        )
    print(f"Downloading {len(target)} file(s) → {dest_dir} …")
    for fname in sorted(target):
        # fname may be a bare filename or a relative/absolute path
        url = base_url.rstrip("/") + "/" + Path(fname).name
        dest = dest_dir / Path(fname).name
        if dest.exists():
            print(f"  (cached) {dest.name}")
            continue
        print(f"  {dest.name} …", end="", flush=True)
        with requests.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(chunk_size=DOWNLOAD_CHUNK_SIZE):
                    fh.write(chunk)
        print(" done")


# ---------------------------------------------------------------------------
# Step 3: Run gplately rotate_grid
# ---------------------------------------------------------------------------


def rotate_grids():
    """Invoke ``gplately rotate_grid`` as a subprocess."""
    cmd = [
        sys.executable,
        "-m",
        "gplately",
        "rotate_grid",
        str(MANTLE_DIR),
        str(ROTATED_DIR),
        "-r",
        "0.2",
        "-j",
        "4",
        "--from-model",
        "zahirovic2022",
        "--to-model",
        "zahirovic2022",
        "--from-anchor",
        "0",
        "--to-anchor",
        "701701",
    ]
    print("\nRunning:", " ".join(cmd))
    subprocess.run(cmd, check=True, text=True)


# ---------------------------------------------------------------------------
# Step 4: Compare rotated grids with reference pmag grids
# ---------------------------------------------------------------------------


def _load_first_var(nc_path):
    """Load the first data variable from a NetCDF file as a numpy array."""
    ds = xr.open_dataset(nc_path)
    var = list(ds.data_vars)[0]
    arr = ds[var].values.astype(float)
    ds.close()
    return arr


def compare_grids():
    """
    Compare each rotated grid against its matching reference pmag grid.

    Matching is done by reconstruction time extracted from filenames.
    Returns True if all comparisons pass, False otherwise.
    """
    rotated_files = sorted(ROTATED_DIR.glob("*.nc"))
    if not rotated_files:
        print("ERROR: No rotated output files found in", ROTATED_DIR)
        return False

    pmag_files = list(PMAG_DIR.glob("*.nc"))

    all_pass = True
    for rot_file in rotated_files:
        t = _time_from_filename(rot_file.name)
        if t is None:
            print(f"  SKIP {rot_file.name}: cannot parse time from filename")
            continue

        # Find the reference pmag file for the same time
        ref_candidates = [f for f in pmag_files if _time_from_filename(f.name) == t]
        if not ref_candidates:
            print(f"  SKIP {rot_file.name}: no reference file found for {t:.0f} Ma")
            continue
        ref_file = ref_candidates[0]

        rot_arr = _load_first_var(rot_file)
        ref_arr = _load_first_var(ref_file)

        # If grids have different shapes, resample the reference to match the
        # rotated grid using xarray's interp (linear, on the lon/lat axes).
        if rot_arr.shape != ref_arr.shape:
            ds_rot = xr.open_dataset(rot_file)
            ds_ref = xr.open_dataset(ref_file)
            rot_var = list(ds_rot.data_vars)[0]
            ref_var = list(ds_ref.data_vars)[0]

            # Identify coordinate names (lat/lon may appear under various names)
            lat_dim = next(
                (c for c in ds_rot[rot_var].dims if "lat" in c.lower()), None
            )
            lon_dim = next(
                (c for c in ds_rot[rot_var].dims if "lon" in c.lower()), None
            )

            if lat_dim and lon_dim:
                ref_resampled = ds_ref[ref_var].interp(
                    {
                        lat_dim: ds_rot[lat_dim],
                        lon_dim: ds_rot[lon_dim],
                    },
                    method="linear",
                )
                ref_arr = ref_resampled.values.astype(float)
            else:
                # Cannot safely align; fall back to flattened comparison
                rot_arr = rot_arr.ravel()
                ref_arr = ref_arr.ravel()
                min_len = min(len(rot_arr), len(ref_arr))
                rot_arr = rot_arr[:min_len]
                ref_arr = ref_arr[:min_len]

            ds_rot.close()
            ds_ref.close()

        rot_v = rot_arr.ravel()
        ref_v = ref_arr.ravel()

        # Mask NaN in either grid
        valid = ~(np.isnan(rot_v) | np.isnan(ref_v))
        rot_v = rot_v[valid]
        ref_v = ref_v[valid]

        if rot_v.size < 10:
            print(f"  SKIP {rot_file.name}: too few valid overlap points")
            continue

        # Pearson correlation
        corr = float(np.corrcoef(rot_v, ref_v)[0, 1])

        # Mean absolute difference as a fraction of the reference value range
        ref_range = float(ref_v.max() - ref_v.min())
        if ref_range == 0:
            # Cannot compute a meaningful fraction; skip the MAD check
            mad_frac = None
        else:
            mad_frac = float(np.mean(np.abs(rot_v - ref_v)) / ref_range)

        corr_ok = corr >= MIN_CORRELATION
        mad_ok = mad_frac is None or mad_frac <= MAX_MEAN_ABS_DIFF_FRACTION
        passed = corr_ok and mad_ok
        status = "PASS" if passed else "FAIL"
        mad_str = "N/A (zero range)" if mad_frac is None else f"{mad_frac:.4f}"
        print(
            f"  [{status}] {t:.0f} Ma  |  corr={corr:.4f}  |  "
            f"mean_abs_diff/range={mad_str}  "
            f"({rot_file.name})"
        )

        if not passed:
            all_pass = False

    return all_pass


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("=== Step 1: Download Optimised Mantle Frame grids ===")
    download_grids(MANTLE_WEBDAV, MANTLE_DIR)

    print("\n=== Step 2: Download Palaeomagnetic Frame grids (reference) ===")
    download_grids(PMAG_WEBDAV, PMAG_DIR)

    print("\n=== Step 3: Rotate mantle-frame grids to pmag reference frame ===")
    rotate_grids()

    print("\n=== Step 4: Compare rotated grids against reference pmag grids ===")
    ok = compare_grids()

    if ok:
        print("\nAll comparisons PASSED.")
    else:
        print("\nSome comparisons FAILED — please review the output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
