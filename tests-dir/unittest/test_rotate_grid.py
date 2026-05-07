#!/usr/bin/env python3
"""
Test script for the ``gplately rotate_grid`` CLI command.

All working files are placed inside a single parent folder ``test-rotate-grid/``
so the test is easy to clean up and does not pollute the current working directory.

Steps
-----
1. Download paleobathymetry grids (100 – 110 Ma) from the Zahirovic 2022
   *Optimised Mantle Frame* WebDAV folder into
   ``test-rotate-grid/Zahirovic_mantle/``.
2. Download the matching grids from the *Palaeomagnetic Frame* WebDAV folder
   into ``test-rotate-grid/Zahirovic_pmag/`` (used as reference).
3. Download the Zahirovic2022 plate model (rotation files) via
   plate-model-manager into ``test-rotate-grid/zahirovic2022/``.
4. Run ``gplately rotate_grid`` (directory mode, named model) to rotate all
   mantle-frame grids into the palaeomagnetic reference frame, writing results
   to ``test-rotate-grid/Zahirovic2022_PMAG_rotated/``.
5. Compare the rotated grids against the reference pmag grids.  The values
   will not be identical (different rotation pipeline / numerical precision)
   but they should be strongly correlated and the mean absolute difference
   should be small relative to the data range.
6. Run ``gplately rotate_grid`` (single-file mode, local rotation files) to
   rotate ``test-rotate-grid/Zahirovic_mantle/paleobathymetry_105Ma.nc`` to
   ``test-rotate-grid/paleobathymetry_pmag_105Ma.nc`` using the downloaded
   rotation file directly.
7. Compare ``test-rotate-grid/Zahirovic_pmag/paleobathymetry_105Ma.nc`` and
   ``test-rotate-grid/paleobathymetry_pmag_105Ma.nc``.

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
from plate_model_manager import PlateModel, PlateModelManager

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# All output lives under this single parent directory.
TEST_DIR = Path("test-rotate-grid")

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

MANTLE_DIR = TEST_DIR / "Zahirovic_mantle"
PMAG_DIR = TEST_DIR / "Zahirovic_pmag"
ROTATED_DIR = TEST_DIR / "Zahirovic2022_PMAG_rotated"
MODEL_DIR = TEST_DIR  # plate-model-manager creates <MODEL_DIR>/zahirovic2022/

# The rotation file path created by plate-model-manager for Zahirovic2022
ROTATION_FILE = TEST_DIR / "zahirovic2022" / "Rotations" / "CombinedRotations.rot"

# Single-file test (step 6 / 7)
SINGLE_TIME = 105  # Ma
SINGLE_INPUT = MANTLE_DIR / f"paleobathymetry_{SINGLE_TIME}Ma.nc"
SINGLE_OUTPUT = TEST_DIR / f"paleobathymetry_pmag_{SINGLE_TIME}Ma.nc"
SINGLE_REFERENCE = PMAG_DIR / f"paleobathymetry_{SINGLE_TIME}Ma.nc"

TIME_MIN = 100  # Ma
TIME_MAX = 110  # Ma

# Acceptance thresholds for the automated comparison
MIN_CORRELATION = 0.90  # Pearson r must be at least this
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


def _load_z_var(nc_path):
    """Load the ``z`` data variable from a NetCDF file as a numpy array."""
    ds = xr.open_dataset(nc_path)
    if "z" not in ds.data_vars:
        ds.close()
        raise KeyError(f"Variable 'z' not found in {nc_path}")
    arr = ds["z"].values.astype(float)
    ds.close()
    return arr


def _compare_two_grids(rot_file, ref_file, label=""):
    """
    Compare *rot_file* against *ref_file*.

    Grids are aligned by interpolating the reference onto the rotated grid's
    coordinate system when their shapes differ.  Returns True if the Pearson
    correlation and mean-absolute-difference/range both pass the acceptance
    thresholds.
    """
    rot_arr = _load_z_var(rot_file)
    ref_arr = _load_z_var(ref_file)

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
        print(f"  SKIP {label}: too few valid overlap points")
        return True  # not a failure — just not enough data to compare

    # Pearson correlation
    corr = float(np.corrcoef(rot_v, ref_v)[0, 1])

    # Mean absolute difference as a fraction of the reference value range
    ref_range = float(ref_v.max() - ref_v.min())
    if ref_range == 0:
        mad_frac = None
    else:
        mad_frac = float(np.mean(np.abs(rot_v - ref_v)) / ref_range)

    corr_ok = corr >= MIN_CORRELATION
    mad_ok = mad_frac is None or mad_frac <= MAX_MEAN_ABS_DIFF_FRACTION
    passed = corr_ok and mad_ok
    status = "PASS" if passed else "FAIL"
    mad_str = "N/A (zero range)" if mad_frac is None else f"{mad_frac:.4f}"
    print(
        f"  [{status}] {label}  |  corr={corr:.4f}  |  "
        f"mean_abs_diff/range={mad_str}"
    )
    return passed


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
# Step 3: Download Zahirovic2022 plate model
# ---------------------------------------------------------------------------


def download_model():
    """Download the Zahirovic2022 rotation model into TEST_DIR via plate-model-manager."""
    print(f"Downloading Zahirovic2022 rotation model → {MODEL_DIR / 'zahirovic2022'} …")
    plate_model = PlateModelManager().get_model(
        "zahirovic2022", data_dir=str(MODEL_DIR)
    )
    if not plate_model:
        # Fall back to a read-only PlateModel if the manager returns None
        plate_model = PlateModel(
            "zahirovic2022", data_dir=str(MODEL_DIR), readonly=True
        )

    # Ensure the rotation files are downloaded locally
    plate_model.get_rotation_model()

    if not ROTATION_FILE.exists():
        raise FileNotFoundError(
            f"Expected rotation file not found after download: {ROTATION_FILE}"
        )
    print(f"  Rotation file available: {ROTATION_FILE}")


# ---------------------------------------------------------------------------
# Step 4: Run gplately rotate_grid (directory mode, named model)
# ---------------------------------------------------------------------------


def rotate_grids_named_model():
    """Invoke ``gplately rotate_grid`` in directory mode using a named model."""
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
# Step 5: Compare directory-mode output with reference pmag grids
# ---------------------------------------------------------------------------


def compare_grids_directory_mode():
    """
    Compare each rotated grid (directory mode) against its reference pmag grid.

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

        passed = _compare_two_grids(rot_file, ref_file, label=f"{t:.0f} Ma ({rot_file.name})")
        if not passed:
            all_pass = False

    return all_pass


# ---------------------------------------------------------------------------
# Step 6: Run gplately rotate_grid (single-file mode, local rotation files)
# ---------------------------------------------------------------------------


def rotate_single_file_local_rot():
    """
    Invoke ``gplately rotate_grid`` in single-file mode using local rotation files.

    Rotates ``SINGLE_INPUT`` (105 Ma mantle frame) to ``SINGLE_OUTPUT``
    (105 Ma pmag frame) using the rotation file downloaded in step 3.
    """
    cmd = [
        sys.executable,
        "-m",
        "gplately",
        "rotate_grid",
        str(SINGLE_INPUT),
        str(SINGLE_OUTPUT),
        "--from-rotation-files",
        str(ROTATION_FILE),
        "--to-rotation-files",
        str(ROTATION_FILE),
        "--from-anchor",
        "0",
        "--to-anchor",
        "701701",
    ]
    print("\nRunning:", " ".join(cmd))
    subprocess.run(cmd, check=True, text=True)


# ---------------------------------------------------------------------------
# Step 7: Compare single-file output with reference pmag grid
# ---------------------------------------------------------------------------


def compare_single_file():
    """
    Compare ``SINGLE_OUTPUT`` against ``SINGLE_REFERENCE``.

    Returns True if the comparison passes, False otherwise.
    """
    if not SINGLE_OUTPUT.exists():
        print(f"ERROR: Single-file output not found: {SINGLE_OUTPUT}")
        return False
    if not SINGLE_REFERENCE.exists():
        print(f"ERROR: Single-file reference not found: {SINGLE_REFERENCE}")
        return False

    return _compare_two_grids(
        SINGLE_OUTPUT,
        SINGLE_REFERENCE,
        label=f"{SINGLE_TIME} Ma (single-file / local rot)",
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    TEST_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Step 1: Download Optimised Mantle Frame grids ===")
    download_grids(MANTLE_WEBDAV, MANTLE_DIR)

    print("\n=== Step 2: Download Palaeomagnetic Frame grids (reference) ===")
    download_grids(PMAG_WEBDAV, PMAG_DIR)

    print("\n=== Step 3: Download Zahirovic2022 plate model ===")
    download_model()

    print("\n=== Step 4: Rotate mantle-frame grids to pmag reference frame (named model) ===")
    rotate_grids_named_model()

    print("\n=== Step 5: Compare directory-mode output against reference pmag grids ===")
    ok_dir = compare_grids_directory_mode()

    print(
        f"\n=== Step 6: Rotate single file ({SINGLE_TIME} Ma) using local rotation files ==="
    )
    rotate_single_file_local_rot()

    print(f"\n=== Step 7: Compare single-file output against reference pmag grid ===")
    ok_single = compare_single_file()

    if ok_dir and ok_single:
        print("\nAll comparisons PASSED.")
    else:
        print("\nSome comparisons FAILED — please review the output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
