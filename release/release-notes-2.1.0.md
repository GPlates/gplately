## GPlately 2.1.0 Release Notes (TBD)

> This file is a living document and work in progress. The checklist below tracks what needs
> to be done before the 2.1.0 release can be announced. Items marked ✅ are complete.

---

### Pre-release Checklist

#### Code / Feature Work

- [x] Raster reconstruction performance tuning ([#419](https://github.com/GPlates/gplately/pull/419))
- [x] Add shallow copy function to `PlateReconstruction` class ([#419](https://github.com/GPlates/gplately/pull/419))
- [x] New function to create a `Raster` object from scattered points ([#419](https://github.com/GPlates/gplately/pull/419))
- [x] Enhance `Raster.copy()` function ([#419](https://github.com/GPlates/gplately/pull/419))
- [x] Enhance how `fill_value` is handled in `Raster` ([#419](https://github.com/GPlates/gplately/pull/419))
- [ ] Improve seed point reconstruction in seafloor gridding ([#413](https://github.com/GPlates/gplately/pull/413)) — milestone: release 2.1 (still open/draft; depends on pyGPlates ≥ 1.1)
- [ ] Deform continental polygons in seafloor gridding ([#414](https://github.com/GPlates/gplately/pull/414)) — milestone: release 2.1 (still open/draft; depends on #413)
- [ ] Support regional raster reconstruction ([#417](https://github.com/GPlates/gplately/issues/417)) — several open design questions need to be answered
- [ ] Create a new `plot` sub-folder and reorganise plotting code ([#409](https://github.com/GPlates/gplately/issues/409)) — decision needed; backward compatibility via re-exports required
- [ ] Fix `PlateModel.__getstate__` `KeyError: 'executor'` during pickling ([#386](https://github.com/GPlates/gplately/issues/386)) — simple `pop(…, None)` fix in `plate_model_manager`; decide whether to fix upstream or work around in gplately
- [ ] Fix runtime error for incomplete pickle support ([#366](https://github.com/GPlates/gplately/issues/366))
- [ ] Document return-tuple order of `Points.reconstruct(time, return_array=True)` ([#398](https://github.com/GPlates/gplately/issues/398))
- [ ] Clarify `get_layer` vs `get_raster` for time-dependent rasters ([#397](https://github.com/GPlates/gplately/issues/397))
- [ ] Update example notebooks in `Notebooks/Examples` ([#420](https://github.com/GPlates/gplately/issues/420) / [#421](https://github.com/GPlates/gplately/pull/421))
- [ ] Fix misspellings across the codebase ([#394](https://github.com/GPlates/gplately/issues/394))
- [ ] Add CONTRIBUTING and CODE_OF_CONDUCT files ([#395](https://github.com/GPlates/gplately/issues/395)) — optional for this release

#### Documentation

- [ ] Write/update the 2.1.0 release notes in this file (new features, bug fixes, breaking changes, performance improvements)
- [ ] Create the 2.1.0 documentation folder in the `gh-pages` branch
- [ ] Check the online documentation for errors or broken links
- [ ] Update the plate models list documentation if any new models have been added
- [ ] Check the [PyPI page](https://pypi.org/project/gplately/) — especially images and links

#### Testing & Validation

- [ ] Run `pytest` against all test cases
- [ ] Make sure all Jupyter notebooks run without errors
- [ ] Make sure all examples run without errors
- [ ] Check pip installation on Windows, macOS and Ubuntu
- [ ] Check conda installation on Windows, macOS and Ubuntu
- [ ] Check all `gplately` CLI commands work as expected
- [ ] Make sure the Docker images work correctly

#### Release Mechanics

- [ ] Set `USING_DEV_VERSION = True` during development (in `gplately/__init__.py`)
- [ ] Set `USING_DEV_VERSION = False` before the official release
- [ ] Make sure DEV warnings are **not** showing, especially the large warning banner
- [ ] Verify the version number is correct (tag / `pyproject.toml` / `__init__.py`)
- [ ] Update `docker/build-docker-image-version.txt` to `v2.1.0`
- [ ] Tag the release commit as `v2.1.0`

---

### New Features (draft — to be expanded)

* Raster reconstruction performance improvements.
* New `PlateReconstruction.shallow_copy()` method.
* New `Raster.from_points()` constructor (create a raster from scattered data points).
* Enhanced `Raster.copy()` and `fill_value` handling.
* *(pending)* Improved seafloor seed-point reconstruction algorithm — naturally deactivates seed points at convergent boundaries without empirical velocity/distance thresholds ([#413](https://github.com/GPlates/gplately/pull/413)).
* *(pending)* Continental polygon deformation support in seafloor gridding ([#414](https://github.com/GPlates/gplately/pull/414)).

### Bug Fixes (draft — to be expanded)

* *(pending)* Fix `PlateModel.__getstate__` `KeyError: 'executor'` during multiprocessing pickling ([#386](https://github.com/GPlates/gplately/issues/386)).
* *(pending)* Fix incomplete pickle support runtime error ([#366](https://github.com/GPlates/gplately/issues/366)).

### Other Changes

See [all changes](https://github.com/GPlates/gplately/compare/v2.0.0...v2.1.0) since GPlately 2.0.0.
