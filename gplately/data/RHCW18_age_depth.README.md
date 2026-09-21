# `RHCW18_age_depth.dat` — provenance and attribution

A two-column (age in Ma, depth in m) age-depth lookup table for the RHCW18 oceanic
lithosphere model, used by `gplately.age_to_basement_depth(model="rhcw18")`. It has no
closed-form expression, which is why it ships as data rather than as code.

## Source

<https://github.com/freddrichards/RHCW18_Plate_Model> — the "preferred parameters"
seafloor-depth-against-age output of that repository.

This copy is byte-identical (SHA-256 `086e1bf68f0d0b5a…`) to the one in EarthByte's
[simple_paleobathymetry](https://github.com/EarthByte/simple_paleobathymetry) workflow, which
GPlately's paleobathymetry support was ported from.

## Citation

The upstream repository asks that the following be cited when the data is used:

> Richards, F.D., M.J. Hoggard, L.R. Cowton & N.J. White (2018). Reassessing the thermal
> structure of oceanic lithosphere with revised global inventories of basement depths and
> heat flow measurements. *Journal of Geophysical Research: Solid Earth*, 123, 9136–9161.

> Richards, F.D., M.J. Hoggard, A. Crosby, S. Ghelichkhan & N.J. White (2020). Structure and
> dynamics of the oceanic lithosphere-asthenosphere system. *Physics of the Earth and
> Planetary Interiors*.

## Licence status — unresolved

**The upstream repository states no licence**, in a `LICENSE` file or in its README; it gives
only the citation request above. There is therefore no licence for GPlately to restate here,
and redistributing the file inside a GPL-2.0 package rests on nothing explicit.

This is recorded rather than resolved. Closing it needs a decision, not a paragraph:

- ask Richards et al. to add a licence upstream, or for permission to redistribute; or
- fetch the table at run time (e.g. through `plate_model_manager`) instead of vendoring it; or
- keep it as-is, having judged the citation request to be the intended terms.

Note that pyBacktrack, which also implements RHCW18, does **not** redistribute this file — it
has its own implementation — so there is no precedent to lean on there.
