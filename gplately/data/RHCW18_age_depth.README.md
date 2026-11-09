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

## Licence status

**The upstream repository states no licence**, in a `LICENSE` file or in its README; it gives
only the citation request above. There is therefore no licence for GPlately to restate.

**Decision taken:** the citation request is treated as the intended terms. The data is
redistributed here on that basis, with the attribution above, and users of GPlately's
`rhcw18` model are asked to cite both papers.

Recorded so it is not re-litigated. It is worth revisiting if upstream later adds explicit
terms, or if a licence is ever needed rather than inferred — in which case the alternatives
are to ask Richards et al. directly, or to fetch the table at run time (e.g. through
`plate_model_manager`) rather than vendoring it. Note that pyBacktrack, which also implements
RHCW18, has its own implementation and does not redistribute this file, so there is no
precedent to lean on either way.
