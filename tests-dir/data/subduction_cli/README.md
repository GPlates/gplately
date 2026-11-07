# Synthetic subduction CLI fixture

These original synthetic files were generated for the GPlately CLI contribution
using public pyGPlates 1.0.0 constructors. They use the repository's GPL-2.0-only
license and contain no third-party plate-model data.

Two closed topological plates share a northward trench at longitude 0 degrees,
from latitude -20 to 20 degrees. Its polarity is Left. The trench and western
plate use plate 101; the eastern/subducting plate is 201. Rotations are relative
to plate 0 about the North Pole: at 0 Ma both are identity; at 10 Ma plate 101
has 10 degrees and plate 201 has -20 degrees. Features are valid from 10 Ma
through 0 Ma.

All boundary geometries reconstruct with plate 101, keeping their shared endpoints
coincident over the test interval. The eastern topological region retains plate
201, supplying the distinct plate velocity. This avoids artificial gaps caused
by moving its closure independently of the shared trench.

This is an integration fixture with nonzero velocities, not a model of Earth or
an independent scientific-accuracy benchmark. Tests compare CLI routing and
serialization with the unchanged public API; modern and legacy algorithm values
are not required to agree.

Feature IDs were generated once and the serialized files are fixed. Tests read
these artifacts; they do not regenerate IDs or download a model.

- topologies.gpml (14970 bytes): SHA256 `ba4425f4f683bdf46de5bb259e11d3d0b57d4da1ff44e150d275fe425f0ce513`
- rotations.rot (148 bytes): SHA256 `2c76cc8137d9b83388dd8f30aacf307d97041fd654accd3e3cf3b4197a1d7676`

Construction API: https://www.gplates.org/docs/pygplates/sample-code/pygplates_create_topological_features
