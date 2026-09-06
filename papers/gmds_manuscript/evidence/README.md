# Focused-paper evidence

The namespaces `study-initializers` and `study-radius` contain byte-preserving
exports from the grip repository. Their original study documentation is
retained as provenance; links within those documents refer to the original
experiment trees, not to a complete mirror here.

`source-manifest.json` records the full source commit and source/export SHA-256
hashes. CSV files are deterministically compressed, without changing their
contents. Normal manuscript rendering uses these local exports and does not
require the grip repository or private planning files.

To verify the exports from the manuscript directory:

```sh
python3 scripts/verify_focused_evidence.py
```

To intentionally refresh them from a specified committed grip checkout:

```sh
python3 scripts/export_focused_evidence.py --grip-source /path/to/grip
```

The importer rejects selected source files that differ from that checkout's
HEAD. Refreshing changes provenance and requires renewed comparison review.
Full fitting scripts and larger input/coordinate archives remain canonical in
the source experiment directories recorded by the manifest.

The two studies have different domains, samples, targets and scale policies.
Their settings must not be pooled as independent sample replications. Existing
verification records document the original analyses; the focused manuscript's
additional checks are recorded separately.

## Additional diagnostic evidence

`coordinate-components/coordinates.csv.gz` contains the 24 highlighted fits
and six corresponding truth arrays, exported from twelve pinned grip Git
blobs. Its manifest records their source hashes and the derived export hash.
`check_coordinate_components.py` reconstructs one global similarity alignment
per fit, component errors without realignment, and the original ambient saddle
edge-length Jacobian ranks. Both `make verify` and `make pdf` run this portable
check. The rank threshold is `max(J.shape) * machine_epsilon * largest_singular_value`;
these are numerical local-rigidity diagnostics, not global-rigidity certificates.
Explicit refresh requires R and the source repository:

```sh
python3 scripts/export_coordinate_components.py --grip-source /path/to/grip
```

`retained-route-validation/` records an optional all-route recalculation for
24 highlighted configurations, independently summing all 28,680 retained
routes per configuration in ordinary R. The script compares cached coordinates
with the portable source arrays and checks both path diagnostics. The manifest
identifies the exact source cache/input files by SHA-256. Rerunning it requires
the original saved fitting caches, plus R and jsonlite; ordinary manuscript
builds check its frozen record against the score table and do not require those
caches or claim to rerun this validation:

```sh
Rscript scripts/check_retained_routes.R /path/to/grip evidence/retained-route-validation
```

`continuation-comparison.csv` is regenerated from paired frozen score rows.
Differences are density continuation minus uniform-only, at identical graph,
initializer and scale policy. Quantiles use pandas' linear interpolation;
changes are in percentage points, with counts using a 1e-8 absolute tolerance
in relative-error units. `numerical-claims.json` also records within-case ranges
across all seven k values for the 84 full-geodesic-initializer cases.
