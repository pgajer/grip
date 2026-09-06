# Manuscript source lineage

- `migration.json`: every file imported from geodesicMDS commit `a09f7991488d0b9f7512b083c42de9f314377bbc` and its original SHA-256.
- `style-revision.json`: the two as-reviewed style source hashes, adopted in a separate grip commit.
- `review-snapshots.json`: preserved circulation payload hashes for five earlier editions.

The audited import is grip commit `bd79725`; the style adoption is `c1c32a1`.
Later commits correct the reviewed source and establish portable packaging.
Original geodesicMDS refs and a verified full Git bundle preserve earlier history;
this import does not claim those commits as grip history. Frozen experiment
provenance in `../evidence/` is unchanged and continues to identify the original
fitting scripts and commits. It must not be rewritten as if the experiments ran
in the migrated directory.

Only `../geodesic_mds.tex` and `../geodesic_mds.bib` are active manuscript masters.
