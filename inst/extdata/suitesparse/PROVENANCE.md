# Bundled SuiteSparse graphs

Release selection: 30_sept_2026. The package catalogue contains 85 graphs:
80 graphs selected in the comparison gallery and six original examples, with
HB/494_bus shared between the two collections. Gallery membership is frozen
in `selection-30_sept_2026.json`; mutable viewer favorites are not build inputs.

Use `suitesparse.graphs()` to browse and `suitesparse.graph("494_bus")` to load.
Only base R is needed to read the installed xz-compressed `graphs.rds`.
No embeddings, distance matrices, or experiment archives are bundled with it.
The entire small bundle is read on first access and cached privately; subsequent
calls return ordinary R objects with copy-on-modify behavior.

## Attribution and modifications

The [SuiteSparse Matrix Collection policy](https://sparse.tamu.edu/about),
checked on 30_sept_2026, distributes matrices and metadata under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), separately from grip's
software license. Credit the original contributors in each graph's
`provenance$source_header`. Those headers are preserved in full, including
source-specific notes and citations, for all 85 records. `provenance$name`
ends in `_unweighted_graph` to distinguish the derived object from its matrix.
Full collection IDs identify sources; component IDs additionally identify the
retained induced graph. The runtime graph is not the original numeric matrix.

Cite Davis, T. A. and Hu, Y. (2011), *The University of Florida Sparse Matrix
Collection*, ACM TOMS 38(1), Article 1,
[DOI](https://doi.org/10.1145/2049662.2049663).
The original six examples were selected from Zheng, Pawar and Goodman (2018),
[Graph Drawing by Stochastic Gradient Descent](https://doi.org/10.1109/TVCG.2018.2859997).
For Les Miserables also cite Knuth (1993), *The Stanford GraphBase*.

Square conversions sum duplicate entries, remove numerical zeros and diagonal
entries, and take the union of the remaining support with its transpose.
Rectangular conversions create distinct row and column vertices, joining each
nonzero row/column pair. Negative values contribute edges. All edges have unit
traversal length; signs, magnitudes and physical units are not distances.
Full graphs retain every source vertex, including isolates. Component graphs
retain the explicitly recorded largest component in original vertex order,
with one-based original indices. They never replace a full graph under its ID.

## Reproduction

From a grip source checkout, run:

```
python3 data-raw/build-suitesparse.py /path/to/suitesparse_embedding_comparison
```

Maintainer dependencies are Python with NumPy/SciPy and R with jsonlite. These
are not runtime dependencies. Obtain the archived acquisition directory from
the maintainer; the selection records relative paths and SHA-256 hashes.
Source URLs are stored in each graph. The builder never downloads or fits
embeddings. It checks graph JSON, source archives, component mappings, original
graphs and the archived six matrices against frozen checksums. It independently
reconstructs numerical support from Matrix Market, verifies every selected
edge and vertex count, then writes a deterministic R version-2 xz object and
checks an exact round trip. The original six metadata are reconstructed by
`data-raw/original_six.R`; their matrices and Les Miserables names remain in
`inst/extdata/zheng-graphs/`.

Only use `--freeze` to deliberately replace the release selection after user
curation. Review and version the resulting manifest before release. Normal
builds ignore the current viewer manifest. Neither operation edits the viewer,
archives, favorites or saved embedding runs.
