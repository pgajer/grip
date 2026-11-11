# Getting-started illustration

`layout.rds` preserves the original `dwt_307` graph, fitted coordinates,
metadata, warnings and recipe used by `vignettes/grip-examples.Rmd`.
The fit uses full SGD, three dimensions, seed 1 and the default 1,000-pass
budget. It is a teaching illustration, not a benchmark or a convergence claim.
Unused preparation matrices are omitted; this is not a complete fit object.

After installing the current source checkout, run from its root:

```sh
Rscript tools/pkg/build-starting-layout.R
```

The generator evaluates the vignette's graph and fitting code and records
source checksums, the installed native library checksum and session details.
The guide validates the graph and rebuilds its figures and scores. The full
fitting recipe remains visible and runnable; the shorter 80-pass fit is still
computed when the guide is built.
