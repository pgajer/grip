# Fixed SGD gallery illustrations

`layouts.rds` contains the six graph inputs, display coordinates, fitting
metadata, warnings and reproduction recipes for `vignettes/gallery.Rmd`.
Every fit uses the original 1,000-epoch budget, seeds 1 through 6, three
dimensions, full SGD and uniform pair weights. The figures and stopping tables
are illustrations, not a benchmark or a claim of convergence.

The fitted objects are reduced to coordinates and metadata; preparation
matrices and other fields unused by the guide are omitted. No coordinates,
graph inputs, stopping records or seeds are reduced or resampled.

From the package source root, install the current source and run:

```sh
Rscript tools/pkg/build-sgd-gallery.R
```

The generator evaluates the graph and fitting chunks in the vignette itself.
It records the supplied source checkout's checksums, the installed native
library checksum, and the R session. Rendering validates the graph inputs
against the saved inputs and rebuilds the figures and tables. Full fitting
recipes are shown in the guide and remain available for optional reproduction.
