# One saddle, four steps

Open `single-saddle-ivue.R` in R and run it section by section:

1. Sample `n` points uniformly in surface area on `z = C * (x^2 - y^2)`.
2. Construct symmetric kNN `graphs` for `ks = 3:20`, with Euclidean edge
   lengths and component-MST repair where needed.
3. Compute three-dimensional metric-MDS coordinates, `mds`, for every graph.
4. Refine each configuration with edge-KK, giving `mds.edge.kk`.

Change `n`, `seed`, `C`, `half.width`, `ks`, and `max.iter` at the top.
For a quick trial, use `n <- 80` and `max.iter <- 10`. The defaults fit all
18 graphs at n = 1,000, so allow time for the MDS and edge-KK loops.

Every list is keyed by the **character representation of k**, not its position:

```r
k <- 15
g <- graphs[[as.character(k)]]
Z <- mds.edge.kk[[as.character(k)]]
ivue::plot3D.graph(list(adj.list = g$adj_list, weight.list = g$weight_list),
                  X = Z, vertices = seq_len(n),
                  values = X[, "x"], scale = colors, edge.col = "gray80")
ivue::plot3D.cont(Z, values = X[, "x"], scale = colors)
```

Rerun only the viewing calls to change `k` or appearance. Colors retain the
original x coordinate, so the same observations have the same colors in every
view. Coordinates are not aligned or rescaled for display; rotating the widget
may be necessary to compare shapes. Graph views use the weighted adjacency lists returned by dgraphs.
When changing `ks`, also choose a displayed `k` that belongs to that range.

Sourcing the entire script **runs all four steps**, leaving the objects in the
R session. To display a widget after sourcing, run a plotting call as above
(or assign it to `w` and call `print(w)`). No files are saved automatically.

## Dependencies

Use `grip >= 0.2.0.9000`, `dgraphs >= 0.2.0`, and the current `ivue` API. No Python
reference solver or separately compiled path-checking code is needed.
If exploring with the local ivue checkout instead of an installed package,
load it before running the script:

```r
pkgload::load_all("~/current_projects/ivue", export_all = FALSE)
```

If grip and dgraphs are in a separate R library, add that directory with
`.libPaths()` before starting. The script does not install packages or change
the library search path.

## Relationship to the sampled-saddle experiment

The surface, area-uniform sampling rule, graph construction, and primary 3D
fitting settings follow the paper. The default seed is 1, giving a new sample.
With `n = 1000`, `seed = 2211005`, `C = 0.8`, and `half.width = 1`, sampling
reproduces the displayed cloud 5 using R's Mersenne-Twister generator.
The paper uses `[-1, 1]^2`; changing `half.width` to 2 explores a larger and
more steeply curved patch.

This is a coordinate-exploration script, **not** the full two-fidelity
experiment. It does not estimate surface geodesics, select k by surface-to-graph
error, compute fidelity scores, or run the additional refinement audit. In
particular, `3:20` is an exploratory range, not the calibrated optimum: the
paper extended the search to `3:80` and displayed k = 73 for cloud 5.

The full sampled-saddle workflow remains unchanged:

- `sample-clouds.R`, `surface-reference.py`, `graph-sweep.R`, `calibrate.R`,
  and `fit-layouts.R` generate the samples, references, graphs, and layouts.
- `export-paper-data.R` writes `../../precomputed/two-fidelity-saddle.rds`.
- `../../scripts/plot-two-fidelity.R`, function `plot_pilot_calibration_layouts()`,
  draws the publication figure; the manuscript supplies its caption.

`check-single-saddle.R` tests the four steps on a small sample, checks the
sampling against the published cloud, and compares the lightweight fitting
calls with the paper's full-path preparation. Run it with the same R libraries
used for the exploratory script.
