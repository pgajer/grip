# Two-fidelity saddle pilot

This standalone five-cloud experiment studies both stages of `X -> G -> Z`.
It does not modify the manuscript or replace the previous two-dimensional
flattening experiment. Generated data, figures, and numeric summaries are in
`../../../build/two-fidelity-pilot/` relative to this directory. Use a new output
directory if changing any scientific protocol setting; saved results are reused
by their recorded cloud, grid, source-count, and neighborhood identifiers.

## Sampling and graph construction

Each cloud contains 1,000 independent observations on
`(x,y,0.8*(x^2-y^2))`, with `(x,y)` in `[-1,1]^2`, without measurement noise.
Rejection sampling gives uniform **surface area**, not uniform parameter area:
accept a uniform-square proposal with probability
`sqrt(1+2.56*(x^2+y^2))/sqrt(6.12)`.
R seeds are `2211000 + replicate`, with Mersenne-Twister, Inversion normals,
and Rejection sampling explicitly selected. Reference-source seeds are
`3211000 + replicate`; 128 vertices are chosen without replacement, independently
of graph construction and layout scores.

Graphs use `dgraphs::create.sknn.graph()` with exact neighbors, symmetric union,
Euclidean edge lengths in the original 3D observations, no pruning, and exact
component-MST repair. The number and lengths of added bridges are saved.
Graph shortest paths for calibration are computed with `igraph::distances()`.
The initial integer sweep `3:20` had boundary minima, as did its extension to
`3:40`; the pilot therefore expanded to `3:80` for every cloud. This adaptive
range finding is part of the pilot, not a predeclared confirmatory study.

## Numerical reference for smooth-surface distances

The reference domain is the **bounded saddle patch**, not the entire unbounded
surface. Every original observation is included as a vertex, together with a
regular parameter grid. Delaunay triangulation is performed in `(x,y)`, then
vertices are lifted to the saddle. Endpoints are not snapped to nearby grid
vertices. Paths can cross triangle interiors and cannot leave the patch.

The mesh solver is `pygeodesic` 0.1.11, an interface to Kirsanov's implementation
of the exact polyhedral geodesic algorithm [@pygeodesic]. "Exact" concerns
the triangulated surface in floating-point arithmetic; it does not make these
distances exact on the smooth saddle.

The primary reference uses an 81-by-81 grid plus the 1,000 observations.
Distances from 128 random source vertices to all observations cover 119,744
distinct unordered pairs, after removing self-pairs and duplicates. All five
clouds also use a 41-by-41 grid with the same 128 sources and a 161-by-161 grid
with 16 of those sources. The most resolution-sensitive cloud receives a further
161-by-161 check on all 128 sources. These comparisons assess numerical
sensitivity rather than provide a certified uniform error bound.

Independent controls solve the smooth-surface geodesic boundary-value equation
for 128 source-target pairs per cloud. For `f(x,y)=0.8*(x^2-y^2)`, its affine
parameterization satisfies

```
u'' = -grad(f) * (u'^T Hess(f) u') / (1 + |grad(f)|^2).
```

Here `grad(f)=(1.6*x,-1.6*y)` and `Hess(f)=diag(1.6,-1.6)`.
`scipy.integrate.solve_bvp()` solves the endpoint problem with tolerance `1e-8`
and at most 10,000 collocation nodes [@scipy-bvp]. Curve length is integrated
on 1,001 evaluation points. Each check requires successful termination and
checks that the evaluated curve stays inside the parameter square. These are
independent numerical controls, not an interval-arithmetic proof. A separate
flat-surface test compares the polyhedral solver with exact Euclidean lengths.

## The two stages and the scale convention

For each cloud, choose the graph minimizing the **unprofiled** surface-to-graph
error on the reference pairs:

```
E_XG = sqrt(sum((d_G - d_X)^2) / sum(d_X^2)).
```

Graph and surface distances retain their common physical units. The full
error-versus-k curve, fitted-scale diagnostic, signed relative bias, fraction
of graph distances shorter than the reference, and 32/64/128-source sensitivity
are retained. The minimizing integer is an **oracle choice against a numerical
reference**. Shallow minima and reference/source-subset sensitivity prevent
interpreting that integer as a uniquely established population optimum.

On the selected graph, all three configurations have three coordinate columns:

1. Original saddle observations (known edge/path-fidelity control).
2. `grip::metric.mds(dim=3)`, classical scaling of graph distances.
3. That same MDS result followed by `grip::edge.kk(dim=3)`.

The primary refinement uses the package's density-based stiffness continuation
with mixing coefficients `0, .25, .5, .75, 1`, at most 200 iterations per stage,
profiled scale, and zero edge-length stabilizer. A separate sensitivity run
continues its result for at most 1,000 additional iterations at mixing
coefficient 1. The additional run never replaces the primary comparison.
There is no direct optimization of the fixed-path GMDS objective in this pilot.

`grip::prepare.geodesic.kk(tie_mode="single")` retains one input route for each
of all 499,500 unordered vertex pairs. All three candidates use these same
routes. Their embedded lengths are sums of Euclidean edge lengths, without
recomputing shortest paths after embedding.

For observed lengths `y` and graph lengths `g`, set `a=sum(y*g)/sum(g^2)`.
The edge and fixed-path relative RMSE statistics are
`sqrt(sum((y-a*g)^2)/sum((a*g)^2))`, with separate scales. MDS Stress-1 uses
endpoint chords and denominator `sum(chords^2)` instead. An independent C++
route sum and a direct R route sample are checked against `score.gmds()`.
The package's near-tie tolerance can create tiny discrepancies against
igraph's strict shortest-path distances; these are checked and recorded.

For an additional **end-to-end** comparison against `D_X`, convert each
configuration into input units using its edge scale (estimated only from
graph edges), then compare its fixed-path lengths with `D_X` without fitting
another scale against the surface. A Procrustes error measures coordinate
shape discrepancy, distinct from all distance-fidelity statistics.

The original 3D coordinates preserve every Euclidean edge and retained graph
path, including MST bridges. This does not make chord stress zero or establish
surface-geodesic accuracy. Edge/path fidelity also need not identify the
original coordinate shape uniquely.

## Outputs and interpretation

`results.md`, numeric CSVs, and five PDF/PNG figures are generated from saved
objects, never by rerunning layouts during plotting. The representative cloud
is the middle cloud when ordered by its minimum surface-to-graph error;
ties are broken by replicate number. Layouts are similarity-aligned only for
display and the separately named Procrustes statistic. All fidelity scores
use unprojected coordinates.

`pilot-summary/` preserves the compact numeric readout and input checksums in
Git. Larger RDS/NPZ datasets, the isolated Python environment, and all PDF/PNG
figures remain in `build/`. `record-results.R` refreshes only this declared
summary snapshot after validation and plotting.

Five clouds characterize feasibility and pilot behavior, not population
precision. Finite iteration-limited fits are retained; termination before a
budget is exhausted is not automatically called convergence. Reference
precision, source-pair sampling, neighborhood selection, optimization budget,
and sampling variation are separate sources of uncertainty.

## Reproduction

From the repository root, with R packages grip 0.2.0, dgraphs 0.2.0, igraph,
and Rcpp installed, run the following. `GRIP_RJOURNAL_PACKAGE_LIBRARY` may point
to an isolated R library. Python dependencies are pinned in `requirements.txt`.

```sh
pilot_src=papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot
pilot_out=papers/grip-software-paper/build/two-fidelity-pilot
python3.12 -m venv "$pilot_out/venv"
"$pilot_out/venv/bin/python" -m pip install -r "$pilot_src/requirements.txt"
Rscript "$pilot_src/sample-clouds.R" "$pilot_out"
"$pilot_out/venv/bin/python" "$pilot_src/surface-reference.py" "$pilot_out" --checks
"$pilot_out/venv/bin/python" "$pilot_src/surface-reference.py" "$pilot_out" --pilot
Rscript "$pilot_src/graph-sweep.R" "$pilot_out" 1,2,3,4,5 80
Rscript "$pilot_src/calibrate.R" "$pilot_out"
Rscript "$pilot_src/fit-layouts.R" "$pilot_out"
"$pilot_out/venv/bin/python" "$pilot_src/surface-reference.py" "$pilot_out" --cloud 5 --sources 128 --resolutions 161
Rscript "$pilot_src/calibrate.R" "$pilot_out"
Rscript "$pilot_src/check-results.R" "$pilot_out"
Rscript "$pilot_src/plot-results.R" "$pilot_out"
Rscript "$pilot_src/record-results.R" "$pilot_out"
```

Saved R objects retain the actual R/package environment; the Python environment
and reference timing are recorded in JSON. Mesh timings measure construction,
exact distance solves, and result-array creation (not later CSV serialization).
Graph timings measure construction and all-pairs distances. Layout timings
separate preparation, MDS, edge-KK, and the additional-budget diagnostic.
Workers run independent reference clouds; worker elapsed time is not CPU time.

The sources below document solver capabilities, not the pilot's mathematical
derivations or empirical conclusions. BibTeX metadata and claim checks are in
`references.bib` and `citation_verification.html`.
