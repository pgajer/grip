# Sparse SGD for graph layout

## Public interface

`metric.mds(approximation = "sparse", backend = "sgd")` dispatches to the
unexported `.sparse.metric.mds()` before dense preparation. `approximation =
"full"` remains the default; existing positional arguments and defaults are
preserved. Both modes default to `dim = 2`. Only omitted sparse arguments take
mode-specific defaults: random initialization, inverse-squared weighting,
30 epochs, and diagnostics disabled. Sparse mode rejects SMACOF, classical
starts, uniform weights, `n_init != 1`, non-NULL prepared objects, enabled dense
diagnostics, and explicitly supplied dense diagnostic controls. Pivot controls
are a strictly named `sparse_control` list (`n_pivots`, `pivots`); nonempty sparse
controls are rejected in full mode. Both results use method `metric_mds` and
record `metadata$engine` and `metadata$approximation`.

## Sparse contract

 Inputs are a connected simple undirected graph, positive finite edge
lengths (or unit lengths), and 2D or 3D output. Use the existing edge-only graph
preparation; never build all-pairs distances, path caches, or classical MDS.
Reject disconnected inputs; do not repair or silently select a component.

The reference is Zheng, Pawar and Goodman (2018), Algorithm 2, adapting Ortmann,
Klimenta and Brandes (2017). Pinned C++ reference: jxz12/s_gd2 revision
52ab0a5bee183b45eed60061cabe8e63f006e8a0, cpp/s_gd2/sparse.cpp. MIT notice in
inst/COPYRIGHTS. Our native implementation supports both 2D and 3D and uses
explicit mt19937_64 sampling, not upstream randomkit.

## Preparation and deterministic choices

Default 200 pivots, capped at n; explicit distinct 1-based `pivots` override
sampling (an explicitly supplied contradictory `n_pivots` is rejected).
First random pivot is uniform. Subsequent choices are proportional to distance
to the nearest selected pivot (not squared distance), excluding selected
vertices. Accumulate in vertex-index order, using a half-open uniform variate.
Graph lengths are divided by their maximum before native heap Dijkstra,
limiting avoidable unit-dependent threshold changes. Comparisons are strict in
these normalized units, without rounding/tolerance; floating-point borderline
ties on general weighted graphs can still depend on arithmetic. Targets are
restored to original units. Opposite pivot-pivot path distances can round
differently; merged pairs take the smaller of the two.
Ties in final nearest-pivot regions belong to the earliest selected pivot.
Explicit pivot order determines the same tie policy. Adjacency and constraints
are sorted by vertex indices. Dijkstra queue ties use vertex index.

Each pivot owns its nearest-pivot region R(p). For every nonadjacent vertex i,
count members j of R(p) with d(p,j) <= d(p,i)/2 (inclusive threshold).
The directed endpoint weight is count/d(p,i)^2; the reverse weight is zero
unless i is also a pivot. Self-pairs are omitted. Merge pivot-pivot pairs once.
All graph edges have weight 1/length^2 at both endpoints and are included once.
Their target is the supplied edge length, following upstream sparse.cpp.
Consequently, a redundant long weighted edge can differ from its shortest-path
distance: this is a local edge-length constraint, not an APSP claim. Reduction
to full shortest-path stress requires edges themselves to be shortest paths
(e.g. unit graphs or Euclidean chord-length graphs). Test a detour triangle.

Retain O(m+nh) constraints and O(nh) temporary pivot distances. Guard native
allocation estimates before allocation and report them, separately from real
peak process memory. No claim that a native allowance limits the entire R
process. Default allowance: 512 MiB, provisional calibration choice.

## Fitting and reporting

Normalize targets by RMS over retained unordered constraints; multiply
endpoint weights by RMS squared. Random initialization uses the unit cube in
these normalized units; supplied coordinates are divided by RMS. No initial
or final fitted scale. Only translate the returned coordinates to zero mean.

Use independent clipped endpoint updates min(eta*w_i,1) and min(eta*w_j,1),
with the shared half-distance correction from Algorithm 2. A zero-weight
endpoint does not move. Collision handling follows grip's native SGD kernel.
Shuffle lexicographically enumerated unique pairs each epoch using the existing
explicit Fisher-Yates policy. A public seed preserves the caller's R RNG;
record native preparation and fitting seeds, selected pivots, and settings.

Return the terminal iterate, not the best checkpoint of a scalar substitute.
Monitor the endpoint-average proxy sum((w_i+w_j)/2 * residual^2) and label it
`sparse_proxy_stress`. This equals full inverse-squared stress in the all-pivot
geodesic-edge limit. It is NOT equation (19) of Zheng's paper and the asymmetric
updates are not claimed to be its gradient. Report its normalized root error
and trace; do not put this number in `raw_stress` or label it full stress.
Default 30 epochs, existing hybrid schedule (.5, .01, switch .4), checkpoint
every epoch, all provisional calibration choices. This intentionally retains
grip's schedule rather than reproducing the paper's weight-dependent annealing.
Return iteration_limit, converged=FALSE. Full stress is only computed by
explicit external evaluation on small examples; sampled evaluation must use
independent pairs on large examples. No automatic dense diagnostics.

## Validation and demonstration

Independent R reference for distances, regions/counts, pair construction,
asymmetric updates and scores. Tests include endpoint immobility, pivot ties,
weighted detours, all-pivot update parity under identical rates/order/scale,
unit invariance, reproducibility and RNG preservation, disconnected inputs,
invalid controls, memory refusal, and absence of all-pairs preparation.
Independent implementation audit before publication.

A bounded saved vignette experiment compares full inverse-squared SGD, sparse
SGD (several pivot counts), and weighted GRIP using repeated seeds on quadform,
helix and graph examples. Include larger sparse-only/GRIP examples when dense
cost is inappropriate; record omissions explicitly. Record preparation and
fitting time, total time, and measured peak process RSS. Preserve failed and
time-limited runs. Plot every demonstrated graph/layout in 3D using ivue;
report full or separately labeled sampled error, never equate the proxy to it.


## Completed pilot and post-run diagnostic

The first bounded 60-fit experiment (three seeds, 100 SGD epochs, 8/32 small
pivots and 32/128/200 large pivots, 90-second workers, 1200-second total limit)
completed all fits. Full SGD was excluded above 512 vertices by design.
The independent audit identified large overall output scales for weighted GRIP.
A subsequent diagnostic fits one scalar to the independent assessment distances
for EACH method, reports both original and profiled relative error, and uses
that scalar for display before orthogonal Procrustes. Saved coordinates and
original errors are untouched. The publication bundle records this amendment;
it was not used to choose settings or rerun fitting. The assessment distances
for the 4096-vertex case are all targets from 32 separately seeded random source
vertices; these sources are shared across methods and independent of pivots.
The profiled score fits scale on that same sample, not a held-out scaling set.
