# Phase 3: radius, graph, initialization, and sampling study

Frozen September 5, 2026, before new radius-study fits. Main manuscript integration
remains separate. This protocol follows the completed five-cloud bounded-saddle
comparison, including its nonuniform initializer advantage and imperfect geometric
recovery. It does not select cases by favorable appearance.

## Design

Surfaces: paraboloid z=x²+y² and saddle z=x²−y², over the disk of radius r.
Sampling: uniform base-disk area and uniform surface area, using the exact radial
CDF with density proportional to rho sqrt(1+4 rho²) and independent uniform angle.
Radii: 1, 2, 4, 8, 16, 32, 64. Primary n=240, three independent samples
(seed 20260905, 20260906, 20260907). Quantiles and angles are coupled across radii
and surfaces. Replicate 1 reuses the audited original distance matrices.
A nested n=480 check at r=64 extends replicate 1 with an independent 240 points;
the original 240 observations retain their indices. Both sampling measures and
surfaces are included. This limited size check does not estimate a large-n limit.

Two graph regimes are kept distinct:

1. Geodesic: neighbor ranking and edge lengths use numerical smooth geodesics.
2. Ambient: neighbor ranking and edge lengths use original Euclidean chords.

Graphs use symmetric kNN union with ties broken by vertex index. The n=240 k grid
is 4, 8, 16, 32, 64, 128, 239. The n=480 check uses 8 (fixed k), 16 (same k/n as
n=240,k=8), 64 (fixed k), 128 (same k/n as n=240,k=64), and 479 (complete graph).
Record components before repair; disconnected cases receive explicitly labeled
MST augmentation under the same regime's edge lengths. Keep all vertices. Record
bridges, degree/edge-length quantiles, and parameter-plane edge separations.

## Fits and controls

All embeddings have three coordinates. Normalize targets and initial coordinates
by RMS full smooth-geodesic distance for each surface/sample/radius case; restore
physical units for saved coordinates and dimensional scores. This avoids the
r² change in numerical units altering the optimizer's step sizes across radius.

For every graph, fit classical MDS and metric stress MDS to strict graph distances,
then refine each with edge-KK. Stress MDS uses three starts (classical plus two
independent Gaussian configurations), 1,000 iterations each, tolerance 1e-8,
choosing the least achieved raw stress. Random starts are matched across k within
a case. Phase 2 already established why multiple starts and local-fit qualifiers
are needed; this study prioritizes sampling replication over six starts per graph.

In the geodesic regime, additionally refine a fixed full-geodesic classical-MDS
configuration, reused across k. This isolates the refinement graph from a change
in MDS input. The primary edge-KK schedule is density-to-uniform
(0,.25,.5,.75,1), 200 iterations per stage, profiled target scale, zero distance
stabilizer. Matched uniform-stiffness-only refinements use 1,000 iterations from
each MDS initializer. Preserve starts, selected objective, trace/stopping records,
and timings. Do not interpret iteration caps as convergence.

Original coordinates are scored on every graph. Ambient edge lengths are exactly
realized by this control; geodesic edge lengths generally are not. Complete graphs
provide an ambient realization control and a geodesic all-pair control. A separate
uniform-stiffness, identity-scale check verifies the raw all-pair edge objective
algebraically against the same chord-stress sum.

## Pilot gate and optimizer sensitivity

First run replicate 1 at r=1,8,64 and k=4,32,239 in both regimes and measures
(72 graphs). Check scores, conditioning, connectivity, timings, and complete-graph
controls before fitting the remaining grid. At r=64,k=8,32,128,239 on replicate 1,
compare unperturbed starts, isotropic perturbations of RMS amplitude 1e-4 in
normalized units, independent random starts, and original-coordinate starts.
For graph-classical and graph-stress primary fits at these settings, separately
continue edge-KK for 2,000 additional uniform steps. These exploratory controls
do not replace primary selected results or imply global optimality.

## Validation and interpretation

Use unchanged routes across candidates on a graph and all unordered vertex pairs.
Use strict graph distances for MDS, retained route lengths for path diagnostics;
record their maximum discrepancy. Compare independently recomputed edge/path
scores with grip. Record raw chord stress and profiled Stress-1 separately.
For graph/reference error, use physical distances with no reference-fitted scale.
For end-to-end path error, calibrate using graph edges only, never the surface
reference. Coordinate Procrustes error uses known correspondence. Save both
smaller singular-value ratios, unnormalized singular values, and r²-normalized
coordinates. Equal physical axis units are mandatory for geometric displays.

New geodesics use the audited direct smooth solvers. Check all-pair chord lower
bounds, explicit surface-path upper bounds, symmetry, and sampled triangle
inequalities. Validate additional difficult/random endpoint pairs independently
with tighter shooting and SciPy collocation/integration at selected large radii.
Sampled validation errors are not certified all-pair bounds. Record input/source
hashes, versions, seeds, and timing; checkpoints must not depend on private notes.

Interpret finite-radius, fixed-sample results as such. In particular, the
one-dimensional limiting saddle tree may have a three-dimensional Euclidean
embedding, and vanishing Gaussian curvature does not establish recovery or
flattening. The uniform edge-to-path bound may explain path agreement, but small
mean edge error does not supply its maximum-relative-error hypothesis. Preserve
negative results, graph effects, and optimizer sensitivity in the final report.
