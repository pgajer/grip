# Targeted follow-up diagnostics for the geodesic MDS report

These checks reproduce the numerical follow-ups suggested by the September 5,
2026 audit. They consume the saved experiments and do not overwrite their
distance matrices, selected fits, or original result tables. No private audit
files are required to build or run them.

```sh
make -C tools/experiments/mds-audit-diagnostics run
make -C tools/experiments/mds-audit-diagnostics verify
```

Outputs are generated under `output/mds-audit-diagnostics/`:

- `optimizer_spread.csv`: best/worst stress and relative spread for the six
  original 3D starts at every radius, surface, and sampling measure (40 rows).
- `additional_paraboloid_starts.csv`: 16 additional fits at radius 64, using
  seed 640905. For each measure: the selected fit as a restart, perturbations
  of sizes 0.005, 0.02, and 0.1 in normalized coordinates, and four new random
  starts. These runs reuse the original optimization settings.
- `planar_instability.csv`: the smallest third-coordinate Hessian eigenvalue
  for each saved saddle 2D fit at radius 64, actual stress decrease at a
  perturbation of size 0.001, and a finite-difference curvature check at 0.00001.
  Perturbation vectors have unit Euclidean norm. The objective is mean pair
  stress after dividing targets and coordinates by RMS input distance.
- `population_spectrum.csv`: leading and second logarithmic-mean kernel
  eigenvalues and competing angular/radial sectors at 256 and 512 quadrature
  nodes, plus the norm of the projected radial logarithmic term.
- `additional_saddle_pairs.csv`: 48 independent length comparisons from the
  actual radius-64 matrices, 24 per sampling measure. Each set contains the
  four shortest, four longest, four pairs with an endpoint nearest zero
  height, and twelve random pairs (seed 20260906, shared RNG across measures).
- `saddle_bvp_attempts.csv`: all 56 collocation attempts, including the eight
  initial attempts exceeding the mesh limit. All pairs completed on retry.
- `refined_saddle_pairs.csv`: the three most discrepant pairs from those
  comparisons, refined at three collocation tolerances through 1e-10.
- `manifest.json`: source and input/output SHA-256 checksums and build time.

The new geodesic comparisons integrate the smooth saddle ODE using SciPy
DOP853 and solve its endpoint problem by collocation. Pairwise relative
disagreement uses the saved pair's own length in the denominator; it is not
normalized by the matrix RMS. Maxima over these selected pairs are not
all-pairs accuracy bounds or interval certificates. A negative transverse
Hessian eigenvalue excludes the particular saved planar fit as a 3D local
minimum, not all conceivable planar configurations. Multiple starts and
quadrature refinement remain numerical evidence, not global optimality or
analytic mode-selection proofs.
