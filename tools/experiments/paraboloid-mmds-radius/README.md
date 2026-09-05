# Euclidean-distance comparison on expanding paraboloid disks

**Corrected scope:** the requested dissimilarities are geodesic distances along
the smooth manifold. The main experiment and answer are now in
[GEODESIC.md](GEODESIC.md), with `make geodesic` as the reproduction target.
This file retains the ambient-Euclidean comparison. Its exact-reconstruction
statements do not apply to the requested geodesic distances.

For Euclidean input dissimilarities, a 3D MDS solution recovers the original
paraboloid exactly. Increasing its radius does not flatten this reconstruction.
After isotropic normalization by its overall size, the cloud approaches a line:
height grows as radius squared, while horizontal extent grows linearly.

This experiment treats `z = x² + y²` literally, with its coefficient fixed at
one in the chosen units. It uses complete Euclidean chord distances. The
geodesic extension below is a mathematical result, not a simulated result.

## Which MDS objective?

Fixed-target, unweighted raw metric stress is

$$S_r(Y)=\sum_{i<j}(\|Y_i-Y_j\|-\delta_{ij}^{(r)})^2.$$

The weighted stress formulation and the use of multiple starting configurations
are described by [Mair, Groenen, and de Leeuw (2022)](https://doi.org/10.18637/jss.v102.i10),
Section 2, Equation (1), and Section 3.1. Here the targets are fixed without
fitting a transformation, and all pair weights equal one. <!-- cite:mair2022 -->

An absolute error of one contributes one regardless of target distance.
For a relative error `e`, the contribution is `delta² e²`. The derivative with
respect to the fitted distance is `2(d - delta)`, and the second derivative is
two. Long targets do not have greater stiffness under equal absolute errors.
Multiplying every target and coordinate by `c` multiplies stress by `c²` and
preserves the minimizing shape up to scale.

Classical MDS instead truncates the positive spectrum of
`B = -J Delta² J / 2`, where the square is entrywise and
`J = I - 11'/n`. This minimizes the Gram-matrix strain `||B - YY'||_F²`
over centered configurations of the requested dimension. The truncation follows
by diagonalizing symmetric `B` and retaining its largest positive eigenvalues.
It is generally different from minimizing raw stress.
The [R documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cmdscale.html)
identifies `cmdscale` as classical MDS and explains its principal-component axes
and rigid-motion ambiguity. <!-- cite:rcmdscale -->

In grip 0.2.0, `R/gmds_layout_interface.R` defined `metric.mds()` using
`grip.classical.mds.embedding()` from `R/gmds_mds_compare.R`, which calls
`stats::cmdscale()` on graph shortest-path distances. From 0.2.0.9000 this
algorithm is named `classical.mds()`; `metric.mds()` instead wraps ratio SMACOF.
The present experiment
reimplements the classical spectral calculation on Euclidean distances and
also fits raw stress. It does not call that graph-specific wrapper or change
package behavior.

## Exact finite-radius result in three dimensions

Write `P_i(r) = (r u_i, r² q_i)`, where `u_i` lies in the unit disk and
`q_i = ||u_i||²`. If `delta_ij = ||P_i-P_j||`, setting `Y=P` gives zero
stress. Since stress is nonnegative, this is a global optimum. Conversely,
any zero-stress fit has the same distances and hence the same centered Gram
matrix `-J Delta² J/2 = P_c P_c'`. It is congruent to `P`, up to translations,
rotations and reflections. Classical MDS in three dimensions recovers this
same Gram matrix exactly.

In particular, replacing the height by `a r² q_i` while keeping the horizontal
coordinates fixed gives pair distances

$$\sqrt{r^2\|u_i-u_j\|^2+a^2r^4(q_i-q_j)^2}.$$

For any pair of unequal radii, decreasing `|a|` below one creates a nonzero
residual. The optimum remains `|a|=1`; the sign accounts for reflection.
The algorithm has no incentive to flatten the finite-radius paraboloid.
A 2D embedding, by contrast, is planar by construction and cannot preserve
all distances of a generic affinely three-dimensional sample.

## The normalized large-radius limit, directly from stress

Fix a finite set of unit-disk points with at least two distinct squared radii.
The rescaled distances are

$$D_{ij}^{(r)}=\frac{\delta_{ij}^{(r)}}{r^2}
=\sqrt{(q_i-q_j)^2+\frac{\|u_i-u_j\|^2}{r^2}}.$$

Therefore

$$0\le D_{ij}^{(r)}-|q_i-q_j|\le\frac{\|u_i-u_j\|}{r}\le\frac2r.$$

This is convergence to the distance matrix of the one-dimensional coordinates
`q_i`. For normalized embedding coordinates `V=Y/r²`,

$$\frac{S_r(Y)}{r^4}=\sum_{i<j}(\|V_i-V_j\|-D_{ij}^{(r)})^2.$$

Let `M=n(n-1)/2`. The line candidate `V_i=(q_i,0,...)` has normalized stress
at most `4M/r²`. Consequently any global minimizer in any fixed dimension
`k >= 1` has normalized stress tending to zero. Every pairwise distance of
such a minimizer converges to `|q_i-q_j|`, because each squared residual is
bounded by their vanishing sum. Its centered Gram matrix thus converges to
`q_c q_c'`, a matrix of rank one. All but the first squared singular value
vanish in normalized coordinates. This also holds for approximate minimizers
whose excess normalized stress tends to zero.

The normalized classical-MDS Gram matrix has the same rank-one limit, so its
positive spectral truncation also approaches a line. These statements concern
scale-normalized geometry. Unscaled coordinates can diverge, and absolute
distance error need not decrease. There is no uniform probability distribution
over the entire infinite-area paraboloid; the experiment concerns a family of
finite-radius distributions and its scaled limit.

For uniform sampling in the base disk the population covariance is

$$\operatorname{Cov}(P)=\operatorname{diag}(r^2/4,r^2/4,r^4/12).$$

Indeed, `q` is uniform on `[0,1]`, its variance is `1/12`, horizontal
variances are `1/4` before scaling, and cross-covariances vanish by symmetry.
For `r >= sqrt(3)`, the first principal component accounts for
`r²/(r²+6)` of the variance, approaching one. At small radii the horizontal
components dominate; the global cloud is then approximately a plane. The
change is therefore not monotonic flattening into a plane.

### If the input distances are continuous surface geodesics

Let `rho_i = r ||u_i||`. Every surface path has length at least
`|rho_i²-rho_j²|`. For an upper bound, travel radially from the larger radius
to the smaller one, then along the shorter circular arc at that radius.
The radial metric is `sqrt(1+4rho²) d rho`, and
`sqrt(1+4rho²) <= 1+2rho`. Hence

$$r^2|q_i-q_j|\le g_{ij}^{(r)}
\le r^2|q_i-q_j|+r(1+\pi).$$

Thus `g_ij/r²` has the same line-distance limit, and the stress argument
applies with `2` replaced by `1+pi`. However, at finite radii the original
3D coordinates generally do not realize geodesic distances with zero stress;
their MDS geometry can differ substantially from the Euclidean experiment.
The bound concerns continuous geodesics. Sparse graph shortest paths require
their own approximation analysis and are not automatically covered.

## Experiment and observed results

The seed is `20260905`; sample size is 240. Radii are
`0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64`. For each sampling measure, the same
uniform radial quantiles and angles are reused across radii. Base-disk sampling
uses `rho = r sqrt(U)`. Surface-area sampling uses the exact inverse CDF

$$\rho=\frac12\sqrt{\left[1+U\{(1+4r^2)^{3/2}-1\}\right]^{2/3}-1}.$$

This follows by integrating the area element
`rho sqrt(1+4rho²) d rho d theta`. In the surface-area family,
`rho²/r² -> U^(2/3)`, so the limiting radial distribution changes but the
normalized distance geometry still approaches a line.

For each radius and sampling measure, the experiment computes classical MDS in
two and three dimensions. Raw stress in two dimensions is minimized using
L-BFGS-B with an analytic gradient, from six starts: classical MDS, the three
coordinate-plane projections, and two random configurations. The best stress
is retained. The largest permitted iteration count is 4000; the objective and
gradient stopping tolerances are `1e-14` and `1e-10`. Targets are divided by
their RMS pair distance for numerical conditioning; this common rescaling
does not change the optimal shape. Optimization is local; six starts do not
certify a global minimum. The global-limit theorem is independent of optimizer
convergence.

The relative distance RMSE is
`sqrt(sum((d-delta)²)/sum(delta²))`; its denominator uses input distances.
The first-PC variance fraction is `s1²/sum(sj²)` for centered embedding
singular values. Values below are deterministic results of this run, not
confidence intervals or averages over repeated independent datasets.

| Radius | 2D raw-stress relative RMSE | 2D first-PC variance fraction |
|---:|---:|---:|
| 0.25 | 0.007528 | 0.565552 |
| 1 | 0.093869 | 0.563869 |
| 2 | 0.197989 | 0.540380 |
| 4 | 0.114124 | 0.798838 |
| 16 | 0.020076 | 0.983305 |
| 64 | 0.002611 | 0.998930 |

These rows use base-disk sampling. At radius 64, surface-area sampling gives
relative RMSE 0.003688 and first-PC variance fraction 0.998431. The full sweep
for both sampling measures is plotted separately in `sampling_comparison.png`.
Three-dimensional classical MDS reconstructs every sample up to rigid motion:
the maximum relative coordinate discrepancy over all 20 cases is
`4.83e-15` (rounded upward). Its coordinates are also global raw-stress
minimizers, so a separate numerical 3D stress optimization is unnecessary.

At radius 64, simply deleting height and optimally rescaling the original
`xy` coordinates gives relative RMSE 0.6415 under base-disk sampling. This
particular flat-disk candidate preserves distances much worse than the
near-line solution. It is not a bound on all planar configurations.

The selected 2D fits illustrate why relative and absolute error must be
distinguished: between radii 4 and 64, relative RMSE decreases from 0.1141 to
0.002611, but absolute pair-distance RMSE increases from 0.8828 to 4.4392.

## Reproduction and checks

From this directory run `make run` followed by `make verify`.
Python >=3.9, NumPy, SciPy and Matplotlib are required. The reference run used
Python 3.9.6, NumPy 1.24.4, SciPy 1.9.1 and Matplotlib 3.7.5.

Generated output is in `output/paraboloid-mmds-radius/` relative to the
repository root. It is ignored by Git and can be regenerated from these
sources. The output bundle includes:

- `embedding_snapshots.png/.pdf`: matched 3D and 2D views at four stated radii,
  scaled by RMS input pair distance, with equal units and common plot limits;
- `radius_diagnostics.png/.pdf`: the complete base-disk sweep;
- `sampling_comparison.png/.pdf`: the complete sweep under both measures;
- `metrics.csv`: all five methods/candidates, both sampling measures and ten radii;
- `optimizer_runs.csv`: all 120 starts, convergence status and final gradients;
- `embeddings.npz`: input samples and classical-3D/raw-stress-2D coordinates;
- `manifest.json`: timestamp, dependency versions, seed, checks and SHA-256 hashes.

The analytic gradient was checked by independent finite differences (L2 error
`1.82e-8`). All 120 optimization runs reported successful termination. The
largest absolute gradient component among selected fits was `4.66e-9`
(rounded upward). Reconstruction and stress-dominance checks run inside the
experiment; `make verify` also checks artifact hashes and citation coverage.
All three PNGs were visually inspected after rendering.

The accompanying `references.bib` and `citation_verification.html` record
claim-level source verification. The paraboloid bounds, covariance calculation,
and limiting argument are derived above; they are not attributed to those
references. This experiment does not modify the R package, its documentation,
or its tests.
