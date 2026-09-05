# Smooth-geodesic MDS of an expanding saddle

The 3D saddle embeddings do **not** approach a plane or a line in this experiment.
They approach configurations with four narrow arms spanning three dimensions.
This holds for both classical MDS and minimization of raw distance stress, and
for both uniform base-disk and uniform surface-area sampling. The intrinsic
large-radius limit is a four-arm tree. A rank argument below proves the
nonplanarity of the classical MDS limit; the raw-stress result is numerical.

## Experiment and results

The surface is `(x,y,x²−y²)`, with `x²+y² ≤ r²`. The experiment uses 240 points,
seed 20260905, and radii 0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64. Radial quantiles
and angles are reused across radii, and match the preceding paraboloid experiment.
Both sampling measures are tested: uniform in the base disk and uniform in
surface area. Since the area density is `sqrt(1+4ρ²)`, the latter uses the exact
radial inverse CDF, with independent uniform angles.

The inputs are numerically computed **smooth-surface geodesic distances**.
There is no nearest-neighbor graph or mesh. Both 2D and 3D embeddings are fitted.
The default run evaluates 573,600 unordered point pairs across its 20 cases.

For centered embedding singular values `s₁ ≥ s₂ ≥ s₃`, a line has `s₂/s₁ → 0`;
a nondegenerate plane has `s₃/s₂ → 0`. The following are the 3D raw-stress results
for uniform base-disk sampling. Relative distance RMSE is
`sqrt(sum((d−Δ)²)/sum(Δ²))`, summed over unordered pairs.

| Radius | Relative distance RMSE | s₂/s₁ | s₃/s₂ |
|---:|---:|---:|---:|
| 0.25 | 0.003398 | 0.8777 | 0.1390 |
| 1 | 0.025017 | 0.8834 | 0.4810 |
| 4 | 0.051282 | 0.8519 | 0.9309 |
| 16 | 0.063947 | 0.9005 | 0.9116 |
| 64 | 0.066702 | 0.9109 | 0.8997 |

At `r=64`, surface-area sampling gives ratios 0.9220 and 0.9177, with relative
distance RMSE 0.061128. Classical 3D MDS with base-disk sampling gives ratios
0.9175 and 0.9103. The corresponding paraboloid raw-stress embedding has
`s₂/s₁ = 0.03743`: it is approaching a line on this normalized scale.

Generated figures and numeric outputs are in `output/saddle-mmds-geodesic/`:

- `saddle_3d_snapshots.{png,pdf}`: original, classical 3D, and raw-stress 3D
  configurations, with equal spatial units and common limits after division by
  RMS input distance. Colors denote the four arms defined below.
- `saddle_diagnostics.{png,pdf}` and `saddle_sampling.{png,pdf}`: stress and shape
  ratios across radius and sampling measures.
- `saddle_2d_snapshots.{png,pdf}`: what happens when a planar output is imposed.
- `saddle_shape.{png,pdf}`: general quadratic fits after rigid alignment to the
  original sample. The coefficient norm is `sqrt(a²+c²+b²/2)` for
  `z=ax²+bxy+cy²+dx+ey+f`. This dimensional quantity can shrink simply because
  all coordinates grow; it is not a planarity criterion. The fit is descriptive,
  not a claim that the MDS output remains exactly a quadratic surface.
- `saddle_vs_paraboloid.{png,pdf}`: comparison with the preceding experiment,
  generated when its metrics file is available.
- `metrics.csv`, `optimizer_runs.csv`, `distance_checks.csv`,
  `geodesic_validation.csv`, `embeddings.npz`, and `manifest.json`: numeric results,
  validation, all coordinates and target matrices, versions, and checksums.

## Objectives and interpretation

Raw metric stress is

\[
 S_r(X)=\sum_{i<j}(\|X_i-X_j\|-\Delta_{ij}(r))^2.
\]

This fixes the disparities and uses unit pair weights. The stress formulation
and the reason to try multiple starts are discussed in [mair2022]. Each fit here
uses six starts with L-BFGS-B and an analytic gradient. All 240 optimizer runs
reported success; the maximum selected gradient component was `4.95e−9` or less
in the normalized objective. Multiple starts do not certify a global minimum.

Classical MDS instead truncates the positive eigenpart of
`B = −½ J Δ² J`, where the square is entrywise and `J` centers observations.
The local `grip::metric.mds()` wrapper calls classical `cmdscale`, so both
interpretations of “metric MDS” are included. MDS and diagnostic functions are
shared with `../paraboloid-mmds-radius/experiment.py`.

The same absolute error contributes the same raw stress at any target distance.
The same relative error contributes proportionally to the squared target
distance. Large distances therefore matter strongly when their relative errors
are comparable, but the stress formula does not independently impose flattening.

## Curvature does not settle the question

For the two graphs, the graph-surface curvature formula gives

\[
 K_{\rm paraboloid}(\rho)=\frac{4}{(1+4\rho^2)^2},\qquad
 K_{\rm saddle}(\rho)=-\frac{4}{(1+4\rho^2)^2}.
\]

Both approach zero in magnitude. The saddle has negative curvature everywhere,
but does not approach a fixed negative curvature. Its large-distance behavior
cannot be inferred just from that sign. The limit considered here expands a
base disk and divides distances by `r²`; it is not a limit of intrinsic balls
with some other radius normalization.

## Why the computed geodesic is the shortest path

The whole saddle is complete (its metric dominates the Euclidean base metric),
simply connected, and negatively curved. Hadamard–Cartan gives a unique
minimizing geodesic between every pair [gorodski2022, Theorem 6.5.2, Corollary
6.5.3 and the paragraph preceding Lemma 6.5.4].

Writing `v=ẋ`, `w=ẏ`, and `ρ²=x²+y²`, its affine geodesic equations are

\[
 \ddot x=-\frac{4x(v^2-w^2)}{1+4\rho^2},\qquad
 \ddot y=\frac{4y(v^2-w^2)}{1+4\rho^2}.
\]

Along a nonconstant geodesic,

\[
 (\rho^2)''=2(v^2+w^2)
 -\frac{8(x^2-y^2)(v^2-w^2)}{1+4\rho^2}
 \geq\frac{2(v^2+w^2)}{1+4\rho^2}>0.
\]

Consequently the path stays inside the largest endpoint base radius. The full
surface geodesic is also the shortest path on the disk-restricted surface.

`smooth_geodesic.py` solves these equations in unit base coordinates, using
adaptive Dormand–Prince integration, analytic initial-velocity sensitivities,
damped Newton shooting, and continuation across radii. It obtains length from
the conserved metric speed. Failed solves are retried with finer continuation;
unresolved pairs stop the run.

Validation uses independent SciPy DOP853 integration and `solve_bvp` collocation
for 60 pairs at five radii. Maximum relative length disagreement was `7.87e−9`;
maximum independent endpoint discrepancy divided by disk radius was `8.11e−7`.
Exact straight rulings and axis-meridian integrals also pass. Every full matrix
passes ambient chord lower bounds, lifted-straight-path upper bounds, the tree
bounds below, and all triangle inequalities, allowing `1e−7` times RMS distance
for rounding. These are numerical accuracy checks, not interval certificates.

## The limiting metric is a four-arm tree

Fix `uᵢ=(aᵢ,bᵢ)` in the unit disk and let

\[
 P_i(r)=(ra_i,rb_i,r^2h_i),\qquad h_i=a_i^2-b_i^2,\qquad t_i=|h_i|.
\]

Assign each point to one of four arms: positive height with positive or negative
`aᵢ`, and negative height with positive or negative `bᵢ`. All zero-height points
map to the common root. Define

\[
 T_{ij}=\begin{cases}
 |t_i-t_j|,&\text{same arm},\\
 t_i+t_j,&\text{different arms}.
 \end{cases}
\]

For the smooth geodesic distance, the following bounds hold uniformly:

\[
 \boxed{r^2T_{ij}\leq\Delta_{ij}(r)\leq r^2T_{ij}+4r.}
\]

For the lower bound, path length dominates total vertical variation. Different
signs require crossing height zero. To connect opposite positive-height arms,
`x` must cross zero, where height is nonpositive; opposite negative-height arms
similarly require crossing a nonnegative height. This costs at least the sum
of the two absolute heights. The same-arm bound is the height difference.

For the upper bound, move each endpoint along its constant-height hyperbola to
the corresponding positive or negative coordinate semiaxis. Each connector
has length at most `r`: on a positive-height hyperbola `|dx/dy|≤1` and
`|y|≤r/√2`, and the negative-height case is symmetric. Connect the semiaxis
points monotonically along one axis when they share an arm, or through the
origin otherwise. The inequality `sqrt(1+4x²)≤1+2|x|` bounds that part by the
required vertical variation plus at most `2r`. All these paths stay in the disk.
The zero-height case uses a straight ruling to the origin.

Thus `Δ(r)/r² → T` uniformly. At `r=64`, replacing the entire smooth distance
matrix by `r²T` has relative distance RMSE 0.001588 for base-disk sampling and
0.001469 for surface-area sampling. For surface-area sampling the coupled unit
coordinates also converge (their radial quantiles tend to `U^(1/3)`), so the
same argument applies with their limiting tree coordinates.

The geometric distinction from the paraboloid is the disconnected height
levels. A paraboloid height level is one circle, whose distances are only `O(r)`.
A positive saddle height has two components separated by distances of order
`r²`. They survive as separate branches in the normalized limit.

For example, the saddle boundary antipodes `(±r,0,r²)` have exact geodesic distance

\[
 2\int_0^r\sqrt{1+4x^2}\,dx
 =r\sqrt{1+4r^2}+\tfrac12\operatorname{asinh}(2r)
 \sim2r^2.
\]

Reflection symmetry and uniqueness force their geodesic onto `y=0`. In contrast,
the paraboloid boundary-antipodal distance is asymptotic to `πr`, as shown in
`../paraboloid-mmds-radius/ANTIPODAL.md`. Saddle antipodes on a zero-height ruling
have distance `2r`, so “antipodal” alone does not specify one saddle asymptotic.
These statements concern input geodesics, not individual fitted MDS distances.

## A nonplanarity theorem for classical MDS

Let `F` be the `n × 4` matrix with row `tᵢ e_arm(i)`. Direct expansion gives

\[
 B_T=-\tfrac12JT^{\circ2}J
     =JF(2I_4-\mathbf1_4\mathbf1_4^T)F^TJ.
\]

The middle matrix has three positive eigenvalues (all 2) and one negative
eigenvalue (−2). If `JF` has rank four, Sylvester's law of inertia shows that
`B_T` has exactly three positive and one negative eigenvalue. The rank condition
holds for generic samples covering all four arms with enough distinct points;
it is checked for both limiting samples used here. Missing arms or special
degenerate samples require separate analysis.

Since `B(r)/r⁴ → B_T`, all three positive eigenvalues used by classical 3D MDS
remain positive in the limit. Its coordinates divided by `r²`, up to rigid
motion, therefore retain three nonzero singular values. **They cannot collapse
to a plane or a line.** For this base-disk sample the limiting classical ratios
are `s₂/s₁ = 0.91839` and `s₃/s₂ = 0.90934`.

## What follows for raw stress

Writing `X=r²Y` yields the exact identity

\[
 r^{-4}S_r(r^2Y)
 =\sum_{i<j}(\|Y_i-Y_j\|-\Delta_{ij}(r)/r^2)^2
 \longrightarrow\sum_{i<j}(\|Y_i-Y_j\|-T_{ij})^2.
\]

For fixed finite samples, this convergence is uniform on bounded configurations.
Centered global minimizers are bounded: comparison with the zero configuration
bounds stress, hence all pair distances and all centered coordinates. Thus any
accumulation point of scaled global minimizers minimizes raw stress for the
tree. This establishes the limiting optimization problem, but does not by
itself prove that every such minimizer spans three dimensions. The classical
inertia proof is not a raw-stress proof. The six-start computations provide
evidence of nonplanarity for raw stress in these samples, not a universal global
optimization theorem. The tree is not a Euclidean distance matrix under the
rank condition above, so its exact stress cannot be zero in any Euclidean
dimension.

## Reproduction

From the repository root:

```sh
make -C tools/experiments/paraboloid-mmds-radius geodesic  # optional comparator
make -C tools/experiments/saddle-mmds-radius run
make -C tools/experiments/saddle-mmds-radius verify
```

Dependencies used: Python 3.9.6, NumPy 1.24.4, SciPy 1.9.1, Numba 0.57.1,
Matplotlib 3.7.5. The manifest records actual versions. Repeated runs reuse
geodesic caches only when solver checksum and point coordinates match. Generated
outputs are intentionally ignored by Git and are recreated by `make run`.
The optional comparison reads the earlier experiment's metrics and records its
checksum; the matched default sample size and seed should be retained for a
paired comparison. `experiment.py --help` exposes sample-size, seed, and output
directory overrides. The verifier audits the standard 240-point reproduction.

## References

- [mair2022] Mair, Groenen, and de Leeuw (2022), *More on Multidimensional Scaling
  and Unfolding in R: smacof Version 2*. [Article and source](https://doi.org/10.18637/jss.v102.i10).
- [gorodski2022] Claudio Gorodski (2022), *An introduction to Riemannian geometry*,
  preliminary version 3, June 27, pp. 124–125.
  [Author's full notes](https://www.ime.usp.br/~gorodski/teaching/mat5771-2024/master02-21-2024.pdf).

Metadata and claim-level evidence are in `references.bib` and
`citation_verification.html`. Curvature calculations, disk convexity, tree
bounds, the classical MDS inertia calculation, and the raw-stress limit above
are derivations supplied with this experiment, not claims attributed to these
sources.
