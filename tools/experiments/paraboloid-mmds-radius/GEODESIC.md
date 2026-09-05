# MDS from smooth-paraboloid geodesic distances

For the boundary-antipodal pair, see [ANTIPODAL.md](ANTIPODAL.md). It resolves
the order-`r` separation that disappears in the order-`r²` line limit below.

Using geodesic distances changes the finite-radius result: 3D MDS does not
recover the original paraboloid exactly. In this experiment, both classical
and raw-stress MDS produce approximately paraboloidal shapes that are shallower
than the original. However, increasing the radius does not progressively
flatten them into a two-dimensional plane. The mathematically established
limit, after isotropic scaling by `r²`, is one-dimensional.

These conclusions concern `z=x²+y²` with coefficient fixed at one, and
shortest distances along the smooth surface over the closed base disk.
No nearest-neighbor graph, triangulation, or mesh shortest path is used.
The earlier ambient-distance experiment is retained in `README.md` as a
comparison, not as the answer for these dissimilarities.

## Direct smooth-geodesic distances

In polar coordinates `(rho,theta)`, the induced metric is

$$ds^2=(1+4\rho^2)d\rho^2+\rho^2d\theta^2.$$

Clairaut's first integral gives `c=rho² dtheta/ds` along a unit-speed
geodesic. The turning radius is `rho=c`. The geodesic equation and an
elementary angular primitive for this paraboloid are given by
[Wheeler (2016), pp. 5–6, Eqs. (9.3)–(9.4)](https://www.reed.edu/physics/faculty/wheeler/documents/Miscellaneous%20Math/Differential%20Geometry/DINI%20SURFACES/Surface%20of%20Revolution%20Geodesics/Geodesic%20Text.pdf).
The same source gives Gaussian curvature `4/(1+4rho²)²` in Eq. (9.2).
<!-- cite:wheeler2016 -->

For this implementation, put `t=sqrt(rho²-c²)`, `A=1+4c²`, and
`v=asinh(2t/sqrt(A))`. The angle and length measured from the turning point
are the following primitives:

$$H(c,t)=2cv+\operatorname{atan2}\!\left(t,c\sqrt{A+4t^2}\right),$$
$$L(c,t)=\tfrac12t\sqrt{A+4t^2}+\tfrac14Av.$$

For example, differentiating the length primitive gives
`dL/dt=sqrt(1+4c²+4t²)`. The angular differential is
`dtheta/drho=c sqrt(1+4rho²)/(rho sqrt(rho²-c²))`. Its integral yields `H`.

Sort the endpoint radii as `a <= b` and take their smaller angular separation
`alpha` in `[0,pi]`. Let `alpha0=H(a,sqrt(b²-a²))`. If `alpha <= alpha0`,
solve the monotone-radial branch

$$H(c,\sqrt{b^2-c^2})-H(c,\sqrt{a^2-c^2})=\alpha.$$

Otherwise solve the turning branch, replacing the difference by a sum.
The geodesic length is the corresponding difference or sum of `L`. The
implementation uses 62 vectorized bisection iterations. It parametrizes
`c=a sin(phi)` and `t_a=a cos(phi)` to avoid cancellation near a turning
point. Pole and same-meridian distances are handled explicitly.

This selects shortest paths, rather than arbitrary winding geodesics. A
path with total angular variation larger than the smaller endpoint separation
can be shortened by compressing its angular coordinate while keeping its
radial profile. In radial arclength coordinates, a geodesic has nonnegative
radial second derivative, so it has no radial maximum and at most one radial
minimum. It stays inside the larger endpoint radius; restricting the surface
to the base disk therefore causes no extra boundary detour.

The angular difference on the monotone branch increases with `c`. On the
turning branch the angular sum is strictly concave in `c`, starts at `pi`
when `c=0`, and ends at `alpha0` when `c=a`. For `alpha<pi` in the turning
range there is exactly one root. Concavity follows from

$$\frac{d}{dc}H(c,\sqrt{\rho^2-c^2})
=2\operatorname{asinh}\!\frac{2\sqrt{\rho^2-c^2}}{\sqrt{1+4c^2}}
-\frac{\sqrt{1+4\rho^2}}{\sqrt{\rho^2-c^2}},$$

whose derivative is negative for `0<c<rho`. At `alpha=pi`, if a positive
root exists it is shorter than the through-apex meridian; otherwise the
meridian is selected. This also follows by integrating `dL/dc=c dH/dc`
along the turning family. These branch-selection arguments and the length
primitive are derived here rather than attributed to Wheeler's discussion
of the full, potentially winding geodesics.

## Numerical experiment

There are 240 points, seed `20260905`, and ten radii:
`0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64`. Both uniform base-disk sampling
and uniform surface-area sampling are run. Their exact sampling formulas
are documented in `README.md`; quantiles and angles are reused across radii.

For each case, classical MDS is computed in 2D and 3D, and raw stress is
optimized in both dimensions. Each stress fit has six starts. The 2D starts
are classical MDS, the three input coordinate-plane projections, and two
random configurations. The 3D starts are classical MDS, the original surface,
two perturbed 2D fits, and two random configurations. Perturbing the third
coordinate avoids trapping a 3D optimizer in a planar stationary subspace.
L-BFGS-B, analytic gradients, common RMS distance normalization and stopping
tolerances match the Euclidean comparison. These are local numerical fits,
not certified global optima.

Relative distance RMSE means
`sqrt(sum((d-g)²)/sum(g²))`, using input geodesics `g` in the denominator.
Let `s1 >= s2 >= s3` denote the singular values of centered 3D coordinates.
The ratio `s2/s1` measures elongation toward a line; `s3/s2` distinguishes
the two smaller directions. A nondegenerate plane has vanishing `s3/s2`
while retaining a nonvanishing `s2/s1`. A line has vanishing `s2/s1`.

The base-disk results for 3D raw-stress MDS are:

| Radius | Relative distance RMSE | s2/s1 | s3/s2 |
|---:|---:|---:|---:|
| 0.25 | 0.005394 | 0.873502 | 0.071654 |
| 1 | 0.058801 | 0.849537 | 0.407892 |
| 2 | 0.073424 | 0.843417 | 0.948380 |
| 4 | 0.051547 | 0.600171 | 0.840463 |
| 16 | 0.010330 | 0.149400 | 0.843284 |
| 64 | 0.001378 | 0.037432 | 0.846729 |

Small disks give nearly planar fits. Intermediate disks give three substantial
directions. Large disks give increasingly elongated fits while the two smaller
directions remain comparable to each other. At `r=64`, the first principal
component accounts for 99.760% of the variance. Surface-area sampling has the
same qualitative pattern: at `r=64`, raw-stress 3D relative RMSE is 0.001931,
`s2/s1=0.045085`, and `s3/s2=0.857030`.

To examine “flatter paraboloids” directly, rigidly align each 3D fit to its
input coordinates without scaling, then regress

$$z=a(x^2+y^2)+b x+c y+d.$$

The linear terms allow a shifted vertex. This is a shape diagnostic in the
original units; it does not assume the MDS points form an exact paraboloid
or estimate an intrinsic curvature of the point cloud. The original has
`a=1`. For base-disk sampling:

| Radius | Classical-MDS coefficient a | Raw-stress-MDS coefficient a |
|---:|---:|---:|
| 0.25 | 0.5543 | 0.4248 |
| 1 | 0.4788 | 0.4849 |
| 4 | 0.5037 | 0.5793 |
| 16 | 0.5045 | 0.5780 |
| 64 | 0.5045 | 0.5611 |

The fits are shallower than the original, but these coefficients do not
decrease toward zero over the tested radii. At `r=64`, the regression R² is
0.9883 for classical MDS and 0.9657 for raw-stress MDS. The coefficient's
infinite-radius limit is not established by this experiment or by the
leading-order distance limit below. `geodesic_shape.png` shows the complete
coefficient and R² sweeps under both sampling measures.

The fact that Gaussian curvature of the original surface tends to zero far
from the vertex is a local statement about distant patches. It does not imply
that enlarging the whole disk yields a nondegenerate planar MDS limit.

## What can be proved from the stress loss?

For fixed unit-disk points `u_i`, put `q_i=||u_i||²` and
`P_i(r)=(r u_i,r²q_i)`. Every surface path has length at least the absolute
height difference `r²|q_i-q_j|`. For an upper bound, move radially to the
smaller endpoint radius, then along its shorter circular arc. Since
`sqrt(1+4rho²) <= 1+2rho`, this route gives

$$r^2|q_i-q_j|\le g_{ij}^{(r)}
\le r^2|q_i-q_j|+(1+\pi)r.$$

Consequently `g_ij/r² -> |q_i-q_j|`, uniformly over pairs. The latter is
the distance matrix of points on a line. Equal-radius points lose their
angular separation on this scale.

For raw stress `S_r(Y)=sum_{i<j}(||Y_i-Y_j||-g_ij)²`, use `V=Y/r²`:

$$r^{-4}S_r(Y)=\sum_{i<j}(\|V_i-V_j\|-g_{ij}^{(r)}/r^2)^2.$$

In any fixed target dimension at least one, the line candidate
`V_i=(q_i,0,...)` has normalized stress at most
`binom(n,2)(1+pi)²/r²`, tending to zero. Every global minimizer does at
least as well. For fixed finite `n`, each normalized pair distance therefore
converges to `|q_i-q_j|`. Double-centering their squared distances gives the
limiting centered Gram matrix `q_c q_c'`. Provided at least two squared radii
differ, this has rank one. Approximate minimizers whose excess normalized
stress tends to zero have the same conclusion.

Classical MDS has the same leading limit because its rescaled input Gram
matrix converges to this rank-one matrix. This proves a normalized line
limit for both objectives; it does not prove convergence to a flat disk,
monotonic behavior at finite radii, or a limit for the fitted coefficient `a`.
Unscaled configurations grow without bound.

## Verification and reproduction

From this directory run `make geodesic`, then `make verify`. The new output
bundle is `output/paraboloid-mmds-geodesic/` relative to the repository root,
separate from the Euclidean comparison. Dependencies are Python >=3.9,
NumPy, SciPy and Matplotlib; actual versions and timestamps are in its manifest.

The distance solver was independently checked on 50 pairs across five radii
against adaptive quadrature of the original angle and length integrals, and
against the Cartesian geodesic ODE
`x''=-4x ||x'||²/(1+4||x||²)`. Maximum errors were `3.05e-15` for relative
length, `1.64e-14` radians for angle, and `3.98e-11` for endpoint error divided
by disk radius (all rounded upward). Pole, duplicate, same-meridian and
opposite-meridian cases are also checked. All 20 full distance matrices
satisfy chord lower bounds, radial-plus-arc upper bounds, and every triangle
inequality within stated floating-point tolerances.

All 240 optimization runs reported successful termination. The largest
absolute gradient component among the 40 selected fits was `6.58e-9`
(rounded upward). Selected 3D fits were checked against the 2D and classical
3D raw-stress values. The five PNG figures were inspected after rendering.
This is one coupled random sample per sampling measure, not a Monte Carlo
study with confidence intervals.

The bundle includes 120 method/case rows in `metrics.csv`, all starts in
`optimizer_runs.csv`, `geodesic_validation.csv`, `distance_bounds.csv`, the
full geodesic matrices and coordinates in `embeddings.npz`, PNG/PDF figures,
and a timestamped SHA-256 manifest. Source files and citation-verification
evidence are versioned; generated outputs are ignored and reproducible.
