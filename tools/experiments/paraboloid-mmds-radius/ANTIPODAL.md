# Antipodal boundary points: geodesic and MDS asymptotics

Let `p_+(r)=(r,0,r²)` and `p_-(r)=(-r,0,r²)`. Their ambient distance is
`2r`. Their smooth geodesic distance satisfies `g_r ~ pi r`. The distance
assigned by MDS depends on the objective, dimension, and surrounding sample.
All MDS distances below use the original input scale, without an additional
normalization of the fitted coordinates.

The main conclusions are:

- Smooth geodesic target: `g_r = pi r - pi³/(96r) + O(r^-3)`.
- Global raw-stress MDS in at least two dimensions, for a fixed finite generic
  interior sample plus this pair: `d_MDS = pi r + O(1)`.
- Classical MDS in three dimensions, for the previous 240-point base-disk
  sample plus this pair: `d_MDS/r -> 2.92190551`.
- Classical MDS in the rotationally symmetric uniform-base-disk population
  limit: the limiting-operator calculation gives `d_MDS/r ≈ 2.86038105`.
  Uniform surface-area sampling gives approximately `2.84825260` instead.

The population constants are numerical eigenoperator evaluations with
quadrature refinement, not exact closed-form constants. The raw-stress result
has a fixed-sample assumption; it must not be interchanged with a population
or growing-sample limit. Adding an entire boundary ring creates additional
equal-radius constraints and is outside that raw-stress theorem.

## Smooth geodesic distance

Use the metric and exact primitives documented and source-verified in
`GEODESIC.md`. A boundary semicircle supplies the upper bound `pi r`.
For a shortest path with minimum radius `m`, radial travel contributes at
least `2(r²-m²)` and angular travel contributes at least `pi m`. Since its
total length is at most `pi r`, we have `m/r -> 1`, and therefore
`pi m <= g_r <= pi r` proves `g_r/r -> pi`.

For the more precise expansion, the path is symmetric around its turning
radius `c`. Put `t=sqrt(r²-c²)`. Its half-angle is `H(c,t)=pi/2` and its
length is `2L(c,t)`, using the primitives in `geodesic.py`. Expanding with
`c=sqrt(r²-t²)` gives

$$H=2t+\frac{-t^3/3+t/4}{r^2}+O(r^{-4}),$$
$$2L=4rt+\frac{-4t^3/3+t/2}{r}+O(r^{-3}).$$

Consequently, with `t0=pi/4`,
`t=t0+(t0³/6-t0/8)/r²+O(r^-4)`, and substitution gives

$$g_r=\pi r-\frac{\pi^3}{96r}+O(r^{-3}),\qquad
c=r-\frac{\pi^2}{32r}+O(r^{-3}).$$

Thus the path dips inward by only order `1/r`. A path through the vertex
has length of order `r²` and is much longer.

## Raw-stress MDS for a fixed generic finite sample

Take a fixed finite set of nonzero unit-disk interior points with distinct
squared radii `q_i`. Adjoin the two boundary points, whose squared radii
are both one. These conditions hold almost surely for a finite uniform
interior sample. Let the target dimension be at least two, use every pair
with equal weight, and consider global minimizers of fixed-target raw stress.

For endpoints with unequal squared radii, their geodesic distance is
`r²|q_i-q_j|+O(1)`. To see the bounded remainder, sort their physical radii
as `a<b`, let `h=log(b/a)>0`, and use the surface path whose angle changes
as `theta(rho)=theta(a)+alpha log(rho/a)/h`. Its length is at most

$$b^2-a^2+\frac h4+\frac{\alpha^2}{4h},$$

by `sqrt(4rho²+B) <= 2rho+B/(4rho)`. Here `h` and `alpha` do not depend
on `r`, so the remainder is bounded for each pair. The constant can be
large when two radial levels are close; it is not uniform in sample size.

Construct a feasible embedding with all interior points at `(0,r²q_i)`
and the boundary points at `(+g_r/2,r²)` and `(-g_r/2,r²)`. This fits the
boundary-pair target exactly. Every unequal-radius pair still has embedding
distance `r²|q_i-q_j|+O(1)`, because the transverse coordinates are `O(r)`.
Its raw stress is thus `O(1)` for fixed sample size. A global minimizer has
stress no larger, and in particular

$$|d_{\rm MDS}(p_+,p_-)-g_r|\le\sqrt{S_r^{\min}}=O(1).$$

It follows that `d_MDS/r -> pi`. This proves the coefficient but does not
assert exact preservation of the pair at finite radius or global convergence
of the numerical optimizer used in the earlier experiment.

## Classical MDS: finite-sample limiting matrix

Classical MDS is the algorithm used by the repository's `metric.mds()`.
For `q_i,q_j>0`, let `alpha_ij` be the smaller angular separation and define
the logarithmic mean

$$\ell(q,p)=\frac{q-p}{\log q-\log p},\qquad \ell(q,q)=q.$$

Expansion of the smooth geodesic primitives gives

$$g_{ij}^{\,2}=r^4(q_i-q_j)^2+r^2C_{ij}+O(1),$$
$$C_{ij}=\tfrac14(q_i-q_j)(\log q_i-\log q_j)
+\ell(q_i,q_j)\alpha_{ij}^{\,2}.$$

At equal radii the formula follows from the boundary-angle version of the
expansion above: `g²=r²q alpha²+O(1)`. For unequal radii, the first-order
distance correction is `|log(q_i/q_j)|/8 + alpha²/(2|log(q_i/q_j)|)`.

Let `q_c` be centered `q`, and let `H` be the orthogonal projection off both
the constant vector and `q_c`. The leading classical-MDS axis has size
`r²` and follows `q_c`. The smaller axes, on the `r` scale, are determined
by the positive eigenvalues of

$$T=-\tfrac12HCH.$$

If `(lambda_1,v_1)` and `(lambda_2,v_2)` are the two largest positive
eigenpairs, the antipodal distance in the 3D embedding satisfies

$$\frac{d_{\rm classical}(p_+,p_-)}r\longrightarrow
\left[\sum_{a=1}^2\lambda_a
\{v_a(p_+)-v_a(p_-)\}^2\right]^{1/2}.$$

The order-`r²` axis contributes zero to this limit because the two points
have the same `q`. Its higher-order contribution to their coordinate
difference is only `O(1)`. The retained transverse eigenvalues in the present
sample have the required separation from the discarded spectrum.

Using seed `20260905` and the same 240 interior points as the preceding
experiment, then adjoining the boundary pair, gives coefficient `2.92190551`.
Direct finite-radius MDS agrees:

| Radius | Geodesic distance / r | Classical 3D distance / r |
|---:|---:|---:|
| 4 | 3.121449 | 2.907121 |
| 16 | 3.140331 | 2.920982 |
| 64 | 3.141514 | 2.921848 |
| 256 | 3.141588 | 2.921902 |
| 1024 | 3.141592 | 2.921905 |

## Rotationally symmetric population calculation

In the uniform-base-disk population, `q` is uniform on `[0,1]` and the angle
is uniform. Fourier expansion of the squared angular distance gives

$$\alpha^2=\pi^2/3+4\sum_{m\ge1}(-1)^m\cos(m\alpha)/m^2.$$

The two leading transverse modes are the angular sine/cosine modes with
`m=1`. Their radial component solves the integral eigenproblem

$$\lambda f(q)=\int_0^1\ell(q,p)f(p)\,dp,\qquad
\int_0^1f(p)^2\,dp=1.$$

For this mode, the embedding radius is `r sqrt(2lambda) f(q)`. Hence the
antipodal coefficient at the boundary is `2 sqrt(2lambda) |f(1)|`.
Gauss–Legendre quadrature with 64, 128, 256 and 512 nodes gives convergent
values, with the last two differing by less than `6e-9`. At 512 nodes,
the coefficient is `2.86038105`.

Mode selection was also checked numerically. At 256 nodes the leading
`m=1` eigenvalue is about 0.511361; the largest remaining radial, even
angular, and higher odd angular eigenvalues are approximately 0.0000194,
0.002891, and 0.056818. Thus no competing mode is close to displacing the
retained pair. These are numerical operator checks, not rigorous quadrature
error bounds.

For uniform surface-area sampling on expanding disks, the limiting density
of `q` is `(3/2)sqrt(q)`. Replacing `dp` by this weighted measure in the
eigenproblem and normalization gives coefficient `2.84825260`.

## Reproduction

Run `make antipodal`, then `make verify`, in this directory. The output bundle
is `output/paraboloid-mmds-antipodal/` relative to the repository root. It
contains the finite-radius table, quadrature refinement table, limiting
matrices, sample coordinates in polar form, and a timestamped hash manifest.
The derivations above are new calculations using the metric and geodesic
primitives already documented in `GEODESIC.md`; no external source is being
credited with the MDS constants or the raw-stress theorem.
