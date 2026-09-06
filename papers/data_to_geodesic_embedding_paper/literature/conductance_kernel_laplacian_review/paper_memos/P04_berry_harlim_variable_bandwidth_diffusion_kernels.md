# Paper Review: Variable Bandwidth Diffusion Kernels

Paper: Berry, Tyrus; Harlim, John. "Variable Bandwidth Diffusion Kernels." Applied and Computational Harmonic Analysis 40(1):68-96, 2016. DOI: 10.1016/j.acha.2015.01.001. arXiv:1406.5064.
Reviewer: Reviewer-P04
Auditor: Auditor-P04
Status: revised after audit
Date: 2026-05-15
Source manifest ID: P04
Canonical reading copy: `/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/literature/conductance_kernel_laplacian_review/sources/pdf/P04_berry_harlim_variable_bandwidth_diffusion_kernels.pdf`

## Whole-Paper Review

### Reader Background Needed

A mathematically capable reader should be comfortable with:

- Graph Laplacians and weighted adjacency kernels: how a kernel matrix becomes a Markov matrix or graph Laplacian.
- Diffusion maps and the Coifman-Lafon `alpha` normalization: `alpha=1` removes sampling-density bias; other `alpha` values retain density drift.
- Kernel density estimation: fixed bandwidth versus local bandwidth, and why sparse tails are hard.
- Riemannian manifolds embedded in Euclidean space: gradients, the Laplace-Beltrami operator, and volume measure.
- Backward Kolmogorov generators: operators of the form `Delta f + b dot grad f`.
- Monte Carlo approximation of integral operators from random point clouds.
- Spectral approximation: eigenvectors of a discrete graph operator as samples of eigenfunctions of a continuum operator.

### What A Non-Expert Should Understand Before Reading This Paper

The paper asks what happens when the Gaussian-like affinity between two data points uses a local scale instead of one global bandwidth. In a fixed-bandwidth graph, two points at a given Euclidean distance receive the same edge weight no matter whether they lie in a dense region or a sparse tail. In a variable-bandwidth graph, the denominator of the kernel includes a local scale `rho(x)rho(y)`, so the same distance can be treated as "near" in sparse regions and "far" in dense regions.

The main message is subtle: changing bandwidth is not merely a numerical stabilization trick. It changes the continuum differential operator unless the normalization parameters are chosen carefully. For `rho=q^beta`, where `q` is the sampling density, the limit is `Delta f + c1 grad f dot grad q/q`; the coefficient `c1` depends on both the diffusion-map normalization `alpha` and the bandwidth exponent `beta`. The finite-sample error also depends on `q` and `rho`, and fixed-bandwidth kernels can have unbounded pointwise errors where `q` approaches zero. This is the paper's strongest connection to local-scale conductance: adaptive bandwidth can preserve effective connectivity in sparse regions while controlling the operator that the graph approximates.

### Problem And Context

The paper generalizes the asymptotic and finite-sample theory for diffusion maps and Laplacian eigenmaps from fixed-bandwidth kernels to self-tuned or variable-bandwidth kernels. The motivation is stated in the Abstract and Section 1 (PDF pp. 1-2): practical kernel methods often use self-tuned kernels, but much of the theory available at the time covered only fixed bandwidths. The authors cite Zelnik-Manor and Perona's self-tuning spectral clustering as a practical example, and Coifman-Lafon diffusion maps plus Singer's finite-sample analysis as the theoretical background.

The key gap is that previous general convergence work by Ting, Huang, and Jordan gave limiting operators for broad graph Laplacian sequences but not rates or discrete error bounds (Section 1, PDF p. 2). Berry and Harlim provide asymptotic expansions and pointwise high-probability errors for discrete variable-bandwidth operators. The most practically important case is non-compact or effectively open data support, where the sampling density is not bounded away from zero.

### Main Method

The method uses a symmetric variable-bandwidth kernel

```text
K_epsilon^S(x,y) = h( ||x-y||^2 / (epsilon rho(x) rho(y)) )
```

from equation (2) on PDF p. 1, then applies a density correction and left Markov normalization analogous to diffusion maps. The bandwidth function is often chosen as `rho=q^beta` or estimated from data using a nearest-neighbor pilot density estimate (Section 3, PDF pp. 5-6).

The theoretical development has two layers:

1. A continuum expansion of the associated integral operators, showing that the symmetric variable-bandwidth operator converges to

```text
L_{alpha,rho} f = Delta f
                 + 2(1-alpha) grad f dot grad q/q
                 + (d+2) grad f dot grad rho/rho
```

from Theorem 1 and equation (6), PDF p. 3.

2. A discrete Monte Carlo error analysis, extending Singer's fixed-bandwidth bound, showing that pointwise errors depend on powers of `q(x_i)`, `rho(x_i)`, `||grad f(x_i)||`, `N`, and `epsilon` (Theorem 1, equation (5), PDF p. 3; Appendix B, PDF pp. 22-26).

### Main Formulas And Operators

The paper's fixed-bandwidth baseline is equation (1), PDF p. 1:

```text
K_epsilon(x,y) = h( ||x-y||^2 / epsilon ).
```

The variable-bandwidth kernel is equation (2), PDF p. 1:

```text
K_epsilon^S(x,y) = h( ||x-y||^2 / (epsilon rho(x)rho(y)) ).
```

For fixed bandwidth, diffusion maps approximate the Kolmogorov operator in equation (3), PDF p. 2:

```text
L_alpha f = Delta f + (2-2alpha) grad f dot grad q/q.
```

Singer's fixed-bandwidth high-probability error, quoted as equation (4), PDF p. 3, is:

```text
L_{epsilon,alpha} f(x_i)
  = (1/epsilon) [ sum_j K_epsilon(x_i,x_j) f(x_j) / sum_j K_epsilon(x_i,x_j) - f(x_i) ]
  = L_alpha f(x_i) + O( epsilon, ||grad f(x_i)|| / sqrt(N epsilon^{1+d/2}) )
```

where the paper writes the denominator as `sqrt(N epsilon^{1/2+d/4})` after the square-root convention in the displayed bound.

Theorem 1, PDF p. 3, defines:

```text
F_i(x_j) = K_epsilon^S(x_i,x_j) f(x_j) / [q_epsilon^S(x_i)^alpha q_epsilon^S(x_j)^alpha]
G_i(x_j) = K_epsilon^S(x_i,x_j)        / [q_epsilon^S(x_i)^alpha q_epsilon^S(x_j)^alpha]
q_epsilon^S(x_i) = sum_l K_epsilon^S(x_i,x_l) / rho(x_i)^d.
```

The corresponding discrete operator in equation (5), PDF p. 3, is:

```text
L_{epsilon,alpha}^S f(x_i)
  = (1 / [epsilon m rho(x_i)^2])
    [ sum_j F_i(x_j) / sum_j G_i(x_j) - f(x_i) ].
```

With high probability,

```text
L_{epsilon,alpha}^S f(x_i)
  = L_{alpha,rho} f(x_i)
    + O( epsilon,
         q(x_i)^{1/2} rho(x_i)^{-d/2} / sqrt(N epsilon^{4+d/2}),
         ||grad f(x_i)|| q(x_i)^{-(1/2-2alpha+2dalpha)}
           rho(x_i)^{-(d/2+1)} / sqrt(N epsilon^{1+d/2}) )
```

using the paper's displayed equation (5) notation:

```text
O( epsilon,
   q(x_i)^{1/2} rho(x_i)^{-d/2} / [sqrt(N) epsilon^{2+d/4}],
   ||grad f(x_i)|| q(x_i)^{-(1/2-2alpha+2dalpha)}
     rho(x_i)^{-(d/2+1)} / [sqrt(N) epsilon^{1/2+d/4}] ).
```

The continuum operator is equation (6), PDF p. 3:

```text
L_{alpha,rho} f =
  Delta f + 2(1-alpha) grad f dot grad q/q
          + (d+2) grad f dot grad rho/rho.
```

Corollary 1, equation (7), PDF p. 3, specializes to `rho=q^beta+O(epsilon)`:

```text
L_{epsilon,alpha,beta}^S f(x_i)
  = L_{alpha,beta} f(x_i)
    + O( epsilon,
         q(x_i)^{(1-d beta)/2} / [sqrt(N) epsilon^{2+d/4}],
         ||grad f(x_i)|| q(x_i)^{-c2} / [sqrt(N) epsilon^{1/2+d/4}] ).
```

The limiting operator is equation (8), PDF p. 3:

```text
L_{alpha,beta} f = Delta f + c1 grad f dot grad q/q,
c1 = 2 - 2alpha + d beta + 2 beta,
c2 = 1/2 - 2alpha + 2dalpha + d beta/2 + beta.
```

For non-compact manifolds where `q` may approach zero, Section 2 (PDF p. 4) says the error powers should satisfy:

```text
(1-d beta)/2 > 0,
c2 < 0.
```

The paper emphasizes `beta<0`, especially `beta=-1/2`, because it enlarges bandwidth in sparse regions and shrinks it in dense regions. For gradient flow estimation with `c1=1`, `beta=-1/2` gives `alpha=-d/4`; for the Laplacian with `c1=0`, `beta=-1/2` gives `alpha=1/2-d/4` (Section 2, PDF p. 4).

Section 3, equation (9), PDF p. 5, gives the stochastic process target:

```text
dx = -c1 grad U(x) dt + sqrt(2) dW_t,
q(x) proportional to exp(-c1 U(x)).
```

Its backward generator is `L_{alpha,beta}` in equation (8).

The nearest-neighbor pilot bandwidth and density estimate in Section 3, equation (10), PDF p. 5, are:

```text
rho_0(x_i) = [ (1/(k0-1)) sum_{j=2}^{k0} ||x_i - x_{I(i,j)}||^2 ]^{1/2},
epsilon_0^{1/2} = (1/N) sum_i rho_0(x_i),
tilde rho_0 = rho_0 / epsilon_0^{1/2},
```

and

```text
q_0(x_i) =
  (2pi)^{-d/2} / [rho_0(x_i)^d N]
  sum_{l=1}^N exp( -||x_i-x_l||^2 / [2 rho_0(x_i)rho_0(x_l)] )
```

equivalently displayed with `epsilon_0` and `tilde rho_0` in equation (10), with error

```text
q_0(x_i) = q(x_i) + O( epsilon_0,
                       sqrt(q(x_i)) / [N^{1/2} epsilon_0^{d/4} tilde rho_0(x_i)^{d/2}] ).
```

The paper then uses `rho=q_0^beta=q^beta+O(epsilon_0)` and notes that balancing the errors gives `epsilon_0=O(N^{-1/(1+d/4)})` (Section 3, PDF p. 5).

The numerical implementation kernel and normalizations, Section 3, PDF pp. 5-6, are:

```text
K_epsilon^S(x_i,x_j) =
  exp( -||x_i-x_j||^2 / [4 epsilon rho(x_i)rho(x_j)] ),

q_epsilon^S(x_i) =
  sum_{j=1}^N K_epsilon^S(x_i,x_j) / rho(x_i)^d,

K_{epsilon,alpha}^S(x_i,x_j) =
  K_epsilon^S(x_i,x_j) / [q_epsilon^S(x_i)^alpha q_epsilon^S(x_j)^alpha],

q_{epsilon,alpha}^S(x_i) =
  sum_{j=1}^N K_{epsilon,alpha}^S(x_i,x_j),

hat K_{epsilon,alpha}^S(x_i,x_j) =
  K_{epsilon,alpha}^S(x_i,x_j) / q_{epsilon,alpha}^S(x_i),

L_{epsilon,alpha,beta}^S(x_i,x_j) =
  [hat K_{epsilon,alpha}^S(x_i,x_j) - delta_ij] / [epsilon rho(x_i)^2].
```

For eigensolvers, Section 3 constructs a symmetric conjugate. With `D_ii=q_{epsilon,alpha}^S(x_i)`, `P_ii=rho(x_i)`, and `K_ij=K_{epsilon,alpha}^S(x_i,x_j)`,

```text
L = P^{-2}(D^{-1}K - I)/epsilon,
S = P D^{1/2},
hat L = (1/epsilon)(S^{-1} K S^{-1} - P^{-2}),
```

with entries

```text
hat L_ij =
  (1/[epsilon rho(x_i)rho(x_j)])
  [ K_{epsilon,alpha}^S(x_i,x_j) /
    sqrt(q_{epsilon,alpha}^S(x_i)q_{epsilon,alpha}^S(x_j))
    - delta_ij ].
```

Eigenvectors are normalized so that `||phi_vec||_{R^N}=sqrt(N)`, matching the Monte Carlo `L^2(q)` norm (Section 3, PDF p. 6). Sparse implementation keeps `k` nearest neighbors and symmetrizes `K` as `(K+K^T)/2` (Section 3, PDF pp. 6-7).

Appendix A gives the continuum operator derivation. It defines the integral operator in equation (A.1), PDF p. 14:

```text
G_epsilon^S f = epsilon^{-d/2} int_M K_epsilon^S(x,y) f(y) dV(y).
```

It also defines right and left formulations:

```text
K_epsilon^R(x,y) = h( ||x-y||^2 / [epsilon rho(y)] )        (A.2, PDF p. 14)
K_epsilon^L(x,y) = h( ||x-y||^2 / [epsilon rho(x)] )        (A.3, PDF p. 15)
```

Lemma 2, equation (A.4), PDF p. 15, extends the fixed-bandwidth expansion to non-compact manifolds under integrability and bounded-density assumptions:

```text
G_epsilon(fq)(x)
  = m0 f(x)q(x)
    + epsilon m2[ omega(x)f(x)q(x) + Delta(fq)(x) ]
    + O(epsilon^2).
```

The left formulation yields `L_epsilon^L f = Delta f + O(epsilon)` in equation (A.10), PDF p. 17. The right formulation yields the drifted operator

```text
L_epsilon^R f =
  Delta f + (d+2) grad rho/rho dot grad f + O(epsilon)
```

after equation (A.11), PDF p. 18. The symmetric formulation has expansion (A.12), PDF p. 19, and normalized limit (A.13):

```text
L_epsilon^S f =
  (1/[epsilon m rho(x)^2]) [G_epsilon^S f(x)/G_epsilon^S 1(x) - f(x)]
  = Delta f + (d+2) grad rho/rho dot grad f + O(epsilon).
```

For non-uniform sampling, Appendix A.5 defines the biased operator `G_{epsilon,q}^S(f)=G_epsilon^S(fq)`, equation (A.15) for density expansion,

```text
G_{epsilon,q}^S(1)
  = m0 rho^d q [1 + epsilon m(tilde omega + L^S q) + O(epsilon^2)],
```

and the right debiasing in equation (A.16):

```text
G_{epsilon,q,alpha}^S(f) =
  G_{epsilon,q}^S( f rho^{d alpha} / G_{epsilon,q}^S(1)^alpha ).
```

The normalized non-uniform result is equation (A.20), PDF p. 22:

```text
L_{epsilon,alpha}^S f(x)
  = Delta f
    + 2(1-alpha) grad f dot grad q/q
    + (d+2) grad f dot grad rho/rho
    + O(epsilon).
```

With `rho=q^beta`, equation (A.21), PDF p. 22, gives:

```text
L_{epsilon,alpha,beta}^S f(x)
  = Delta f + c1 grad f dot grad q/q + O(epsilon).
```

Appendix B defines the discrete density normalization factor in (B.1), PDF p. 23:

```text
q_epsilon^S(x_j) = sum_l K_epsilon^S(x_j,x_l) / rho(x_j)^d.
```

The discrete ratio approximating the continuous operator is equation (B.2), PDF p. 23:

```text
L_{epsilon,alpha}^S f(x_i)
  approx (1/[epsilon m rho(x_i)^2])
         [sum_j F_i(x_j) / sum_j G_i(x_j) - f(x_i)].
```

The density renormalization error is controlled by (B.7), PDF p. 25:

```text
q(x_j)^{1/2} rho(x_j)^{-d/2} / [N^{1/2} epsilon^{2+d/4}] = O(1).
```

The statistical bias error is controlled by (B.11), PDF p. 26:

```text
||grad f(x_i)|| q(x_i)^{-(1/2-2alpha+2dalpha)}
rho(x_i)^{-(d/2+1)}
/ [sqrt(N) epsilon^{1/2+d/4}] = O(1).
```

### Figures And Experiments

Figure 1 (Section 4, PDF p. 7) compares fixed bandwidth (`alpha=1/2, beta=0`) and variable bandwidth (`alpha=-1/4, beta=-1/2`) for the fourth eigenfunction of the 1D Ornstein-Uhlenbeck generator using 2000 deterministic quantile-spaced points from a standard normal density. Fixed bandwidth is visibly sensitive to `epsilon`; variable bandwidth matches the analytic curve across a wider range.

Figure 2 (Section 4, PDF p. 7) repeats the 1D Ornstein-Uhlenbeck experiment with 20000 quantile-spaced points. The variable-bandwidth MSE has a broad low-error region, while the fixed-bandwidth MSE remains poor over the tested range; the right panel shows the best fixed-bandwidth eigenfunction still failing in the tails.

Figure 3 (Section 4, PDF p. 8) tests a proposed fixed-bandwidth workaround: remove outliers. Even with deterministic samples and removal of `sqrt(N)` low-probability points, the fixed-bandwidth approximation converges slowly; the authors estimate that matching the variable-bandwidth MSE of 0.002 from 1000 points would require about `6*10^7` fixed-bandwidth points under the observed power law.

Figure 4 (Section 4, PDF p. 9) uses 10 randomly sampled 1D Ornstein-Uhlenbeck data sets with `N=20000`. Fixed bandwidth is tuned over 65 epsilon values for each data set but still does not recover the analytic fourth eigenfunction well. Variable bandwidth is substantially better, though one random data set gives a poor estimate, consistent with the theorem's pointwise high-probability nature.

Figure 5 (Section 4, PDF p. 10) studies the 2D Ornstein-Uhlenbeck process with standard Gaussian invariant measure and analytic fourth eigenfunction `phi_4(x,y)=xy`. Variable bandwidth produces a visibly reasonable estimate; fixed bandwidth does not, despite empirical epsilon tuning.

Figure 6 (Section 5, PDF p. 11) studies the Laplacian on a unit circle sampled with density `q(theta)=(1/(4pi))(2+cos(theta))`. Variable bandwidth (`alpha=1/4, beta=-1/2`) yields lower MSE over a much wider epsilon range than fixed bandwidth (`alpha=1, beta=0`), both for deterministic inverse-CDF samples and perturbed samples. The right panel shows the automatic epsilon tuning statistic `S(epsilon) proportional to epsilon^a`, whose maximum slope is near `d/2=1/2`.

Figure 7 (Section 5, PDF p. 12) studies a non-uniformly sampled unit sphere. Variable bandwidth (`alpha=0, beta=-1/2`) is less sensitive to epsilon and works with the automatic tuning method. Fixed bandwidth (`alpha=1, beta=0`) can work only with careful tuning; the automatic epsilon is too small and disconnects sparse regions from the rest of the manifold.

Figure A.8 (Appendix A.4, PDF p. 20) numerically verifies that left-bandwidth kernels converge to `Delta f`, while right and symmetric bandwidth formulations include the drift `(d+2) grad rho/rho dot grad f`. It uses `rho(theta)=exp(cos(theta))` on a circle and a torus.

Figure A.9 (Appendix A.5, PDF p. 23) numerically verifies the `rho=q^beta` expansion for non-uniform sampling on the circle. With `beta=-1/2`, `alpha=1/4` recovers the Laplacian (`c1=0`), and `alpha=-1/4` recovers the backward Kolmogorov operator with `c1=1`.

### Theoretical Claims

Theorem 1 (Section 2, PDF p. 3) is the central claim: for a smooth function on a manifold without boundary and a density \(q\) that is bounded above and in \(L^1 \cap C^3\), the discrete variable-bandwidth operator converges pointwise with high probability to \(L_{\alpha,\rho}\) and has three error components: asymptotic bias \(O(\epsilon)\), density-renormalization sampling error, and ratio-estimator sampling error.

Corollary 1 (Section 2, PDF pp. 3-4) specializes the theory to `rho=q^beta+O(epsilon)`. This is the bridge from local bandwidth to sampling-density-adaptive operators. It identifies `c1`, the density drift in the limiting operator, and `c2`, the exponent controlling sparse-density blow-up in the finite-sample error.

Appendix A proves the continuum expansions. Its main conceptual result is that symmetric variable bandwidth acts like right bandwidth after a local coordinate scaling, which is why the term `(d+2) grad rho/rho dot grad f` appears.

Appendix B proves the finite-sample rate by adapting Singer's Chernoff-bound analysis. It shows that the density estimate error and the ratio-estimator error behave differently, and that fixed bandwidth places sampling density in the denominator of the operator-estimation error.

### Limitations And Scope

The results are pointwise for smooth functions, not a full spectral convergence theorem (Section 2, PDF p. 5; Section 6, PDF p. 13). The numerical experiments compare eigenfunctions, but the authors explicitly say extending spectral convergence to variable-bandwidth kernels is beyond the paper's scope.

The theory excludes manifolds with boundary, especially non-compact manifolds with non-compact boundaries (Section 2, PDF pp. 4-5; Section 6, PDF p. 13). Fixed-bandwidth diffusion maps on compact manifolds have implicit Neumann boundary conditions, but the authors do not extend that proof here.

The recommended `rho=q^beta` construction requires intrinsic dimension `d`, both for choosing `alpha,beta` and for the `rho^d` density normalization (Section 6, PDF p. 13). The paper suggests automatic epsilon tuning may help estimate dimension but does not establish a robust dimension-estimation method.

The kernels are isotropic scalar-bandwidth kernels. The authors mention anisotropic generalizations but leave them outside the paper (Section 6, PDF p. 13).

### Historical / Methodological Importance

This paper is a key bridge between practical self-tuned graph kernels and diffusion-map continuum theory. It explains why variable bandwidth is more than a heuristic for clustering: with the correct normalization it can approximate Laplacians or gradient-flow generators, and with the wrong normalization it approximates a different drift-diffusion operator. It also makes a strong finite-sample argument that fixed bandwidth can fail on non-compact or heavy-tailed data because sparse-density points enter the graph as high-error locations.

For a conductance-kernel review, P04 is important because it links local affinity strength, local length scale, density normalization, and limiting differential operators in one explicit formula family. It is also the conceptual parent of later kNN self-tuned graph Laplacian theory.

## Conductance / Kernel Extraction

### Conductance, Affinity, Or Kernel Formula(s)

`[explicit]` Fixed affinity, equation (1), PDF p. 1:

```text
K_epsilon(x,y) = h( ||x-y||^2 / epsilon ).
```

`[explicit]` Symmetric variable-bandwidth affinity, equation (2), PDF p. 1:

```text
K_epsilon^S(x,y) = h( ||x-y||^2 / [epsilon rho(x)rho(y)] ).
```

`[explicit]` Gaussian implementation, Section 3, PDF p. 5:

```text
K_epsilon^S(x_i,x_j) =
  exp( -||x_i-x_j||^2 / [4 epsilon rho(x_i)rho(x_j)] ).
```

`[derived]` Local conductance interpretation for nearby points: if `rho(y) approx rho(x)`, then a distance `delta=||x-y||` gets approximate affinity

```text
K(x,y) approx exp( -delta^2 / [4 epsilon rho(x)^2] ),
```

so the local `exp(-1)` length scale is about `2 sqrt(epsilon) rho(x)`. If `rho=q^{-1/2}`, sparse regions with small `q` have larger local length scale and therefore larger edge affinity at the same physical distance. This is derived from the Section 3 Gaussian formula, not stated as a conductance theorem by the authors.

`[explicit]` Right and left bandwidth kernels in Appendix A, equations (A.2)-(A.3), PDF pp. 14-15:

```text
K_epsilon^R(x,y) = h( ||x-y||^2 / [epsilon rho(y)] ),
K_epsilon^L(x,y) = h( ||x-y||^2 / [epsilon rho(x)] ).
```

### Graph, Laplacian, Or Diffusion Operator

`[explicit]` The graph/diffusion operator is the row-normalized, density-debiased kernel operator

```text
L_{epsilon,alpha,beta}^S(x_i,x_j) =
  [hat K_{epsilon,alpha}^S(x_i,x_j) - delta_ij] / [epsilon rho(x_i)^2],
```

from Section 3, PDF p. 6.

`[explicit]` Its symmetric conjugate is

```text
hat L = (1/epsilon)(S^{-1} K S^{-1} - P^{-2}),
S = P D^{1/2},
P_ii=rho(x_i),
D_ii=q_{epsilon,alpha}^S(x_i).
```

This is used for eigenvalue computation (Section 3, PDF p. 6).

`[explicit]` The continuum limit for arbitrary `rho` is `L_{alpha,rho}` in equation (6), PDF p. 3.

`[explicit]` The continuum limit for `rho=q^beta` is `L_{alpha,beta}` in equation (8), PDF p. 3.

### Task

The paper is about graph-based operator estimation, diffusion geometry, and spectral/eigenfunction recovery. Its numerical tasks are:

- estimating the generator of Ornstein-Uhlenbeck gradient-flow processes on `R` and `R^2`;
- estimating Laplacian eigenfunctions on non-uniformly sampled compact manifolds, specifically the circle and sphere;
- evaluating epsilon sensitivity and automatic bandwidth selection;
- validating continuum expansions by pointwise operator application in the appendices.

It is not a semi-supervised learning paper, but its results directly inform graph regression and graph signal smoothing because those methods depend on the graph operator induced by the chosen edge weights.

### Explicit Author Motivations

`[explicit]` The Abstract, PDF p. 1, says variable-bandwidth kernels are common in applications, but existing theory mostly covers fixed bandwidths.

`[explicit]` Section 1, PDF p. 2, says the main contribution is to extend Coifman-Lafon asymptotics and Singer's discrete analysis to variable-bandwidth kernels and to give rigorous error bounds.

`[explicit]` Section 1, PDF p. 2, says fixed-bandwidth errors become unbounded as sampling density approaches zero and that choosing bandwidth inversely proportional to density can control sparse-region errors.

`[explicit]` Section 2, PDF p. 4, says `beta<0` is expected to work best because it increases bandwidth in sparse regions and decreases bandwidth in dense regions.

`[explicit]` Section 5, PDF pp. 10-12, motivates variable bandwidth on compact manifolds by reduced epsilon sensitivity under non-uniform sampling.

### Derived Or Implied Motivations

`[derived]` For local-scale conductance, the paper implies that a density-adaptive graph can keep sparse-tail points connected without globally increasing epsilon. A global epsilon large enough for sparse regions would over-connect dense regions; a global epsilon small enough for dense regions can disconnect sparse regions.

`[derived]` The paper also implies that local bandwidth and diffusion-map density normalization should be considered jointly. Changing only edge conductance without changing normalization changes the limiting operator.

`[derived]` The finite-sample bounds imply that adaptive bandwidth can be necessary for convergence in operator-estimation settings on open/non-compact supports, not merely beneficial.

### Effect On Eigenfunctions / Diffusion / Smoothing

`[explicit]` The operator drift changes with `rho`: equation (6) adds `(d+2) grad f dot grad rho/rho`. Thus eigenfunctions and diffusion geometry can change if bandwidth is altered without compensating `alpha,beta` choices.

`[explicit]` For `rho=q^beta`, equation (8) says the density drift is governed by `c1`. Setting `c1=0` recovers the Laplacian and removes density bias; setting `c1=1` recovers the generator for the sampled gradient flow in the authors' examples.

`[explicit]` Figures 1-7 show that variable bandwidth improves eigenfunction recovery or epsilon stability in all tested non-uniform examples, with especially large benefits on non-compact Gaussian examples and sparse regions of the sphere.

`[derived]` For graph smoothing, the practical effect is a spatially varying smoothing length. In sparse regions, `rho=q^{-1/2}` makes the smoother borrow from farther neighbors; in dense regions, it remains local. The continuum effect is not just "more smoothing": it is a specific diffusion-plus-drift operator controlled by `alpha,beta,d`.

### Relationship To Adaptive-Scale Graph Construction

P04 is directly about adaptive-scale graph construction. It formalizes the same intuition as self-tuned spectral clustering, but in a diffusion-map/operator-estimation framework with continuum limits and finite-sample error terms.

For conductance comparators, the most portable formula is the Section 3 Gaussian:

```text
w_ij = exp( -||x_i-x_j||^2 / [4 epsilon rho_i rho_j] ).
```

If `rho_i=q_i^beta`, this becomes a density-adaptive conductance. If `beta=-1/2`, then sparse regions have larger local scale. The kernel itself is symmetric, but subsequent Markov normalizations and the final `rho_i^{-2}` generator scaling make the operator non-symmetric unless one uses the symmetric conjugate for eigensolvers.

### What The Paper Does Not Claim

The paper does not claim that every self-tuned kernel recovers the Laplacian. It shows the opposite: variable bandwidth adds drift unless `alpha,beta` are chosen to cancel density effects.

The paper does not prove spectral convergence for variable-bandwidth kernels, despite using eigenfunction experiments.

The paper does not handle non-compact manifolds with boundary.

The paper does not give a ready-made method for intrinsic dimension estimation; it only suggests an epsilon-tuning slope heuristic may help.

The paper does not define conductance for gflow or SIMODS directly. It defines kernel affinities and operators whose edge weights may be used as conductances in a graph construction.

### Relevance To gflow / SIMODS

P04 should be kept distinct from current `gflow` overlap-density smoothing semantics. The current overlap-density smoother should not be retroactively described as Berry-Harlim variable-bandwidth diffusion unless it actually uses the paper's `rho`, `q_epsilon^S`, `alpha`, and `rho^{-2}` generator normalization.

The paper is most relevant to planned length/kernel-conductance comparators:

- A local-length comparator can set `rho_i` from a pilot density estimate or kNN scale and use `w_ij=exp(-d_ij^2/(4 epsilon rho_i rho_j))`.
- A kernel-conductance comparator can compare fixed bandwidth, kNN/self-tuned bandwidth, and `rho=q^beta` bandwidth while tracking the induced operator normalization separately.
- A density-adaptive smoothing comparator must choose whether it wants Laplacian-like behavior (`c1=0`) or density-biased gradient-flow behavior (`c1>0`), because P04 shows that this choice is mathematical, not cosmetic.

For SIMODS narrative, P04 supports the idea that local scale is an operator-design parameter. It also warns that sparse-region robustness from larger bandwidth comes with a drift term and dimension-dependent normalization.

## Figure Handling

### Copied Paper Figures Used

Reproduced cropped figure panels for Figures 1--7 and Appendix Figures A.8--A.9 are managed by `/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/literature/conductance_kernel_laplacian_review/paper_figure_screenshots.yml` and embedded in the generated HTML memo next to the primary figure descriptions. These are internal-review cropped figure panels from the canonical reading copy, not manuscript-ready reused figures.

### Original Explanatory Figures Proposed Or Created

Created:

`/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/literature/conductance_kernel_laplacian_review/figures/P04_local_bandwidth_conductance_curve.png`

This original figure illustrates the Section 3 Gaussian implementation formula with a 1D Gaussian toy density and `rho=q^{-1/2}`. It shows relative density, local length scale `2 sqrt(epsilon) rho(x)`, and the approximate edge affinity for a fixed small displacement. It is meant to clarify the derived local-conductance implication: sparse regions receive larger bandwidth and therefore higher affinity for the same physical separation.

![Original internal explanatory figure: P04 local-bandwidth conductance curve](/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/literature/conductance_kernel_laplacian_review/figures/P04_local_bandwidth_conductance_curve.png)

## Evidence Table

| Claim | Label | Source reference | Notes |
| --- | --- | --- | --- |
| Fixed kernels use `K_epsilon(x,y)=h(||x-y||^2/epsilon)`. | explicit | Eq. (1), Section 1, PDF p. 1 | Baseline fixed bandwidth. |
| Variable kernels use `K_epsilon^S(x,y)=h(||x-y||^2/(epsilon rho(x)rho(y)))`. | explicit | Eq. (2), Section 1, PDF p. 1 | Symmetric local scale. |
| Diffusion-map fixed-bandwidth limit is `Delta f + (2-2alpha) grad f dot grad q/q`. | explicit | Eq. (3), Section 2, PDF p. 2 | Recovers Laplacian at `alpha=1`. |
| Theorem 1 gives pointwise high-probability errors for the discrete variable-bandwidth operator. | explicit | Theorem 1, Eq. (5), PDF p. 3 | Three error components. |
| Arbitrary `rho` adds `(d+2) grad f dot grad rho/rho` to the limiting operator. | explicit | Eq. (6), PDF p. 3; Eq. (A.20), PDF p. 22 | Core operator effect. |
| `rho=q^beta` yields `L_{alpha,beta}=Delta f+c1 grad f dot grad q/q`. | explicit | Corollary 1, Eq. (8), PDF p. 3 | `c1=2-2alpha+d beta+2 beta`. |
| Sparse-region errors are controlled by `c2`; fixed bandwidth can blow up as `q->0`. | explicit | Corollary 1 discussion, PDF p. 4 | For `beta=0`, `alpha>0`, `c2>0`. |
| `beta<0` increases bandwidth in sparse regions and decreases it in dense regions. | explicit | Section 2, PDF p. 4 | Author's stated intuition. |
| The practical Gaussian affinity is `exp(-||x_i-x_j||^2/(4 epsilon rho_i rho_j))`. | explicit | Section 3, PDF p. 5 | Conductance-ready formula. |
| For nearby points, local length scale is approximately `2 sqrt(epsilon) rho(x)`. | derived | Section 3 Gaussian formula, PDF p. 5 | Derived by setting exponent to `-1` and `rho(y) approx rho(x)`. |
| The pilot density estimate uses kNN scale `rho_0` and produces `q_0=q+O(...)`. | explicit | Eq. (10), Section 3, PDF p. 5 | Used to set `rho=q_0^beta`. |
| Sparse matrices are built with kNN truncation and symmetrized by `(K+K^T)/2`. | explicit | Section 3, PDF pp. 6-7 | Implementation detail. |
| 1D and 2D Ornstein-Uhlenbeck experiments show fixed bandwidth failing in sparse tails. | explicit | Figures 1-5, Section 4, PDF pp. 7-10 | Empirical validation of sparse-density error story. |
| Circle and sphere experiments show reduced epsilon sensitivity for variable bandwidth on compact non-uniform data. | explicit | Figures 6-7, Section 5, PDF pp. 11-12 | Compact-manifold evidence. |
| Variable bandwidth is necessary for convergence in some operator-estimation problems with densities not bounded away from zero. | explicit | Conclusion, PDF p. 13 | Author claim. |
| The paper does not prove spectral convergence. | explicit | Section 2, PDF p. 5; Section 6, PDF p. 13 | Only pointwise convergence is proved. |
| P04 informs planned SIMODS length/kernel-conductance comparators but is not identical to current overlap-density smoothing. | derived | Whole paper plus local project distinction | The paper gives kernel/operator formulas, not current gflow semantics. |

## Open Questions For Auditor

1. Should the memo treat `rho=q^beta` as a candidate SIMODS comparator directly, or should planned comparators use a kNN-scale `rho_i` closer to Zelnik-Manor/Perona and Cheng/Wu?
2. For local conductance, should the audit emphasize the raw symmetric affinity `w_ij` or the full generator normalization with `rho_i^{-2}` and `alpha` debiasing?
3. Does the project want the Laplacian target (`c1=0`) or a density-biased/gradient-flow target (`c1>0`) for planned smoothing comparators?
4. Should dimension `d` be fixed from known SIMODS geometry or estimated from data when using `beta != 0`?
5. Should future memos standardize a conductance notation that separates edge affinity, row-stochastic transition probability, and graph generator scaling?

## Revision Notes

### Post-Audit Revision, 2026-05-15

- Auditor-P04 flagged Eq. (4) denominator wording. The final synthesis should
  preserve the paper's displayed high-probability denominator convention and
  avoid translating it into a different square-root form unless the algebra is
  shown explicitly.
- The paper's \(\Delta\) sign convention should be stated wherever Eq. (6) or
  the generator \(L_{\alpha,\rho}\) is summarized, because sign conventions
  differ across graph Laplacian and PDE communities.
- SIMODS/gflow clarification for synthesis: P04 supports variable-bandwidth
  kernel-conductance comparators and density-normalized operator variants. It
  does not describe the current `fit.rdgraph.regression()` overlap-density /
  Riemannian-complex conductance.
- Any phrase such as "conceptual parent" for later adaptive smoothers should
  be labeled `contextual`, not attributed directly to Berry and Harlim.

- 2026-05-15: Drafted by Reviewer-P04 from the canonical P04 PDF. Rendered and visually inspected pages containing Figures 1-7 and Appendix Figures A.8-A.9. Created one original explanatory conductance/local-bandwidth figure.
