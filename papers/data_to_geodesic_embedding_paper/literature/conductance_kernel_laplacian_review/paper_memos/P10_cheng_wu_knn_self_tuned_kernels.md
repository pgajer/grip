# Paper Review: kNN Self-Tuned Graph Laplacian Convergence

Paper: Cheng, Xiuyuan; Wu, Hau-Tieng. "Convergence of Graph Laplacian with kNN Self-tuned Kernels." arXiv:2011.01479v2 [math.ST], 2021.
Reviewer: Reviewer-P10
Auditor: Auditor-P10
Status: revised after audit
Date: 2026-05-15
Source manifest ID: P10
Canonical reading copy: `literature/conductance_kernel_laplacian_review/sources/pdf/P10_cheng_wu_knn_self_tuned_kernels.pdf`

## Whole-Paper Review

### Reader Background Needed

A mathematically capable reader should be comfortable with:

- Weighted graphs, affinity matrices, degree matrices, unnormalized graph Laplacians `D - W`, and random-walk graph Laplacians `I - D^{-1}W`.
- Kernel graph construction from point clouds, especially Gaussian or radial kernels of the form `k0(||x_i-x_j||^2/epsilon)`.
- k-nearest-neighbor distances and their use as local scale or density estimators.
- Riemannian manifolds embedded in Euclidean space, including volume measure, gradients, the Laplace-Beltrami operator, and weighted Laplacians.
- Diffusion maps and density normalization: how kernel normalization changes the continuum drift term.
- Dirichlet forms such as `int p |grad f|^2 dV` and graph analogues such as `f^T(D-W)f`.
- High-probability finite-sample error bounds, bias-variance tradeoffs, and joint limits where `N -> infinity`, `epsilon -> 0`, and `k -> infinity`.

### What A Non-Expert Should Understand Before Reading This Paper

The paper studies a practical graph-building trick that is common in spectral clustering: instead of using one global kernel bandwidth for all points, each data point receives a local bandwidth equal to a kNN distance. Sparse regions get larger local bandwidths, dense regions get smaller local bandwidths, and edge weights are computed using the product of the two endpoint scales.

The hard question is not whether this works heuristically, but what differential operator the resulting graph Laplacian approximates as the sample size grows. Fixed-bandwidth graph Laplacian convergence was already well studied. For kNN self-tuned bandwidths, the scale function is random, nonsmooth, and has derivatives that do not converge to the derivatives of the population scale. Cheng and Wu's contribution is to prove convergence anyway, with rates, by first proving uniform relative `C0` consistency of the kNN bandwidth estimate and then conditioning graph-Laplacian analysis on that event.

The most important practical message is that local bandwidths reduce variance problems in low-density regions. For fixed-bandwidth kernels, pointwise variance terms worsen where the sampling density `p(x)` is small. For the kNN self-tuned kernels in this paper, the corresponding variance factor is favorable in low-density regions because the bandwidth expands there.

### Problem And Context

The paper begins from the standard kernelized affinity

```text
W_ij = k0(||x_i - x_j||^2 / epsilon)
```

in equation (1), PDF p. 1. Given `W`, the standard graph Laplacians are `D - W` and `I - D^{-1}W`, with `D_ii = sum_j W_ij` (Section 1, PDF pp. 1-2). The authors recall that, for samples from a smooth density `p` on a smooth compact `d`-dimensional manifold `M`, graph Laplacians can converge to intrinsic manifold operators such as the Laplace-Beltrami operator `Delta_M` or the weighted Laplacian

```text
Delta_p = Delta_M + (grad_M p / p) dot grad_M
```

from equation (2), PDF p. 2.

The specific gap is self-tuned bandwidth theory. Zelnik-Manor and Perona's self-tuned spectral clustering used

```text
W_ij = k0(||x_i - x_j||^2 / (R_i R_j))
```

where `R_i` is a kNN distance; this is equation (3), PDF p. 2. Cheng and Wu state that convergence results for this type of graph Laplacian were incomplete, especially when the bandwidth function is estimated from data rather than assumed known (Section 1, PDF pp. 2-3; Section 1.2, PDF pp. 5-6).

### Main Method

The paper introduces a one-parameter family of self-tuned affinities:

```text
W_ij^(alpha) =
  k0(||x_i - x_j||^2 / (epsilon rho_hat(x_i)rho_hat(x_j)))
  / (rho_hat(x_i)^alpha rho_hat(x_j)^alpha).
```

This is equation (4), PDF p. 2. The local bandwidth function `rho_hat` is estimated from a stand-alone sample `Y`, while the graph itself is built on sample `X`. This split reduces dependence in the proof and gives a clean two-step analysis (Section 1.1, PDF p. 5).

The theoretical pipeline is:

1. Estimate the kNN bandwidth function `rho_hat` and prove uniform relative convergence to `rho_bar = p^{-1/d}`. This is Section 2, with the main result in Theorem 2.3, PDF p. 8.
2. Condition on the good event that `sup_M |rho_hat-rho_bar|/rho_bar < epsilon_rho`.
3. Prove convergence of the kernelized Dirichlet form, graph Dirichlet form, and graph Laplacian operators. These are Proposition 3.2 and Theorems 3.3, 3.5, 3.6, and 3.7, PDF pp. 10-14.
4. Compare with fixed-bandwidth, density-normalized kernels in Theorems 3.4 and 3.8, PDF pp. 12 and 14.

Algorithm 1, PDF pp. 4-5, gives the practical construction using the raw kNN distances `R_hat_i`:

```text
W_ij = k0(||x_i - x_j||^2 / (sigma_0^2 R_hat_i R_hat_j))
       / (R_hat_i^alpha R_hat_j^alpha).
```

The algorithm computes either `L_un` through `eig(D-W)` or the modified random-walk operator through the generalized eigenproblem `eig-gen(D-W, D D_Rhat^2)`.

### Main Formulas And Operators

[explicit] Fixed-bandwidth baseline: equation (1), PDF p. 1,

```text
W_ij = k0(||x_i - x_j||^2 / epsilon).
```

[explicit] Original self-tuned affinity: equation (3), PDF p. 2,

```text
W_ij = k0(||x_i - x_j||^2 / (R_hat_i R_hat_j)).
```

[explicit] Cheng-Wu self-tuned family: equation (4), PDF p. 2,

```text
W_ij^(alpha) =
  k0(||x_i - x_j||^2 / (epsilon rho_hat(x_i)rho_hat(x_j)))
  / (rho_hat(x_i)^alpha rho_hat(x_j)^alpha).
```

[explicit] kNN bandwidth estimator: equation (6), PDF p. 7,

```text
rho_hat(x) = R_hat(x) * ( (1/m0[h]) * (k/N_y) )^{-1/d},
R_hat(x) = inf{ r > 0 : sum_{j=1}^{N_y} 1_{||y_j-x||<r} >= k }.
```

Here `R_hat(x)` is the distance from `x` to its k-th nearest neighbor in `Y`, and `m0[h] = int_R^d h(|u|^2) du`; for `h=1_[0,1)`, `m0[h]` is the unit `d`-ball volume (Section 2.1, PDF p. 7).

[explicit] Population bandwidth:

```text
rho_bar = p^{-1/d}.
```

This appears in the Abstract, Table 1 on PDF p. 3, Section 2, and Theorem 2.3.

[explicit] Limiting operator for a smooth positive bandwidth `rho`: equation (9), PDF p. 9,

```text
L_rho^(alpha) =
  Delta + 2 (grad p / p) dot grad
        + (d - 2 alpha + 2) (grad rho / rho) dot grad.
```

[explicit] Limiting operator after substituting `rho_bar = p^{-1/d}`: equation (10), PDF p. 9,

```text
L^(alpha) =
  Delta + (1 + 2(alpha-1)/d) (grad p / p) dot grad.
```

[explicit] Graph Dirichlet form: equation (11), PDF p. 10,

```text
E_N(f,f) =
  (2/(epsilon m2 N^2)) epsilon^{-d/2} f^T(D-W)f
  = (1/(epsilon m2 N^2)) sum_{i,j} epsilon^{-d/2} W_ij (f_i-f_j)^2.
```

[explicit] Kernelized Dirichlet-form kernel: equations (12)-(13), PDF p. 10,

```text
K_hat(x,y) =
  epsilon^{-d/2} k0(||x-y||^2 / (epsilon rho_hat(x)rho_hat(y)))
  / (rho_hat(x)^alpha rho_hat(y)^alpha),

E^(alpha)(f,f) =
  (1/(epsilon m2)) int_M int_M
    (f(x)-f(y))^2 K_hat(x,y) p(x)p(y) dV(x)dV(y).
```

[explicit] Self-tuned unnormalized graph Laplacian operator: equation (17), PDF p. 12,

```text
L_un^(alpha) f(x) =
  (2 epsilon^{-d/2-1}/(m2 rho_hat(x)^alpha))
  (1/N) sum_j k0(||x-x_j||^2/(epsilon rho_hat(x)rho_hat(x_j)))
       (f(x_j)-f(x))/rho_hat(x_j)^alpha.
```

[explicit] Modified random-walk graph Laplacian operator: equation (18), PDF p. 13,

```text
L_rw0^(alpha) f(x) =
  (2/(epsilon m2 rho_hat(x)^2))
  [ sum_j k0(...) f(x_j)/rho_hat(x_j)^alpha
    / sum_j k0(...)/rho_hat(x_j)^alpha
    - f(x) ].
```

The paper notes that this differs from the usual random-walk matrix by multiplying another diagonal matrix `D_rho_hat^{-2}` up to constants and sign (Section 3.4, PDF p. 13).

### Figures And Experiments

The paper contains both theory-supporting synthetic experiments and an MNIST embedding demonstration.

- Figure 1, PDF p. 15, compares kNN estimation of `rho_bar=p^{-1/d}` to fixed-bandwidth KDE estimation of `p` on a circle. The top row shows `rho_hat` and its relative error across `k_y`; the bottom row shows KDE `p_hat` and relative error across `epsilon`. The figure supports Theorem 2.3 and Remark 2.2: kNN `rho_hat` has more spatially uniform relative error, while KDE relative error worsens where `p` is small.
- Figure 2, PDF p. 15, defines the 1D synthetic manifold in `R^4`, its nonuniform density `p`, the test function `f`, and the target `Delta_p f`. This supplies the testbed for Figures 3-6.
- Figure 3, PDF p. 16, plots relative errors for the Dirichlet form and pointwise `L_N f` errors versus `epsilon`, with `N_y=4000` and `k_y=32,256`. It shows pointwise `L_N f` errors growing at small `epsilon` as variance dominates, while the Dirichlet form is less sensitive. This visually matches Theorems 3.3 and 3.5.
- Figure 4, PDF p. 17, shows single-realization `L_N f` estimates using estimated `rho_hat` and population `rho_bar` for two `k_y` values and two `epsilon` values. It makes the bias-variance tradeoff concrete: smaller `epsilon` oscillates more; larger `epsilon` biases more; larger `k_y` smooths `rho_hat` but can bias low-density regions.
- Figure 5, PDF p. 17, repeats the visualization using `rho_hat_X` estimated from the graph sample `X` itself rather than stand-alone `Y`. This addresses the practical self-tuned setup even though the main theory uses a stand-alone `Y`.
- Figure 6, PDF p. 18, compares Dirichlet-form and `Err_infty` errors using `rho_bar`, `rho_hat_X`, and `rho_hat_Y` over multiple `N_y,k_y` choices. Section 4.3 concludes that a large stand-alone `Y` can reduce bandwidth-estimation error and oscillation, but splitting a limited dataset may worsen performance by reducing the graph sample size.
- Figure 7, PDF p. 18, shows recovery of Laplace-Beltrami eigenfunctions on `S^1`. The paper tests `alpha=1-d/2` and a mixed `rho_hat,p_hat` normalization from equation (21). The random-walk variants look better visually than unnormalized variants in this example; the authors defer fuller random-walk self-tuned study.
- Figure 8, PDF p. 19, compares 3D eigenvector embeddings of 1000 MNIST images under fixed-bandwidth diffusion-map normalization and two self-tuned kernels from equation (22). Self-tuned kernels remain informative across a range of `sigma_0`.
- Figure 9, PDF p. 19, identifies outliers in the fixed-bandwidth embedding as samples with large `R_hat_i`, supporting the claim that fixed bandwidth is unstable for low-density points.
- Figure 10, PDF p. 20, illustrates the hyperplane/polygon partition used to prove Lemma 2.1: on each polygon, the kNN point is fixed and `R_hat(x)=||x-y_p||`.
- Figure 11, PDF p. 21, illustrates the proof of Proposition 2.2 by plotting empirical `mu_hat(x,r)`, population `mu(x,r)`, and the local approximation `m0[h]p(x)r^d`, with `R_hat`, `R_bar`, and bracketing radii `R_-/R_+`.

Table 1, PDF p. 3, is a useful notation table. Table 2, PDF p. 6, is a roadmap linking bandwidth estimation to the convergence theorems for Dirichlet forms, pointwise operators, and weak forms.

### Theoretical Claims

[explicit] Under Assumption 2.1, `M` is a smooth compact boundaryless `d`-manifold embedded in `R^D`, and `p` is smooth and bounded above and below (PDF p. 7). Under Assumption 3.1, `k0` is nonnegative, sufficiently smooth, and exponentially decaying through four derivatives (PDF p. 9).

[explicit] Lemma 2.1, PDF p. 7, proves that \(\widehat R\) is globally Lipschitz with Lipschitz constant at most 1 and is \(C^\infty\) off a finite union of hyperplanes.

[explicit] Proposition 2.2, PDF pp. 7-8, gives a fixed-point relative error:

```text
|rho_hat(x)-rho_bar(x)|/rho_bar(x)
  = O^[p]((k/N)^(2/d)) + O(sqrt(log N/k))
```

up to constants stated in the theorem. The theorem's displayed bound is

```text
O^[p]((k/N)^(2/d)) + sqrt(3 s log N)/(d sqrt(k)).
```

[explicit] Theorem 2.3, PDF p. 8, upgrades this to uniform high-probability relative error:

```text
sup_x |rho_hat(x)-rho_bar(x)|/rho_bar(x)
  = O^[p]((k/N)^(2/d)) + sqrt(39 log N)/(d sqrt(k))
```

with probability at least `1 - N^{-10}` for large `N`.

[explicit] Remark 2.1, PDF p. 8, balances the bias `(k/N)^(2/d)` and variance `sqrt(log N/k)` terms, giving `k ~ N^(1/(1+d/4))` up to logs and `epsilon_rho ~ N^(-1/(2+d/2))`.

[explicit] Remark 2.2, PDF p. 8, contrasts kNN `rho_hat` with fixed-bandwidth KDE. KDE relative error has a variance factor proportional to `p(x)^(-1/2)`, while kNN `rho_hat` has a uniform variance term independent of `p(x)`.

[explicit] Section 2.3, PDF pp. 8-9, proves a major technical warning: although `rho_hat` is `C0` consistent, its derivative diverges. At differentiability points, `|grad_bar R_hat(x)|=1`, so

```text
|grad_bar rho_hat(x)| = ((1/m0[h]) k/N_y)^(-1/d),
```

which diverges like `(k/N_y)^(-1/d)` when `k/N_y -> 0`. Thus `grad rho_hat` is not pointwise consistent for `grad rho_bar`.

[explicit] Proposition 3.2, PDF p. 10, proves kernelized Dirichlet-form convergence:

```text
E^(alpha)(f,f) =
  E_{p_alpha}(f,f)(1 + O^[alpha](epsilon_rho)) + O^[f,p](epsilon),
p_alpha = p^(1 + 2(alpha-1)/d).
```

[explicit] Theorem 3.3, PDF p. 11, proves graph Dirichlet-form convergence under `epsilon=o(1)`, `epsilon^(d/2) N = Omega(log N)`, and `epsilon_rho=o(1)`:

```text
E_N(f,f) =
  E_{p_alpha}(f,f)(1 + O^[alpha](epsilon_rho))
  + O^[f,p](epsilon)
  + O^[1]( sqrt( (log N)/(N epsilon^(d/2))
                 int_M |grad f|^4 p^(1+4(alpha-1)/d) ) ).
```

[explicit] Section 3.3, PDF p. 11, records the important special cases: `alpha=1` gives `p_alpha=p` and recovers `Delta_p`; `alpha=0` gives `p_0=p^(1-2/d)` and corresponds to the original self-tuned graph's modified density; `alpha=1-d/2` gives constant `p_alpha` and recovers `Delta_M`.

[explicit] Theorem 3.5, PDF p. 13, gives pointwise convergence for the modified random-walk operator under `epsilon=o(1)`, `epsilon^(d/2+1)N=Omega(log N)`, and `epsilon_rho=o(epsilon)`:

```text
L_rw0^(alpha) f(x) =
  L^(alpha) f(x)
  + O^[f,p](epsilon, epsilon_rho/epsilon)
  + O^[1](||grad f||_infty p(x)^(1/d)
          sqrt(log N/(N epsilon^(d/2+1)))).
```

[explicit] Theorem 3.6, PDF p. 14, gives the analogous pointwise result for the unnormalized operator:

```text
L_un^(alpha) f(x) =
  p(x)^(2(alpha-1)/d) L^(alpha) f(x)
  + O^[f,p](epsilon, epsilon_rho/epsilon)
  + O^[1](||grad f||_infty p(x)^((2alpha-1)/d)
          sqrt(log N/(N epsilon^(d/2+1)))).
```

[explicit] Theorem 3.7, PDF p. 14, gives weak convergence for `L_un^(alpha)` with error `O(epsilon, epsilon_rho)` plus a variance term of order `sqrt(log N/(N epsilon))`. It removes the `epsilon_rho/epsilon` term but applies only in weak form and only to the unnormalized operator.

[explicit] Theorem 3.8, PDF p. 14, gives the fixed-bandwidth random-walk comparison. Its variance factor is `p(x)^(-1/2)`, whereas Theorem 3.5 has the self-tuned factor `p(x)^(1/d)`. The authors explicitly interpret this as the theoretical advantage of self-tuned bandwidths in low-density regions (end of Section 3.4, PDF p. 14).

### Limitations And Scope

[explicit] The main theory assumes compact, smooth, boundaryless manifolds and smooth density bounded away from zero and infinity (Assumption 2.1, PDF p. 7). This is strong relative to many real datasets with boundaries, singularities, mixed dimensions, or outliers.

[explicit] The graph-theory proof uses a stand-alone `Y` sample to estimate `rho_hat`; using `X` itself is discussed and tested, but not the main theorem setting (Section 1.1, PDF p. 5; Section 4.3, PDF pp. 17-18).

[explicit] The pointwise convergence theorem requires `epsilon_rho=o(epsilon)`, while Dirichlet-form convergence only requires `epsilon_rho=o(1)` (Theorems 3.3 and 3.5, PDF pp. 11 and 13). This is a real distinction for practical operator estimation.

[explicit] Section 2.3, PDF pp. 8-9, shows that the kNN bandwidth's derivatives do not converge. The paper therefore cannot simply apply smooth variable-bandwidth theory to the random kNN scale.

[explicit] The authors defer fuller study of the random-walk graph Laplacian with self-tuned kernel after the Laplace-Beltrami eigenfunction experiment (Section 4.4, PDF p. 18).

[contextual] The paper proves operator and Dirichlet-form convergence, not spectral convergence of eigenvectors/eigenvalues with rates. It uses eigenvector experiments, but the main theorem set is not a spectral convergence theorem.

[contextual] The paper is about kernel weights on sample clouds, not graph conductance in the electrical-network sense. For our review, the affinity can be interpreted as a conductance-like edge weight, but that is our synthesis rather than the paper's terminology.

### Historical / Methodological Importance

This paper is a direct theoretical bridge between the practical self-tuning kernel of Zelnik-Manor and Perona and the variable-bandwidth diffusion-kernel theory represented by Berry-Harlim-style smooth bandwidth analysis. Its key methodological move is to handle the actual kNN-estimated scale rather than assuming a smooth density-derived bandwidth is already available.

For adaptive graph construction, the paper clarifies three issues that are easy to conflate:

- kNN local scale estimation gives `rho_hat -> p^{-1/d}` in relative `C0` error, with better low-density variance behavior than fixed-bandwidth KDE.
- The normalization exponent `alpha` determines which weighted Laplacian the graph approximates.
- The kNN scale is nonsmooth, so pointwise graph-Laplacian convergence is harder and has stricter rate requirements than Dirichlet-form convergence.

## Conductance / Kernel Extraction

### Conductance, Affinity, Or Kernel Formula(s)

[explicit] Standard fixed-bandwidth kernel, equation (1), PDF p. 1:

```text
W_ij = k0(||x_i-x_j||^2 / epsilon).
```

[explicit] Original self-tuned kernel, equation (3), PDF p. 2:

```text
W_ij = k0(||x_i-x_j||^2 / (R_hat_i R_hat_j)).
```

[explicit] Cheng-Wu self-tuned family, equation (4), PDF p. 2:

```text
W_ij^(alpha) =
  k0(||x_i-x_j||^2 / (epsilon rho_hat(x_i)rho_hat(x_j)))
  / (rho_hat(x_i)^alpha rho_hat(x_j)^alpha).
```

[explicit] Algorithm 1 implementation kernel, equation (5), PDF p. 5:

```text
W_ij =
  k0(||x_i-x_j||^2 / (sigma_0^2 R_hat_i R_hat_j))
  / (R_hat_i^alpha R_hat_j^alpha).
```

The paper explains that `sigma_0 R_hat_i = sqrt(epsilon) rho_hat(x_i)` and

```text
epsilon = sigma_0^2 * ((1/m0[h]) * (k/N_y))^(2/d).
```

Thus the practical algorithm does not need to know `d`, even though the theory uses the normalized `rho_hat` (Section 1.1, PDF p. 5).

[explicit] Mixed normalization for Laplace-Beltrami recovery without knowing `d`, equation (21), PDF p. 18:

```text
W_ij =
  k0(||x_i-x_j||^2 / (epsilon rho_hat(x_i)rho_hat(x_j)))
  / ((rho_hat p_hat^(1/2))(x_i) (rho_hat p_hat^(1/2))(x_j))
  = W_ij^(1) / (p_hat(x_i)^(1/2) p_hat(x_j)^(1/2)).
```

[explicit] MNIST self-tuned kernels, equation (22), PDF p. 19:

```text
W_ij^(1) =
  k0(||x_i-x_j||^2 / (sigma_0^2 R_hat_i R_hat_j))
  / (sigma_0^2 R_hat_i R_hat_j),

W_ij^0 =
  k0(||x_i-x_j||^2 / (sigma_0^2 R_hat_i R_hat_j))
  / (sigma_0^2 R_hat_i R_hat_j sqrt(mu_hat_i mu_hat_j)).
```

### Graph, Laplacian, Or Diffusion Operator

[explicit] Basic graph Laplacians are `D-W` and `I-D^{-1}W`, with `D_ii=sum_j W_ij` (Section 1, PDF p. 1).

[explicit] The practical algorithm computes eigenpairs using either `eig(D-W)` for `L_un` or `eig-gen(D-W, D D_Rhat^2)` for the modified random-walk operator (Algorithm 1, PDF p. 5).

[explicit] The normalized graph Dirichlet form is equation (11), PDF p. 10:

```text
E_N(f,f) =
  (1/(epsilon m2 N^2)) sum_{i,j}
    epsilon^{-d/2} W_ij (f_i-f_j)^2.
```

[explicit] The limiting operator family is equation (10), PDF p. 9:

```text
L^(alpha) =
  Delta + (1 + 2(alpha-1)/d) (grad p/p) dot grad.
```

[explicit] The associated weighted density for the Dirichlet form is

```text
p_alpha = p^(1 + 2(alpha-1)/d).
```

This is stated in Proposition 3.2 and Theorem 3.3, PDF pp. 10-11.

### Task

[explicit] The paper's tasks are graph-based geometric data analysis, unsupervised learning, spectral clustering, dimension-reduced embedding, and diffusion/operator approximation. These applications are named in Section 1, PDF pp. 1-2.

[explicit] The numerical tasks are kNN bandwidth estimation (Section 4.1), Dirichlet-form and `L_N f` estimation (Section 4.2), stand-alone `Y` comparison (Section 4.3), Laplace-Beltrami eigenfunction recovery (Section 4.4), and MNIST eigenvector embedding (Section 4.5).

### Explicit Author Motivations

[explicit] The paper motivates self-tuning because global bandwidth choice is challenging in high-dimensional or uneven-density data, and graph-Laplacian performance can be sensitive to that choice (Section 1, PDF p. 2).

[explicit] The authors state that existing convergence results for self-tuned kernels were incomplete when the bandwidth function is estimated from data (Section 1, PDF p. 2; Section 1.2, PDF pp. 5-6).

[explicit] The paper aims to show convergence of graph Laplacian operators and graph Dirichlet forms with rates for kNN self-tuned kernels (Abstract and Section 1, PDF pp. 1, 3-4).

[explicit] A major motivation is low-density robustness: self-tuned bandwidths enlarge in sparse regions, reducing variance relative to fixed bandwidths. This is stated in the Abstract, discussed after Theorem 3.8 on PDF p. 14, and demonstrated in Figures 1 and 8-9.

### Derived Or Implied Motivations

[derived] The stand-alone `Y` sample is partly a proof device and partly a practical memory/computation option: if more samples are available than can be used to build the graph, using the extra samples for bandwidth estimation may improve `rho_hat` without enlarging `W`. This inference follows from Section 1.1, PDF p. 5, and Section 4.3, PDF pp. 17-18.

[derived] The `alpha` parameter is a conductance-normalization knob. Increasing or decreasing `alpha` does not just rescale weights; it changes the continuum density in the Dirichlet form from `p` to `p_alpha=p^(1+2(alpha-1)/d)`. This follows from Proposition 3.2 and Theorem 3.3, PDF pp. 10-11.

[derived] The paper suggests that evaluating a local-scale graph smoother requires checking both the edge-length scale and the operator normalization. A graph can be locally adaptive but still converge to a density-biased operator if `alpha` is not chosen for the intended limit.

### Effect On Eigenfunctions / Diffusion / Smoothing

[explicit] For `alpha=1`, `L^(alpha)=Delta_p`, so the graph approximates the weighted Laplacian associated with the sampling measure `p dV` (Abstract; Section 3.3, PDF p. 11).

[explicit] For `alpha=1-d/2`, `p_alpha` is constant and the Dirichlet form corresponds to the Laplace-Beltrami operator `Delta_M` (Section 3.3, PDF p. 11; Section 4.4, PDF p. 18).

[explicit] For `alpha=0`, the original self-tuned graph approximates a weighted Laplacian with modified density `p^(1-2/d)` (Section 3.3, PDF p. 11).

[explicit] Figure 7, PDF p. 18, shows that the first four nontrivial graph eigenvectors approximate sine/cosine Laplace-Beltrami eigenfunctions on `S^1`, with random-walk variants visually better in that example.

[explicit] The modified random-walk pointwise variance term in Theorem 3.5 has factor `p(x)^(1/d)`, while the fixed-bandwidth comparison in Theorem 3.8 has `p(x)^(-1/2)` (PDF p. 14). This means the self-tuned operator should be less noisy in sparse regions.

[derived] For smoothing, the Dirichlet-form result is more directly relevant to energy penalties such as `f^T L f`; it has an `O(epsilon_rho)` bandwidth-estimation contribution rather than the stricter `O(epsilon_rho/epsilon)` pointwise operator contribution. This follows from Theorems 3.3 and 3.5.

### Relationship To Adaptive-Scale Graph Construction

[explicit] The local scale is the kNN radius `R_hat_i`, and the edge scale is the geometric product `R_hat_i R_hat_j` or normalized version `rho_hat(x_i)rho_hat(x_j)` (equations (3)-(6), PDF pp. 2, 5, 7).

[explicit] The population local scale is `rho_bar=p^{-1/d}` (Abstract; Table 1, PDF p. 3; Theorem 2.3, PDF p. 8). Thus denser regions have smaller bandwidths and sparser regions have larger bandwidths.

[explicit] The paper connects directly to original self-tuning spectral clustering and nearest-neighbor density estimation in Section 1.2, PDF pp. 5-6.

[contextual] Relative to smooth variable-bandwidth papers, this paper handles the empirically common kNN scale. It is not merely a special case of smooth variable bandwidth, because Section 2.3 shows derivative inconsistency.

[derived] For a conductance comparator, equation (4) can be read as a local-scale conductance rule: edge weight is a decreasing function of squared distance divided by local endpoint scales, with an additional endpoint amplitude normalization. The distance denominator and the `alpha` normalization should be treated as separate design choices.

### What The Paper Does Not Claim

[explicit] The paper does not claim that kNN `rho_hat` is `C1` consistent. It proves the opposite in Section 2.3, PDF pp. 8-9.

[explicit] The main theory does not require knowing `d` in Algorithm 1, but the theoretical statement and normalization of `rho_hat` use `d` (Section 1.1, PDF p. 5; equation (6), PDF p. 7).

[explicit] The main proofs do not analyze the fully dependent case where the same finite sample `X` is both graph sample and bandwidth-estimation sample, although Section 4.3 studies it empirically.

[contextual] The paper does not claim that all self-tuned kernels recover `Delta_M` or `Delta_p`. The limit depends on `alpha`, and recovering `Delta_M` requires `alpha=1-d/2` or the mixed normalization in equation (21).

[contextual] The paper does not develop SIMODS/gflow, Mapper/Riemannian complexes, overlap-density smoothing, or electrical network conductance models.

### Relevance To gflow / SIMODS

[contextual] Current `gflow` overlap-density smoothing should remain distinct from Cheng-Wu self-tuned graph kernels. The current `fit.rdgraph.regression()` semantics are based on overlap-density/Riemannian-complex graph structure, not on Euclidean point-cloud kNN radii or the self-tuned kernel in equation (4).

[derived] P10 is highly relevant for planned length/kernel-conductance comparators. It gives a principled candidate edge-weight family:

```text
conductance_ij ~
  k0(length_ij^2 / (epsilon local_scale_i local_scale_j))
  / (local_scale_i^alpha local_scale_j^alpha),
```

where the population interpretation of `local_scale` is density-adaptive length. In a SIMODS comparator, `local_scale_i` could be defined from kNN distances, local cell density, or graph-neighborhood lengths, but that would be a new comparator, not a reinterpretation of the current overlap-density smoother.

[derived] The paper's separation between Dirichlet-form convergence and pointwise operator convergence is useful for SIMODS. If the planned comparator is used as a smoothing penalty, the Dirichlet-form theory is the closer analogy. If it is used as a pointwise Laplacian operator, the stricter pointwise conditions and `epsilon_rho/epsilon` term become more relevant.

[derived] The density-effect result suggests a testable SIMODS hypothesis: fixed global length kernels may over-penalize or disconnect sparse regions, while local-scale kernels may stabilize smoothing across nonuniform sampling. This should be tested as a comparator against the existing overlap-density smoother rather than blended into its current semantics.

## Figure Handling

### Copied Paper Figures Used

Reproduced cropped figure panels for Figures 1--11 are managed by `paper_figure_screenshots.yml` and embedded in the generated HTML memo next to the primary figure descriptions. These are internal-review cropped figure panels from the canonical reading copy, not manuscript-ready reused figures.

### Original Explanatory Figures Proposed Or Created

No original explanatory figure was created for P10. A useful synthesis figure, if desired later, would compare fixed bandwidth, kNN local scale, and `alpha` endpoint normalization as three separate knobs in a conductance/kernel comparator.

## Evidence Table

| Claim | Label | Source reference | Notes |
| --- | --- | --- | --- |
| The paper studies graph Laplacian convergence for kNN self-tuned kernels with estimated bandwidths. | explicit | Abstract, PDF p. 1; Section 1, PDF pp. 2-4 | Central paper claim. |
| Standard fixed-bandwidth affinity is `W_ij=k0(||x_i-x_j||^2/epsilon)`. | explicit | Eq. (1), PDF p. 1 | Baseline kernel. |
| Original self-tuned affinity uses kNN distances in the denominator product `R_hat_i R_hat_j`. | explicit | Eq. (3), PDF p. 2 | Cites Zelnik-Manor and Perona. |
| Cheng-Wu's family divides by `rho_hat_i^alpha rho_hat_j^alpha`. | explicit | Eq. (4), PDF p. 2 | `alpha` selects limiting operator. |
| The practical algorithm uses raw kNN distances `R_hat_i` and parameter `sigma_0`. | explicit | Algorithm 1 and Eq. (5), PDF pp. 4-5 | Does not require knowing `d`. |
| The kNN bandwidth estimator is normalized from the k-th neighbor distance and targets `rho_bar=p^{-1/d}`. | explicit | Eq. (6), PDF p. 7; Theorem 2.3, PDF p. 8 | Theoretical normalization uses `d`. |
| `R_hat` is Lipschitz and piecewise smooth off hyperplanes. | explicit | Lemma 2.1, PDF p. 7; Fig. 10, PDF p. 20 | Used for uniform convergence proof. |
| Uniform relative error of `rho_hat` is bounded by bias `(k/N)^(2/d)` plus variance `sqrt(log N/k)`. | explicit | Theorem 2.3, PDF p. 8 | High-probability `C0` result. |
| Optimal kNN scaling is `k ~ N^(1/(1+d/4))` up to logs. | explicit | Remark 2.1, PDF p. 8 | Balances bias and variance. |
| kNN `rho_hat` has better low-density relative variance behavior than fixed-bandwidth KDE. | explicit | Remark 2.2, PDF p. 8; Fig. 1, PDF p. 15 | KDE relative variance includes `p(x)^(-1/2)`. |
| `rho_hat` is not `C1` consistent; its derivative diverges. | explicit | Section 2.3, PDF pp. 8-9 | Key technical obstacle. |
| The limiting random-walk operator is `L^(alpha)=Delta+(1+2(alpha-1)/d) grad p/p dot grad`. | explicit | Eq. (10), PDF p. 9 | After substituting `rho_bar=p^{-1/d}`. |
| The graph Dirichlet form converges to the differential Dirichlet form for `p_alpha=p^(1+2(alpha-1)/d)`. | explicit | Proposition 3.2 and Theorem 3.3, PDF pp. 10-11 | Relevant to smoothing penalties. |
| `alpha=1` recovers `Delta_p`; `alpha=1-d/2` recovers `Delta_M`; `alpha=0` gives modified density `p^(1-2/d)`. | explicit | Section 3.3, PDF p. 11 | Important normalization map. |
| Pointwise convergence requires `epsilon_rho=o(epsilon)`. | explicit | Theorem 3.5, PDF p. 13; Theorem 3.6, PDF p. 14 | Stricter than Dirichlet-form convergence. |
| Weak convergence of `L_un` removes the `epsilon_rho/epsilon` term. | explicit | Theorem 3.7, PDF p. 14 | Only weak form and unnormalized operator. |
| Fixed-bandwidth pointwise variance worsens in low-density regions; self-tuned variance improves there. | explicit | Theorem 3.8 and following paragraph, PDF p. 14 | Fixed factor `p^(-1/2)` versus self-tuned `p^(1/d)` in Theorem 3.5. |
| Stand-alone `Y` can improve bandwidth estimation when many extra samples are available, but splitting limited data may hurt. | explicit | Section 4.3, PDF pp. 17-18; Figs. 5-6 | Practical caveat. |
| MNIST fixed-bandwidth embeddings show outliers with large kNN distances; self-tuned kernels are more stable. | explicit | Section 4.5, PDF pp. 19-20; Figs. 8-9 | Empirical support for low-density robustness. |
| P10's affinity can serve as a planned SIMODS length/kernel-conductance comparator. | derived | Eq. (4), PDF p. 2; SIMODS review context | Comparator only; not current `gflow` overlap-density semantics. |
| Current `fit.rdgraph.regression()` overlap-density smoothing is distinct from this paper's kNN self-tuned point-cloud kernel. | contextual | SIMODS/gflow review context; P10 methods throughout | Important scope boundary for H005. |

## Open Questions For Auditor

1. Should the synthesis treat P10 primarily as "kNN self-tuned kernel theory" or as a direct sequel to the P04 smooth variable-bandwidth diffusion-kernel memo?
2. For SIMODS, should the proposed comparator use the raw Algorithm 1 form with `R_hat_i` and `sigma_0`, or the theory-normalized `rho_hat` form with explicit dependence on estimated intrinsic dimension?
3. Should the final review emphasize Dirichlet-form convergence as the closer analogue to smoothing penalties, while marking pointwise operator convergence as a stricter and separate use case?
4. If a length/kernel-conductance comparator is implemented, should `alpha=1` be the default because it recovers `Delta_p`, or should the comparator expose `alpha` explicitly because `alpha` controls density bias?
5. Should the MNIST outlier result in Figures 8-9 be cited as empirical evidence for low-density stabilization, or kept secondary because it is an embedding demonstration rather than a controlled convergence experiment?

## Revision Notes

### Post-Audit Revision, 2026-05-15

- Auditor-P10 requested that the final synthesis distinguish kNN bandwidth
  estimation from graph support. P10 studies kNN/self-tuned scales and graph
  Laplacian convergence; do not imply that bandwidth estimation and support
  construction are the same design choice.
- SIMODS/gflow clarification for synthesis: current
  `fit.rdgraph.regression()` uses overlap-density/Riemannian-complex
  conductance. P10 supports planned kNN self-tuned kernel-conductance
  comparators, not the existing gflow operator semantics.
- Add an operator-normalization caveat for special \(\alpha\) choices. The
  meaning of \(\alpha=1\) or other values depends on the paper's operator and
  normalization, so avoid transferring a default between P03/P04/P10 without
  restating the construction.
- Treat the MNIST/outlier result as secondary empirical evidence, not as a
  convergence proof.

- 2026-05-15: Initial Reviewer-P10 memo drafted from the 60-page canonical PDF. Full text, theorem statements, proof figures, and experimental figures/tables were reviewed.
