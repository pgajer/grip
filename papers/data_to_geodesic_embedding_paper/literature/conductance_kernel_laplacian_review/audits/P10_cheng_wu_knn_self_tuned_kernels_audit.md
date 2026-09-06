# Paper Audit: Cheng and Wu, kNN Self-Tuned Graph Laplacian Convergence

Paper: Cheng and Wu, "Convergence of Graph Laplacian with kNN Self-tuned Kernels"
Reviewer memo: `paper_memos/P10_cheng_wu_knn_self_tuned_kernels.md`
Auditor: Auditor-P10
Status: draft
Date: 2026-05-15

## Audit Findings

### Formula And Reference Checks

The memo is broadly faithful to the canonical 60-page PDF. It correctly records the fixed-bandwidth baseline in Eq. (1), the Zelnik-Manor/Perona self-tuned affinity in Eq. (3), Cheng-Wu's alpha-family in Eq. (4), Algorithm 1 / Eq. (5), the normalized kNN bandwidth estimator in Eq. (6), the limiting operators in Eqs. (9)-(10), the graph and kernelized Dirichlet forms in Eqs. (11)-(13), and the two pointwise graph-Laplacian operators in Eqs. (17)-(18).

The kNN local-scale description is correct. The paper's theoretical `rho_hat` is the k-th neighbor distance from the stand-alone sample `Y`, rescaled by `((1/m0[h]) k/N_y)^(-1/d)`, and targets `rho_bar = p^(-1/d)`. The memo also correctly notes the practical Algorithm 1 form using raw `R_hat_i` and `sigma_0`, which avoids requiring `d` in implementation.

The convergence assumptions and rates are accurately summarized. Theorem 2.3 requires `k=o(N)` and `k=Omega(log N)` and gives uniform relative `C0` error with bias `(k/N)^(2/d)` plus variance `sqrt(log N/k)`. Theorems 3.3, 3.5, 3.6, and 3.7 are reported with the correct qualitative conditions: Dirichlet-form convergence only needs `epsilon_rho=o(1)`, pointwise convergence needs `epsilon_rho=o(epsilon)`, and weak convergence of `L_un` removes the `epsilon_rho/epsilon` penalty.

One small precision note should be carried into synthesis: statements like "`alpha=1-d/2` recovers `Delta_M`" are correct for the limiting `L^(alpha)` / modified random-walk operator and for the associated Dirichlet form, but Theorem 3.6 says the unnormalized pointwise operator converges to `p(x)^(2(alpha-1)/d) L^(alpha) f(x)`. The memo does include this formula, so this is not a formula error, but later prose should avoid saying every graph operator recovers `Delta_M` without naming the normalization.

The low-density variance comparison is correct. The memo's `p(x)^(1/d)` versus `p(x)^(-1/2)` comparison matches Theorems 3.5 and 3.8 and the authors' paragraph at the end of Section 3.4. The Figure 1 comparison between kNN `rho_hat` relative error and fixed-bandwidth KDE relative error is also correctly tied to Remark 2.2.

### Missing Items

No major paper-local item is missing. The memo covers the abstract, motivation, assumptions, theorem roadmap, the kNN estimator, derivative divergence of `rho_hat`, Dirichlet-form and pointwise operator results, fixed-bandwidth comparisons, all main figures, proof figures, and the MNIST demonstration.

Two synthesis-facing clarifications should be added or preserved explicitly:

- P10 uses kNN distances to estimate bandwidths; this is not automatically a kNN graph-support rule. Edge support is determined by the kernel support, dense/sparse implementation, or any later sparsification. This matters for the audit-template distinction between graph construction support and conductance weighting.
- The current `fit.rdgraph.regression()` precomputed-graph path treats `weight.list` as positive edge lengths for local ordering/truncation, while the smoothing conductance is overlap-density/Riemannian-complex based, approximately `c_e^rho = 1 / max(rho1(e), 1e-10)` in the current default regime. P10 supports a planned local-scale kernel/conductance comparator, not a reinterpretation of that current smoother.

### Overclaims Or Ambiguities

The memo avoids major overclaims. It correctly states that P10 proves pointwise operator convergence and Dirichlet-form convergence, not full spectral convergence with rates. The eigenfunction and MNIST sections are treated as experiments and demonstrations rather than theorem-backed spectral guarantees.

The conductance language is appropriately cautious. P10 defines affinities, graph Laplacians, and Dirichlet forms; translating `W_ij` into conductance-like weights for SIMODS is a derived comparator design. That translation should keep separate the distance denominator `epsilon local_scale_i local_scale_j`, the amplitude normalization `local_scale_i^alpha local_scale_j^alpha`, and the downstream choice of unnormalized versus modified random-walk operator.

For the reviewer's open questions, the audit answers are:

1. Treat P10 primarily as kNN self-tuned kernel theory. It is a natural bridge from P01 and P04, but it is not merely a direct sequel to smooth variable-bandwidth theory because Section 2.3 proves the kNN scale is not `C1` consistent.
2. For a practical SIMODS comparator, the Algorithm 1 raw `R_hat_i, sigma_0` form is the cleanest implementation analogue. The theory-normalized `rho_hat` form should be retained for analysis and for explaining density/operator limits. If a continuum target such as `Delta_M` is desired, either estimate/specify intrinsic `d` or use a mixed normalization like Eq. (21), with all choices reported.
3. Yes: emphasize Dirichlet-form convergence as the closer analogue for smoothing penalties such as `f^T L f`. Pointwise operator convergence is stricter and should remain a separate use case.
4. Expose `alpha` explicitly in any comparator. `alpha=1` is a defensible default if the target is `Delta_p`, but `alpha` controls density bias; `alpha=1-d/2` or Eq. (21)-style normalization is relevant for `Delta_M`.
5. Cite the MNIST result only as secondary empirical evidence for low-density stabilization. The controlled convergence evidence is in the theorem comparison and synthetic experiments.

### Evidence Label Corrections

The evidence labels are mostly sound. The main equations, assumptions, theorem statements, figure descriptions, and paper-scope limitations marked `explicit` are supported by the PDF.

Suggested minor label discipline:

- Keep "P10 is a theoretical bridge between P01 and P04" as contextual synthesis, not as a direct author claim.
- Keep "P10 affinity can serve as a SIMODS conductance comparator" as `derived`; the paper says affinity/weights, not electrical conductance or SIMODS conductance.
- Keep absence claims about SIMODS, Mapper/Riemannian complexes, and overlap-density smoothing as contextual audit judgments.

### Figure Handling Checks

The memo reports no copied or screenshot paper figures, and I found no P10 figure artifact added to the review repository. This is compliant.

The figure references are accurate enough for internal review: Figures 1-2 introduce kNN/KDE estimation and the synthetic data; Figures 3-6 support Dirichlet-form versus pointwise error and stand-alone `Y`; Figure 7 is the Laplace-Beltrami eigenfunction demonstration; Figures 8-9 are the MNIST stability/outlier demonstration; Figures 10-11 are proof illustrations for Lemma 2.1 and Proposition 2.2. If a later synthesis figure is created, it should be labeled as original/derived and should not imply it is reproduced from Cheng and Wu.

### gflow / SIMODS Relevance Checks

The memo satisfies the core SIMODS separation requirement, with the clarifications above recommended before final synthesis:

- current `fit.rdgraph.regression()` overlap-density smoothing: kept separate from P10's Euclidean/manifold kNN self-tuned kernel;
- inverse-length conductance: a planned comparator, not P10 and not current `rdgraph`;
- Gaussian/RBF conductance: P10 supports a local-scale kernel-weight comparator, with `k0` such as a Gaussian;
- local/self-tuned kernels: central to P10 and accurately described;
- graph construction support versus conductance weighting: should be stated more explicitly, because kNN estimates bandwidths and does not by itself define graph support;
- row-normalized diffusion versus symmetric graph Laplacian smoothing: the memo correctly distinguishes the modified random-walk operator, unnormalized operator, and Dirichlet form.

Current `fit.rdgraph.regression()` overlap-density smoothing must remain distinct from planned length/kernel-conductance comparators. P10 can justify testing a local-scale kernel comparator against the current overlap-density smoother; it should not be cited as evidence that the current smoother already implements kNN self-tuned kernel conductance.

## Required Revisions

1. Add one sentence distinguishing kNN bandwidth estimation from kNN/adaptive graph support: P10's kNN radius sets local scale, while edge support comes from the kernel support or implementation sparsification.
2. Add one sentence in the SIMODS relevance section naming current `fit.rdgraph.regression()` semantics: `weight.list` supplies positive edge lengths for local ordering/truncation, and current smoothing conductance is overlap-density/Riemannian-complex based, not direct inverse-length or P10 local-scale kernel conductance.
3. Add a small operator-normalization caveat wherever special `alpha` values are summarized: `alpha=1-d/2` recovers `Delta_M` for `L^(alpha)` / modified random-walk and the Dirichlet-form target, while the unnormalized pointwise operator includes the density prefactor from Theorem 3.6.
4. Keep the MNIST Figures 8-9 claim secondary: empirical stability/outlier evidence, not controlled convergence or spectral-convergence proof.

## Verdict

minor revisions
