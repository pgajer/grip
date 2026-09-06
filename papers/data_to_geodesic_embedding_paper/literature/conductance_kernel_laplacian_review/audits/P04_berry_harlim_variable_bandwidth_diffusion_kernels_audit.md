# Paper Audit: Berry and Harlim, Variable Bandwidth Diffusion Kernels

Paper: Berry and Harlim, "Variable Bandwidth Diffusion Kernels"
Reviewer memo: `paper_memos/P04_berry_harlim_variable_bandwidth_diffusion_kernels.md`
Auditor: Auditor-P04
Status: draft
Date: 2026-05-15

## Audit Findings

### Formula And Reference Checks

The memo is broadly faithful to the canonical 26-page PDF. It correctly reports the fixed kernel in Eq. (1), the symmetric variable-bandwidth kernel in Eq. (2), the fixed-bandwidth diffusion-map limit in Eq. (3), Theorem 1's discrete operator and continuum limit in Eqs. (5)-(6), Corollary 1's `rho = q^beta + O(epsilon)` specialization in Eqs. (7)-(8), the Section 3 Gaussian implementation, and the Appendix A/B continuum and finite-sample derivations.

One notation correction is required in the fixed-bandwidth baseline around Eq. (4). The memo's equivalent form

```text
||grad f(x_i)|| / sqrt(N epsilon^{1+d/2})
```

is fine, because it equals `||grad f(x_i)|| / [sqrt(N) epsilon^{1/2+d/4}]`. But the following sentence should not say that the paper writes `sqrt(N epsilon^{1/2+d/4})`. The PDF writes the denominator as `sqrt(N) epsilon^{1/2+d/4}`.

The variable-bandwidth error terms are otherwise correct. The memo's two equivalent presentations of Eq. (5),

```text
q(x_i)^{1/2} rho(x_i)^{-d/2} / [sqrt(N) epsilon^{2+d/4}]
||grad f(x_i)|| q(x_i)^{-(1/2-2alpha+2dalpha)}
  rho(x_i)^{-(d/2+1)} / [sqrt(N) epsilon^{1/2+d/4}]
```

match Theorem 1 and Appendix B. The Corollary 1 coefficients `c1 = 2 - 2alpha + d beta + 2 beta` and `c2 = 1/2 - 2alpha + 2dalpha + d beta/2 + beta` are also correct.

The memo should add a small sign-convention note. The paper states that its `Delta` is the negative of the usual Laplace-Beltrami operator for convenience, and Section 3 describes it as having negative eigenvalues. The formulas are quoted correctly, but a reader could otherwise miss that `Delta` is being used with the paper's sign convention.

The Section 3 stochastic-process description is accurate: Eq. (9) is `dx = -c1 grad U(x) dt + sqrt(2) dW_t`, with invariant density proportional to `exp(-c1 U)`, and the backward generator is the operator in Eq. (8) under the chosen `alpha,beta`.

The appendix references are accurate. The memo correctly distinguishes the appendix-only left and right formulations from the main symmetric kernel; it also correctly records that the left formulation yields a Laplacian limit while right and symmetric formulations add the `(d+2) grad rho/rho dot grad f` drift.

### Missing Items

No major paper-local item is missing. The memo covers the abstract and motivation, main theorem, corollary, numerical implementation, sparse kNN truncation and symmetrization, all main figures, Appendix Figures A.8-A.9, finite-sample error interpretation, limitations, and the dimension requirement.

Two small additions would improve downstream synthesis:

- State explicitly that `d` is the intrinsic manifold dimension, not the ambient dimension and not automatically the graph dimension. The PDF uses `M subset R^n` with intrinsic dimension `d`; this `d` enters `rho^d`, `c1`, `c2`, and the tuning-slope heuristic.
- When discussing SIMODS comparators, name the current supported `fit.rdgraph.regression()` semantics: supplied `weight.list` values are positive graph edge lengths for local ordering/truncation, while the current smoothing conductance is overlap-density/Riemannian-complex based, approximately `c_e^rho = 1 / max(rho1(e), 1e-10)` in the current default regime. This is separate from direct inverse-length conductance and from P04 Gaussian/RBF variable-bandwidth affinities.

### Overclaims Or Ambiguities

The memo mostly avoids overclaiming. It correctly says that the paper proves pointwise convergence for smooth functions, not full spectral convergence, even though the experiments compare eigenvectors/eigenfunctions.

The phrase "conceptual parent of later kNN self-tuned graph Laplacian theory" is reasonable as synthesis for this literature review, but it is not a claim made inside P04. It should be marked or mentally treated as contextual cross-paper synthesis rather than paper evidence.

The conductance language is handled carefully. The memo correctly says P04 defines kernel affinities and normalized operators, not gflow/SIMODS conductances. For the open question about raw affinities versus generator-normalized operators, the audit answer is: keep both, but never collapse them. The raw symmetric affinity

```text
w_ij = exp(-d_ij^2 / [4 epsilon rho_i rho_j])
```

is the conductance-like edge-weight candidate for a planned comparator. The paper's actual operator, however, also includes `q_epsilon^S` alpha-debiasing, row normalization, and final generator scaling by `1/(epsilon rho_i^2)`. Claims about continuum limits, Laplacian recovery, or gradient-flow generators require the full normalized operator, not just the raw edge affinity.

The memo's statement that variable bandwidth can preserve effective connectivity in sparse regions is acceptable as a derived local-scale interpretation, especially because it is tied to the Gaussian formula and not presented as an author-defined conductance theorem.

### Evidence Label Corrections

Evidence labels are mostly sound.

Suggested minor label polish:

- Keep "P04 informs planned SIMODS length/kernel-conductance comparators" as `derived` or `contextual`, never `explicit`.
- Treat "conceptual parent of later kNN self-tuned graph Laplacian theory" as contextual synthesis if it remains in the memo.
- Absence claims such as "the paper does not define conductance for gflow or SIMODS directly" are best labeled `contextual` or `derived`, because they are audit judgments about scope and terminology.

The major formulas, figure descriptions, and theorem/corollary claims marked `explicit` are supported by the PDF.

### Figure Handling Checks

The memo reports no copied or screenshot paper figures, and I found only the original explanatory PNG:

```text
literature/conductance_kernel_laplacian_review/figures/P04_local_bandwidth_conductance_curve.png
```

The figure is a 1134 x 1260 PNG and is described as an original toy illustration of the Section 3 Gaussian formula with `rho = q^{-1/2}`. That is acceptable. The memo should keep it labeled as derived/original and not imply that the figure appears in Berry and Harlim.

The summaries of Figures 1-7 and Appendix Figures A.8-A.9 match the canonical PDF. The figure-page references are usable as internal review references.

### gflow / SIMODS Relevance Checks

The memo satisfies the main SIMODS separation requirement. It explicitly says current `gflow` overlap-density smoothing should not be retroactively described as Berry-Harlim variable-bandwidth diffusion unless it uses the paper's `rho`, `q_epsilon^S`, `alpha`, row normalization, and `rho^{-2}` generator scaling.

Required distinctions for synthesis:

- current `fit.rdgraph.regression()` overlap-density smoothing: separate from P04 and should be described with current project semantics, not P04 formulas;
- inverse-length conductance: a planned comparator, not P04 unless separately defined from lengths;
- Gaussian/RBF conductance: P04 supports a variable-bandwidth Gaussian affinity comparator;
- local/self-tuned kernels: central to P04, but P04's `rho = q^beta` theory is not identical to kNN self-tuning unless the pilot scale/density estimate and normalization are specified;
- graph construction support versus conductance weighting: P04's kNN truncation is an implementation sparsification step, while `K_epsilon^S` and its normalizations define weights/operators;
- row-normalized diffusion versus symmetric graph Laplacian smoothing: P04 constructs a row-normalized Markov/generator operator and then uses a symmetric conjugate for eigensolvers, not a plain unnormalized symmetric Laplacian.

For the open question about `d`, the audit answer is: do not silently choose it. For a P04-style comparator, either fix `d` from known SIMODS geometry or estimate it and report sensitivity, because `d` affects `rho^d`, `c1`, `c2`, and the choice of `alpha` for Laplacian-like versus gradient-flow-like behavior.

Current `fit.rdgraph.regression()` overlap-density smoothing should remain distinct from planned length/kernel-conductance comparators. P04 can justify designing a future density-adaptive Gaussian comparator; it should not be cited as evidence that the current overlap-density/Riemannian-complex smoother is already a variable-bandwidth diffusion-kernel method.

## Required Revisions

1. Fix the Eq. (4) denominator wording: the paper's denominator is `sqrt(N) epsilon^{1/2+d/4}`, equivalently `sqrt(N epsilon^{1+d/2})`, not `sqrt(N epsilon^{1/2+d/4})`.
2. Add one sentence noting the paper's `Delta` sign convention, since the PDF defines its Laplacian as the negative of the usual Laplace-Beltrami operator and later refers to negative eigenvalues.
3. Strengthen the SIMODS relevance section with one compact sentence naming current `fit.rdgraph.regression()` overlap-density/Riemannian-complex conductance semantics and keeping them separate from direct inverse-length and P04 Gaussian/RBF variable-bandwidth comparators.
4. If the "conceptual parent" phrasing remains, mark it as contextual synthesis rather than a direct P04 claim.

## Verdict

minor revisions
