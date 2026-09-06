# Paper Audit: Zhu, Ghahramani, and Lafferty, Gaussian Fields and Harmonic Functions

Paper: Zhu, Xiaojin; Ghahramani, Zoubin; Lafferty, John. "Semi-Supervised Learning Using Gaussian Fields and Harmonic Functions." ICML 2003.
Reviewer memo: `paper_memos/P05_zhu_ghahramani_lafferty_gaussian_fields.md`
Auditor: Auditor-P05
Status: draft
Date: 2026-05-15

## Audit Findings

### Formula And Reference Checks

The memo is broadly faithful to the canonical 8-page PDF. It correctly reports the graph setup, binary boundary labels, Gaussian feature-space affinity in Eq. (1), quadratic energy in Eq. (2), Gaussian-field density under clamped labels, harmonic mean-value equation in Eq. (3), block solve in Eq. (5), random-walk hitting-probability view, electrical-network view, CMN rule in Eq. (9), dongle/external-classifier formula in Eq. (10), entropy objective and gradients in Eqs. (11)-(14), CMN probability in Eq. (15), and text-document weight in Eq. (16).

One formula correction is required. Section 3.2 prints the full graph heat kernel as

```text
K_t = e^{-t Delta}
```

and the restricted Dirichlet heat kernel as

```text
K'_t = e^{-t Delta_uu}.
```

The memo currently says the rendered PDF uses `K_t = exp(t Delta)` for the full graph. That is incorrect. The key inverse formula in Eq. (6), `G = integral_0^infty e^{-t Delta_uu} dt = (D_uu - W_uu)^-1`, is otherwise correctly stated.

Add a small factor-convention caution for the normalized-cut comparator. The paper's Eq. (8) prints

```text
R(f) = f^T Delta f / f^T D f
     = sum_ij w_ij (f(i)-f(j))^2 / sum_i d_i f(i)^2.
```

With symmetric `W` and the usual double-sum convention, `f^T Delta f = (1/2) sum_ij w_ij(f_i-f_j)^2`. The missing factor is irrelevant to the Rayleigh quotient/eigenproblem but should be flagged for H005 formula reuse, especially because Eq. (2) correctly includes the `1/2` in the energy.

The block linear system is correct: with labeled nodes first, `Delta_uu = D_uu - W_uu`, and the Dirichlet solution is `f_u = (D_uu - W_uu)^-1 W_ul f_l`, equivalently `(I-P_uu)^-1 P_ul f_l`. The memo's sign and block orientation are sound.

### Missing Items

No major paper-local item is missing for H005. The memo covers the central method, all equivalence interpretations, main algorithmic extensions, figures, experiments, limitations, and relevance to conductance/kernel choices.

Two small additions would improve final synthesis use:

- State explicitly that the paper's heat-kernel discussion uses the positive combinatorial Laplacian `Delta = D-W` together with `e^{-t Delta}`.
- In the SIMODS section, add one compact sentence distinguishing the unnormalized Dirichlet Laplacian solve from its row-normalized Markov equivalent. Both appear in the paper; neither is current `fit.rdgraph.regression()` overlap-density smoothing.

### Overclaims Or Ambiguities

The memo is appropriately conservative about gflow. It does not claim that Zhu et al. use overlap-density/Riemannian-complex conductance, inverse-length conductance, or current SIMODS operators.

The conductance language is source-supported but should be used with care. The paper explicitly says in Section 3.1 to imagine graph edges as resistors with conductance `W`. Therefore it is acceptable to call `w_ij` conductance in the electrical-network/Dirichlet-energy interpretation. In formula tables and graph-construction prose, prefer "weight/affinity `w_ij`, interpretable as conductance" so the paper's broader terminology is preserved.

The multi-label discussion is correctly limited. The PDF states that CMN extends naturally to the general multi-label case and that the harmonic solution remains efficiently computable in the multi-label case, but the displayed derivations and formulas are mostly binary. The final report should include this as a short scope note, not expand the binary formulas into unsupported detailed multi-class algebra.

The memo's historical importance and gflow transfer statements are synthesis claims and are already mostly labeled `contextual` or `derived`. Keep them that way.

### Evidence Label Corrections

Most evidence labels are sound.

Required label/wording corrections:

- Correct the heat-kernel row and any associated evidence entry: the printed full graph formula is `K_t=e^{-t Delta}`, not `e^{t Delta}`.
- Add a factor-convention note to the normalized-cut/Eq. (8) row if Eq. (8) is used as a formula source.
- Treat "weight/kernel selection controls graph interpolation/regression" as `derived`, as the memo already does; the paper explicitly shows dependence on `W` but does not use gflow/SIMODS regression language.
- Keep absence claims about manifold convergence, gflow overlap-density conductance, and modern comparator families as `contextual` audit/synthesis judgments.

### Figure Handling Checks

The memo reports no copied or screenshot paper figures, and I found no copied source figure in the P05 figure handling. This passes the audit requirement.

The original explanatory SVG exists at:

```text
figures/P05_three_node_harmonic_interpolation.svg
```

It is an original three-node specialization of Eq. (3), not copied from Zhu et al. The figure is acceptable as an internal explanatory aid. If it is used in the final report, keep it labeled as an original derived schematic and not as a paper figure.

The paper figure summaries match the PDF at audit granularity: Figure 2 shows synthetic three-band/two-spiral harmonic minimization; Figure 3 contains digits "1" vs. "2", all 10 digits, and odd-vs-even with VP combination; Figure 4 contains the three 20-newsgroups binary tasks; Figure 5 shows the entropy/smoothing toy parameter-learning example; Figure 6 shows learned digit feature scales.

### gflow / SIMODS Relevance Checks

The memo satisfies the main SIMODS separation requirement. It distinguishes:

- current `fit.rdgraph.regression()` overlap-density/Riemannian-complex smoothing from Zhu et al.'s Gaussian and text affinities;
- inverse-length conductance as a planned comparator, not a P05 formula;
- Gaussian/RBF conductance as source-supported only through Eq. (1)-style feature affinities;
- local/self-tuned kernels as absent from P05, aside from global feature-wise `sigma_d` learning;
- graph construction support from conductance weighting, especially for the sparse 10-nearest-neighbor text graph;
- row-normalized Markov propagation `P=D^-1 W` from the symmetric/combinatorial Laplacian Dirichlet solve `Delta=D-W`.

For H005 synthesis, use P05 mainly as a finite-graph Dirichlet interpolation and conductance-as-electrical-weight reference. Do not cite it as evidence for current overlap-density smoothing or for local adaptive-bandwidth graph construction.

### Open Questions Resolved

1. Sign convention: use the PDF's printed `K_t=e^{-t Delta}` and `K'_t=e^{-t Delta_uu}` with `Delta=D-W`.
2. Conductance wording: call `w_ij` a weight/affinity generally; call it conductance when invoking the electrical-network interpretation, which is explicit in Section 3.1.
3. Three-node SVG: acceptable as an original derived schematic; include only if the final report wants a small explanatory figure in house style.
4. Multi-label handling: mention the paper's brief extension claims, but keep detailed formulas binary unless another source is added.
5. gflow separation: the memo's contextual paragraph is appropriate; preserve the separation between current overlap-density smoothing and planned length/kernel-conductance comparators.

## Required Revisions

1. Correct the Section 3.2 heat-kernel sign statement: the full graph heat kernel is printed as `K_t=e^{-t Delta}`, not `e^{t Delta}`.
2. Add a short factor-convention note for Eq. (8): under the standard symmetric double-sum identity, `f^T Delta f` carries a `1/2`, while the PDF prints the unhalved numerator in the Rayleigh quotient.
3. In the final synthesis-facing wording, use "weight/affinity, interpretable as conductance" outside the electrical-network paragraph, and reserve plain "conductance" for the Section 3.1 interpretation or H005-derived comparator language.
4. Add one compact SIMODS clarification that P05's row-normalized Markov view and unnormalized Laplacian Dirichlet solve are paper-equivalent views, not the current `fit.rdgraph.regression()` overlap-density smoother.

## Verdict

minor revisions
