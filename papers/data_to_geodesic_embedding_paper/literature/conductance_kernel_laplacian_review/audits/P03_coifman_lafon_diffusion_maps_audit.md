# Paper Audit: Coifman & Lafon, Diffusion Maps

Paper: Coifman, Ronald R. and Lafon, Stephane. "Diffusion Maps." Applied and Computational Harmonic Analysis 21(1):5--30, 2006. DOI: 10.1016/j.acha.2006.04.006.
Reviewer memo: `literature/conductance_kernel_laplacian_review/paper_memos/P03_coifman_lafon_diffusion_maps.md`
Auditor: Auditor-P03
Status: audit_ready
Date: 2026-05-15

## Audit Findings

### Formula And Reference Checks

- Kernel and Markov normalization are accurate. The memo correctly records the symmetric nonnegative base kernel `k(x,y)`, degree `d(x)=integral k(x,y) dmu(y)`, row-normalized transition kernel `p(x,y)=k(x,y)/d(x)`, and conservation `integral p(x,y)dmu(y)=1` from Section 2.1.
- The memo's note on the Section 2.1 displayed operator is appropriate. The PDF extraction shows `Pf(x)=integral a(x,y)f(y)dmu(y)`, but the surrounding text has just defined `p` as the Markov transition kernel; treating this as `p(x,y)` in the memo is the mathematically consistent reading.
- Stationarity, reversibility, and eigenpair statements are accurate for the finite connected/ergodic setting and the Appendix A compact/symmetrizable setting: `pi(y)=d(y)/sum_z d(z)`, detailed balance, `P psi_l=lambda_l psi_l`, and `1=lambda_0>|lambda_1|>=...` as stated in the paper.
- Diffusion distance is accurately stated as the weighted `L^2(X,dmu/pi)` distance between `t`-step transition densities, and the eigen-expansion with `lambda_l^{2t}(psi_l(x)-psi_l(y))^2` is correct. The truncation rule `|lambda_l|^t > delta |lambda_1|^t` and diffusion map coordinates `lambda_l^t psi_l(x)` match Section 2.4.
- Appendix A symmetrization is consistent with the PDF formula `a(x,y)=sqrt(pi(x)/pi(y))p(x,y)=k(x,y)/(sqrt(pi(x))sqrt(pi(y)))`. If reused in synthesis, remember this is a symmetrizing kernel for spectral analysis, not the row-stochastic transition kernel itself.
- Section 3.1 alpha normalization is correctly captured: `q_epsilon(x)=integral k_epsilon(x,y)q(y)dy`, `k_epsilon^(alpha)(x,y)=k_epsilon(x,y)/(q_epsilon(x)^alpha q_epsilon(y)^alpha)`, then a second degree normalization produces `p_{epsilon,alpha}` and `P_{epsilon,alpha}`.
- Theorem 2 is quoted accurately at the level used in the memo. For synthesis prose, it may be worth noting that the displayed limit is the operator on `f`, while the paper also describes the conjugate Schrödinger operator on `phi=f q^{1-alpha}`.
- The alpha special cases are correctly assigned: `alpha=0` is the classical normalized graph Laplacian on isotropic weights with maximal density influence; `alpha=1/2` connects to Fokker-Planck/Langevin dynamics; `alpha=1` recovers the Laplace-Beltrami/Neumann heat-kernel limit.
- Directed kernel formula in Section 4 is accurately summarized as an application-specific anisotropic kernel favoring diffusion along level sets of `f`.

### Missing Items

- No major missing formula or section was found for H005. The memo covers Sections 2--6, Appendix A, and the Section 3/Appendix B asymptotic role at sufficient depth for the conductance/kernel review.
- The finite-sample rates from Section 5 are summarized rather than fully transcribed. That is acceptable for this memo because P03 is mainly being used for normalization/eigenbasis/comparator logic; detailed convergence-rate auditing can be left to P08/P09 or a later synthesis pass.
- The memo does not create a separate original figure file. This is acceptable: the inline Mermaid pipeline is sufficient for a review memo unless H005 later requests a rendered cross-paper synthesis figure.

### Overclaims Or Ambiguities

- Minor: the Figure 3 bullet is labeled `[explicit]` but includes the phrase "main figure for graph low-pass smoothing." The figure explicitly shows decay of numerical rank under powers of `P`; "graph low-pass smoothing" is a derived graph-signal interpretation. Split that sentence or mark the low-pass clause as `[derived]` if the reviewer revises.
- Minor: the "generic conductance-like edge weight" language is acceptable because the memo immediately says the paper does not use electrical conductance terminology. Do not let this become a direct claim that Coifman-Lafon defines electrical-network conductances.
- The paper's robustness-to-noise language in Section 2.4 is correctly reported, but synthesis should keep it qualitative unless paired with Section 5's scale condition `sqrt(epsilon)` larger than perturbation size.
- The Section 3.4 sentence "Finally, when alpha = 0" should be treated as an apparent typo in the canonical PDF, not as a formal erratum unless an external erratum is found. The heading, Proposition 3, limit `L_{epsilon,1}`, subsequent sentence "By setting alpha = 1", Figure 4, and Appendix B all support `alpha=1`.

### Evidence Label Corrections

- Keep the low-pass/filtering language as `[derived]`, not `[explicit]`, except where quoting the paper's own wording about low-frequency content of leading eigenfunctions and increasing oscillation deeper in the spectrum.
- The gflow/SIMODS mapping is correctly labeled `[derived]`; it is project interpretation, not a claim made by Coifman and Lafon.
- The Section 3.4 typo note is correctly labeled `[uncertain]`; the evidence is internal consistency of the paper, not an author correction.

### Figure Handling Checks

- No copied paper figures were included in the memo. This passes the figure-handling requirement.
- Figure references 1--6 match the canonical PDF pages and section roles. The memo's Figure 4 description correctly distinguishes `alpha=0` graph-Laplacian embeddings from `alpha=1` Laplace-Beltrami embeddings.
- The inline Mermaid pipeline is original explanatory material, not a reproduced paper figure.

### gflow / SIMODS Relevance Checks

- The memo correctly distinguishes current `fit.rdgraph.regression()` overlap-density smoothing from Coifman-Lafon diffusion maps. This distinction should remain mandatory in synthesis.
- Current supported SIMODS/gflow semantics, per the project operator note, are: supplied `weight.list` values are positive edge lengths for local neighborhood ordering/truncation; the current smoothing conductance is overlap-density/Riemannian-complex based, approximately `c_e^rho = 1/max(rho1(e), 1e-10)` in the default regime; it is not a direct inverse-length conductance, Gaussian/RBF affinity, row-stochastic transition matrix, or diffusion-map kernel.
- The most direct Coifman-Lafon comparator for SIMODS is a row-stochastic diffusion smoother `P^t` built from a chosen symmetric nonnegative length/kernel affinity. A symmetric normalized Laplacian version is a closely related implementation/eigensolver form via reversible symmetrization, but it should be named separately from row-normalized diffusion. Diffusion-coordinate truncation is best treated as an embedding/low-rank representation comparator unless explicitly used as a smoother.
- A density-sensitivity comparator is scientifically useful: compare `alpha=0` and `alpha=1` on the same Euclidean/kernel graph family when enough information exists to estimate `q_epsilon`. This should be a planned comparator, not a reinterpretation of current overlap-density smoothing.
- Local/self-tuned kernels belong to P01/P04/P10-style extensions. P03 supplies the normalization and diffusion-map baseline, not the local bandwidth rule.

## Required Revisions

Minor revisions recommended before final synthesis use:

1. Split or relabel the Figure 3 sentence so rank decay is `[explicit]` and "low-pass smoothing" is `[derived]`.
2. In the SIMODS relevance section, add one sentence using the current project wording: `fit.rdgraph.regression()` uses edge lengths for neighborhood ordering/truncation and an overlap-density/Riemannian-complex conductance, not direct length or Gaussian-kernel conductance.
3. Preserve the Section 3.4 `alpha=0` issue as an internal "apparent PDF typo" note; do not cite it as a formal erratum.

## Verdict

minor revisions
