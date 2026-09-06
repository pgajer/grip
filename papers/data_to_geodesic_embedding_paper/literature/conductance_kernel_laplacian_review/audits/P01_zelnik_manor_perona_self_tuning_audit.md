# Paper Audit: Self-Tuning Spectral Clustering

Paper: Zelnik-Manor, Lihi and Pietro Perona. "Self-Tuning Spectral Clustering." NeurIPS, 2004.
Reviewer memo: `paper_memos/P01_zelnik_manor_perona_self_tuning.md`
Auditor: Auditor-P01
Status: draft
Date: 2026-05-15

## Audit Findings

### Formula And Reference Checks

The memo accurately records the paper's central self-tuned affinity:
\[
\hat A_{ij}=\exp\left(-\frac{d^2(s_i,s_j)}{\sigma_i\sigma_j}\right),
\]
with \(\hat A_{ii}=0\) in the algorithm. The product \(\sigma_i\sigma_j\), not a sum, average, squared local bandwidth, or one-sided scale, is the correct denominator in Eq. (1).

The local-scale definition is also correct: Eq. (2) sets \(\sigma_i=d(s_i,s_K)\), where \(s_K\) is the \(K\)th neighbor of \(s_i\). The memo correctly states that all reported experiments use \(K=7\), while the paper gives no operational recipe beyond saying \(K\) is scale-independent and related to embedding dimension.

The graph operator is correctly described as the symmetric normalized affinity \(L=D^{-1/2}\hat A D^{-1/2}\), with \(D_{ii}=\sum_j \hat A_{ij}\). The memo properly avoids claiming that the paper analyzes the combinatorial Laplacian \(D-\hat A\), the random-walk operator, or a continuum graph-Laplacian limit.

The rotational-alignment section is broadly accurate. The paper defines \(Z=XR\), \(M_i=\max_j Z_{ij}\), and prints Eq. (3) as
\[
J=\sum_i\sum_j \frac{Z_{ij}^2}{M_i^2}.
\]
The memo's concern about Figure 5 is warranted: an exactly one-sparse row would contribute 1 to the printed raw sum, so the Figure 5 values near 0 cannot be the unshifted, unaveraged printed \(J\). Appendix A's derivative line includes a visible "-1" term and Figure 5 is labeled as alignment cost/quality, so the plotted quantity appears to be a shifted and likely normalized excess cost such as a per-row deviation from the one-sparse ideal. The paper does not spell this out. The memo should carry this forward as an ambiguity, not as a resolved implementation detail.

Figure references are accurate: Figure 1 demonstrates global-\(\sigma\) sensitivity; Figure 2 is the strongest kernel/conductance illustration; Figure 3 reports automatic group-number results; Figure 4 supports the critique of eigenvalue-gap selection; Figure 5 plots the alignment criterion for the top-row Figure 3 data sets; Figure 6 shows automatic image segmentation.

### Missing Items

No major missing paper content was found. The memo covers the kernel, local scale, normalized operator, rotation-based cluster-number selection, final non-maximum assignment, limitations, and figure set.

One small synthesis-facing addition would help: explicitly state that P01 is a local-bandwidth weighting method, not a cKNN/adaptive-support graph construction method. It estimates \(\sigma_i\) from nearest-neighbor distances, but then uses those scales in a dense affinity matrix unless a later implementation sparsifies it. This is important for the final concept map's distinction between adaptive-radius/cKNN graph support and conductance weighting.

### Overclaims Or Ambiguities

The memo is appropriately conservative about continuum theory, electrical conductance language, and gflow transfer. I found no substantive overclaim that P01 proves smoothing consistency, Laplace-Beltrami convergence, or regression behavior.

The only ambiguity requiring attention is the Eq. (3)/Figure 5 cost scale. The memo already marks this as uncertain in the figures section and open questions. For the final synthesis, do not cite Figure 5's near-zero y-axis values as the literal printed \(J\); cite only the ranking/selection rule unless the implementation convention is separately verified.

### Evidence Label Corrections

Evidence labels are mostly correct. The conductance translation \(c_{ij}=\hat A_{ij}\) is correctly labeled derived, because the paper uses "affinity" language rather than physical conductance. Claims about current `fit.rdgraph.regression()` are correctly contextual/derived rather than paper-explicit.

Suggested label refinement: classify the statement that P01 is relevant to adaptive/cKNN graph construction as derived/contextual unless it is phrased narrowly as "nearest-neighbor local bandwidth selection." The paper uses a nearest-neighbor distance to set bandwidths; it does not introduce a cKNN graph-support rule.

### Figure Handling Checks

The memo reports that no copied paper figures were added. That is compliant. If copied figures are later used, Figure 2 is the best candidate and must be labeled as reproduced from Zelnik-Manor and Perona (2004), Figure 2, internal review only unless permission or replacement is obtained.

The original toy figure `figures/P01_self_tuned_kernel_toy.png` is sufficient for internal explanation of Eq. (1) versus a global Gaussian baseline. It is original, not a copied paper figure, and its caption makes clear that edge thickness encodes Gaussian affinity. For the final report, the toy figure can stand in for the paper's Figure 2 unless the report specifically needs a historical visual from the source paper.

### gflow / SIMODS Relevance Checks

The memo correctly distinguishes current `fit.rdgraph.regression()` overlap-density/Riemannian-complex smoothing from a planned P01-style comparator. It should remain explicit that the current smoother uses overlap-density semantics, while P01 supplies an adaptive Gaussian edge-weight formula that could inspire a new length/kernel-conductance comparator.

The audit template's distinctions are satisfied with one caveat:

- current `fit.rdgraph.regression()` overlap-density smoothing: distinguished correctly;
- inverse-length conductance: not part of P01 and not conflated;
- Gaussian/RBF conductance: correctly treated as affinity-to-conductance translation;
- local/self-tuned kernels: central and accurately described;
- graph construction support versus conductance weighting: needs the minor clarification above, because nearest-neighbor scale estimation is not itself cKNN support construction;
- row-normalized diffusion versus symmetric graph Laplacian smoothing: the memo correctly treats P01's \(L=D^{-1/2}\hat A D^{-1/2}\) as a symmetric normalized affinity and does not rebrand it as the current gflow smoother.

## Required Revisions

Minor revisions recommended before synthesis:

1. Move the Eq. (3)/Figure 5 ambiguity from "Open Questions" into the main formula or figure discussion as a standing caution: Figure 5 appears to plot a shifted/normalized alignment cost, while Eq. (3) is printed as the raw sum.
2. Add one sentence in the adaptive-scale section distinguishing nearest-neighbor bandwidth selection from cKNN/adaptive graph-support construction.

No reread is required.

## Verdict

minor revisions
