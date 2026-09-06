# H005 Execution Plan: Conductance And Kernel Choices For Graph Laplacians

Built 2026-05-15 08:16:46 EDT

## Scope

This is a long standalone literature review for the SIMODS/gflow graph
low-pass smoothing benchmark. It is not a citation inventory. The review must
explain how conductance and kernel choices
\[
  c_{ij}=\phi(\ell_{ij};\theta)
\]
shape graph Laplacian eigenfunctions, diffusion, low-pass smoothing,
semi-supervised learning, manifold learning, and graph-signal filtering.

The final report must distinguish the current `fit.rdgraph.regression()`
overlap-density/Riemannian-complex smoother,
\[
  c_e^\rho = 1/\max(\rho_1(e),10^{-10}),
\]
from planned length/kernel conductance comparators such as
\[
  c_e^{\mathrm{len}} = 1/(\ell_e+\epsilon).
\]

## Tier 1 Papers

Each Tier 1 paper has one dedicated reviewer-auditor pair.

| ID | Paper | Pair |
| --- | --- | --- |
| P01 | Zelnik-Manor & Perona, *Self-Tuning Spectral Clustering* | Reviewer-P01 / Auditor-P01 |
| P02 | Belkin & Niyogi, *Laplacian Eigenmaps and Spectral Techniques for Embedding and Clustering* | Reviewer-P02 / Auditor-P02 |
| P03 | Coifman & Lafon, *Diffusion Maps* | Reviewer-P03 / Auditor-P03 |
| P04 | Berry & Harlim, *Variable Bandwidth Diffusion Kernels* | Reviewer-P04 / Auditor-P04 |
| P05 | Zhu, Ghahramani & Lafferty, *Gaussian Fields and Harmonic Functions* | Reviewer-P05 / Auditor-P05 |
| P06 | Belkin, Niyogi & Sindhwani, *Manifold Regularization* | Reviewer-P06 / Auditor-P06 |
| P07 | Shuman et al., *The Emerging Field of Signal Processing on Graphs* | Reviewer-P07 / Auditor-P07 |
| P08 | Hein, Audibert & von Luxburg, *Graph Laplacians and their Convergence on Random Neighborhood Graphs* | Reviewer-P08 / Auditor-P08 |
| P09 | Belkin & Niyogi, *Towards a Theoretical Foundation for Laplacian-Based Manifold Methods* | Reviewer-P09 / Auditor-P09 |
| P10 | Cheng & Wu, *Convergence of Graph Laplacian with kNN Self-tuned Kernels* | Reviewer-P10 / Auditor-P10 |

## Tier 2 Candidates

Tier 2 papers will be used to fill gaps after Tier 1 reviews:

- von Luxburg, Belkin & Bousquet, *Consistency of Spectral Clustering*.
- Singer, *From graph to manifold Laplacian: The convergence rate*.
- Smola & Kondor, *Kernels and Regularization on Graphs*.
- Hammond, Vandergheynst & Gribonval, *Wavelets on Graphs via Spectral Graph Theory*.
- Chung, *Spectral Graph Theory* selected sections.

## Required Memo Template Modifications

Compared with the H005 template, every Tier 1 memo must add:

- `Reader Background Needed`.
- `What a non-expert should understand before reading this paper.`
- `Copied Paper Figures Used`, with citation, figure number, page/section, and
  reason for inclusion.
- `Original Explanatory Figures Proposed Or Created`, for toy examples,
  conductance-curve plots, toy self-tuned kernels, three-node Laplacian
  examples, or pipeline diagrams.
- Explicit evidence labels for all substantive claims:
  `explicit`, `derived`, `contextual`, or `uncertain`.

## Output Paths

Base directory:

```text
/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/literature/conductance_kernel_laplacian_review
```

Planned outputs:

- LaTeX: `conductance_kernel_laplacian_review.tex`
- PDF: `conductance_kernel_laplacian_review.pdf`
- HTML: `conductance_kernel_laplacian_review.html`
- bibliography: `references.bib`
- source manifest: `source_manifest.yml`
- reviewer memos: `paper_memos/P##_*.md`
- audit memos: `audits/P##_audit.md`
- tables: `tables/*.tex`
- generated explanatory figures: `figures/P##_*`

## Build Method

PDF:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error conductance_kernel_laplacian_review.tex
```

HTML:

```bash
pandoc conductance_kernel_laplacian_review.tex \
  --from latex --to html5 --standalone --mathjax --citeproc \
  --bibliography=references.bib \
  -o conductance_kernel_laplacian_review.html
```

## Final Report Concept Map

The report must include a short concept map connecting:

- inverse length conductance;
- heat kernels;
- self-tuned/local-scale kernels;
- adaptive-radius/cKNN graph construction;
- diffusion/Markov normalization;
- graph low-pass smoothing.

## Phase Order

1. Acquire Tier 1 PDFs/HTML where available and record them in
   `source_manifest.yml`.
2. Create review and audit templates.
3. Launch ten dedicated Tier 1 reviewers with disjoint memo paths.
4. Launch corresponding auditors after each review memo exists.
5. Resolve audit findings and mark accepted/revision-needed memos.
6. Fill Tier 2 gaps with focused single-reviewer memos.
7. Synthesize conductance-family, task, and theory tables.
8. Write final LaTeX report using `references.bib`.
9. Build PDF and HTML.
10. Update resource dashboard only after the final report is accepted.

## Availability Notes

All ten Tier 1 PDFs were cached locally on 2026-05-15. The manifest records
original URL, DOI/arXiv ID when available, local path, access date, canonical
reading-copy status, page count, and SHA-256.

P05's older CMU `aladdin` host did not resolve during acquisition; the AAAI
ICML PDF was cached as the canonical reading copy.

## Risks And Approval Questions

- Tier 1 now contains ten full papers; this should produce a long review and
  may need multiple synthesis passes.
- Copied figures must remain internal review aids unless the final manuscript
  secures permission or replaces them with original explanatory figures.
- Auditors should reject memos that treat contextual SIMODS interpretations as
  author claims.
- The final report should cite stale or unavailable sources only as provenance,
  not as live evidence.
