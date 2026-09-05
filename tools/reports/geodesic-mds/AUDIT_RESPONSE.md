# Response to the September 5, 2026 audit

The revision preserves the central theorems, original distance matrices,
selected embeddings, and original numerical tables. It adds targeted numerical
checks and narrows claims where the evidence is numerical. The report now has
39 pages, 15 figures, and eight tables.

## Main comments

1. **Raw-stress convergence wording.** The saddle experiment README, report
   opening, Section 5.1, and sampling-figure caption now describe four-arm
   behavior as finite-radius numerical evidence. Sections 5.6–5.7 distinguish
   convergence of the rescaled objective and the proved classical-MDS result
   from the unresolved universal claim about raw-stress minimizers.

2. **Scope of geodesic validation.** Appendix C.1 and the saddle README identify
   the original maximum as a statistic of 60 validation comparisons, not an
   all-pairs bound. A reproducible additional check covers 48 pairs from the
   actual radius-64 matrices. Refinement of the three most discrepant pairs
   gives a maximum relative disagreement of approximately 8.483e-8. The text
   records pair selection, solver retries, and refinement; distinguishes
   independent trajectory energy drift from internal shooting residuals; and
   uses “numerical tolerance” where floating-point rounding alone is not
   established as the explanation.

3. **Population mode selection.** Section 4.4 now qualifies the leading-mode
   identification immediately as numerical. It displays the radial logarithmic
   kernel and shows why projecting off the constant and height directions
   annihilates it. The competing Fourier sectors and their signs are explicit.
   New Table 5 reports the second logarithmic-mean eigenvalue and competing
   sectors; checks at 256 and 512 quadrature nodes agree on their ordering.
   Section 4.3 also derives the unequal-radius correction directly from the
   angular and length integrands. The companion antipodal note is aligned.

4. **Gorodski citation locator.** Both citation-verification HTML files and the
   saddle README now locate the minimizing statement in the proof of Lemma
   6.5.4, printed page 125 (one-based PDF page 131). The underlying claim and
   cited source are unchanged.

## Broader observations

- **Intrinsic dimension versus Euclidean span.** The opening foregrounds the
  one-dimensional saddle tree whose classical embedding spans three Euclidean
  dimensions. The direct saddle/paraboloid comparison is now Figure 1 in the
  main text; all 14 other figures remain in the atlas.
- **Optimizer spread.** New Section 5.7 and Table 7 report all six original
  radius-64 starts. Worst/best stress ratios are approximately 1.965 and 2.365
  for the two paraboloid measures, whereas saddle relative spreads are below
  4.4e-12. Sixteen additional paraboloid fits improve the selected stress only
  at relative order 1e-9. The saved selected fits are retained.
- **Planar instability.** Section 5.7 derives the transverse Hessian and reports
  negative eigenvalues for both saved planar saddle fits, corroborated by
  finite differences and actual stress-decreasing perturbations. This excludes
  those configurations as three-dimensional local minima; it does not exclude
  every possible planar global minimizer.
- **Tree rank condition.** Section 5.5 gives a readily checked sufficient
  condition and proof: represent every branch by a nonroot point and include
  either two distinct positive lengths on one branch or the root.
- **Antipodal bound constants.** Section 4.2 explains why the finite-sample
  bounded-error statement can have large constants and need not describe
  moderate radii well. It records the small radial gap and growing candidate
  stresses relevant to this sample.
- **Quadratic regression.** The existing shape coefficients remain descriptive
  summaries of numerical fits. The text reinforces this interpretation using
  the observed optimizer spread. An anisotropic quadratic-regression study
  remains a possible extension; no revised conclusion depends on it.

## Reproducibility and validation

The new formal experiment lives at
`tools/experiments/mds-audit-diagnostics/`. It reads the existing experiment
outputs and requires no private audit files. Its saved outputs include
optimizer spreads at every radius, the 16 extra fits, transverse Hessian checks,
population-sector estimates, all 48 pair comparisons and solver attempts, nine
refinements, and checksums.

The report build verifies the original experiment bundles and the new diagnostic
bundle before generating tables. It also checks citations, cross-references,
all 15 figure captions, overfull boxes, and source/output checksums. The revised
PDF was rendered and visually reviewed. The source ZIP includes this response,
the new CSV evidence, and everything needed to compile the report independently.
The full original matrix and optimization sweeps were not rerun for this prose
and targeted-diagnostics revision.
