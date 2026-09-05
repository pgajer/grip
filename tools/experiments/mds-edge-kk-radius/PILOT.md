# Pilot gate: passed with a scale-policy amendment

The initial profiled-scale pilot exposed uniform contraction of configurations.
The minimum saved edge scale was 1.50e-7. The loss identity E(cZ)=c²E(Z)
explains why this cannot be treated as shape convergence. The complete initial
fit audit is in summary/initial-pilot-scale-audit.csv. The source, protocol,
and local checkpoints were retained rather than silently overwritten.

The amended 72-graph pilot completed 972 scored candidates and 216 MDS starts.
It added fixed-target-scale controls and computed package diagnostics after
edge calibration. Maximum independently recomputed package-score difference was
8.44e-15; maximum explicit path-sum discrepancy was 8.88e-16 in normalized
units. All 24 complete-graph objective checks passed. Twelve graphs required
labeled MST augmentation. Fixed-scale candidates had edge-calibration factors
between 0.99555 and 1.01757. These are successful numerical controls, not a
claim of global optimizer convergence.

Before expansion, 16 large-radius pilot graphs (k=32,239) received 256 additional
scored candidates covering perturbations, random/original starts, and additional
budgets under both scale policies. For stress-initialized branches, perturbation
changed s3/s2 by at most 0.00018 and continuation by at most 0.0047. Random starts
could change it by more than 0.4. Longer profiled runs contracted further, to an
edge scale around 5.7e-10; fixed-scale controls remained stable. These results
justify reporting achieved MDS-initialized fits and optimizer variation
separately. They do not justify claims of unique recovery or global optimality.

All 88 surface/sample/radius inputs were ready before expansion. Independent
validation used 256 actual r=64 pairs spanning all 16 high-radius cases,
including both sample sizes. Tighter shooting and independent collocation or
quadrature agreed within the predeclared 1e-5 relative tolerance. The complete
per-pair errors are in summary/geodesic-validation.csv; their maximum is a
sampled discrepancy, never an all-pair accuracy bound.

Decision: proceed with the full frozen radius/k/replicate grid and the nested
size check, preserving both scale policies and all negative findings. No
package default or manuscript result is changed by this decision.
