# Two-fidelity saddle pilot: results

Generated 2026-09-02 20:23:33 EDT.

Five independent surface-area-uniform clouds, each n=1000. All fitted layouts are 3D. No manuscript files were modified.

## Graph calibration

Loss values below are percentages. The numerical surface reference uses 128 random sources and 119,744 distinct pairs per cloud. X->G uses physical units, without scale profiling.

| Cloud | Selected k | X->G error | MDS path error | MDS + edge-KK path error | MDS edge error | Refined edge error | MDS Stress-1 | Refined Stress-1 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 71 | 0.2167 | 3.5063 | 0.2572 | 7.2064 | 0.3995 | 3.4590 | 4.6401 |
| 2 | 70 | 0.2047 | 3.3907 | 0.2376 | 7.1332 | 0.3829 | 3.3402 | 4.4828 |
| 3 | 67 | 0.2125 | 3.6564 | 0.2920 | 7.3990 | 0.3817 | 3.5757 | 4.9131 |
| 4 | 72 | 0.2032 | 3.4507 | 0.2155 | 7.0877 | 0.3681 | 3.4603 | 4.7152 |
| 5 | 73 | 0.2111 | 3.4495 | 0.2447 | 7.1005 | 0.3717 | 3.4387 | 4.9086 |

## Accounting and numerical checks

- 390 graph cases; five MDS fits, five primary edge-KK fits, and five separate additional-budget fits. No failed or excluded clouds.
- Primary edge-KK: 25/25 stages reached their iteration cap; additional-budget runs: 5/5. Earlier termination is not assumed to mean certified convergence.
- Original-coordinate edge error <= 5.04e-17; path error <= 0 (dimensionless).
- Maximum disagreement with package scoring: 0. Independent R route-sum check: 4.44e-16.
- Smooth-surface BVP spot-check error: 0.002935% to 0.003141% relative RMSE across clouds; each cloud has 128 controls.
- Mesh 81 versus mesh 161 on 16 sources per cloud: 0.002167% to 0.002417% relative difference.

- The extra 128-source fine-mesh check on cloud 5 retains k=73 (primary k=73); the full-source reference difference is 0.002266%.

## Pilot findings

The initial k=3:20 and expanded k=3:40 sweeps had boundary minima. The final k=3:80 sweep gives interior reference minima at 71, 70, 67, 72, 73. Their surface-to-graph errors span 0.203% to 0.217%. Nearby k values can be nearly indistinguishable, and smaller reference-source subsets sometimes select a different integer.

Primary edge-KK reduces fixed-path error in 5/5 clouds and edge error in 5/5; it increases MDS Stress-1 in 5/5. This is a change in which geometric quantity is preserved, not uniform improvement under all criteria.

Procrustes shape discrepancy changes from a median 29.2% under MDS to 30.7% after edge-KK. Small graph-path error therefore should not be presented as recovery of the original saddle coordinates. These finite-budget results do not show that exact recovery is impossible: the original coordinates provide an exact zero-loss realization.

At k=3, repaired-graph surface errors span 62.2% to 94.7% even though the original coordinates preserve all graph paths exactly. Figure 5 makes that separation explicit.

## Interpretation limits

The per-cloud k is an oracle choice using a numerical reference for the known surface. It is not an automatic rule available for unknown real-data geometry. Calibration uses a random source subset, while graph-to-embedding scoring uses all unordered pairs. Exact mesh distances are approximations to smooth-surface distances; mesh comparisons, smooth BVP checks, and source-subset sensitivity are recorded separately. Positive fitted errors are achieved values, not proven minima. Five clouds establish pilot behavior, not population precision.

The original 3D observations realize all Euclidean graph edges and every retained path exactly in mathematics, including MST bridges. This does not imply zero chord stress, exact surface-geodesic approximation, or shape recovery by the fitted methods. Endpoint chords and surface geodesics are distinct quantities.

## Files

- graph-calibration.pdf/png: complete error-versus-k curves and enlarged minima.
- representative-layouts.pdf/png: median-selected three-panel 3D comparison.
- paired-layout-scores.pdf/png: all five clouds and additional-budget diagnostics.
- two-stage-errors.pdf/png: surface-to-graph and end-to-end embedded-path errors.
- why-both-fidelities-matter.pdf/png: the same original coordinates preserve a poor graph and a calibrated graph exactly.
- calibration-curves.csv, reference-validation.csv, layout-scores.csv, additional-budget-scores.csv, layout-timings.csv, fit-status.csv: complete numeric summaries.
- Source README.md: protocol, equations, reproduction commands, dependency versions, and interpretation.
