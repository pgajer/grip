# Weighted GRIP ND Parity Checkpoint, 2026-05-19

## Summary

The dimension-general weighted backend is now signature-compatible with
`grip.layout.weighted.legacy()` at the public R wrapper level and implements the
previously missing legacy-weighted public paths:

- coarse active-set repulsion controls;
- metric-neighbor cap;
- final anchor force, final move scaling, and alternate final KK-repulse mode;
- disconnected no-anchor random insertion behavior;
- post-layout LGKK polish;
- interleaved compiled multiscale LGKK stages.

The default `grip.layout.weighted()` wrapper routes to the ND backend. The
legacy backend remains available as `grip.layout.weighted.legacy()`.

## Checks Completed

Focused and gated checks completed after the LGKK port and formatting audit:

```sh
Rscript -e 'devtools::test(filter = "layout-weighted-nd$")'
```

Result: 91 pass.

```sh
GRIP_RUN_GFLOW_LGKK_PARITY_TESTS=true \
GRIP_ENFORCE_GFLOW_LGKK_PARITY=true \
Rscript -e 'devtools::test(filter = "gflow-parity")'
```

Result: 12 pass, 7 gated skips.

```sh
Rscript -e 'devtools::test(filter = "layout-weighted")'
```

Result: 153 pass, 10 gated skips.

```sh
GRIP_RUN_GFLOW_FINAL_ANCHOR_PARITY_TESTS=true \
GRIP_ENFORCE_GFLOW_FINAL_ANCHOR_PARITY=true \
Rscript -e 'devtools::test(filter = "gflow-parity")'
```

Result: 6 pass.

```sh
GRIP_RUN_TRACE_PARITY_TESTS=true \
Rscript -e 'devtools::test(filter = "layout-weighted-nd-trace-parity")'
```

Result: 13 pass, 1 gated skip.

```sh
GRIP_RUN_FINAL_ANCHOR_TRACE_PARITY_TESTS=true \
Rscript -e 'devtools::test(filter = "layout-weighted-nd-trace-parity")'
```

Result: 6 pass, 1 gated skip.

```sh
GRIP_RUN_GFLOW_FULL_PARITY_TESTS=true \
GRIP_ENFORCE_GFLOW_FULL_PARITY=true \
Rscript -e 'devtools::test(filter = "gflow-parity")'
```

Result: 13 pass, 7 gated skips.

```sh
GRIP_RUN_GFLOW_STRESS_PARITY_TESTS=true \
GRIP_ENFORCE_GFLOW_STRESS_PARITY=true \
Rscript -e 'devtools::test(filter = "gflow-parity")'
```

Result: 13 pass, 7 gated skips.

```sh
GRIP_RUN_GFLOW_TRACE_ALL_TESTS=true \
GRIP_ENFORCE_GFLOW_TRACE_ALL_PARITY=true \
Rscript -e 'devtools::test(filter = "gflow-parity")'
```

Result: 10 pass, 7 gated skips.

```sh
GRIP_RUN_GFLOW_REFINEMENT_STEP_TRACE_TESTS=true \
Rscript -e 'devtools::test(filter = "gflow-parity")'
```

Result: skipped because no divergent gflow trace frame was found at tolerance
`1e-08`. This means there was no current focal divergence for the
refinement-step diagnostic to expand.

```sh
make check-fast
```

Result: completed with the known 2 warnings and 1 note:

- non-portable compiler warning flags;
- undeclared `geometry` usage in tests;
- CRAN new-submission note.

```sh
make check
```

Result: completed with 2 warnings and 3 notes:

- non-portable compiler warning flags;
- undeclared `geometry` usage in tests;
- CRAN new-submission note;
- local environment note: unable to verify current time;
- HTML Tidy note during manual validation.

## Gflow Fixture Status

The previous external gflow parse blocker in
`R/transported_graph_hessian.R` was cleared in the gflow repository. After that
fix, the stress parity and cross-dimensional trace-all gates ran successfully.

The refinement-step gflow diagnostic still skipped, but for a different reason:
the trace comparison found no divergent frame at tolerance `1e-08`, leaving no
focal divergent frame for the refinement-step expander to inspect. This is not a
parity failure.

## Signature Audit

`formals(grip.layout.weighted.legacy())` and
`formals(grip.layout.weighted())` were compared after the LGKK port:

- no legacy-only arguments;
- no ND/default-wrapper-only arguments;
- no default-value differences.

## Remaining Notes

The untracked directories below are unrelated to this weighted-GRIP parity work
and were intentionally left untouched:

```text
dev/design/misf-edge-kk/
dev/misf-edge-kk/
```
