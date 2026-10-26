# L-infinity preparation validation

`prepare.sparse.mds(method = "linf")` was checked against the preserved
2026-10-25 max-normalization experiment for the comb-V3V4 rebuilt-reference
population. The data and original results remain external to this package.

## Comparison

Both implementations used maximum absolute coordinate differences on the same
max-normalized profiles, exact 10-nearest-neighbor symmetric-union local edges,
and regional endpoint weights. No PCA projection was used. The Python experiment
uses NumPy's random generator; the package uses its existing C++ generator.
The same seed therefore does not promise the same pivots across implementations.
For the historical comparison only, the native preparation received the saved
pivot indices. The public interface continued to select its own pivots.

| Check | Pilot | Full population |
| --- | ---: | ---: |
| Distinct profiles | 1,024 | 24,412 |
| Features | 1,159 | 1,159 |
| Preserved pivots | 200 | 400 |
| Local edges, exactly reproduced | 7,990 | 188,465 |
| Preserved constraints, exactly reproduced | 190,310 | 9,867,465 |

Every canonical pair, target distance, and both endpoint counts matched exactly.
The full population represents 24,465 source samples after consolidation of
identical profiles. The public pilot preparation's internally selected pivots
and constraints were also compared against the independent Python calculation;
all constraints matched exactly for those selected pivots.

The pilot was refit with SGD, seed 17, one supplied initial configuration,
200 epochs, the Zheng schedule with epsilon 0.1, and inverse-squared target
weights. Its maximum coordinate difference from the saved fit was
5.55e-16. RMS relative error on all 523,776 pilot pairs was 21.0960795%;
the sparse fitting proxy was 18418.5856817. The proxy is not the all-pair
relative-error percentage. Completing 200 epochs does not certify convergence.
No new full-population embedding was fitted, and no original result was changed.

The optimized public preparations took 0.727 seconds for the pilot and
192.113 seconds for the full population on the validation machine. These are
single-run timings, not comparative benchmarks. Each preparation was given an
8 GiB workspace allowance; this is not a process-memory cap. An earlier attempt
was stopped after discovering reused debug object files. Its output was retained
separately and its timing was not used. All native objects were then rebuilt
with optimization enabled before the completed validation.

## Reproduce

First build the local source with optimization enabled. Cleaning native objects
is necessary when documentation generation has previously compiled debug objects:

```r
pkgbuild::clean_dll("/path/to/grip")
pkgbuild::compile_dll("/path/to/grip", debug = FALSE, force = TRUE)
```

Run the comparison using a Python environment with NumPy and SciPy:

```sh
python tools/validation/linf-maxnorm.py /path/to/comb_v3v4_maxnorm /new/output/directory
```

The experiment directory must contain its preserved `scripts/geometry.py` and
`results/2026-10-25-v1` inputs and outputs. The validator checks input and saved
constraint hashes, records the current package source hashes and R session,
and refuses to overwrite an existing output directory. It limits numerical
libraries to one thread. Private matrices and generated results must remain
outside version control.

The standalone package tests in `test-linf-sparse-mds.R` cover dense distance
and regional-weight oracles, exact neighbor ties, randomized and farthest
selection, scale invariance, RNG preservation, SGD integration, invalid inputs,
and workspace refusal. They require none of the cohort data.

## Package checks

- Focused sparse-MDS checks: 1,679 assertions passed, with no failures, warnings,
  or skips.
- Complete source-tree suite with the local `dgraphs` dependency: 7,034
  assertions passed, with no failures, warnings, or skips.
- `make check` (`R CMD check --as-cran`): zero errors, zero warnings, two notes.
  The notes concern the development version/tarball size and the locally outdated
  HTML Tidy validator. Examples, tests, vignette rebuilds, and the PDF manual
  passed. Its isolated dependency library skipped 24 optional tests; the separate
  source-tree suite above exercised all tests without skips.
- Repository source-only checks and the public API/documentation audit passed.
