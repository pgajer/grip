# README preview fixture

`saddle-preview.rds` contains only the coordinates, aligned configurations,
and triangle connectivity used by the public README saddle animations.
It is a frozen display fixture; it contains no manuscript text, review notes,
study summaries, or fitting caches.

Run `make readme-saddle`, `make readme-saddle-animation`, or
`make saddle-overlay-animation` from the package root to rebuild the previews.
These commands require no manuscript repository and do not rerun fits.

The [manifest](saddle-preview.json) records the fixture checksum and known display
structure. Fitting settings and solver provenance are unavailable in this
fixture; its historical method labels do not identify the current metric-MDS
implementation. The website recipe distinguishes replay from fitting.
