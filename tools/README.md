# Package tools

`tools/pkg/` contains package documentation and release utilities, including
self-contained inputs for the public README previews. Repository checks live
at this directory's top level.

Manuscript experiments, figure builders, reports, and submission utilities
belong inside their owning private manuscript workspace outside this repository.
The public package must build and pass checks without those workspaces.

`make repo-hygiene` checks public source placement, development-file hygiene,
Rd documentation for exported functions, and private-data/reference exclusion.
`make test-release-content` exercises the release-content check in isolated
fixtures, including absent documentation and private-path failures. Historical
manuscript catalogs are not inputs to either target.
