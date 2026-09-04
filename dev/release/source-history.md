# Repository source history

## September 3, 2026 cleanup

The repository history was rewritten to remove internal working records that
were not intended for distribution. The cleanup also removed retired, unused
OpenGL trackball source and revised the introduction of a historical API
catalog without changing its technical entries.

This changed affected commit identifiers and the Git object behind the
`v0.2.0` tag. The cleaned tag resolves to
[`520902ad3f1b2aeabd287a379f5c08729c7b2c5d`](https://github.com/pgajer/grip/tree/520902ad3f1b2aeabd287a379f5c08729c7b2c5d).
All files shared by that release snapshot and the published CRAN source
archive were verified unchanged. This was a repository maintenance operation,
not a new package release or a change to the published CRAN archive.

The published `grip_0.2.0.tar.gz` archive has SHA-256
`5fbaecee890cd1f402e6776c661f5934325d37ac93a7d39120b54cf7acaf31a8`.
The software paper's reproducibility instructions use the cleaned Git
identifier and record this archive checksum separately.

For clones made before the cleanup, preserve any uncommitted work and use a
fresh clone. Port necessary unpublished changes onto the cleaned history;
do not merge or push branches based on the old history, because doing so can
restore the removed records.
