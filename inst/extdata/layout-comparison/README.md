# Fixed GRIP layout comparison

`comparison.rds` preserves the `494_bus` graph and the complete comparison
used by `vignettes/grip-real-data.Rmd`: three GRIP presets, each with seeds
1 through 3, in three dimensions. All nine layouts, individual scores,
summary records, settings and warnings are retained. This is a teaching
comparison, not a tuned recommendation for power networks. Elapsed-time
fields describe the saved runs, not the current document build.

After installing the current source checkout, run from its root:

```sh
Rscript tools/pkg/build-layout-comparison.R
```

The generator evaluates the graph and comparison code in the guide and
records source checksums, the installed native library checksum and session
information. The guide validates the graph input and rebuilds every table
and figure from the saved comparison. The fitting recipe remains visible
and runnable for optional reproduction.
