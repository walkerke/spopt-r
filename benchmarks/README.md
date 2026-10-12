# Benchmarks

Reproducible performance and equivalence evidence for the regionalization
functions. Not part of the built package (`^benchmarks$` is in `.Rbuildignore`).

- `data/meck.gpkg`, `data/vancouver.gpkg`: the two real tract datasets used by
  the equivalence harness (Mecklenburg County NC, ACS 2018-2022, 293 tracts with
  median income, % bachelor's, median rent and median age; Vancouver, Burnaby and
  New Westminster BC, 2021 Census, 185 tracts with ten occupation shares).
  `fetch-tracts.R` regenerates them with tidycensus and cancensus (API keys
  needed); the committed copies let the harness run without either.
- `fixtures.R`: `make_grid_fixture(n, p, seed, islands)` builds a hexagonal grid
  of `n` polygons with `p` spatially smooth standardized attributes and a
  population column; `islands` deletes vertical bands to disconnect the grid.
- `stage-timing.R [libpath] [label] [sizes]`: times each public call with the
  neighbour list supplied (so `poly2nb` is reported once, separately), writes
  `stage-timing-<label>.csv`.
- `equivalence.R capture|compare <libpath> <dir>`: saves or checks complete
  outputs (partition up to relabelling, objective, region count, nb objects)
  for 29 cases over `nc`, Mecklenburg and Vancouver tracts and synthetic grids,
  including an integer-valued fixture with exact ties and a three-island
  fixture. `baseline/` holds outputs captured from the pre-optimization build
  (branch `fix/azp-maxp-search`), except the two island cases, re-captured
  after the `bridge_islands` connectivity fix.

The harness covers all regionalization functions. Cases for the rewritten AZP /
max-p engine and for `redcap()` are skipped automatically when the installed
build does not have them, so it runs on any subset of the regionalization pull
requests; their baseline files are kept here so the cases run once those land.
Conversely, on a build that lacks this pull request's SKATER, Ward and bridging
changes, three cases differ by design: `gint_ward` (the previous Ward broke
ties in hash-map order, so it does not reproduce itself) and the two `gisl_*`
cases (the previous `bridge_islands` left that fixture disconnected).

Install a build into its own library and point the scripts at it:

```sh
R CMD INSTALL --library=/tmp/rlib-x --no-docs .
Rscript benchmarks/equivalence.R compare /tmp/rlib-x benchmarks/baseline
Rscript benchmarks/stage-timing.R /tmp/rlib-x x "500,2500,10000"
```

Results: `before-*.csv` / `after-*.csv` and the table in `results.md`; summarized in NEWS.md. Machine:
Apple Silicon laptop, R 4.5, one process; absolute times vary across
machines, scaling does not.
