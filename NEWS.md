# spopt (development version)

* `spider_lines()` gains a `facilities` argument to draw links for selected facilities only. Filtering happens before any geometry is built, so one store's trade area from a large `huff()` result takes milliseconds instead of seconds.



## Bug fixes

* `bridge_islands = TRUE` could return a graph that was still disconnected:
  each component was linked once to its closest unit in any other component,
  so two islands closer to each other than to the mainland linked to each
  other only. Bridging now repeats until one component remains. `skater()` on
  such a forest previously returned one region more than requested.

## Performance (regionalization)

Measured on synthetic hexagonal grids with four smooth attributes
(`benchmarks/stage-timing.R`); partitions are unchanged on every captured
fixture (`benchmarks/equivalence.R`).

* `skater()`: cuts are scored from per-subtree sums in one pass per region and
  each region's best cut is cached, O(k n p) instead of O(k n^2 p).
  2,500 units: 2.2 s -> 0.09 s; 10,000 units: 31 s -> 0.8 s.
* `ward_spatial()`: candidate merges are kept in a priority queue with lazy
  invalidation instead of rescanning every adjacent pair per merge. Results are
  now deterministic (tie-breaking previously followed hash-map order).
  2,500 units: 0.39 s -> 0.09 s; 10,000 units: 7.3 s -> 0.8 s.
* Island bridging finds each component's closest cross-component pair with one
  distance matrix instead of a nested R loop (2,235 units, three islands:
  2.3 s -> 0.1 s excluding `poly2nb`).

# spopt 0.1.3

* New `spider_lines()` turns location-allocation results from `p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, and `huff()` into an sf layer of demand-to-facility lines. It supports split `cflp()` allocations, `huff()` probabilities, and road-following lines through a user-supplied `route_fun`. `r5r_route_fun()` builds a `route_fun` from an r5r network.
* `p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, and `huff()` now add a `.cost` column to the demand output: the cost to the assigned facility (or primary store), taken from the original cost matrix before any NA/Inf replacement.
* `p_median()`, `cflp()`, and `mclp()` now record `weight_col` in their result metadata.
* `distance_matrix()` now reduces any non-POINT geometry (including MULTIPOINT) to its centroid, so each feature gets one row or column. Previously, MULTIPOINT inputs produced a mis-sized matrix in projected CRSs.

# spopt 0.1.2

* Initial CRAN release.
