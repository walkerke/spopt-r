# spopt (development version)

* `spider_lines()` gains a `facilities` argument to draw links for selected facilities only. Filtering happens before any geometry is built, so one store's trade area from a large `huff()` result takes milliseconds instead of seconds.

# spopt 0.1.3

* New `spider_lines()` turns location-allocation results from `p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, and `huff()` into an sf layer of demand-to-facility lines. It supports split `cflp()` allocations, `huff()` probabilities, and road-following lines through a user-supplied `route_fun`. `r5r_route_fun()` builds a `route_fun` from an r5r network.
* `p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, and `huff()` now add a `.cost` column to the demand output: the cost to the assigned facility (or primary store), taken from the original cost matrix before any NA/Inf replacement.
* `p_median()`, `cflp()`, and `mclp()` now record `weight_col` in their result metadata.
* `distance_matrix()` now reduces any non-POINT geometry (including MULTIPOINT) to its centroid, so each feature gets one row or column. Previously, MULTIPOINT inputs produced a mis-sized matrix in projected CRSs.

# spopt 0.1.2

* Initial CRAN release.
