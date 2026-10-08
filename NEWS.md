# spopt (development version)

* New `spider_lines()` turns location-allocation results from `p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, and `huff()` into an sf layer of demand-to-facility lines. It supports split `cflp()` allocations, `huff()` probabilities, and road-following lines through a user-supplied `route_fun`.
* `p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, and `huff()` now add a `.cost` column to the demand output: the cost to the assigned facility (or primary store), taken from the original cost matrix before any NA/Inf replacement.
* `p_median()`, `cflp()`, and `mclp()` now record `weight_col` in their result metadata.

# spopt 0.1.2

* Initial CRAN release.
