# Changelog

## spopt 0.1.3

- New
  [`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)
  turns location-allocation results from
  [`p_median()`](https://walker-data.com/spopt-r/reference/p_median.md),
  [`p_center()`](https://walker-data.com/spopt-r/reference/p_center.md),
  [`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md),
  [`lscp()`](https://walker-data.com/spopt-r/reference/lscp.md),
  [`mclp()`](https://walker-data.com/spopt-r/reference/mclp.md), and
  [`huff()`](https://walker-data.com/spopt-r/reference/huff.md) into an
  sf layer of demand-to-facility lines. It supports split
  [`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md)
  allocations,
  [`huff()`](https://walker-data.com/spopt-r/reference/huff.md)
  probabilities, and road-following lines through a user-supplied
  `route_fun`.
  [`r5r_route_fun()`](https://walker-data.com/spopt-r/reference/r5r_route_fun.md)
  builds a `route_fun` from an r5r network.
- [`p_median()`](https://walker-data.com/spopt-r/reference/p_median.md),
  [`p_center()`](https://walker-data.com/spopt-r/reference/p_center.md),
  [`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md),
  [`lscp()`](https://walker-data.com/spopt-r/reference/lscp.md),
  [`mclp()`](https://walker-data.com/spopt-r/reference/mclp.md), and
  [`huff()`](https://walker-data.com/spopt-r/reference/huff.md) now add
  a `.cost` column to the demand output: the cost to the assigned
  facility (or primary store), taken from the original cost matrix
  before any NA/Inf replacement.
- [`p_median()`](https://walker-data.com/spopt-r/reference/p_median.md),
  [`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md), and
  [`mclp()`](https://walker-data.com/spopt-r/reference/mclp.md) now
  record `weight_col` in their result metadata.
- [`distance_matrix()`](https://walker-data.com/spopt-r/reference/distance_matrix.md)
  now reduces any non-POINT geometry (including MULTIPOINT) to its
  centroid, so each feature gets one row or column. Previously,
  MULTIPOINT inputs produced a mis-sized matrix in projected CRSs.

## spopt 0.1.2

CRAN release: 2026-04-22

- Initial CRAN release.
