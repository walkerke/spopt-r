# Spider lines for location-allocation results

Builds "spider lines" (also called allocation lines or desire lines):
one line from each demand location to the facility that serves it. Works
with results from
[`p_median()`](https://walker-data.com/spopt-r/reference/p_median.md),
[`p_center()`](https://walker-data.com/spopt-r/reference/p_center.md),
[`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md),
[`lscp()`](https://walker-data.com/spopt-r/reference/lscp.md),
[`mclp()`](https://walker-data.com/spopt-r/reference/mclp.md), and
[`huff()`](https://walker-data.com/spopt-r/reference/huff.md). Lines are
straight by default; supply `route_fun` to draw lines that follow a road
network instead.

## Usage

``` r
spider_lines(
  result,
  allocations = c("primary", "all"),
  min_share = 0,
  cost_matrix = NULL,
  anchor = c("centroid", "point_on_surface"),
  route_fun = NULL,
  facilities = NULL
)
```

## Arguments

- result:

  A result from
  [`p_median()`](https://walker-data.com/spopt-r/reference/p_median.md),
  [`p_center()`](https://walker-data.com/spopt-r/reference/p_center.md),
  [`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md),
  [`lscp()`](https://walker-data.com/spopt-r/reference/lscp.md),
  [`mclp()`](https://walker-data.com/spopt-r/reference/mclp.md), or
  [`huff()`](https://walker-data.com/spopt-r/reference/huff.md), with
  demand and facility rows in the order the solver returned them.

- allocations:

  Which links to draw. `"primary"` (default) draws at most one line per
  demand row: the assigned facility (`.facility`) or the
  highest-probability store (`.primary_store`). `"all"` also draws split
  allocations from
  [`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md) and
  every store with non-zero probability from
  [`huff()`](https://walker-data.com/spopt-r/reference/huff.md); for the
  other solvers it is the same as `"primary"`.

- min_share:

  Drop links whose `share` is less than or equal to this value. Applies
  to primary and additional links alike; shares are never renormalized.
  Must be in `[0, 1)`. Most useful for
  [`huff()`](https://walker-data.com/spopt-r/reference/huff.md) results
  with `allocations = "all"`.

- cost_matrix:

  Optional demand x facility cost matrix, normally the one passed to the
  solver. Not needed for primary links, which take their cost from the
  solver's `.cost` column. If supplied, it takes precedence and provides
  `cost` for every link, including additional
  [`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md) and
  [`huff()`](https://walker-data.com/spopt-r/reference/huff.md) links.

- anchor:

  How polygon (or multipoint) geometries are reduced to line endpoints:
  `"centroid"` (default, matching
  [`distance_matrix()`](https://walker-data.com/spopt-r/reference/distance_matrix.md))
  or `"point_on_surface"` (guaranteed to fall inside the polygon).

- route_fun:

  Optional function that returns network-following geometries. It is
  called once as `route_fun(from, to)`, where `from` and `to` are sf
  POINT data frames in EPSG:4326 with the same number of rows and an
  `id` column; row *i* of `from` should be routed to row *i* of `to`.
  Identical coordinate pairs are sent only once. It must return an sf or
  sfc object with one LINESTRING or MULTILINESTRING per row, in input
  order, with a known CRS. Return an empty geometry for a pair that
  cannot be routed; those links fall back to straight lines.
  [`r5r_route_fun()`](https://walker-data.com/spopt-r/reference/r5r_route_fun.md)
  builds such a function from an r5r network.

- facilities:

  Optional integer row indices of `result$facilities` (or
  `result$stores`) to draw links for. Other facilities' links are
  skipped before any geometry is built, so drawing one store's trade
  area from a large
  [`huff()`](https://walker-data.com/spopt-r/reference/huff.md) result
  stays fast.

## Value

An sf object in the CRS of `result$demand`, one row per link, with
columns:

- `demand_index`: row in `result$demand`

- `facility_index`: row in `result$facilities` (or `result$stores`)

- `share`: 1 for whole assignments; the allocation fraction for
  [`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md); the
  choice probability for
  [`huff()`](https://walker-data.com/spopt-r/reference/huff.md)

- `weight`: demand weight multiplied by `share`, when the solver
  recorded a weight column

- `cost`: link cost from `cost_matrix`, or from the solver's `.cost`
  column for primary links; omitted when unavailable

- `length`: length of the drawn line

- `routed`: only when `route_fun` is used; `FALSE` where the router
  failed and a straight line was drawn instead

Geometry is LINESTRING, or MULTILINESTRING when `route_fun` is used.
Counts of unassigned demand, filtered links, and routing failures are
stored in the "spopt" attribute.

## Details

Indices and matrices in spopt results are positional, so do not reorder
or filter `result$demand` or `result$facilities` before calling
`spider_lines()`. Join attributes onto the output with `demand_index`
and `facility_index` instead.

Demand that the solver left unassigned (uncovered points in
[`lscp()`](https://walker-data.com/spopt-r/reference/lscp.md) or
[`mclp()`](https://walker-data.com/spopt-r/reference/mclp.md), or an
infeasible
[`p_center()`](https://walker-data.com/spopt-r/reference/p_center.md))
produces no lines.

For data in longitude/latitude, lines are drawn as two-vertex segments.
Use
[`sf::st_segmentize()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
on the output if you need great-circle arcs.

## See also

[`r5r_route_fun()`](https://walker-data.com/spopt-r/reference/r5r_route_fun.md)
for road-following lines with r5r;
[`p_median()`](https://walker-data.com/spopt-r/reference/p_median.md),
[`mclp()`](https://walker-data.com/spopt-r/reference/mclp.md),
[`cflp()`](https://walker-data.com/spopt-r/reference/cflp.md),
[`huff()`](https://walker-data.com/spopt-r/reference/huff.md)

## Examples

``` r
# \donttest{
library(sf)

demand <- st_as_sf(data.frame(
  x = runif(100), y = runif(100), population = rpois(100, 500)
), coords = c("x", "y"))
facilities <- st_as_sf(data.frame(x = runif(20), y = runif(20)),
                       coords = c("x", "y"))

result <- p_median(demand, facilities, n_facilities = 5,
                   weight_col = "population")
lines <- spider_lines(result)

plot(st_geometry(lines), col = lines$facility_index)
plot(st_geometry(result$facilities[result$facilities$.selected, ]),
     pch = 19, cex = 1.5, add = TRUE)

# }

if (FALSE) { # \dontrun{
# Lines that follow the road network, using a local OSRM server
osrm_lines <- function(from, to) {
  geoms <- lapply(seq_len(nrow(from)), function(i) {
    r <- tryCatch(osrm::osrmRoute(src = from[i, ], dst = to[i, ]),
                  error = function(e) NULL)
    if (is.null(r)) sf::st_linestring() else sf::st_geometry(r)[[1]]
  })
  sf::st_sfc(geoms, crs = 4326)
}

routed <- spider_lines(result, route_fun = osrm_lines)
} # }
```
