# Build a Corridor Graph for Cached Routing

Pre-build the routing graph from a cost surface so that multiple
[`route_corridor`](https://walker-data.com/spopt-r/reference/route_corridor.md)
calls can skip graph construction. The returned object is a snapshot of
the cost surface at build time; subsequent edits to the raster do not
affect the graph.

## Usage

``` r
corridor_graph(cost_surface, neighbours = 8L, resolution_factor = 1L)
```

## Arguments

- cost_surface:

  A terra SpatRaster (single band). Same requirements as
  [`route_corridor`](https://walker-data.com/spopt-r/reference/route_corridor.md).

- neighbours:

  Integer. Cell connectivity: 4, 8 (default), or 16. Fixed at build
  time.

- resolution_factor:

  Integer, default 1L. If \> 1, the surface is aggregated before graph
  construction. Fixed at build time.

## Value

An opaque `spopt_corridor_graph` object for use with
[`route_corridor`](https://walker-data.com/spopt-r/reference/route_corridor.md).

## Details

The object retains both the CSR graph (in Rust) and an independent copy
of the raster (for coordinate mapping). The printed `graph_storage`
reflects only the CSR arrays, not the raster copy.

## Examples

``` r
# \donttest{
library(terra)
#> terra 1.9.46
r <- rast(nrows = 500, ncols = 500, xmin = 0, xmax = 500000,
          ymin = 0, ymax = 500000, crs = "EPSG:32614")
values(r) <- runif(ncell(r), 0.5, 2.0)

g <- corridor_graph(r, neighbours = 8L)
print(g)
#> Corridor graph
#>   Grid: 500 x 500 (250,000 cells)
#>   Cell size: 1000.0 x 1000.0
#>   Neighbours: 8 (1,994,004 edges)
#>   Build time: 0.005s | Graph storage: ~33.9 MB
path <- route_corridor(g, c(50000, 50000), c(450000, 450000))
# }
```
