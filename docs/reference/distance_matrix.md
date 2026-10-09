# Compute distance matrix between sf objects

Computes pairwise distances between geometries. For geographic (longlat)
coordinate systems, uses great circle distances in meters via
[`sf::st_distance()`](https://r-spatial.github.io/sf/reference/geos_measures.html).
For projected coordinate systems, uses fast Euclidean distance in the
CRS units (typically meters).

## Usage

``` r
distance_matrix(
  x,
  y = NULL,
  type = c("euclidean", "manhattan"),
  use_centroids = NULL
)
```

## Arguments

- x:

  An sf object (demand points for facility location, or areas for
  regionalization).

- y:

  An sf object (facility locations). If NULL, computes distances within
  x.

- type:

  Distance type: "euclidean" (default) or "manhattan". Note that for
  geographic CRS, only "euclidean" (great circle) distance is available.

- use_centroids:

  Logical. If TRUE, use geometry centroids. Defaults to TRUE when any
  geometry is not a POINT (polygons, multipoints, lines), so each
  feature contributes exactly one row or column.

## Value

A numeric matrix of distances. Rows correspond to x, columns to y. For
geographic CRS, distances are in meters. For projected CRS, distances
are in the CRS units (usually meters).

## Examples

``` r
# \donttest{
library(sf)
demand <- st_as_sf(data.frame(x = runif(10), y = runif(10)), coords = c("x", "y"))
facilities <- st_as_sf(data.frame(x = runif(5), y = runif(5)), coords = c("x", "y"))
d <- distance_matrix(demand, facilities)
# }
```
