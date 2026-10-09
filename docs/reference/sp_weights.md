# Create spatial weights from an sf object

Constructs spatial weights (neighborhood structure) from sf geometries.
Wraps spdep functions with a convenient interface.

## Usage

``` r
sp_weights(
  data,
  type = c("queen", "rook", "knn", "distance"),
  k = NULL,
  d = NULL,
  ...
)
```

## Arguments

- data:

  An sf object with polygon or point geometries.

- type:

  Type of weights. One of:

  - "queen" (default): Polygons sharing any boundary point are neighbors

  - "rook": Polygons sharing an edge are neighbors

  - "knn": K-nearest neighbors based on centroid distance

  - "distance": All units within a distance threshold are neighbors

- k:

  Number of nearest neighbors. Required when `type = "knn"`.

- d:

  Distance threshold. Required when `type = "distance"`. Units match the
  CRS of the data (e.g., meters for projected CRS).

- ...:

  Additional arguments passed to spdep functions.

## Value

A neighbors list object (class "nb") compatible with spdep.

## Details

**Choosing a weight type:**

- Use **queen/rook** for polygon data where physical adjacency matters

- Use **knn** when you need guaranteed connectivity (no isolates) or for
  point data

- Use **distance** for point data or when interaction depends on
  proximity

**KNN weights** always produce a connected graph (if k \>= 1), making
them useful for datasets with islands or disconnected polygons.

## Examples

``` r
# \donttest{
library(sf)
nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#> Reading layer `nc' from data source 
#>   `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/sf/shape/nc.shp' 
#>   using driver `ESRI Shapefile'
#> Simple feature collection with 100 features and 14 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -84.32385 ymin: 33.88199 xmax: -75.45698 ymax: 36.58965
#> Geodetic CRS:  NAD27

# Queen contiguity (default)
w_queen <- sp_weights(nc, type = "queen")

# K-nearest neighbors (guarantees connectivity)
w_knn <- sp_weights(nc, type = "knn", k = 6)

# Distance-based (e.g., 50km for projected data)
nc_proj <- st_transform(nc, 32119)  # NC State Plane
w_dist <- sp_weights(nc_proj, type = "distance", d = 50000)
# }
```
