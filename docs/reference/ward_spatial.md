# Ward Spatial Clustering

Performs spatially-constrained hierarchical clustering using Ward's
minimum variance method. Only spatially contiguous areas can be merged,
ensuring all resulting regions are spatially connected.

## Usage

``` r
ward_spatial(
  data,
  attrs = NULL,
  n_regions,
  weights = "queen",
  bridge_islands = FALSE,
  scale = TRUE,
  verbose = FALSE
)
```

## Arguments

- data:

  An sf object with polygon or point geometries.

- attrs:

  Character vector of column names to use for clustering (e.g.,
  `c("var1", "var2")`). If NULL, uses all numeric columns.

- n_regions:

  Integer. Number of regions (clusters) to create.

- weights:

  Spatial weights specification. Can be:

  - `"queen"` (default): Polygons sharing any boundary point are
    neighbors

  - `"rook"`: Polygons sharing an edge are neighbors

  - An `nb` object from spdep or created with
    [`sp_weights()`](https://walker-data.com/spopt-r/reference/sp_weights.md)

  - A list for other weight types: `list(type = "knn", k = 6)` for
    k-nearest neighbors, or `list(type = "distance", d = 5000)` for
    distance-based weights

- bridge_islands:

  Logical. If TRUE, automatically connect disconnected components (e.g.,
  islands) using nearest-neighbor edges. If FALSE (default), the
  function will error when the spatial weights graph is disconnected.

- scale:

  Logical. If TRUE (default), standardize attributes before clustering.

- verbose:

  Logical. Print progress messages.

## Value

An sf object with a `.region` column containing cluster assignments.
Metadata is stored in the "spopt" attribute.

## Details

This function implements spatially-constrained agglomerative
hierarchical clustering using Ward's minimum variance criterion. Unlike
standard Ward clustering, this version enforces spatial contiguity by
only allowing clusters that share a border to be merged.

The algorithm:

1.  Starts with each observation as its own cluster

2.  At each step, finds the pair of **adjacent** clusters with minimum
    Ward distance (increase in total within-cluster variance)

3.  Merges them into a single cluster

4.  Repeats until the desired number of regions is reached

The result guarantees that all regions are spatially contiguous.

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

# Cluster into 8 spatially-contiguous regions
result <- ward_spatial(nc, attrs = c("SID74", "SID79"), n_regions = 8)
plot(result[".region"])

# }
```
