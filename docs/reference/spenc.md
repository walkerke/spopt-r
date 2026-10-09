# Spatially-Encouraged Spectral Clustering (SPENC)

Performs spectral clustering with spatial constraints by combining
spatial connectivity with attribute similarity using kernel methods.
This approach is useful for clustering with highly non-convex clusters
or irregular topologies in geographic contexts.

## Usage

``` r
spenc(
  data,
  attrs = NULL,
  n_regions,
  weights = "queen",
  bridge_islands = FALSE,
  gamma = 1,
  scale = TRUE,
  seed = NULL,
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

- gamma:

  Numeric. RBF kernel parameter controlling attribute similarity decay.
  Larger values = faster decay = more local similarity. Default is 1.
  Can also be "auto" to estimate from data.

- scale:

  Logical. If TRUE (default), standardize attributes before clustering.

- seed:

  Optional integer for reproducibility.

- verbose:

  Logical. Print progress messages.

## Value

An sf object with a `.region` column containing cluster assignments.
Metadata is stored in the "spopt" attribute, including:

- algorithm: "spenc"

- n_regions: Number of regions created

- objective: Within-cluster sum of squared distances in embedding space

- gamma: The gamma parameter used

- solve_time: Time to solve in seconds

## Details

SPENC (Wolf, 2021) extends spectral clustering to incorporate spatial
constraints. The algorithm:

1.  Computes attribute affinity using an RBF (Gaussian) kernel

2.  Multiplies element-wise with spatial weights (only neighbors have
    affinity)

3.  Computes the normalized Laplacian of the combined affinity matrix

4.  Extracts the k smallest eigenvectors as a spectral embedding

5.  Applies k-means clustering to the embedding

Key advantages:

- Can find non-convex cluster shapes

- Respects spatial connectivity

- Balances attribute similarity with spatial proximity

The gamma parameter controls how quickly attribute similarity decays
with distance in attribute space. Larger values create more localized
clusters.

## References

Wolf, L. J. (2021). Spatially-encouraged spectral clustering: a
technique for blending map typologies and regionalization. International
Journal of Geographical Information Science, 35(11), 2356-2373.
[doi:10.1080/13658816.2021.1934475](https://doi.org/10.1080/13658816.2021.1934475)

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

# Basic SPENC with 8 regions
result <- spenc(nc, attrs = c("SID74", "SID79"), n_regions = 8)

# Adjust gamma for different cluster tightness
result <- spenc(nc, attrs = c("SID74", "SID79"), n_regions = 8, gamma = 0.5)

# View results
plot(result[".region"])

# }
```
