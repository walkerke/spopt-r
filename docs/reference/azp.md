# Automatic Zoning Procedure (AZP)

Performs regionalization using the Automatic Zoning Procedure algorithm.
AZP uses local search to minimize within-region heterogeneity while
maintaining spatial contiguity. Three variants are available: basic
(greedy), tabu search, and simulated annealing.

## Usage

``` r
azp(
  data,
  attrs = NULL,
  n_regions,
  weights = "queen",
  bridge_islands = FALSE,
  method = c("tabu", "basic", "sa"),
  max_iterations = 100L,
  tabu_length = 10L,
  cooling_rate = 0.99,
  initial_temperature = 0,
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

- method:

  Character. Optimization method: "basic" (greedy local search), "tabu"
  (tabu search), or "sa" (simulated annealing). Default is "tabu".

- max_iterations:

  Integer. Maximum number of iterations (default 100).

- tabu_length:

  Integer. Length of tabu list for tabu method (default 10).

- cooling_rate:

  Numeric. Cooling rate for SA method, between 0 and 1 (default 0.99).

- initial_temperature:

  Numeric. Initial temperature for SA method. If 0 (default),
  automatically set based on initial objective.

- scale:

  Logical. If TRUE (default), standardize attributes before clustering.

- seed:

  Optional integer for reproducibility.

- verbose:

  Logical. Print progress messages.

## Value

An sf object with a `.region` column containing cluster assignments.
Metadata is stored in the "spopt" attribute, including:

- algorithm: "azp"

- method: The optimization method used

- n_regions: Number of regions created

- objective: Total within-region sum of squared deviations

- solve_time: Time to solve in seconds

## Details

The Automatic Zoning Procedure (AZP) was introduced by Openshaw (1977)
and refined by Openshaw & Rao (1995). It is a local search algorithm
that:

1.  Starts with an initial random partition into n_regions

2.  Iteratively moves border areas between regions to reduce
    heterogeneity

3.  Maintains spatial contiguity throughout

4.  Terminates when no improving moves are found

Three variants are available:

- **basic**: Greedy local search that only accepts improving moves

- **tabu**: Tabu search that can accept non-improving moves to escape
  local optima, with a tabu list preventing cycling

- **sa**: Simulated annealing that accepts worse moves with decreasing
  probability as temperature cools

## References

Openshaw, S. (1977). A geographical solution to scale and aggregation
problems in region-building, partitioning and spatial modelling.
Transactions of the Institute of British Geographers, 2(4), 459-472.

Openshaw, S., & Rao, L. (1995). Algorithms for reengineering 1991 Census
geography. Environment and Planning A, 27(3), 425-446.

## Examples

``` r
# \donttest{
library(sf)
#> Linking to GEOS 3.13.0, GDAL 3.8.5, PROJ 9.5.1; sf_use_s2() is TRUE
nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#> Reading layer `nc' from data source 
#>   `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/sf/shape/nc.shp' 
#>   using driver `ESRI Shapefile'
#> Simple feature collection with 100 features and 14 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -84.32385 ymin: 33.88199 xmax: -75.45698 ymax: 36.58965
#> Geodetic CRS:  NAD27

# Basic AZP with 8 regions
result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 8)

# Tabu search (often finds better solutions)
result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 8,
              method = "tabu", tabu_length = 15)

# Simulated annealing
result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 8,
              method = "sa", cooling_rate = 0.95)

# View results
plot(result[".region"])

# }
```
