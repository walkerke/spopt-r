# Huff Model for Market Share Analysis

Computes probability surfaces to predict market share and sales
potential based on distance decay and store attractiveness. The Huff
model is widely used in retail site selection to estimate the
probability that a consumer at a given location will choose a particular
store.

## Usage

``` r
huff(
  demand,
  stores,
  attractiveness_col,
  attractiveness_exponent = 1,
  distance_exponent = -1.5,
  sales_potential_col = NULL,
  cost_matrix = NULL,
  distance_metric = "euclidean"
)
```

## Arguments

- demand:

  An sf object representing demand points or areas. Can be customer
  locations, census block groups, grid cells, etc.

- stores:

  An sf object representing store/facility locations.

- attractiveness_col:

  Character vector. Column name(s) in `stores` containing attractiveness
  values (e.g., square footage, parking spaces). Multiple columns can be
  specified for composite attractiveness.

- attractiveness_exponent:

  Numeric vector. Exponent(s) for attractiveness (default 1). Must be
  same length as `attractiveness_col` or length 1 (recycled). Higher
  values increase the importance of that variable.

- distance_exponent:

  Numeric. Distance decay exponent (default -1.5).

  Should be negative; more negative = faster decay with distance.

- sales_potential_col:

  Optional character. Column name in `demand` containing sales potential
  values (e.g., disposable income, population). If NULL, each demand
  point is weighted equally.

- cost_matrix:

  Optional. Pre-computed distance/cost matrix (demand x stores). If
  NULL, Euclidean distance is computed from geometries.

- distance_metric:

  Distance metric if cost_matrix is NULL: "euclidean" (default) or
  "manhattan".

## Value

A list with:

- `$demand`: Original demand sf with added columns:

  - `.primary_store`: Row index (1-based) of the highest-probability
    store

  - `.cost`: Cost to the primary store, from the original cost matrix

  - `.entropy`: Competition measure (higher = more competition)

  - `.prob_<store_id>`: Probability columns for each store

- `$stores`: Original stores sf with added columns:

  - `.market_share`: Proportion of total market captured

  - `.expected_sales`: Expected sales (sum of prob × potential)

- `$probability_matrix`: Full probability matrix (n_demand × n_stores)

Metadata in "spopt" attribute includes parameters used.

## Details

The Huff model calculates the probability that a consumer at location i
will choose store j using:

\$\$P\_{ij} = \frac{A_j \times D\_{ij}^\beta}{\sum_k A_k \times
D\_{ik}^\beta}\$\$

Where:

- \\(A_j\\) is the composite attractiveness of store j

- \\(D\_{ij}\\) is the distance from i to j

- \\(\beta\\) is the distance decay exponent (default -1.5)

When multiple attractiveness variables are specified, the composite
attractiveness is computed as:

\$\$A_j = \prod_m V\_{jm}^{\alpha_m}\$\$

Where \\(V\_{jm}\\) is the value of attractiveness variable m for store
j, and \\(\alpha_m\\) is the corresponding exponent.

The distance exponent is typically negative because probability
decreases with distance. Common values range from -1 to -3.

## Outputs

**Market Share**: The weighted average probability across all demand
points, representing the proportion of total market potential captured
by each store.

**Expected Sales**: The sum of (probability × sales_potential) for each
store, representing the expected sales volume.

**Entropy**: A measure of local competition. Higher entropy indicates
more competitive areas where multiple stores have similar probabilities.

## References

Huff, D. L. (1963). A Probabilistic Analysis of Shopping Center Trade
Areas. Land Economics, 39(1), 81-90.
[doi:10.2307/3144521](https://doi.org/10.2307/3144521)

Huff, D. L. (1964). Defining and Estimating a Trading Area. Journal of
Marketing, 28(3), 34-38.
[doi:10.1177/002224296402800307](https://doi.org/10.1177/002224296402800307)

## Examples

``` r
# \donttest{
library(sf)

# Create demand grid with spending potential
demand <- st_as_sf(expand.grid(x = 1:10, y = 1:10), coords = c("x", "y"))
demand$spending <- runif(100, 1000, 5000)

# Existing stores with varying sizes (attractiveness)
stores <- st_as_sf(data.frame(
  id = c("Store_A", "Store_B", "Store_C"),
  sqft = c(50000, 25000, 75000),
  parking = c(200, 100, 300),
  x = c(2, 8, 5), y = c(2, 8, 5)
), coords = c("x", "y"))

# Single attractiveness variable
result <- huff(demand, stores,
               attractiveness_col = "sqft",
               distance_exponent = -2,
               sales_potential_col = "spending")

# Multiple attractiveness variables with different exponents
# Composite: A = sqft^1.0 * parking^0.5
result_multi <- huff(demand, stores,
                     attractiveness_col = c("sqft", "parking"),
                     attractiveness_exponent = c(1.0, 0.5),
                     distance_exponent = -2,
                     sales_potential_col = "spending")

# View market shares
result_multi$stores[, c("id", "sqft", "parking", ".market_share", ".expected_sales")]
#> Simple feature collection with 3 features and 5 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 2 ymin: 2 xmax: 8 ymax: 8
#> CRS:           NA
#>        id  sqft parking .market_share .expected_sales    geometry
#> 1 Store_A 50000     200     0.2401693        70463.39 POINT (2 2)
#> 2 Store_B 25000     100     0.1624989        47675.65 POINT (8 8)
#> 3 Store_C 75000     300     0.5973318       175251.51 POINT (5 5)

# Evaluate a new candidate store
candidate <- st_as_sf(data.frame(
  id = "New_Store", sqft = 40000, parking = 250, x = 3, y = 7
), coords = c("x", "y"))

all_stores <- rbind(stores, candidate)
result_with_candidate <- huff(demand, all_stores,
                              attractiveness_col = c("sqft", "parking"),
                              attractiveness_exponent = c(1.0, 0.5),
                              distance_exponent = -2,
                              sales_potential_col = "spending")

# Compare market shares with and without candidate
result_with_candidate$stores[, c("id", ".market_share")]
#> Simple feature collection with 4 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 2 ymin: 2 xmax: 8 ymax: 8
#> CRS:           NA
#>          id .market_share    geometry
#> 1   Store_A     0.1965657 POINT (2 2)
#> 2   Store_B     0.1343459 POINT (8 8)
#> 3   Store_C     0.4432445 POINT (5 5)
#> 4 New_Store     0.2258438 POINT (3 7)
# }
```
