# Capacitated Facility Location Problem (CFLP)

Solves the Capacitated Facility Location Problem: minimize total
weighted distance from demand points to facilities, subject to capacity
constraints at each facility. Unlike standard p-median, facilities have
limited capacity and demand may need to be split across multiple
facilities.

## Usage

``` r
cflp(
  demand,
  facilities,
  n_facilities,
  weight_col,
  capacity_col,
  facility_cost_col = NULL,
  cost_matrix = NULL,
  distance_metric = "euclidean",
  verbose = FALSE
)
```

## Arguments

- demand:

  An sf object representing demand points.

- facilities:

  An sf object representing candidate facility locations.

- n_facilities:

  Integer. Number of facilities to locate. Set to 0 if using
  `facility_cost_col` to determine optimal number.

- weight_col:

  Character. Column name in `demand` containing demand weights (e.g.,
  population, customers, volume).

- capacity_col:

  Character. Column name in `facilities` containing capacity of each
  facility.

- facility_cost_col:

  Optional character. Column name in `facilities` containing fixed cost
  to open each facility. If provided and `n_facilities = 0`, the solver
  determines the optimal number of facilities to minimize total cost.

- cost_matrix:

  Optional. Pre-computed distance/cost matrix (demand x facilities).

- distance_metric:

  Distance metric: "euclidean" (default) or "manhattan".

- verbose:

  Logical. Print solver progress.

## Value

A list with two sf objects:

- `$demand`: Original demand sf with `.facility` column (primary
  assignment) and `.split` column (TRUE if demand is split across
  facilities), and `.cost` column (cost to the primary facility, from
  the original cost matrix)

- `$facilities`: Original facilities sf with `.selected`, `.n_assigned`,
  and `.utilization` columns

Metadata is stored in the "spopt" attribute, including:

- `objective`: Total cost (transportation + facility costs if
  applicable)

- `mean_distance`: Mean weighted distance

- `n_split_demand`: Number of demand points split across facilities

- `allocation_matrix`: Full allocation matrix (n_demand x n_facilities)

## Details

The CFLP extends the p-median problem by adding capacity constraints.
Each facility \\j\\ has a maximum capacity \\Q_j\\, and the total demand
assigned to it cannot exceed this capacity.

When demand exceeds available capacity at the nearest facility, the
solver may split demand across multiple facilities. The `.split` column
indicates which demand points have been split, and the
`allocation_matrix` in metadata shows the exact fractions.

Two modes of operation:

1.  **Fixed number**: Set `n_facilities` to select exactly that many
    facilities

2.  **Cost-based**: Set `n_facilities = 0` and provide
    `facility_cost_col` to let the solver determine the optimal number
    based on fixed + variable costs

## References

Daskin, M. S. (2013). Network and discrete location: Models, algorithms,
and applications (2nd ed.). John Wiley & Sons.
[doi:10.1002/9781118537015](https://doi.org/10.1002/9781118537015)

Sridharan, R. (1995). The capacitated plant location problem. European
Journal of Operational Research, 87(2), 203-213.
[doi:10.1016/0377-2217(95)00042-O](https://doi.org/10.1016/0377-2217%2895%2900042-O)

## Examples

``` r
# \donttest{
library(sf)

# Demand points with population
demand <- st_as_sf(data.frame(
  x = runif(100), y = runif(100), population = rpois(100, 500)
), coords = c("x", "y"))

# Facilities with varying capacities
facilities <- st_as_sf(data.frame(
  x = runif(15), y = runif(15),
  capacity = c(rep(5000, 5), rep(10000, 5), rep(20000, 5)),
  fixed_cost = c(rep(100, 5), rep(200, 5), rep(400, 5))
), coords = c("x", "y"))

# Fixed number of facilities
result <- cflp(demand, facilities, n_facilities = 5,
               weight_col = "population", capacity_col = "capacity")

# Check utilization
result$facilities[result$facilities$.selected, c("capacity", ".utilization")]
#> Simple feature collection with 5 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 0.01087246 ymin: 0.1990008 xmax: 0.7607611 ymax: 0.7502479
#> CRS:           NA
#>    capacity .utilization                     geometry
#> 7     10000      1.00000  POINT (0.6332316 0.7502479)
#> 11    20000      0.27815  POINT (0.7607611 0.5438312)
#> 12    20000      0.35185 POINT (0.01087246 0.1990008)
#> 14    20000      0.76560  POINT (0.2638418 0.7502398)
#> 15    20000      0.59740  POINT (0.5938366 0.2709778)

# Cost-based (optimal number of facilities)
result <- cflp(demand, facilities, n_facilities = 0,
               weight_col = "population", capacity_col = "capacity",
               facility_cost_col = "fixed_cost")
attr(result, "spopt")$n_selected
#> [1] 10
# }
```
