# Maximum Coverage Location Problem (MCLP)

Solves the Maximum Coverage Location Problem: maximize total weighted
demand covered by locating exactly p facilities.

## Usage

``` r
mclp(
  demand,
  facilities,
  service_radius,
  n_facilities,
  weight_col,
  cost_matrix = NULL,
  distance_metric = "euclidean",
  fixed_col = NULL,
  verbose = FALSE
)
```

## Arguments

- demand:

  An sf object representing demand points.

- facilities:

  An sf object representing candidate facility locations.

- service_radius:

  Numeric. Maximum distance for coverage.

- n_facilities:

  Integer. Number of facilities to locate (p).

- weight_col:

  Character. Column name in `demand` containing demand weights.

- cost_matrix:

  Optional. Pre-computed distance matrix.

- distance_metric:

  Distance metric: "euclidean" (default) or "manhattan".

- fixed_col:

  Optional column name in `facilities` indicating which facilities are
  pre-selected. The column should be logical (`TRUE` for fixed) or
  character (`"required"`/`"candidate"`). Fixed facilities are always
  selected; the solver optimizes the remaining slots. Useful for
  expansion planning where some facilities already exist.

- verbose:

  Logical. Print solver progress.

## Value

A list with two sf objects:

- `$demand`: Original demand sf with `.covered`, `.facility`, and
  `.cost` columns (`.cost` is the cost to the assigned facility, from
  the original cost matrix; `NA` if uncovered)

- `$facilities`: Original facilities sf with `.selected` column

Metadata is stored in the "spopt" attribute.

## Details

The MCLP maximizes the total weighted demand covered by locating exactly
p facilities. Unlike
[`lscp()`](https://walker-data.com/spopt-r/reference/lscp.md), which
requires full coverage, MCLP accepts partial coverage and optimizes for
the best possible outcome given a fixed budget (number of facilities).

The integer programming formulation is: \$\$\max \sum_i w_i z_i\$\$
Subject to: \$\$\sum_j y_j = p\$\$ \$\$z_i \leq \sum_j a\_{ij} y_j \quad
\forall i\$\$ \$\$y_j, z_i \in \\{0,1\\}\$\$

Where \\(w_i\\) is the weight (demand) at location i, \\(y_j = 1\\) if
facility j is selected, \\(z_i = 1\\) if demand i is covered, and
\\(a\_{ij} = 1\\) if facility j can cover demand i.

## Use Cases

MCLP is appropriate when you have a fixed budget or capacity constraint:

- **Healthcare access**: Locating p clinics to maximize the population
  within a 30-minute drive, given limited funding

- **Retail site selection**: Choosing p store locations to maximize the
  number of potential customers within a trade area

- **Emergency services**: Placing p ambulance stations to maximize the
  population reachable within an 8-minute response time

- **Conservation**: Selecting p reserve sites to maximize the number of
  species or habitat area protected

- **Telecommunications**: Locating p cell towers to maximize population
  coverage when full coverage is not economically feasible

For situations where complete coverage is required, use
[`lscp()`](https://walker-data.com/spopt-r/reference/lscp.md) to find
the minimum number of facilities needed.

## References

Church, R., & ReVelle, C. (1974). The Maximal Covering Location Problem.
Papers in Regional Science, 32(1), 101-118.
[doi:10.1007/BF01942293](https://doi.org/10.1007/BF01942293)

## See also

[`lscp()`](https://walker-data.com/spopt-r/reference/lscp.md) for
finding the minimum facilities needed for complete coverage

## Examples

``` r
# \donttest{
library(sf)

# Create demand with weights
demand <- st_as_sf(data.frame(
  x = runif(50), y = runif(50), population = rpois(50, 100)
), coords = c("x", "y"))
facilities <- st_as_sf(data.frame(x = runif(10), y = runif(10)), coords = c("x", "y"))

# Maximize population coverage with 3 facilities
result <- mclp(demand, facilities, service_radius = 0.3,
               n_facilities = 3, weight_col = "population")

attr(result, "spopt")$coverage_pct
#> [1] 78.15953
# }
```
