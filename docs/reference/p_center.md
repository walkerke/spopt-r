# P-Center Problem

Solves the P-Center problem: minimize the maximum distance from any
demand point to its nearest facility by locating exactly p facilities.
This is an equity-focused (minimax) objective that ensures no demand
point is too far from service.

## Usage

``` r
p_center(
  demand,
  facilities,
  n_facilities,
  cost_matrix = NULL,
  distance_metric = "euclidean",
  method = c("binary_search", "mip"),
  fixed_col = NULL,
  verbose = FALSE
)
```

## Arguments

- demand:

  An sf object representing demand points.

- facilities:

  An sf object representing candidate facility locations.

- n_facilities:

  Integer. Number of facilities to locate (p).

- cost_matrix:

  Optional. Pre-computed distance matrix.

- distance_metric:

  Distance metric: "euclidean" (default) or "manhattan".

- method:

  Algorithm to use: "binary_search" (default, faster) or "mip" (direct
  mixed-integer programming formulation).

- fixed_col:

  Optional column name in `facilities` indicating which facilities are
  pre-selected. The column should be logical (`TRUE` for fixed) or
  character (`"required"`/`"candidate"`). Fixed facilities are always
  selected; the solver optimizes the remaining slots.

- verbose:

  Logical. Print solver progress.

## Value

A list with two sf objects:

- `$demand`: Original demand sf with `.facility` column and `.cost`
  column (cost to the assigned facility, from the original cost matrix;
  `NA` if unassigned)

- `$facilities`: Original facilities sf with `.selected` column

Metadata includes `max_distance` (the objective value).

## Details

The p-center problem minimizes the maximum distance between any demand
point and its assigned facility. This "minimax" objective ensures
equitable access by focusing on the worst-served location rather than
average performance.

Two algorithms are available:

- `"binary_search"` (default): Binary search over distances with set
  covering subproblems. This converts the difficult minimax objective
  into simpler feasibility problems and is typically much faster.

- `"mip"`: Direct mixed-integer programming formulation with the minimax
  objective. Can be slower but may be preferred for small problems or
  when exact optimality certificates are needed.

The direct MIP formulation is: \$\$\min W\$\$ Subject to: \$\$\sum_j y_j
= p\$\$ \$\$\sum_j x\_{ij} = 1 \quad \forall i\$\$ \$\$x\_{ij} \leq y_j
\quad \forall i,j\$\$ \$\$\sum_j d\_{ij} x\_{ij} \leq W \quad \forall
i\$\$ \$\$x\_{ij}, y_j \in \\{0,1\\}\$\$

Where W is the maximum distance to minimize, \\(d\_{ij}\\) is the
distance from demand i to facility j, \\(x\_{ij} = 1\\) if demand i is
assigned to facility j, and \\(y_j = 1\\) if facility j is selected.

## Use Cases

P-center is appropriate when equity and worst-case performance matter:

- **Emergency services**: Fire stations or ambulance depots where
  response time standards must be met for all residents

- **Equity-focused planning**: Ensuring no community is underserved,
  even if it increases average travel distance

- **Critical infrastructure**: Backup facilities or emergency shelters
  where everyone must be within reach

- **Service level guarantees**: When contracts or regulations specify
  maximum acceptable distance or response time

For efficiency-focused objectives that minimize total travel, consider
[`p_median()`](https://walker-data.com/spopt-r/reference/p_median.md)
instead.

## References

Hakimi, S. L. (1965). Optimum Distribution of Switching Centers in a
Communication Network and Some Related Graph Theoretic Problems.
Operations Research, 13(3), 462-475.
[doi:10.1287/opre.13.3.462](https://doi.org/10.1287/opre.13.3.462)

## See also

[`p_median()`](https://walker-data.com/spopt-r/reference/p_median.md)
for minimizing total weighted distance (efficiency objective)

## Examples

``` r
# \donttest{
library(sf)

demand <- st_as_sf(data.frame(x = runif(50), y = runif(50)), coords = c("x", "y"))
facilities <- st_as_sf(data.frame(x = runif(15), y = runif(15)), coords = c("x", "y"))

# Minimize maximum distance with 4 facilities
result <- p_center(demand, facilities, n_facilities = 4)

# Maximum distance any demand point must travel
attr(result, "spopt")$max_distance
#> [1] 0.3515138
# }
```
