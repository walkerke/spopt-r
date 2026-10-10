# P-Median Problem

Solves the P-Median problem: minimize total weighted distance from
demand points to their assigned facilities by locating exactly p
facilities. This is an efficiency-focused objective that minimizes
overall travel burden.

## Usage

``` r
p_median(
  demand,
  facilities,
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
  selected; the solver optimizes the remaining slots.

- verbose:

  Logical. Print solver progress.

## Value

A list with two sf objects:

- `$demand`: Original demand sf with `.facility` column (assigned
  facility) and `.cost` column (cost to the assigned facility, from the
  original cost matrix)

- `$facilities`: Original facilities sf with `.selected` and
  `.n_assigned` columns

Metadata is stored in the "spopt" attribute.

## Details

The p-median problem minimizes the total weighted distance (or travel
cost) between demand points and their nearest assigned facility. It is
the most widely used location model for efficiency-oriented facility
siting.

The integer programming formulation is: \$\$\min \sum_i \sum_j w_i
d\_{ij} x\_{ij}\$\$ Subject to: \$\$\sum_j y_j = p\$\$ \$\$\sum_j
x\_{ij} = 1 \quad \forall i\$\$ \$\$x\_{ij} \leq y_j \quad \forall
i,j\$\$ \$\$x\_{ij}, y_j \in \\{0,1\\}\$\$

Where \\(w_i\\) is the demand weight at location i, \\(d\_{ij}\\) is the
distance from demand i to facility j, \\(x\_{ij} = 1\\) if demand i is
assigned to facility j, and \\(y_j = 1\\) if facility j is selected.

## Use Cases

P-median is appropriate when minimizing total travel cost or distance:

- **Public facilities**: Schools, libraries, or community centers where
  the goal is to minimize total student/patron travel

- **Warehouses and distribution**: Locating distribution centers to
  minimize total shipping costs to customers

- **Healthcare**: Positioning clinics to minimize aggregate patient
  travel time across a population

- **Service depots**: Locating maintenance facilities to minimize total
  technician travel to service calls

For equity-focused objectives where no demand point should be too far,
consider
[`p_center()`](https://walker-data.com/spopt-r/reference/p_center.md)
instead.

## References

Hakimi, S. L. (1964). Optimum Locations of Switching Centers and the
Absolute Centers and Medians of a Graph. Operations Research, 12(3),
450-459.
[doi:10.1287/opre.12.3.450](https://doi.org/10.1287/opre.12.3.450)

## See also

[`p_center()`](https://walker-data.com/spopt-r/reference/p_center.md)
for minimizing maximum distance (equity objective)

## Examples

``` r
# \donttest{
library(sf)

demand <- st_as_sf(data.frame(
  x = runif(100), y = runif(100), population = rpois(100, 500)
), coords = c("x", "y"))
facilities <- st_as_sf(data.frame(x = runif(20), y = runif(20)), coords = c("x", "y"))

# Locate 5 facilities minimizing total weighted distance
result <- p_median(demand, facilities, n_facilities = 5, weight_col = "population")

# Mean distance to assigned facility
attr(result, "spopt")$mean_distance
#> [1] 0.1743598
# }
```
