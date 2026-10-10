# P-Dispersion Problem

Solves the P-Dispersion problem: maximize the minimum distance between
any two selected facilities. This "maximin" objective ensures facilities
are spread out as much as possible.

## Usage

``` r
p_dispersion(
  facilities,
  n_facilities,
  cost_matrix = NULL,
  distance_metric = "euclidean",
  verbose = FALSE
)
```

## Arguments

- facilities:

  An sf object representing candidate facility locations. Note: This
  problem does not use demand points.

- n_facilities:

  Integer. Number of facilities to locate (p).

- cost_matrix:

  Optional. Pre-computed inter-facility distance matrix.

- distance_metric:

  Distance metric: "euclidean" (default) or "manhattan".

- verbose:

  Logical. Print solver progress.

## Value

An sf object (the facilities input) with a `.selected` column. Metadata
includes `min_distance` (the objective value).

## Details

The p-dispersion problem selects p facilities from a set of candidates
such that the minimum pairwise distance between any two selected
facilities is maximized. Unlike p-median or p-center, this problem does
not consider demand points—it focuses solely on spreading facilities
apart.

The mixed integer programming formulation uses a Big-M approach:
\$\$\max D\$\$ Subject to: \$\$\sum_j y_j = p\$\$ \$\$D \leq d\_{ij} +
M(2 - y_i - y_j) \quad \forall i \< j\$\$ \$\$y_j \in \\{0,1\\}, \quad D
\geq 0\$\$

Where D is the minimum separation distance to maximize, \\(d\_{ij}\\) is
the distance between facilities i and j, \\(y_j = 1\\) if facility j is
selected, and M is a large constant. When both facilities i and j are
selected (\\(y_i = y_j = 1\\)), the constraint reduces to \\(D \leq
d\_{ij}\\), ensuring D is at most the distance between any pair of
selected facilities.

## Use Cases

P-dispersion is appropriate when facilities should be spread apart:

- **Obnoxious facilities**: Hazardous waste sites, prisons, or other
  undesirable facilities that should be separated from each other

- **Franchise territories**: Retail locations where stores should not
  cannibalize each other's market

- **Redundant systems**: Backup servers or emergency caches that should
  be geographically distributed for resilience

- **Monitoring networks**: Air quality sensors or seismic monitors that
  should cover distinct areas without overlap

- **Spatial sampling**: Selecting representative sample locations that
  are well-distributed across a study area

## References

Kuby, M. J. (1987). Programming Models for Facility Dispersion: The
p-Dispersion and Maxisum Dispersion Problems. Geographical Analysis,
19(4), 315-329.
[doi:10.1111/j.1538-4632.1987.tb00133.x](https://doi.org/10.1111/j.1538-4632.1987.tb00133.x)

## Examples

``` r
# \donttest{
library(sf)

facilities <- st_as_sf(data.frame(x = runif(20), y = runif(20)), coords = c("x", "y"))

# Select 5 facilities maximally dispersed
result <- p_dispersion(facilities, n_facilities = 5)

# Minimum distance between any two selected facilities
attr(result, "spopt")$min_distance
#> [1] 0.4916683
# }
```
