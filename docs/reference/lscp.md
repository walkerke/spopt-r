# Location Set Covering Problem (LSCP)

Solves the Location Set Covering Problem: find the minimum number of
facilities needed to cover all demand points within a given service
radius.

## Usage

``` r
lscp(
  demand,
  facilities,
  service_radius,
  cost_matrix = NULL,
  distance_metric = "euclidean",
  verbose = FALSE
)
```

## Arguments

- demand:

  An sf object representing demand points (or polygons, using
  centroids).

- facilities:

  An sf object representing candidate facility locations.

- service_radius:

  Numeric. Maximum distance for a facility to cover a demand point.

- cost_matrix:

  Optional. Pre-computed distance matrix (demand x facilities). If NULL,
  computed from geometries.

- distance_metric:

  Distance metric: "euclidean" (default) or "manhattan".

- verbose:

  Logical. Print solver progress.

## Value

A list with two sf objects:

- `$demand`: Original demand sf with `.covered` column (logical),
  `.facility` (nearest selected facility), and `.cost` (cost to that
  facility, from the original cost matrix; `NA` if uncovered)

- `$facilities`: Original facilities sf with `.selected` column
  (logical)

Metadata is stored in the "spopt" attribute.

## Details

The LSCP minimizes the number of facilities required to ensure that
every demand point is within the service radius of at least one
facility. This is a mandatory coverage model where full coverage is
required.

The integer programming formulation is: \$\$\min \sum_j y_j\$\$ Subject
to: \$\$\sum_j a\_{ij} y_j \geq 1 \quad \forall i\$\$ \$\$y_j \in
\\0,1\\\$\$

Where \\y_j = 1\\ if facility j is selected, and \\a\_{ij} = 1\\ if
facility j can cover demand point i (distance \\\leq\\ service radius).

## Use Cases

LSCP is appropriate when complete coverage is mandatory:

- **Emergency services**: Fire stations, ambulance depots, or hospitals
  where every resident must be reachable within a response time standard

- **Public services**: Schools, polling places, or post offices where
  universal access is required by law or policy

- **Infrastructure**: Cell towers or utility substations where gaps in
  coverage are unacceptable

- **Retail/logistics**: Warehouse locations to ensure all customers can
  receive same-day or next-day delivery

For situations where complete coverage is not required or not feasible
within budget constraints, consider
[`mclp()`](https://walker-data.com/spopt-r/reference/mclp.md) instead.

## References

Toregas, C., Swain, R., ReVelle, C., & Bergman, L. (1971). The Location
of Emergency Service Facilities. Operations Research, 19(6), 1363-1373.
[doi:10.1287/opre.19.6.1363](https://doi.org/10.1287/opre.19.6.1363)

## See also

[`mclp()`](https://walker-data.com/spopt-r/reference/mclp.md) for
maximizing coverage with a fixed number of facilities

## Examples

``` r
# \donttest{
library(sf)

# Create demand and facility points
demand <- st_as_sf(data.frame(x = runif(50), y = runif(50)), coords = c("x", "y"))
facilities <- st_as_sf(data.frame(x = runif(10), y = runif(10)), coords = c("x", "y"))

# Find minimum facilities to cover all demand within 0.3 units
result <- lscp(demand, facilities, service_radius = 0.3)

# View selected facilities
result$facilities[result$facilities$.selected, ]
#> Simple feature collection with 6 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 0.04811234 ymin: 0.2235654 xmax: 0.9033195 ymax: 0.9433377
#> CRS:           NA
#>                        geometry .selected .n_assigned
#> 1   POINT (0.6124158 0.2804039)      TRUE          12
#> 3   POINT (0.7246128 0.5457883)      TRUE           5
#> 6    POINT (0.2448466 0.871377)      TRUE           6
#> 7   POINT (0.9033195 0.2235654)      TRUE           5
#> 9  POINT (0.04811234 0.3990994)      TRUE           9
#> 10  POINT (0.6471519 0.9433377)      TRUE           8
# }
```
