# Flow Refueling Location Model (FRLM)

Solves the Flow Refueling Location Model to optimally place refueling
facilities along network paths. This model maximizes the volume of
origin-destination flows that can be served given vehicle range
constraints.

## Usage

``` r
frlm(
  flows,
  candidates,
  network = NULL,
  vehicle_range,
  n_facilities,
  method = c("greedy"),
  verbose = FALSE
)
```

## Arguments

- flows:

  A data frame or sf object containing flow information with columns:

  - `origin`: Origin identifier

  - `destination`: Destination identifier

  - `volume`: Flow volume (e.g., number of trips)

- candidates:

  An sf object with candidate facility locations (points).

- network:

  Optional. A distance matrix between candidates. If NULL (default),
  Euclidean distances are computed from candidate geometries. For
  network distances, compute externally using packages like r5r or dodgr
  and pass the resulting matrix here.

- vehicle_range:

  Numeric. Maximum vehicle range (same units as network distances).

- n_facilities:

  Integer. Number of facilities to place.

- method:

  Character. Optimization method: "greedy" (default and currently only
  option).

- verbose:

  Logical. Print progress messages.

## Value

A list with class "spopt_frlm" containing:

- `facilities`: The candidates sf object with a `.selected` column

- `selected_indices`: 1-based indices of selected facilities

- `coverage`: Coverage statistics

Metadata is stored in the "spopt" attribute.

## Details

The Flow Refueling Location Model (Kuby & Lim, 2005) addresses the
problem of locating refueling stations for range-limited vehicles (e.g.,
electric vehicles, hydrogen fuel cell vehicles) along travel paths.

A flow (origin-destination path) is "covered" if a vehicle can complete
the **round trip** with refueling stops at the selected facilities. The
model assumes:

- Vehicles start at the origin with **half a tank** (can travel R/2)

- At each open station, vehicles refuel to full (can travel R)

- The round trip must be completable without running out of fuel

For a flow to be covered, three conditions must be met:

1.  First open station must be within R/2 from origin (half-tank start)

2.  Each subsequent open station must be within R of the previous

3.  Last open station must be within R/2 of destination (to allow
    return)

This implementation uses a greedy heuristic that iteratively selects the
facility providing the greatest marginal increase in covered flow
volume.

## Input Format

For simple cases, you can provide:

- `flows`: Data frame with origin, destination, volume

- `candidates`: sf points for potential facility locations

- `network`: Pre-computed distance matrix (optional)

## References

Kuby, M., & Lim, S. (2005). The flow-refueling location problem for
alternative-fuel vehicles. Socio-Economic Planning Sciences, 39(2),
125-145.
[doi:10.1016/j.seps.2004.03.001](https://doi.org/10.1016/j.seps.2004.03.001)

Capar, I., & Kuby, M. (2012). An efficient formulation of the flow
refueling location model for alternative-fuel stations. IIE
Transactions, 44(8), 622-636.
[doi:10.1080/0740817X.2011.635175](https://doi.org/10.1080/0740817X.2011.635175)

## Examples

``` r
# \donttest{
# Simple example with distance matrix
library(sf)

# Create candidate locations
candidates <- st_as_sf(data.frame(
  id = 1:10,
  x = runif(10, 0, 100),
  y = runif(10, 0, 100)
), coords = c("x", "y"))

# Create flows (using candidate indices as origins/destinations)
flows <- data.frame(
  origin = c(1, 1, 3, 5),
  destination = c(8, 10, 7, 9),
  volume = c(100, 200, 150, 300)
)

# Solve with vehicle range of 50 units
result <- frlm(flows, candidates, vehicle_range = 50, n_facilities = 3)

# View selected facilities
result$facilities[result$facilities$.selected, ]
#> Simple feature collection with 1 feature and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 30.21872 ymin: 27.714 xmax: 30.21872 ymax: 27.714
#> CRS:           NA
#>    id                geometry .selected
#> 10 10 POINT (30.21872 27.714)      TRUE
# }
```
