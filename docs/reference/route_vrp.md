# Vehicle Routing Problem (VRP)

Solves the Capacitated Vehicle Routing Problem (CVRP): find minimum-cost
routes for a fleet of vehicles, each with a capacity limit, to serve all
customers from a central depot. Uses Clarke-Wright savings heuristic for
construction with 2-opt, relocate, and swap local search improvement.

## Usage

``` r
route_vrp(
  locations,
  depot = 1L,
  demand_col,
  vehicle_capacity,
  n_vehicles = NULL,
  cost_matrix = NULL,
  distance_metric = "euclidean",
  method = "2-opt",
  service_time = NULL,
  max_route_time = NULL,
  balance = NULL,
  earliest = NULL,
  latest = NULL
)
```

## Arguments

- locations:

  An sf object representing all locations (depot + customers).

- depot:

  Integer. Row index of the depot in `locations` (1-based). Defaults to
  1.

- demand_col:

  Character. Column name in `locations` containing the demand at each
  stop. The depot's demand is ignored (set to 0 internally).

- vehicle_capacity:

  Numeric. Maximum capacity per vehicle.

- n_vehicles:

  Integer or NULL. Maximum number of vehicles. If NULL (default), uses
  as many as needed to satisfy capacity constraints.

- cost_matrix:

  Optional. Pre-computed square cost/distance matrix (n x n). If NULL,
  computed from `locations` using `distance_metric`.

- distance_metric:

  Distance metric when computing from geometry: "euclidean" (default) or
  "manhattan".

- method:

  Algorithm: "2-opt" (default, savings + local search) or "savings"
  (Clarke-Wright construction only).

- service_time:

  Optional service/dwell time at each stop. Supply either a numeric
  vector of length `nrow(locations)` or a column name in `locations`.
  The depot's service time is forced to 0.

- max_route_time:

  Optional maximum total time per route (travel time + service time).
  Routes that would exceed this limit are split.

- balance:

  Optional balancing mode. Currently supports `"time"` to minimize the
  longest route's total time after cost optimization. This runs a
  post-optimization phase that may slightly increase total cost (bounded
  to a small percentage). Ignored when `method = "savings"`.

- earliest:

  Optional earliest arrival/service time at each stop. Supply either a
  numeric vector of length `nrow(locations)` or a column name in
  `locations`. Must be supplied together with `latest`.

- latest:

  Optional latest arrival/service time at each stop. Must be supplied
  together with `earliest`. A vehicle may arrive early and wait, but
  cannot begin service after this time.

## Value

An sf object (the input `locations`) with added columns:

- `.vehicle`: Vehicle assignment (1, 2, ...). Depot is 0.

- `.visit_order`: Visit sequence within each vehicle's route.

Metadata is stored in the "spopt" attribute, including `n_vehicles`,
`total_cost`, `total_time`, per-vehicle `vehicle_costs`,
`vehicle_times`, `vehicle_loads`, and `vehicle_stops`.

## Details

The CVRP extends the TSP to multiple vehicles with capacity constraints.
The solver uses a two-phase approach:

1.  **Construction** (Clarke-Wright savings): Starts with one route per
    customer, then iteratively merges routes that produce the greatest
    distance savings while respecting capacity limits.

2.  **Improvement** (if `method = "2-opt"`): Applies intra-route 2-opt
    (reversing subtours), inter-route relocate (moving a customer
    between routes), and inter-route swap (exchanging customers between
    routes).

## Units and cost matrices

The `cost_matrix` can contain travel times, distances, or any
generalized cost. The solver minimizes the sum of matrix values along
each route. When using time-based features (`service_time`,
`max_route_time`, `earliest`/`latest`, `balance = "time"`), the cost
matrix should be in time-compatible units (e.g., minutes). Otherwise the
time-based outputs and constraints will be dimensionally inconsistent.

In the per-vehicle summary, **Cost** is the raw matrix objective (travel
only), while **Time** adds service duration and any waiting induced by
time windows. When no service time or windows are used, Time equals Cost
and is not shown separately.

## Use Cases

- **Oilfield logistics**: Route vacuum trucks to well pads with fluid
  volume constraints

- **Delivery routing**: Multiple delivery vehicles with weight/volume
  limits

- **Service dispatch**: Assign and sequence jobs across a fleet of
  technicians or crews

- **Waste collection**: Route collection vehicles with load capacity

## See also

[`route_tsp()`](https://walker-data.com/spopt-r/reference/route_tsp.md)
for single-vehicle routing

## Examples

``` r
# \donttest{
library(sf)

# Depot + 20 customers with demands
locations <- st_as_sf(
  data.frame(
    id = 1:21, x = runif(21), y = runif(21),
    demand = c(0, rpois(20, 10))
  ),
  coords = c("x", "y")
)

result <- route_vrp(locations, depot = 1, demand_col = "demand", vehicle_capacity = 40)

# How many vehicles needed?
attr(result, "spopt")$n_vehicles
#> [1] 7

# Per-vehicle costs
attr(result, "spopt")$vehicle_costs
#> [1] 1.94 1.20 1.43 1.64 1.32 2.15 0.54
# }
```
