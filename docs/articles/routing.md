# Route optimization with spopt

*Route optimization* asks: given a set of stops that need to be visited,
what’s the best sequence? This is the Traveling Salesman Problem (TSP)
when you have one vehicle, and the Vehicle Routing Problem (VRP) when
you have a fleet. These problems show up everywhere - delivery routing,
field service scheduling, sales territory coverage, equipment inspection
rounds.

spopt’s routing functions take a travel-time matrix (from r5r, OSRM,
mapboxapi, or any routing engine) and find the optimal stop sequence.
The solvers handle real-world constraints: open-ended routes, fixed
start/end points, vehicle capacities, and time windows.

## Setting up the problem

Route optimization has three components:

1.  **Stops**: Locations that need to be visited
2.  **Depot**: Where vehicles start (and possibly return)
3.  **Cost matrix**: Travel time or distance between every pair of stops

The cost matrix is the key input. Euclidean distance works for quick
prototyping, but real routing needs real travel times - a stop across a
river might be 2 miles as the crow flies but 15 minutes by car. That’s
where an external routing engine comes in.

Let’s work through a realistic scenario: a delivery station in west Fort
Worth dispatching drivers to 25 residential addresses across Fort Worth,
Benbrook, and Aledo.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`spopt`](https://walker-data.com/spopt-r/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyverse`](https://tidyverse.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sf`](https://r-spatial.github.io/sf/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`mapgl`](https://walker-data.com/mapgl/)`)`\
\
`# Load the bundled delivery data`\
[`data`](https://rdrr.io/r/utils/data.html)`(``delivery_data``)`\
\
`stops`` ``<-`` ``delivery_data``$``stops`\
`ttm`` ``<-`` ``delivery_data``$``matrix`\
\
`# Take a look at what we have`\
`stops`` ``|>`\
`  `[`st_drop_geometry`](https://r-spatial.github.io/sf/reference/st_geometry.html)`(``)`` ``|>`\
`  `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``id``, ``address``, ``packages``)`

    # A tibble: 26 × 3
       id    address                                    packages
       <chr> <chr>                                         <int>
     1 depot 3700 San Jacinto Dr, Fort Worth, TX 76116         0
     2 D01   6301 Camp Bowie Blvd, Fort Worth, TX 76116        3
     3 D02   5800 Lovell Ave, Fort Worth, TX 76107             3
     4 D03   4200 Birchman Ave, Fort Worth, TX 76107           4
     5 D04   6100 Wester Ave, Fort Worth, TX 76116             3
     6 D05   3516 Corto Ave, Fort Worth, TX 76109              5
     7 D06   7400 Oakmont Blvd, Fort Worth, TX 76132           5
     8 D07   7200 Oakmont Blvd, Fort Worth, TX 76132           3
     9 D08   8517 Calmont Ave, Fort Worth, TX 76116            5
    10 D09   3905 Las Vegas Trail, Fort Worth, TX 76116        4
    # ℹ 16 more rows

The `delivery_data` dataset includes 26 locations: an Amazon delivery
station (the depot) and 25 residential stops across western Tarrant and
Parker counties. Each stop has a package count, and `ttm` is a 26x26
travel-time matrix in minutes computed with r5r from an OpenStreetMap
road network.

\
`# Depot is the first row`\
`depot`` ``<-`` ``stops`` ``|>`` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``id`` ``==`` ``"depot"``)`\
`deliveries`` ``<-`` ``stops`` ``|>`` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``id`` ``!=`` ``"depot"``)`\
\
[`maplibre`](https://walker-data.com/mapgl/reference/maplibre.html)`(``style ``=`` `[`openfreemap_style`](https://walker-data.com/mapgl/reference/openfreemap_style.html)`(``"bright"``)``, bounds ``=`` ``stops``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"deliveries"``,`\
`    source ``=`` ``deliveries``,`\
`    circle_color ``=`` ``"steelblue"``,`\
`    circle_radius ``=`` ``6``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``2``,`\
`    tooltip ``=`` ``"address"`\
`  ``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"depot"``,`\
`    source ``=`` ``depot``,`\
`    circle_color ``=`` ``"#e41a1c"``,`\
`    circle_radius ``=`` ``10``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``3``,`\
`    tooltip ``=`` ``"address"`\
`  ``)`

The red marker is the depot (Amazon DDA9 on San Jacinto Dr). Blue
markers are delivery stops spread from central Fort Worth west into
Aledo and Willow Park.

## Generating the travel-time matrix

The bundled matrix was built with r5r using the same workflow described
in the [travel-time matrices
vignette](https://walker-data.com/spopt-r/articles/travel-time-matrices.md).
The key difference for routing is that you need a **square** matrix -
every stop to every other stop - rather than the rectangular
demand-to-facility matrix used in facility location.

\
`# How the matrix was generated (requires Java 21 + OSM data)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`r5r`](https://github.com/ipeaGIT/r5r)`)`\
[`options`](https://rdrr.io/r/base/options.html)`(``java.parameters ``=`` ``"-Xmx4G"``)`\
`rJavaEnv``::`[`java_quick_install`](https://www.ekotov.pro/rJavaEnv/reference/java_quick_install.html)`(``version ``=`` ``21``)`\
\
`r5r_core`` ``<-`` `[`build_network`](https://ipeagit.github.io/r5r/reference/build_network.html)`(``data_path ``=`` ``"path/to/osm/directory"``)`\
\
`# Prepare points for r5r`\
`r5r_pts`` ``<-`` ``stops`` ``|>`\
`  `[`st_coordinates`](https://r-spatial.github.io/sf/reference/st_coordinates.html)`(``)`` ``|>`\
`  `[`as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)`(``)`` ``|>`\
`  `[`rename`](https://dplyr.tidyverse.org/reference/rename.html)`(``lon ``=`` ``X``, lat ``=`` ``Y``)`` ``|>`\
`  `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``id ``=`` ``stops``$``id``)`\
\
`# Square matrix: all stops to all stops`\
`ttm_long`` ``<-`` `[`travel_time_matrix`](https://ipeagit.github.io/r5r/reference/travel_time_matrix.html)`(`\
`  ``r5r_core``,`\
`  origins ``=`` ``r5r_pts``,`\
`  destinations ``=`` ``r5r_pts``,`\
`  mode ``=`` ``"CAR"``,`\
`  departure_datetime ``=`` `[`as.POSIXct`](https://rdrr.io/r/base/as.POSIXlt.html)`(``"2025-03-15 08:00:00"``)``,`\
`  max_trip_duration ``=`` ``120`\
`)`\
\
`# Reshape to matrix`\
`ttm`` ``<-`` ``ttm_long`` ``|>`\
`  `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``from_id``, ``to_id``, ``travel_time_p50``)`` ``|>`\
`  `[`pivot_wider`](https://tidyr.tidyverse.org/reference/pivot_wider.html)`(``names_from ``=`` ``to_id``, values_from ``=`` ``travel_time_p50``)`` ``|>`\
`  `[`arrange`](https://dplyr.tidyverse.org/reference/arrange.html)`(`[`match`](https://rspatial.github.io/terra/reference/match.html)`(``from_id``, ``stops``$``id``)``)`` ``|>`\
`  `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``-``from_id``)`` ``|>`\
`  `[`as.matrix`](https://rspatial.github.io/terra/reference/coerce.html)`(``)`

Any routing engine that produces pairwise travel times works here; r5r
is a great choice, but there are many others with R bindings, including
OSRM, Valhalla (open source), Mapbox, and Google (commercial), among
others.

## Single-driver route: TSP

The simplest case: one driver, all 25 stops, return to the depot. The
[`route_tsp()`](https://walker-data.com/spopt-r/reference/route_tsp.md)
function finds the shortest sequence.

\
`result`` ``<-`` `[`route_tsp`](https://walker-data.com/spopt-r/reference/route_tsp.md)`(``stops``, start ``=`` ``1``, cost_matrix ``=`` ``ttm``)`\
\
`meta`` ``<-`` `[`attr`](https://rdrr.io/r/base/attr.html)`(``result``, ``"spopt"``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Optimized route: %.0f minutes\n"``, ``meta``$``total_cost``)``)`

    Optimized route: 135 minutes

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Nearest-neighbor baseline: %.0f minutes\n"``, ``meta``$``nn_cost``)``)`

    Nearest-neighbor baseline: 150 minutes

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Improvement: %.1f%%\n"``, ``meta``$``improvement_pct``)``)`

    Improvement: 10.0%

The solver starts with a nearest-neighbor heuristic (always go to the
closest unvisited stop), then improves it with 2-opt and or-opt local
search. The improvement over nearest-neighbor is typical - greedy
heuristics get you close, but the local search squeezes out the
remaining inefficiency.

Let’s look at the visit order:

\
`result`` ``|>`\
`  `[`st_drop_geometry`](https://r-spatial.github.io/sf/reference/st_geometry.html)`(``)`` ``|>`\
`  `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``.visit_order`` ``<=`` ``10``)`` ``|>`\
`  `[`arrange`](https://dplyr.tidyverse.org/reference/arrange.html)`(``.visit_order``)`` ``|>`\
`  `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``id``, ``address``, ``.visit_order``)`

    # A tibble: 10 × 3
       id    address                                    .visit_order
       <chr> <chr>                                             <int>
     1 depot 3700 San Jacinto Dr, Fort Worth, TX 76116             1
     2 D15   1201 Sproles Dr, Benbrook, TX 76126                   2
     3 D13   9200 Westpark Dr, Benbrook, TX 76126                  3
     4 D11   1100 Mercedes St, Benbrook, TX 76126                  4
     5 D12   1000 Winscott Rd, Benbrook, TX 76126                  5
     6 D07   7200 Oakmont Blvd, Fort Worth, TX 76132               6
     7 D06   7400 Oakmont Blvd, Fort Worth, TX 76132               7
     8 D23   4801 Old Granbury Rd, Fort Worth, TX 76133            8
     9 D04   6100 Wester Ave, Fort Worth, TX 76116                 9
    10 D25   6500 Woodway Dr, Fort Worth, TX 76133                10

The optimizer clusters geographically - it handles the Benbrook stops
together rather than bouncing back and forth. As the bundled
`delivery_data` includes pre-computed route geometries from r5r, we can
draw the actual driving paths:

\
`tsp_route`` ``<-`` ``delivery_data``$``tsp_route`\
\
`# Number stops by visit order`\
`result_ordered`` ``<-`` ``result`` ``|>`\
`  `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``.visit_order``)``)`` ``|>`\
`  `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``label ``=`` `[`as.character`](https://rdrr.io/r/base/character.html)`(``.visit_order``)``)`\
\
[`maplibre`](https://walker-data.com/mapgl/reference/maplibre.html)`(``style ``=`` `[`openfreemap_style`](https://walker-data.com/mapgl/reference/openfreemap_style.html)`(``"bright"``)``, bounds ``=`` ``stops``)`` ``|>`\
`  `[`add_line_layer`](https://walker-data.com/mapgl/reference/add_line_layer.html)`(`\
`    id ``=`` ``"route"``,`\
`    source ``=`` ``tsp_route``,`\
`    line_color ``=`` ``"#2563eb"``,`\
`    line_width ``=`` ``3``,`\
`    line_opacity ``=`` ``0.8`\
`  ``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"stops"``,`\
`    source ``=`` ``result_ordered`` ``|>`` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``id`` ``!=`` ``"depot"``)``,`\
`    circle_color ``=`` ``"#2563eb"``,`\
`    circle_radius ``=`` ``7``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``2``,`\
`    tooltip ``=`` ``"address"`\
`  ``)`` ``|>`\
`  `[`add_symbol_layer`](https://walker-data.com/mapgl/reference/add_symbol_layer.html)`(`\
`    id ``=`` ``"labels"``,`\
`    source ``=`` ``result_ordered`` ``|>`` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``id`` ``!=`` ``"depot"``)``,`\
`    text_field ``=`` `[`get_column`](https://walker-data.com/mapgl/reference/get_column.html)`(``"label"``)``,`\
`    text_size ``=`` ``11``,`\
`    text_color ``=`` ``"white"``,`\
`    text_halo_color ``=`` ``"#1e40af"``,`\
`    text_halo_width ``=`` ``1.5`\
`  ``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"depot"``,`\
`    source ``=`` ``depot``,`\
`    circle_color ``=`` ``"#dc2626"``,`\
`    circle_radius ``=`` ``10``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``3``,`\
`    tooltip ``=`` ``"address"`\
`  ``)`

The route follows actual roads between stops. Notice the optimized route
saves the far western stops (Aledo, Willow Park) for last before
returning to the depot.

### Open routes

Drivers don’t always return to the depot. A field technician might start
at the office and end at their last appointment. Set `end = NULL` for an
open route:

\
`result_open`` ``<-`` `[`route_tsp`](https://walker-data.com/spopt-r/reference/route_tsp.md)`(``stops``, start ``=`` ``1``, end ``=`` ``NULL``, cost_matrix ``=`` ``ttm``)`\
`meta_open`` ``<-`` `[`attr`](https://rdrr.io/r/base/attr.html)`(``result_open``, ``"spopt"``)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Closed route: %.0f minutes\n"``, ``meta``$``total_cost``)``)`

    Closed route: 135 minutes

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Open route:   %.0f minutes\n"``, ``meta_open``$``total_cost``)``)`

    Open route:   120 minutes

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Saved by not returning: %.0f minutes\n"``,`\
`            ``meta``$``total_cost`` ``-`` ``meta_open``$``total_cost``)``)`

    Saved by not returning: 15 minutes

Dropping the return leg saves real time, especially when the last stop
is far from the depot.

### Fixed start and end

Sometimes the start and end points are different - a courier picks up
from a warehouse and must end at a specific drop-off location. Use
`start` and `end` to fix both endpoints:

\
`# Start at depot (1), end at the farthest delivery (index of a Parker County stop)`\
`parker_stop`` ``<-`` `[`which`](https://rdrr.io/r/base/which.html)`(`[`grepl`](https://rdrr.io/r/base/grep.html)`(``"Willow Park|Hudson Oaks|Aledo"``, ``stops``$``address``)``)``[``1``]`\
\
`result_path`` ``<-`` `[`route_tsp`](https://walker-data.com/spopt-r/reference/route_tsp.md)`(``stops``, start ``=`` ``1``, end ``=`` ``parker_stop``, cost_matrix ``=`` ``ttm``)`\
`meta_path`` ``<-`` `[`attr`](https://rdrr.io/r/base/attr.html)`(``result_path``, ``"spopt"``)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Route type: %s\n"``, ``meta_path``$``route_type``)``)`

    Route type: path

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Total time: %.0f minutes\n"``, ``meta_path``$``total_cost``)``)`

    Total time: 124 minutes

## Fleet routing: VRP

Real delivery operations have multiple drivers, and each van has a
capacity limit. The
[`route_vrp()`](https://walker-data.com/spopt-r/reference/route_vrp.md)
function solves the capacitated vehicle routing problem: assign stops to
vehicles and sequence each vehicle’s route, all while respecting
capacity constraints.

Our stops have package counts ranging from 3 to 5. With a van capacity
of 35 packages:

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Total packages: %d\n"``, `[`sum`](https://rdrr.io/r/base/sum.html)`(``stops``$``packages``)``)``)`

    Total packages: 98

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Van capacity: 35\n"``)``)`

    Van capacity: 35

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Minimum vans needed: %d\n"``, `[`ceiling`](https://rdrr.io/r/base/Round.html)`(`[`sum`](https://rdrr.io/r/base/sum.html)`(``stops``$``packages``)`` ``/`` ``35``)``)``)`

    Minimum vans needed: 3

\
`result_vrp`` ``<-`` `[`route_vrp`](https://walker-data.com/spopt-r/reference/route_vrp.md)`(`\
`  ``stops``,`\
`  depot ``=`` ``1``,`\
`  demand_col ``=`` ``"packages"``,`\
`  vehicle_capacity ``=`` ``35``,`\
`  cost_matrix ``=`` ``ttm`\
`)`\
\
`meta_vrp`` ``<-`` `[`attr`](https://rdrr.io/r/base/attr.html)`(``result_vrp``, ``"spopt"``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Vehicles used: %d\n"``, ``meta_vrp``$``n_vehicles``)``)`

    Vehicles used: 3

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Total drive time: %.0f minutes\n"``, ``meta_vrp``$``total_cost``)``)`

    Total drive time: 167 minutes

The solver uses a Clarke-Wright savings heuristic for initial route
construction, then improves with intra-route 2-opt and or-opt
(resequencing stops within each route) and inter-route relocate and swap
(moving stops between vehicles).

\
[`summary`](https://rspatial.github.io/terra/reference/summary.html)`(``result_vrp``)`

    VRP routes: 26 locations, 3 vehicles (depot: 1)
      Method: 2-opt | Total cost: 167.0 | Improvement: 5.7%
      Capacity: 35
      Solve time: 0.000s

    Per-vehicle summary:
     Vehicle Stops Load Cost
           1     8   33   44
           2     9   34   59
           3     8   31   64

### Visualizing vehicle routes

The bundled data includes road geometries for each vehicle’s route.
Let’s see how the fleet covers the delivery area:

\
`vrp_route`` ``<-`` ``delivery_data``$``vrp_route`` ``|>`\
`  `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``vehicle ``=`` `[`as.character`](https://rdrr.io/r/base/character.html)`(``vehicle``)``)`\
\
`vrp_stops`` ``<-`` ``result_vrp`` ``|>`\
`  `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``.vehicle`` ``>`` ``0``)`` ``|>`\
`  `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``.vehicle ``=`` `[`as.character`](https://rdrr.io/r/base/character.html)`(``.vehicle``)``)`\
\
`vehicle_colors`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#e41a1c"``, ``"#377eb8"``, ``"#4daf4a"``)`\
`vehicle_ids`` ``<-`` `[`as.character`](https://rdrr.io/r/base/character.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(``meta_vrp``$``n_vehicles``)``)`\
\
[`maplibre`](https://walker-data.com/mapgl/reference/maplibre.html)`(``style ``=`` `[`openfreemap_style`](https://walker-data.com/mapgl/reference/openfreemap_style.html)`(``"bright"``)``, bounds ``=`` ``stops``)`` ``|>`\
`  `[`add_line_layer`](https://walker-data.com/mapgl/reference/add_line_layer.html)`(`\
`    id ``=`` ``"routes"``,`\
`    source ``=`` ``vrp_route``,`\
`    line_color ``=`` `[`match_expr`](https://walker-data.com/mapgl/reference/match_expr.html)`(`\
`      column ``=`` ``"vehicle"``,`\
`      values ``=`` ``vehicle_ids``,`\
`      stops ``=`` ``vehicle_colors``[`[`seq_along`](https://rdrr.io/r/base/seq.html)`(``vehicle_ids``)``]`\
`    ``)``,`\
`    line_width ``=`` ``3``,`\
`    line_opacity ``=`` ``0.8`\
`  ``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"stops"``,`\
`    source ``=`` ``vrp_stops``,`\
`    circle_color ``=`` `[`match_expr`](https://walker-data.com/mapgl/reference/match_expr.html)`(`\
`      column ``=`` ``".vehicle"``,`\
`      values ``=`` ``vehicle_ids``,`\
`      stops ``=`` ``vehicle_colors``[`[`seq_along`](https://rdrr.io/r/base/seq.html)`(``vehicle_ids``)``]`\
`    ``)``,`\
`    circle_radius ``=`` ``7``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``2``,`\
`    tooltip ``=`` ``"address"`\
`  ``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"depot"``,`\
`    source ``=`` ``depot``,`\
`    circle_color ``=`` ``"black"``,`\
`    circle_radius ``=`` ``10``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``3``,`\
`    tooltip ``=`` ``"address"`\
`  ``)`

Each color is a different van’s route, with all routes ending at the
depot (the black marker). You can see how the three routes are
partitioned geographically; the green driver handles the Parker County /
western stops, whereas the blue driver largely handles the southern
stops and the red driver the northern stops.

### Constraining the fleet size

If you have exactly 2 vans available, set `n_vehicles`:

\
`result_2vans`` ``<-`` `[`route_vrp`](https://walker-data.com/spopt-r/reference/route_vrp.md)`(`\
`  ``stops``,`\
`  depot ``=`` ``1``,`\
`  demand_col ``=`` ``"packages"``,`\
`  vehicle_capacity ``=`` ``50``,`\
`  n_vehicles ``=`` ``2``,`\
`  cost_matrix ``=`` ``ttm`\
`)`\
\
`meta_2`` ``<-`` `[`attr`](https://rdrr.io/r/base/attr.html)`(``result_2vans``, ``"spopt"``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Vehicles: %d\n"``, ``meta_2``$``n_vehicles``)``)`

    Vehicles: 2

\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Total time: %.0f min (vs %.0f min with %d vans)\n"``,`\
`            ``meta_2``$``total_cost``, ``meta_vrp``$``total_cost``, ``meta_vrp``$``n_vehicles``)``)`

    Total time: 192 min (vs 167 min with 3 vans)

Fewer vehicles means longer routes per driver. The solver finds the best
assignment given the constraint, but you’ll want to check that the
resulting route times are realistic for a single shift.

### Shift limits and service time

If drivers have a 90-minute shift limit, set `max_route_time`. The
solver will split routes so that no vehicle exceeds the time budget. You
can also account for time spent at each stop with `service_time` –
loading/unloading, signatures, etc.

\
`n`` ``<-`` `[`nrow`](https://rspatial.github.io/terra/reference/dimensions.html)`(``stops``)`\
`result_shift`` ``<-`` `[`route_vrp`](https://walker-data.com/spopt-r/reference/route_vrp.md)`(`\
`  ``stops``,`\
`  depot ``=`` ``1``,`\
`  demand_col ``=`` ``"packages"``,`\
`  vehicle_capacity ``=`` ``50``,`\
`  cost_matrix ``=`` ``ttm``,`\
`  service_time ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``3``, ``n``)``,   ``# 3 minutes at each stop`\
`  max_route_time ``=`` ``90``         ``# 90-minute shift limit`\
`)`\
\
[`summary`](https://rspatial.github.io/terra/reference/summary.html)`(``result_shift``)`

    VRP routes: 26 locations, 3 vehicles (depot: 1)
      Method: 2-opt | Total cost: 162.0 | Improvement: 6.4%
      Capacity: 50 | Max route time: 90
      Solve time: 0.000s

    Per-vehicle summary:
     Vehicle Stops Load Cost Time
           1     9   38   54   81
           2     9   34   50   77
           3     7   26   58   79

      Cost = matrix objective (travel only); Time = travel + service + waiting

`max_route_time` is a hard constraint: total time (travel + service) on
every route stays within the limit. The objective is still minimizing
total travel time, so the solver won’t use extra vehicles unless the
time constraint forces it.

### Balancing route times

The cost optimizer finds the cheapest set of routes, but that can
produce uneven workloads – one driver finishes in 30 minutes while
another is out for 80. Setting `balance = "time"` runs a
post-optimization phase that redistributes stops to reduce the longest
route time, at the cost of a small bounded cost increase.

\
`result_balanced`` ``<-`` `[`route_vrp`](https://walker-data.com/spopt-r/reference/route_vrp.md)`(`\
`  ``stops``,`\
`  depot ``=`` ``1``,`\
`  demand_col ``=`` ``"packages"``,`\
`  vehicle_capacity ``=`` ``35``,`\
`  cost_matrix ``=`` ``ttm``,`\
`  service_time ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``3``, ``n``)``,`\
`  balance ``=`` ``"time"`\
`)`\
\
[`summary`](https://rspatial.github.io/terra/reference/summary.html)`(``result_balanced``)`

    VRP routes: 26 locations, 3 vehicles (depot: 1)
      Method: 2-opt | Total cost: 167.0 | Improvement: 5.7%
      Capacity: 35 | Balance: time (1 move)
      Solve time: 0.000s

    Per-vehicle summary:
     Vehicle Stops Load Cost Time
           1     8   34   45   69
           2     9   34   59   86
           3     8   30   63   87

      Cost = matrix objective (travel only); Time = travel + service + waiting

The balancing phase only runs when `method = "2-opt"` (the default). It
targets the longest route and tries moving or swapping stops to shorten
it. Capacity and time constraints are still respected.

### Time windows

If customers have availability windows, set `earliest` and `latest`. The
solver respects these when constructing and improving routes. A vehicle
may arrive early and wait, but it cannot begin service after the
`latest` time.

\
`# Simulate availability windows: depot open all day, customers available`\
`# within a 30-minute window starting at staggered times`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`window_open`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, `[`runif`](https://rdrr.io/r/stats/Uniform.html)`(``n`` ``-`` ``1``, ``0``, ``40``)``)``   ``# depot at 0`\
`window_close`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``200``, ``window_open``[``-``1``]`` ``+`` ``30``)`` ``# 30-minute windows`\
\
`result_tw`` ``<-`` `[`route_vrp`](https://walker-data.com/spopt-r/reference/route_vrp.md)`(`\
`  ``stops``,`\
`  depot ``=`` ``1``,`\
`  demand_col ``=`` ``"packages"``,`\
`  vehicle_capacity ``=`` ``50``,`\
`  cost_matrix ``=`` ``ttm``,`\
`  service_time ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``3``, ``n``)``,`\
`  earliest ``=`` ``window_open``,`\
`  latest ``=`` ``window_close`\
`)`\
\
[`summary`](https://rspatial.github.io/terra/reference/summary.html)`(``result_tw``)`

    VRP routes: 26 locations, 6 vehicles (depot: 1)
      Method: 2-opt | Total cost: 229.0 | Improvement: 6.2%
      Capacity: 50 | Time windows: active
      Solve time: 0.000s

    Per-vehicle summary:
     Vehicle Stops Load Cost  Time
           1     4   17   42 54.00
           2     5   19   33 68.28
           3     4   16   25 53.76
           4     4   16   41 68.13
           5     3   10   41 59.16
           6     5   20   47 74.22

      Cost = matrix objective (travel only); Time = travel + service + waiting

Arrival and departure times are added to the output as `.arrival_time`
and `.departure_time` columns. Time windows interact with all other
constraints: capacity limits, `max_route_time` (which includes waiting
time), and fleet size. Note that `balance = "time"` is disabled when
time windows are active.

## Generating route geometries for mapping

The maps above use pre-bundled road geometries. Routing engines return
travel-time matrices and driving geometries through separate API calls:
the matrix gives you the numbers for optimization, and a
directions/itinerary call gives you the road geometry for visualization.
This two-step workflow is standard across all routing engines (r5r,
OSRM, mapboxapi, Google).

To generate your own route lines: extract the stop sequence from the
solver output, then request driving directions for each consecutive pair
of stops.

Here’s how to do it with r5r:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`r5r`](https://github.com/ipeaGIT/r5r)`)`\
\
`# Prepare stop coordinates for r5r`\
`r5r_pts`` ``<-`` ``stops`` ``|>`\
`  `[`st_coordinates`](https://r-spatial.github.io/sf/reference/st_coordinates.html)`(``)`` ``|>`\
`  `[`as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)`(``)`` ``|>`\
`  `[`rename`](https://dplyr.tidyverse.org/reference/rename.html)`(``lon ``=`` ``X``, lat ``=`` ``Y``)`` ``|>`\
`  `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``id ``=`` ``stops``$``id``)`\
\
`# Get road geometries for a sequence of stops`\
`get_route_legs`` ``<-`` ``function``(``tour_indices``, ``r5r_core``)`` ``{`\
`  `[`map`](https://purrr.tidyverse.org/reference/map.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`length`](https://rdrr.io/r/base/length.html)`(``tour_indices``)`` ``-`` ``1``)``, \``(``i``)`` ``{`\
`    `[`detailed_itineraries`](https://ipeagit.github.io/r5r/reference/detailed_itineraries.html)`(`\
`      ``r5r_core``,`\
`      origins ``=`` ``r5r_pts`` ``|>`` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``id`` ``==`` ``stops``$``id``[``tour_indices``[``i``]``]``)``,`\
`      destinations ``=`` ``r5r_pts`` ``|>`` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``id`` ``==`` ``stops``$``id``[``tour_indices``[``i`` ``+`` ``1``]``]``)``,`\
`      mode ``=`` ``"CAR"``,`\
`      departure_datetime ``=`` `[`as.POSIXct`](https://rdrr.io/r/base/as.POSIXlt.html)`(``"2025-03-15 08:00:00"``)``,`\
`      shortest_path ``=`` ``TRUE`\
`    ``)`\
`  ``}``)`` ``|>`\
`    `[`bind_rows`](https://dplyr.tidyverse.org/reference/bind_rows.html)`(``)`\
`}`\
\
`# TSP: use the tour from metadata`\
`tsp_meta`` ``<-`` `[`attr`](https://rdrr.io/r/base/attr.html)`(``result``, ``"spopt"``)`\
`tsp_route`` ``<-`` ``get_route_legs``(``tsp_meta``$``tour``, ``r5r_core``)`\
\
`# VRP: build each vehicle's tour (depot -> stops -> depot)`\
`vrp_route`` ``<-`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``meta_vrp``$``n_vehicles``)`` ``|>`\
`  `[`map`](https://purrr.tidyverse.org/reference/map.html)`(``\``(``v``)`` ``{`\
`    ``vehicle_stops`` ``<-`` ``result_vrp`` ``|>`\
`      `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``.vehicle`` ``==`` ``v``)`` ``|>`\
`      `[`arrange`](https://dplyr.tidyverse.org/reference/arrange.html)`(``.visit_order``)`\
\
`    ``vehicle_tour`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, `[`match`](https://rspatial.github.io/terra/reference/match.html)`(``vehicle_stops``$``id``, ``stops``$``id``)``, ``1``)`\
`    ``get_route_legs``(``vehicle_tour``, ``r5r_core``)`` ``|>`\
`      `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``vehicle ``=`` ``v``)`\
`  ``}``)`` ``|>`\
`  `[`bind_rows`](https://dplyr.tidyverse.org/reference/bind_rows.html)`(``)`

The same pattern works with other routing engines. With OSRM
([`osrm::osrmRoute()`](https://riatelab.r-universe.dev/osrm/reference/osrmRoute.html)),
you’d request each leg as a pair of coordinates. With mapboxapi
(`mb_directions()`), you’d pass origin/destination pairs. Pay attention
to the stop sequence: spopt gives you the optimal order, and the routing
engine gives you the road geometry to draw.

## TSP vs VRP: when to use which

| Scenario | Function | Key parameters |
|----|----|----|
| One driver, return to base | `route_tsp(start = 1)` | Default closed tour |
| One driver, end anywhere | `route_tsp(start = 1, end = NULL)` | Open route |
| One driver, fixed endpoint | `route_tsp(start = 1, end = 5)` | Path between two points |
| Multiple drivers, capacity limits | `route_vrp(demand_col, vehicle_capacity)` | Fleet optimization |
| Fixed fleet size | `route_vrp(n_vehicles = 3)` | Constrained fleet |
| Shift time limits | `route_vrp(max_route_time = 90)` | Per-route duration cap |
| Customer availability | `route_vrp(earliest, latest)` | Time window constraints |

For most delivery and field service problems, start with
[`route_tsp()`](https://walker-data.com/spopt-r/reference/route_tsp.md)
to understand the single-vehicle baseline, then move to
[`route_vrp()`](https://walker-data.com/spopt-r/reference/route_vrp.md)
when you need to split work across a fleet.

## Next steps

- [Travel-Time Cost
  Matrices](https://walker-data.com/spopt-r/articles/travel-time-matrices.md) -
  Generate travel-time matrices with r5r
- [Facility
  Location](https://walker-data.com/spopt-r/articles/facility-location.md) -
  Optimize depot and warehouse placement
- [Getting
  Started](https://walker-data.com/spopt-r/articles/getting-started.md) -
  Overview of all spopt capabilities
