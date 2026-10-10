# Travel-time cost matrices

By default, spopt’s facility location functions calculate Euclidean
(straight-line) distances between points. While this is fast and often
sufficient for initial analysis, real-world accessibility depends on
road networks, traffic patterns, and travel mode. A location that’s
“close” as the crow flies might be far away by car if separated by a
river, highway, or geographic barrier.

This vignette demonstrates how to generate travel-time matrices using
external routing engines and integrate them with spopt for more
realistic facility location analysis.

## Why travel times matter

Consider placing emergency medical services in a metropolitan area.
Straight-line distance might suggest one configuration, but accounting
for actual drive times - influenced by road networks, speed limits, and
traffic - could yield very different optimal locations.

The difference can be substantial. In areas with limited river
crossings, highway barriers, or mountainous terrain, Euclidean distance
can significantly underestimate actual travel costs. For high-stakes
decisions like emergency service placement or critical infrastructure
siting, using real travel times is essential.

## Routing options in R

Several R packages can generate travel-time matrices, including but not
limited to:

- **r5r**: Fast routing using R5 engine; supports driving, walking,
  cycling, and public transit; handles large-scale analyses
- **dodgr**: Street network routing using OpenStreetMap data; good for
  walking and cycling
- **osrm**: Interface to OSRM routing engine; fast car routing with
  traffic-like weights
- **mapboxapi**: Interface to Mapbox’s routing services; includes
  traffic-aware routing
- **googleway**: Interface to Google’s Directions API

In this vignette, we’ll use r5r to generate a travel-time matrix, then
demonstrate how to use it with spopt’s facility location algorithms.

## Generating a travel-time matrix with r5r

r5r requires Java 21 and OpenStreetMap data. Here’s the workflow used to
generate a travel-time matrix for Tarrant County, Texas:

\
`# Install r5r and set up Java`\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"r5r"``)`\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"rJavaEnv"``)`\
`rJavaEnv``::`[`java_quick_install`](https://www.ekotov.pro/rJavaEnv/reference/java_quick_install.html)`(``version ``=`` ``21``)`\
\
`# Set Java memory and load libraries`\
[`options`](https://rdrr.io/r/base/options.html)`(``java.parameters ``=`` ``"-Xmx8G"``)`` ``# Allocating 8GB RAM`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`r5r`](https://github.com/ipeaGIT/r5r)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sf`](https://r-spatial.github.io/sf/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyverse`](https://tidyverse.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidycensus`](https://walker-data.com/tidycensus/)`)`\
\
`# Download and clip OSM data to study area`\
`# First, fetch OSM data from https://download.geofabrik.de/north-america/us-south-latest.osm.pbf (3.7GB)`\
`# Then, use osmium to clip the US South extract to Tarrant County:`\
`# osmium extract -b -97.55,32.55,-97.0,33.05 us-south.osm.pbf -o tarrant.osm.pbf`\
\
`# Build routing network`\
`data_path`` ``<-`` ``"path/to/osm/directory"`\
`r5r_core`` ``<-`` `[`build_network`](https://ipeagit.github.io/r5r/reference/build_network.html)`(``data_path ``=`` ``data_path``)`\
\
`# Get tract centroids as demand points`\
`tarrant_tracts`` ``<-`` `[`get_acs`](https://walker-data.com/tidycensus/reference/get_acs.html)`(`\
` geography ``=`` ``"tract"``,`\
` variables ``=`` ``"B01003_001"``,`\
` state ``=`` ``"TX"``,`\
` county ``=`` ``"Tarrant"``,`\
` geometry ``=`` ``TRUE``,`\
` year ``=`` ``2023`\
`)`` ``|>`\
` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``estimate`` ``>`` ``0``)`` ``|>`\
` `[`rename`](https://dplyr.tidyverse.org/reference/rename.html)`(``population ``=`` ``estimate``)`\
\
`demand_pts`` ``<-`` ``tarrant_tracts`` ``|>`\
` `[`st_centroid`](https://r-spatial.github.io/sf/reference/geos_unary.html)`(``)`` ``|>`\
` `[`st_transform`](https://r-spatial.github.io/sf/reference/st_transform.html)`(``4326``)`` ``|>`\
` `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``id ``=`` `[`row_number`](https://dplyr.tidyverse.org/reference/row_number.html)`(``)``)`\
\
`# Sample 30 candidate facility locations`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1983``)`\
`candidate_pts`` ``<-`` `[`st_sample`](https://r-spatial.github.io/sf/reference/st_sample.html)`(`[`st_union`](https://r-spatial.github.io/sf/reference/geos_combine.html)`(``tarrant_tracts``)``, ``30``)`` ``|>`\
` `[`st_as_sf`](https://r-spatial.github.io/sf/reference/st_as_sf.html)`(``)`` ``|>`\
` `[`st_transform`](https://r-spatial.github.io/sf/reference/st_transform.html)`(``4326``)`` ``|>`\
` `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``id ``=`` `[`row_number`](https://dplyr.tidyverse.org/reference/row_number.html)`(``)``)`\
\
`# Prepare points for r5r (requires id, lon, lat columns)`\
`demand_r5r`` ``<-`` ``demand_pts`` ``|>`\
` `[`st_coordinates`](https://r-spatial.github.io/sf/reference/st_coordinates.html)`(``)`` ``|>`\
` `[`as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)`(``)`` ``|>`\
` `[`rename`](https://dplyr.tidyverse.org/reference/rename.html)`(``lon ``=`` ``X``, lat ``=`` ``Y``)`` ``|>`\
` `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``id ``=`` `[`as.character`](https://rdrr.io/r/base/character.html)`(``demand_pts``$``id``)``)`\
\
`candidates_r5r`` ``<-`` ``candidate_pts`` ``|>`\
` `[`st_coordinates`](https://r-spatial.github.io/sf/reference/st_coordinates.html)`(``)`` ``|>`\
` `[`as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)`(``)`` ``|>`\
` `[`rename`](https://dplyr.tidyverse.org/reference/rename.html)`(``lon ``=`` ``X``, lat ``=`` ``Y``)`` ``|>`\
` `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``id ``=`` `[`as.character`](https://rdrr.io/r/base/character.html)`(``candidate_pts``$``id``)``)`\
\
`# Calculate travel-time matrix`\
`ttm`` ``<-`` `[`travel_time_matrix`](https://ipeagit.github.io/r5r/reference/travel_time_matrix.html)`(`\
` ``r5r_core``,`\
` origins ``=`` ``demand_r5r``,`\
` destinations ``=`` ``candidates_r5r``,`\
` mode ``=`` ``"CAR"``,`\
` departure_datetime ``=`` `[`as.POSIXct`](https://rdrr.io/r/base/as.POSIXlt.html)`(``"2025-03-15 08:00:00"``)``,`\
` max_trip_duration ``=`` ``120``,`\
` progress ``=`` ``TRUE`\
`)`\
\
`# Reshape to matrix format`\
`ttm_matrix`` ``<-`` ``ttm`` ``|>`\
` `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``from_id``, ``to_id``, ``travel_time_p50``)`` ``|>`\
` `[`pivot_wider`](https://tidyr.tidyverse.org/reference/pivot_wider.html)`(``names_from ``=`` ``to_id``, values_from ``=`` ``travel_time_p50``)`` ``|>`\
` `[`arrange`](https://dplyr.tidyverse.org/reference/arrange.html)`(`[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(``from_id``)``)`` ``|>`\
` `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``-``from_id``)`` ``|>`\
` `[`as.matrix`](https://rspatial.github.io/terra/reference/coerce.html)`(``)`\
\
`ttm_matrix``[`[`is.na`](https://rdrr.io/r/base/NA.html)`(``ttm_matrix``)``]`` ``<-`` ``Inf`\
\
[`stop_r5`](https://ipeagit.github.io/r5r/reference/stop_r5.html)`(``r5r_core``)`

## Using the bundled travel-time data

spopt includes a pre-computed travel-time matrix for Tarrant County that
we’ll use for the examples below. This data was generated using the
workflow above.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`spopt`](https://walker-data.com/spopt-r/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sf`](https://r-spatial.github.io/sf/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyverse`](https://tidyverse.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`mapgl`](https://walker-data.com/mapgl/)`)`\
\
`# Load the bundled data`\
[`data`](https://rdrr.io/r/utils/data.html)`(``tarrant_travel_times``)`\
\
`# Extract components`\
`tracts`` ``<-`` ``tarrant_travel_times``$``tracts`\
`demand`` ``<-`` ``tarrant_travel_times``$``demand`\
`candidates`` ``<-`` ``tarrant_travel_times``$``candidates`\
`ttm`` ``<-`` ``tarrant_travel_times``$``matrix`\
\
`# Check dimensions`\
[`dim`](https://rdrr.io/r/base/dim.html)`(``ttm``)`

    [1] 448  30

The matrix has 448 rows (demand points) and 30 columns (candidate
facilities), with travel times in minutes.

## P-Median: Euclidean vs travel time

Let’s compare P-Median solutions using Euclidean distance versus actual
travel times:

\
`# Solution using travel-time matrix`\
`result_tt`` ``<-`` `[`p_median`](https://walker-data.com/spopt-r/reference/p_median.md)`(`\
` demand ``=`` ``demand``,`\
` facilities ``=`` ``candidates``,`\
` n_facilities ``=`` ``5``,`\
` weight_col ``=`` ``"population"``,`\
` cost_matrix ``=`` ``ttm`\
`)`\
\
`# Solution using Euclidean distance`\
`result_euc`` ``<-`` `[`p_median`](https://walker-data.com/spopt-r/reference/p_median.md)`(`\
` demand ``=`` ``demand``,`\
` facilities ``=`` ``candidates``,`\
` n_facilities ``=`` ``5``,`\
` weight_col ``=`` ``"population"`\
`)`\
\
`# Get selected facility IDs from each solution`\
`selected_tt_ids`` ``<-`` ``result_tt``$``facilities`` ``|>`\
` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``.selected``)`` ``|>`\
` `[`pull`](https://dplyr.tidyverse.org/reference/pull.html)`(``id``)`\
\
`selected_euc_ids`` ``<-`` ``result_euc``$``facilities`` ``|>`\
` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``.selected``)`` ``|>`\
` `[`pull`](https://dplyr.tidyverse.org/reference/pull.html)`(``id``)`\
\
`# Categorize facilities`\
`candidates_compared`` ``<-`` ``candidates`` ``|>`\
` `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(`\
`   selected_tt ``=`` ``id`` `[`%in%`](https://rspatial.github.io/terra/reference/match.html)` ``selected_tt_ids``,`\
`   selected_euc ``=`` ``id`` `[`%in%`](https://rspatial.github.io/terra/reference/match.html)` ``selected_euc_ids``,`\
`   category ``=`` `[`case_when`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html)`(`\
`     ``selected_tt`` ``&`` ``selected_euc`` ``~`` ``"Both methods"``,`\
`     ``selected_tt`` ``~`` ``"Travel time only"``,`\
`     ``selected_euc`` ``~`` ``"Euclidean only"``,`\
`     ``TRUE`` ``~`` ``"Not selected"`\
`   ``)`\
` ``)`\
\
`# Count by category`\
`candidates_compared`` ``|>`\
` `[`st_drop_geometry`](https://r-spatial.github.io/sf/reference/st_geometry.html)`(``)`` ``|>`\
` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``category`` ``!=`` ``"Not selected"``)`` ``|>`\
` `[`count`](https://dplyr.tidyverse.org/reference/count.html)`(``category``)`

              category n
    1     Both methods 1
    2   Euclidean only 4
    3 Travel time only 4

The two methods select quite different facility locations. Let’s
visualize the comparison:

\
`# Filter to only selected facilities`\
`selected_facilities`` ``<-`` ``candidates_compared`` ``|>`\
` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``category`` ``!=`` ``"Not selected"``)`\
\
[`maplibre`](https://walker-data.com/mapgl/reference/maplibre.html)`(``bounds ``=`` ``tracts``)`` ``|>`\
` `[`add_fill_layer`](https://walker-data.com/mapgl/reference/add_fill_layer.html)`(`\
`   id ``=`` ``"tracts"``,`\
`   source ``=`` ``tracts``,`\
`   fill_color ``=`` ``"lightgray"``,`\
`   fill_opacity ``=`` ``0.3`\
` ``)`` ``|>`\
` `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`   id ``=`` ``"facilities"``,`\
`   source ``=`` ``selected_facilities``,`\
`   circle_radius ``=`` ``10``,`\
`   circle_color ``=`` `[`match_expr`](https://walker-data.com/mapgl/reference/match_expr.html)`(`\
`     column ``=`` ``"category"``,`\
`     values ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Both methods"``, ``"Travel time only"``, ``"Euclidean only"``)``,`\
`     stops ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#9b59b6"``, ``"#e74c3c"``, ``"#3498db"``)`\
`   ``)``,`\
`   circle_stroke_color ``=`` ``"white"``,`\
`   circle_stroke_width ``=`` ``2`\
` ``)`

Purple markers indicate facilities selected by both methods, red markers
are selected only when using travel times, and blue markers are selected
only with Euclidean distance.

## Why the solutions differ

The travel-time solution accounts for the actual road network.
Facilities shift toward locations with better highway access, even if
they’re farther in straight-line distance. Let’s examine the objective
values:

\
`# Compare objective values`\
`tt_obj`` ``<-`` `[`attr`](https://rdrr.io/r/base/attr.html)`(``result_tt``, ``"spopt"``)``$``objective`\
`euc_obj`` ``<-`` `[`attr`](https://rdrr.io/r/base/attr.html)`(``result_euc``, ``"spopt"``)``$``objective`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Travel-time solution objective:"``, `[`round`](https://rspatial.github.io/terra/reference/math-generics.html)`(``tt_obj``, ``0``)``, ``"person-minutes\n"``)`

    Travel-time solution objective: 29305339 person-minutes

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Euclidean solution objective:"``, `[`round`](https://rspatial.github.io/terra/reference/math-generics.html)`(``euc_obj``, ``0``)``, ``"person-meters\n"``)`

    Euclidean solution objective: 16669528322 person-meters

Note that the objectives use different units (minutes vs meters), so
they’re not directly comparable. What matters is that each solution is
optimal for its respective cost measure.

## Visualizing service areas

We can also visualize how demand points are assigned to facilities under
each solution:

\
`# Get demand assignments from travel-time solution`\
`demand_tt`` ``<-`` ``result_tt``$``demand`` ``|>`\
` `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``facility_id ``=`` `[`as.character`](https://rdrr.io/r/base/character.html)`(``.facility``)``)`\
\
`# Get the selected facilities`\
`facilities_tt`` ``<-`` ``result_tt``$``facilities`` ``|>`\
` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``.selected``)`` ``|>`\
` `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``facility_id ``=`` `[`as.character`](https://rdrr.io/r/base/character.html)`(``id``)``)`\
\
[`maplibre`](https://walker-data.com/mapgl/reference/maplibre.html)`(``bounds ``=`` ``tracts``)`` ``|>`\
` `[`add_fill_layer`](https://walker-data.com/mapgl/reference/add_fill_layer.html)`(`\
`   id ``=`` ``"tracts"``,`\
`   source ``=`` ``tracts``,`\
`   fill_color ``=`` ``"lightgray"``,`\
`   fill_opacity ``=`` ``0.2`\
` ``)`` ``|>`\
` `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`   id ``=`` ``"demand"``,`\
`   source ``=`` ``demand_tt``,`\
`   circle_radius ``=`` ``4``,`\
`   circle_opacity ``=`` ``0.7``,`\
`   circle_color ``=`` `[`match_expr`](https://walker-data.com/mapgl/reference/match_expr.html)`(`\
`     column ``=`` ``"facility_id"``,`\
`     values ``=`` ``facilities_tt``$``facility_id``,`\
`     stops ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#e41a1c"``, ``"#377eb8"``, ``"#4daf4a"``, ``"#984ea3"``, ``"#ff7f00"``)`\
`   ``)`\
` ``)`` ``|>`\
` `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`   id ``=`` ``"facilities"``,`\
`   source ``=`` ``facilities_tt``,`\
`   circle_radius ``=`` ``12``,`\
`   circle_color ``=`` `[`match_expr`](https://walker-data.com/mapgl/reference/match_expr.html)`(`\
`     column ``=`` ``"facility_id"``,`\
`     values ``=`` ``facilities_tt``$``facility_id``,`\
`     stops ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#e41a1c"``, ``"#377eb8"``, ``"#4daf4a"``, ``"#984ea3"``, ``"#ff7f00"``)`\
`   ``)``,`\
`   circle_stroke_color ``=`` ``"white"``,`\
`   circle_stroke_width ``=`` ``3`\
` ``)`

Each demand point is colored by its assigned facility. Notice how the
service areas follow the road network structure rather than forming
simple circular regions; look in particular along major highways.

### Spider lines with travel times

[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)
draws a line from each demand point to its assigned facility. When the
solver used a travel-time matrix, each demand point’s `.cost` column
holds the drive time to its assigned facility.
[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)
picks this up automatically as a `cost` column, so the lines can be
styled by drive time:

\
`tt_lines`` ``<-`` `[`spider_lines`](https://walker-data.com/spopt-r/reference/spider_lines.md)`(``result_tt``)`\
\
[`summary`](https://rspatial.github.io/terra/reference/summary.html)`(``tt_lines``$``cost``)`

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
       2.00   10.00   13.00   13.62   17.00   36.00 

\
[`maplibre`](https://walker-data.com/mapgl/reference/maplibre.html)`(``bounds ``=`` ``tracts``)`` ``|>`\
`  `[`add_fill_layer`](https://walker-data.com/mapgl/reference/add_fill_layer.html)`(`\
`    id ``=`` ``"tracts"``,`\
`    source ``=`` ``tracts``,`\
`    fill_color ``=`` ``"lightgray"``,`\
`    fill_opacity ``=`` ``0.2`\
`  ``)`` ``|>`\
`  `[`add_line_layer`](https://walker-data.com/mapgl/reference/add_line_layer.html)`(`\
`    id ``=`` ``"allocations"``,`\
`    source ``=`` ``tt_lines``,`\
`    line_color ``=`` ``interpolate``(`\
`      column ``=`` ``"cost"``,`\
`      values ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``5``, ``15``, ``30``)``,`\
`      stops ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#fee08b"``, ``"#f46d43"``, ``"#a50026"``)`\
`    ``)``,`\
`    line_width ``=`` ``1.5`\
`  ``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"facilities"``,`\
`    source ``=`` ``facilities_tt``,`\
`    circle_radius ``=`` ``10``,`\
`    circle_color ``=`` ``"black"``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``2`\
`  ``)`` ``|>`\
`  ``add_legend``(`\
`    ``"Drive time (minutes)"``,`\
`    values ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``5``, ``15``, ``30``)``,`\
`    colors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#fee08b"``, ``"#f46d43"``, ``"#a50026"``)`\
`  ``)`

The lines are straight, but their color reflects drive time over the
road network, so the longest trips stand out.

To draw lines that follow the roads themselves, pass a routing function
to the `route_fun` argument of
[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md).
[`r5r_route_fun()`](https://walker-data.com/spopt-r/reference/r5r_route_fun.md)
builds one from the same r5r network used to compute the travel-time
matrix above:

\
`routed_lines`` ``<-`` `[`spider_lines`](https://walker-data.com/spopt-r/reference/spider_lines.md)`(`\
`  ``result_tt``,`\
`  route_fun ``=`` `[`r5r_route_fun`](https://walker-data.com/spopt-r/reference/r5r_route_fun.md)`(`\
`    ``r5r_core``,`\
`    mode ``=`` ``"CAR"``,`\
`    departure_datetime ``=`` `[`as.POSIXct`](https://rdrr.io/r/base/as.POSIXlt.html)`(``"2026-10-07 08:00:00"``)`\
`  ``)`\
`)`

Each spider line then traces the driving route between a tract and its
assigned facility. r5r routes every pair in one batched call. Any pair
it can’t route falls back to a straight line, flagged by
`routed == FALSE`. For other routers such as OSRM, see the example in
[`?spider_lines`](https://walker-data.com/spopt-r/reference/spider_lines.md).

## Best practices

When working with travel-time matrices:

- **Match row/column order**: Ensure your cost matrix rows correspond to
  demand points and columns to facilities in the same order as your sf
  objects.

- **Handle unreachable pairs**: Set unreachable origin-destination pairs
  to `Inf` rather than `NA`.

- **Cache your matrices**: Travel-time computation is expensive. Save
  matrices with
  [`readr::write_rds()`](https://readr.tidyverse.org/reference/read_rds.html)
  for reuse.

- **Optimize your source data**: To speed up your matrix calculations,
  clip your OSM source data as close as possible to your area of
  analysis with `osmium` before building your routing network. If using
  a third-party API, pay close attention to costs. For many-to-many
  matrices required for solving optimization problems, Google / Mapbox
  APIs can get expensive quickly.

## Next steps

- [Facility
  Location](https://walker-data.com/spopt-r/articles/facility-location.md) -
  Overview of all location algorithms
- [Huff Model](https://walker-data.com/spopt-r/articles/huff-model.md) -
  Travel times also work with the Huff model
- [Regionalization](https://walker-data.com/spopt-r/articles/regionalization.md) -
  Build custom regions
