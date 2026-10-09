# Road-following spider lines with r5r

Creates a routing function for the `route_fun` argument of
[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md),
using
[`r5r::detailed_itineraries()`](https://ipeagit.github.io/r5r/reference/detailed_itineraries.html)
so that spider lines follow the street network instead of running
straight. Use the same r5r network you used to build your travel-time
matrix.

## Usage

``` r
r5r_route_fun(r5r_network, mode = "CAR", departure_datetime = Sys.time(), ...)
```

## Arguments

- r5r_network:

  An r5r network built with
  [`r5r::build_network()`](https://ipeagit.github.io/r5r/reference/build_network.html).

- mode:

  Transport mode passed to
  [`r5r::detailed_itineraries()`](https://ipeagit.github.io/r5r/reference/detailed_itineraries.html).
  Defaults to `"CAR"`.

- departure_datetime:

  Departure time passed to
  [`r5r::detailed_itineraries()`](https://ipeagit.github.io/r5r/reference/detailed_itineraries.html).
  Defaults to the current time.

- ...:

  Further arguments passed to
  [`r5r::detailed_itineraries()`](https://ipeagit.github.io/r5r/reference/detailed_itineraries.html),
  such as `max_trip_duration`.

## Value

A function taking `(from, to)` that returns one MULTILINESTRING per
origin-destination pair (EPSG:4326), with empty geometries for pairs r5r
could not route. Pass it to
[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)
as `route_fun`.

## Details

r5r returns each itinerary as one or more legs. The returned function
keeps the first itinerary option for each pair, orders its legs, and
combines them into a single geometry. Every pair is routed in one
batched call, and
[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)
sends each distinct pair only once.

## See also

[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(r5r)
options(java.parameters = "-Xmx4G")

net <- build_network(data_path = "path/to/osm")
result <- p_median(demand, facilities, n_facilities = 5,
                   weight_col = "population", cost_matrix = ttm)

lines <- spider_lines(
  result,
  route_fun = r5r_route_fun(net, mode = "CAR",
                            departure_datetime = as.POSIXct("2026-10-07 08:00"))
)
} # }
```
