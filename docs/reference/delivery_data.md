# Delivery routing data for Fort Worth, TX

Pre-computed delivery stop locations and travel-time matrix for a route
optimization scenario in western Tarrant and Parker counties, Texas. The
depot is an Amazon delivery station on San Jacinto Dr in Fort Worth,
with 25 residential delivery stops in Fort Worth, Benbrook, Aledo, and
Willow Park.

## Usage

``` r
delivery_data
```

## Format

A list with four elements:

- stops:

  An sf object with 26 rows (1 depot + 25 delivery stops) and columns:

  id

  : Stop identifier ("depot", "D01", ..., "D25")

  address

  : Street address

  packages

  : Number of packages for delivery (0 for depot)

- matrix:

  A 26x26 numeric matrix of driving travel times in minutes, computed
  with r5r using OpenStreetMap road network data. Row and column names
  correspond to stop IDs.

- tsp_route:

  An sf LINESTRING object with road-snapped route geometries for the
  optimized single-vehicle tour, one row per leg.

- vrp_route:

  An sf LINESTRING object with road-snapped route geometries for the
  optimized multi-vehicle solution, with a `vehicle` column indicating
  which van each leg belongs to.

## Source

Travel times computed with r5r v2.3 using OpenStreetMap data for Tarrant
and Parker counties, Texas. Addresses geocoded with Nominatim via
tidygeocoder.
