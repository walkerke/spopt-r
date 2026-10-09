# Corridor routing with spopt

The [route optimization
vignette](https://walker-data.com/spopt-r/articles/routing.md) covers
problems where you’re sequencing stops along an existing road network.
*Corridor routing* is a different class of problem. There is no road
network to follow. Instead, you’re finding the best path across raw
terrain. The key input is a *cost surface*: a raster grid where each
cell’s value represents how expensive it is to cross that piece of land.
The solver finds the path from origin to destination that minimizes
total accumulated cost.

This kind of analysis shows up in pipeline routing, transmission line
siting, trail planning, and wildlife connectivity modeling. The analyst
builds a cost surface that reflects the relevant obstacles and
preferences, then the algorithm finds the cheapest way through.

This vignette walks through a realistic pipeline routing scenario in
rural West Virginia. We’ll build a cost surface from publicly available
terrain and land cover data, find the optimal corridor between two real
infrastructure points, generate diverse alternative corridors, and route
a multi-origin gathering system using a cached graph. All of the well
and disposal facility locations used here are real, sourced from the
[West Virginia DEP Office of Oil and Gas](https://tagis.dep.wv.gov/oog/)
well database.

## The problem

Operators developing Marcellus Shale wells in West Virginia produce
significant volumes of water as a byproduct of natural gas extraction.
This produced water needs to be transported to saltwater disposal (SWD)
facilities, typically via gathering pipelines. As development expands,
new pipeline routes must be planned across hilly, heavily forested
Appalachian terrain, balancing construction cost, environmental
constraints, and the practicalities of right-of-way acquisition from
landowners.

In this example, we’ll route a produced water pipeline from the Antero
Resources Farley Unit pad in Doddridge County to the Ritchie Hunter 2
saltwater disposal well in Ritchie County, a distance of roughly 25 km
across terrain that includes forested ridges, river valleys, and
scattered rural development.

## Building a cost surface

The cost surface encodes our knowledge about the landscape into a form
the routing algorithm can use. Each cell holds a *friction* value: the
cost per unit distance to cross that cell. Higher friction means more
expensive terrain. Cells set to `NA` are treated as impassable barriers.

The resulting friction values are unitless. They represent *relative*
construction difficulty, not dollar costs. A cell with friction 4.0
costs four times as much to cross as a cell with friction 1.0,
regardless of what either costs in absolute terms. The algorithm
optimizes over relative friction, and translating the results to dollar
estimates requires a separate post-processing step with construction
cost benchmarks.

The weights used below are illustrative. In a real siting study, these
would be calibrated against actual construction cost data, regulatory
requirements, and stakeholder input. The general structure, however, is
representative of how pipeline routing cost surfaces are built in
practice: combine terrain difficulty with land cover constraints, and
exclude areas where construction is infeasible.

Our cost surface combines two publicly available layers: terrain slope
derived from a USGS digital elevation model, and land cover from the
National Land Cover Database (NLCD).

``` r

library(spopt)
library(terra)
library(sf)
library(tidyverse)
library(mapgl)
```

``` r

dem  <- rast("/vsicurl/https://walker-data.com/maps/data/wv_dem.tif")
nlcd <- rast("/vsicurl/https://walker-data.com/maps/data/wv_nlcd.tif")

map_bounds  <- c(-81.25, 39.05, -80.65, 39.45)
```

``` r

maplibre(
  style = openfreemap_style("positron"),
  bounds = map_bounds
) |>
  add_image_source(
    id = "dem", data = dem,
    colors = hcl.colors(100, "Earth")
  ) |>
  add_raster_layer(
    id = "dem-layer", source = "dem",
    raster_opacity = 0.8, raster_resampling = "nearest"
  )
maplibre(
  style = openfreemap_style("positron"),
  bounds = map_bounds
) |>
  add_image_source(
    id = "nlcd", data = nlcd
  ) |>
  add_raster_layer(
    id = "nlcd-layer", source = "nlcd",
    raster_opacity = 0.8, raster_resampling = "nearest"
  )
```

Elevation (DEM)

Land Cover (NLCD)

The friction model has two components. The first is a slope-based cost:
trenching and pipe-laying become progressively more difficult on steeper
terrain. We’ll first compute a slope raster with
[`terrain()`](https://rspatial.github.io/terra/reference/terrain.html),
then assign a base friction of 1.0 for flat ground, increasing by 0.3
per degree of slope.

``` r

slope <- terrain(dem, v = "slope", unit = "degrees")

friction <- 1.0 + slope * 0.3
```

The second component is a set of land cover multipliers. Different NLCD
classes carry different construction and permitting costs. Open
pastureland is the cheapest to cross. Forested land requires clearing
and restoration. Developed areas involve expensive right-of-way
negotiations. Water bodies and wetlands are treated as impassable due to
permitting constraints under Section 404 of the Clean Water Act.

``` r

# Land cover weight table
#
# NLCD Class              Multiplier  Rationale
# ─────────────────────── ────────── ──────────────────────────────────
# 81  Pasture/hay          1.0x       Easiest: flat, minimal clearing
# 82  Cultivated crops     1.2x       Crop damage compensation
# 52  Shrub/scrub          1.3x       Light clearing required
# 21  Developed, open      1.5x       Lawns, parks — minor ROW cost
# 41  Deciduous forest     2.0x       Timber clearing + restoration
# 42  Evergreen forest     2.0x       Same
# 43  Mixed forest         2.0x       Same
# 22  Developed, low       3.0x       Residential — ROW negotiation
# 23  Developed, medium    5.0x       Commercial/suburban — expensive ROW
# 24  Developed, high      10.0x      Urban core — avoid if possible
# 11  Open water           NA         Impassable (directional drill needed)
# 90  Woody wetlands       NA         Impassable (Section 404 permitting)
# 95  Herbaceous wetlands  NA         Impassable (Section 404 permitting)

lc_weights <- c(
  "81" = 1.0, "82" = 1.2, "52" = 1.3, "21" = 1.5,
  "41" = 2.0, "42" = 2.0, "43" = 2.0,
  "22" = 3.0, "23" = 5.0, "24" = 10.0
)
lc_exclude <- c(11, 90, 95)

nlcd_vals <- values(nlcd, mat = FALSE)
lc_mult <- classify(
  nlcd,
  rcl = cbind(as.integer(names(lc_weights)), unname(lc_weights)),
  others = 1.0
)

cost_surface <- friction * lc_mult
cost_surface[nlcd_vals %in% lc_exclude] <- NA
```

The resulting surface reflects the terrain of this part of West
Virginia. Valleys and pastureland have low friction and are cheap to
build through. Forested ridges have moderate friction due to clearing
costs. Water features and wetlands are impassable and must be routed
around.

## Defining the endpoints

The origin and destination for our routing problem are real
infrastructure locations from the [WV DEP Office of Oil and Gas well
database](https://tagis.dep.wv.gov/oog/), queried in March 2026. The
origin is the Antero Resources Farley Unit pad, a pair of Marcellus
horizontal wells completed in May 2025 in Doddridge County (API:
4701706970, 4701706971). The destination is the Ritchie Hunter 2
saltwater disposal well, an active brine disposal facility in Ritchie
County (API: 4708510142). Coordinates were taken from the latitude and
longitude fields in the well database.

``` r

farley_ll  <- c(-80.832905, 39.178088)
swd_ll     <- c(-81.096363, 39.256618)

# Corridor routing requires projected coordinates matching the cost surface CRS
to_utm <- function(ll) {
  st_sfc(st_point(ll), crs = 4326) |>
    st_transform(crs(cost_surface)) |>
    st_coordinates() |>
    as.numeric()
}

farley  <- to_utm(farley_ll)
swd <- to_utm(swd_ll)
```

``` r

# Mapping helpers (reused across sections)
sites <- tibble(
  label = c("Farley Unit", "Ritchie Hunter SWD"),
  type  = c("well_pad", "disposal")
) |>
  st_as_sf(geometry = st_sfc(st_point(farley_ll), st_point(swd_ll), crs = 4326))

cost_colors <- hcl.colors(100, "YlOrRd", rev = TRUE)
```

## Finding the least-cost corridor

With the cost surface and endpoints defined,
[`route_corridor()`](https://walker-data.com/spopt-r/reference/route_corridor.md)
finds the path that minimizes total accumulated friction.

``` r

path <- route_corridor(cost_surface, farley, swd)
path
```

    Least-cost corridor
      Method: dijkstra
      Total cost: 101739.5
      Path distance: 31230
      Cells traversed: 883
      Sinuosity: 1.282
      Solve time: 0.251 s

The returned object is an sf LINESTRING with several useful columns.
`total_cost` is the accumulated friction along the path, a relative
measure that can be compared across corridors. `path_dist` is the actual
path length in CRS units (meters). `sinuosity` is the ratio of path
length to straight-line distance; a value of 1.0 would mean a perfectly
straight route, and values above 1.0 indicate how much the path deviates
to find easier terrain.

``` r

maplibre(style = openfreemap_style("positron"), bounds = map_bounds) |>
  add_image_source(id = "cost", data = cost_surface, colors = cost_colors) |>
  add_raster_layer(
    id = "cost-layer", source = "cost",
    raster_opacity = 0.7, raster_resampling = "nearest"
  ) |>
  add_line_layer(
    id = "corridor", source = st_transform(path, 4326),
    line_color = "#2563eb", line_width = 3, line_opacity = 0.9
  ) |>
  add_circle_layer(
    id = "sites", source = sites,
    circle_color = match_expr("type",
      values = c("well_pad", "disposal"),
      stops = c("#1e1e1e", "#dc2626")
    ),
    circle_radius = 8, circle_stroke_color = "white",
    circle_stroke_width = 2, tooltip = "label"
  )
```

The corridor threads through lower-cost terrain, favoring valleys and
pastureland while routing around forested ridges and water features. The
curved shape of the path is characteristic of raster-based least-cost
path analysis. The algorithm moves between adjacent grid cells, so the
result reflects optimal traversal of a discrete grid. In practice, an
engineering team would use this corridor as a starting zone and lay out
straight pipe segments with planned deflection points within it.

## Generating alternative corridors

The mathematically optimal route is frequently not the one that gets
built. A landowner along the corridor may refuse to grant an easement. A
field survey may discover a wetland that wasn’t captured in the NLCD
data. A road crossing may turn out to be more expensive than the model
assumed. For these reasons, planners need to evaluate multiple
geographically distinct alternatives.

[`route_k_corridors()`](https://walker-data.com/spopt-r/reference/route_k_corridors.md)
generates a set of spatially diverse corridor alternatives using a
heuristic penalty-based approach. After finding the optimal path, it
penalizes the cost surface in a buffer zone around that path, then
solves again. Each successive corridor is pushed into different
geography by the cumulative penalties from prior corridors. The result
is a set of alternatives that use different valleys, different ridge
crossings, and different stretches of terrain. See
[`?route_k_corridors`](https://walker-data.com/spopt-r/reference/route_k_corridors.md)
for references and methodological details.

``` r

alternatives <- route_k_corridors(
  cost_surface, farley, swd,
  k = 5,
  penalty_factor = 2.0
)

alternatives
```

    k-Diverse Corridor Routing (spopt)
      Corridors found: 5 of 5 requested
      Penalty: 2.0x within 1217.8 of each prior path
      Routing time: 1.580s (solve: 1.392s, graph build: 0.189s)

                             Cost    Distance  Sinuosity     Spacing  Overlap
      Optimal             101,740       31230      1.282           -        -
      Alternative 1       115,133       47507      1.951      6559.0     7.0%
      Alternative 2       125,941       36104      1.482      2424.5    17.1%
      Alternative 3       115,480       34647      1.423      2082.7    33.7%
      Alternative 4       147,929       40777      1.674      2379.6    17.1%

A few things are worth noting in the output. All costs are computed on
the *original* (unpenalized) surface, so they are directly comparable.
The cost difference between the optimal corridor and any alternative
reflects the real terrain trade-off. The `mean_spacing` and
`pct_overlap` columns are measured against the optimal corridor,
providing a quantitative measure of how geographically different each
alternative is. The alternatives are not guaranteed to be ordered by
cost; later alternatives sometimes discover cheaper corridors through
terrain that earlier iterations didn’t explore.

``` r

alt_colors <- c("#dc2626", "#f97316", "#eab308", "#22c55e", "#3b82f6")

alternatives_4326 <- alternatives |>
  mutate(
    color = alt_colors[alternative],
    alt_label = if_else(alternative == 1L,
      sprintf("Optimal (cost: %s)", format(round(total_cost), big.mark = ",")),
      sprintf("Alternative %d (cost: %s)", alternative - 1L,
              format(round(total_cost), big.mark = ",")))
  ) |>
  st_transform(4326)

maplibre(style = openfreemap_style("positron"), bounds = map_bounds) |>
  add_image_source(id = "cost", data = cost_surface, colors = cost_colors) |>
  add_raster_layer(id = "cost-layer", source = "cost", raster_opacity = 0.7) |>
  add_line_layer(
    id = "corridors", source = alternatives_4326,
    line_color = get_column("color"),
    line_width = 3, line_opacity = 0.9
  ) |>
  add_circle_layer(
    id = "sites", source = sites,
    circle_color = match_expr("type",
      values = c("well_pad", "disposal"),
      stops = c("#1e1e1e", "#dc2626")
    ),
    circle_radius = 8, circle_stroke_color = "white",
    circle_stroke_width = 2, tooltip = "label"
  ) |>
  add_categorical_legend(
    "Corridors",
    values = alternatives_4326$alt_label,
    colors = alt_colors[seq_len(nrow(alternatives))],
    position = "bottom-left", patch_shape = "line"
  )
```

The corridors fan out across the landscape, each taking advantage of
different terrain features. The penalty-based approach ensures that the
alternatives are not minor detours from the optimal path but instead
represent different routing strategies through this part of the
Appalachian Basin.

## Gathering system routing

A single disposal facility often serves multiple well pads. Antero
operates several Marcellus pad sites spread across Doddridge County, and
each one needs a pipeline connection to the Ritchie Hunter SWD facility.
Rather than rebuilding the routing graph from scratch for each
origin-destination pair,
[`corridor_graph()`](https://walker-data.com/spopt-r/reference/corridor_graph.md)
pre-builds the graph once so that all three pads can be routed
efficiently.

``` r

# Three Antero pad sites spread across Doddridge County
pads <- tribble(
  ~pad_name,  ~lon,        ~lat,
  "Farley",   -80.832905,  39.178088,
  "Holtz",    -80.781658,  39.303559,
  "Michels",  -80.805644,  39.333569
)

pads_utm <- pads |>
  mutate(utm = map2(lon, lat, ~ to_utm(c(.x, .y))))
```

``` r

g <- corridor_graph(cost_surface, neighbours = 8L)
g
```

    Corridor graph
      Grid: 1502 x 1463 (2,197,426 cells)
      Cell size: 29.6 x 29.6
      Neighbours: 8 (17,473,816 edges)
      Build time: 0.033s | Graph storage: ~297.2 MB

The graph object is a snapshot of the cost surface at build time. Once
built, any number of origin-destination pairs can be routed on it
without reprocessing the raster.

``` r

pad_colors <- c("#2563eb", "#16a34a", "#f97316")

pad_paths <- pads_utm |>
  mutate(
    corridor = map(utm, ~ route_corridor(g, .x, swd)),
    color = pad_colors[row_number()]
  )

# Combine corridors into a single sf for mapping
gathering <- pad_paths |>
  mutate(corridor = map2(corridor, pad_name, ~ {
    .x$pad_name <- .y
    class(.x) <- setdiff(class(.x), "spopt_corridor")
    attr(.x, "spopt") <- NULL
    .x
  })) |>
  pull(corridor) |>
  bind_rows() |>
  left_join(tibble(pad_name = pads$pad_name, color = pad_colors), by = "pad_name") |>
  st_transform(4326)

# Pad locations for mapping
pads_sf <- pads |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
```

``` r

swd_pt <- sites |> filter(type == "disposal")

maplibre(style = openfreemap_style("positron"), bounds = map_bounds) |>
  add_image_source(id = "cost", data = cost_surface, colors = cost_colors) |>
  add_raster_layer(id = "cost-layer", source = "cost", raster_opacity = 0.7) |>
  add_line_layer(
    id = "gathering", source = gathering,
    line_color = get_column("color"),
    line_width = 2.5, line_opacity = 0.9
  ) |>
  add_circle_layer(
    id = "pads", source = pads_sf,
    circle_color = "#1e1e1e", circle_radius = 7,
    circle_stroke_color = "white", circle_stroke_width = 2,
    tooltip = "pad_name"
  ) |>
  add_circle_layer(
    id = "swd", source = swd_pt,
    circle_color = "#dc2626", circle_radius = 9,
    circle_stroke_color = "white", circle_stroke_width = 2,
    tooltip = "label"
  ) |>
  add_categorical_legend(
    "Gathering System",
    values = c(pads$pad_name, "Ritchie Hunter SWD"),
    colors = c(pad_colors, "#dc2626"),
    position = "bottom-left",
    patch_shape = "circle"
  )
```

Each pad takes a different route to the disposal facility, reflecting
its position relative to the terrain. The Farley pad in the south has
the longest route, crossing several ridge-valley systems. The Holtz and
Michels pads are located further north and are close to one another; you
can see how the paths merge, showing how this tool can inform pipeline
*networks*, not just individual pipelines. The cost surface only needs
to be processed once; routing each additional pad is just another solve
on the same graph.

## Corridor characterization

Once a corridor has been identified, a natural next step is
understanding what the pipeline would actually cross. The `cell_indices`
stored in the corridor’s metadata make it possible to extract terrain
and land cover attributes for every cell along the path.

``` r

meta <- attr(path, "spopt")

path_profile <- tibble(
  elevation_m = terra::extract(dem, meta$cell_indices)[, 1],
  slope_deg   = terra::extract(slope, meta$cell_indices)[, 1],
  nlcd_class  = as.character(terra::extract(nlcd, meta$cell_indices)[, 1])
)

cat(sprintf("Corridor: Farley Unit -> Ritchie Hunter SWD\n"))
```

    Corridor: Farley Unit -> Ritchie Hunter SWD

``` r

cat(sprintf("  Length: %.1f km (straight-line: %.1f km)\n",
    path$path_dist / 1000, path$straight_line_dist / 1000))
```

      Length: 31.2 km (straight-line: 24.4 km)

``` r

cat(sprintf("  Sinuosity: %.3f\n", path$sinuosity))
```

      Sinuosity: 1.282

``` r

cat(sprintf("  Elevation: %d - %d m\n",
    min(path_profile$elevation_m), max(path_profile$elevation_m)))
```

      Elevation: 216 - 335 m

``` r

cat(sprintf("  Slope: mean %.1f deg, max %.1f deg\n",
    mean(path_profile$slope_deg), max(path_profile$slope_deg)))
```

      Slope: mean 4.0 deg, max 18.5 deg

The land cover breakdown gives a rough picture of where construction
complexity lies. Forest crossings require clearing and post-construction
restoration. Developed areas require right-of-way negotiation with
landowners. Pastureland is generally the simplest and cheapest to cross.

``` r

lc_labels <- tribble(
  ~nlcd_class, ~label,
  "21", "Developed, open space",
  "22", "Developed, low intensity",
  "23", "Developed, medium",
  "41", "Deciduous forest",
  "42", "Evergreen forest",
  "43", "Mixed forest",
  "52", "Shrub/scrub",
  "71", "Grassland",
  "81", "Pasture/hay",
  "82", "Cultivated crops"
)

path_profile |>
  count(nlcd_class, name = "cells") |>
  left_join(lc_labels, by = "nlcd_class") |>
  mutate(
    label = coalesce(label, paste("Class", nlcd_class)),
    pct = round(100 * cells / sum(cells), 1)
  ) |>
  arrange(desc(cells)) |>
  select(label, cells, pct)
```

    # A tibble: 8 × 3
      label                             cells   pct
      <chr>                             <int> <dbl>
    1 Class Pasture/Hay                   340  38.5
    2 Class Deciduous Forest              278  31.5
    3 Class Developed, Open Space         214  24.2
    4 Class Mixed Forest                   32   3.6
    5 Class Grassland/Herbaceous           10   1.1
    6 Class Developed, Low Intensity        7   0.8
    7 Class Developed, High Intensity       1   0.1
    8 Class Developed, Medium Intensity     1   0.1

This kind of per-corridor characterization bridges the gap between the
routing output and the downstream analyses that matter to project teams:
construction cost estimation, ROW acquisition planning, and
environmental permitting. spopt identifies the corridors; those
subsequent steps consume the corridor data to inform project decisions.

## Next steps

- [Route
  optimization](https://walker-data.com/spopt-r/articles/routing.md) -
  TSP and VRP with travel-time matrices
- [Facility
  Location](https://walker-data.com/spopt-r/articles/facility-location.md) -
  Optimize depot and warehouse placement
- [Getting
  Started](https://walker-data.com/spopt-r/articles/getting-started.md) -
  Overview of all spopt capabilities
