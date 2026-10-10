# Getting started with spopt

spopt provides R-native implementations of spatial optimization
algorithms for regionalization, facility location, and market analysis.
The package is inspired by Python’s PySAL spopt library, my all-time
favorite Python package. Its aim is to bring spopt’s powerful algorithms
to R users with an sf-first API and a Rust backend for performance.

### Installation

The easiest way to install spopt is from my r-universe repository, which
provides pre-built binaries:

\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"spopt"``, repos ``=`` ``"https://walkerke.r-universe.dev"``)`

Once installed, load the package along with sf for spatial data handling
and tidyverse for data manipulation:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`spopt`](https://walker-data.com/spopt-r/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sf`](https://r-spatial.github.io/sf/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyverse`](https://tidyverse.tidyverse.org)`)`

### What can spopt do?

spopt includes four families of spatial optimization algorithms:

**Regionalization**: Build spatially-contiguous regions from smaller
geographies. This is useful when you need to aggregate Census blocks
into larger areas, create balanced sales territories, or design compact
political districts. Algorithms include SKATER, AZP, Max-P, SPENC, and
spatially-constrained Ward clustering.

**Facility location**: Find optimal locations for facilities given
demand points and candidate sites. Whether you’re siting fire stations
to minimize response times, placing retail stores to maximize coverage,
or locating EV charging stations along highway corridors, spopt has
algorithms for these problems. Options include P-Median, P-Center, MCLP,
LSCP, CFLP, P-Dispersion, and FRLM.

**Route optimization**: Find the best stop sequence for delivery
drivers, field technicians, or any visit-based workflow. Single-vehicle
routing with open routes and fixed endpoints
([`route_tsp()`](https://walker-data.com/spopt-r/reference/route_tsp.md)),
or multi-vehicle dispatch with capacity constraints
([`route_vrp()`](https://walker-data.com/spopt-r/reference/route_vrp.md)).
Bring a travel-time matrix from r5r, OSRM, or any routing engine; spopt
handles the sequencing.

**Market analysis**: Model consumer behavior and market competition with
the Huff model. This classic retail gravity model predicts market share
and expected sales based on store attractiveness and distance to
consumers.

### A quick example

Let’s run a quick facility location analysis to see spopt in action.
We’ll find optimal locations for 5 facilities to serve Census tracts in
Tarrant County, Texas (home of Fort Worth).

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidycensus`](https://walker-data.com/tidycensus/)`)`\
\
`# Get population data for Tarrant County tracts`\
`tarrant`` ``<-`` `[`get_acs`](https://walker-data.com/tidycensus/reference/get_acs.html)`(`\
`  geography ``=`` ``"tract"``,`\
`  variables ``=`` ``"B01003_001"``,`\
`  state ``=`` ``"TX"``,`\
`  county ``=`` ``"Tarrant"``,`\
`  geometry ``=`` ``TRUE``,`\
`  year ``=`` ``2023`\
`)`\
\
`# Use tract centroids as both demand points and candidate facility sites`\
`tarrant_pts`` ``<-`` ``tarrant`` ``|>`\
`  `[`st_centroid`](https://r-spatial.github.io/sf/reference/geos_unary.html)`(``)`` ``|>`\
`  `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``estimate``)``)`\
\
`# Solve the P-Median problem: minimize total weighted distance`\
`result`` ``<-`` `[`p_median`](https://walker-data.com/spopt-r/reference/p_median.md)`(`\
`  demand ``=`` ``tarrant_pts``,`\
`  facilities ``=`` ``tarrant_pts``,`\
`  n_facilities ``=`` ``5``,`\
`  weight_col ``=`` ``"estimate"`\
`)`\
\
`# View selected facility locations`\
`result``$``facilities`` ``|>`` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``.selected``)`

    Simple feature collection with 5 features and 8 fields
    Geometry type: POINT
    Dimension:     XY
    Bounding box:  xmin: -97.40206 ymin: 32.66538 xmax: -97.12351 ymax: 32.89751
    Geodetic CRS:  NAD83
            GEOID                                        NAME   variable estimate
    1 48439111516 Census Tract 1115.16; Tarrant County; Texas B01003_001     7137
    2 48439113637 Census Tract 1136.37; Tarrant County; Texas B01003_001     4732
    3 48439110402 Census Tract 1104.02; Tarrant County; Texas B01003_001     5387
    4 48439105800    Census Tract 1058; Tarrant County; Texas B01003_001     4430
    5 48439113946 Census Tract 1139.46; Tarrant County; Texas B01003_001     6509
       moe                   geometry .selected .fixed .n_assigned
    1 1007 POINT (-97.12453 32.66618)      TRUE  FALSE         112
    2  703 POINT (-97.12351 32.84844)      TRUE  FALSE          90
    3   16 POINT (-97.40206 32.80275)      TRUE  FALSE          69
    4  698 POINT (-97.34284 32.66538)      TRUE  FALSE         107
    5 3880 POINT (-97.28117 32.89751)      TRUE  FALSE          71

The
[`p_median()`](https://walker-data.com/spopt-r/reference/p_median.md)
function returns a list with allocated demand points and selected
facilities. Each demand point is assigned to its nearest selected
facility, and the solution minimizes the total population-weighted
distance.

### SF-first design

All spopt functions accept and return sf objects. This means you can:

- Pass sf polygons or points directly to functions
- Use any coordinate reference system (functions handle transformations
  internally)
- Pipe results directly into visualization with ggplot2, mapview, or
  mapgl
- Integrate seamlessly with tidycensus, tigris, and other spatial R
  packages

### Next steps

- [Regionalization](https://walker-data.com/spopt-r/articles/regionalization.md) -
  Build spatially-contiguous regions with SKATER, AZP, Max-P, and more
- [Facility
  Location](https://walker-data.com/spopt-r/articles/facility-location.md) -
  Solve location-allocation problems with P-Median, MCLP, and other
  algorithms
- [Route
  Optimization](https://walker-data.com/spopt-r/articles/routing.md) -
  Optimize delivery routes and fleet dispatch with TSP and VRP
- [Huff Model](https://walker-data.com/spopt-r/articles/huff-model.md) -
  Model market share and retail competition
- [Travel-Time Cost
  Matrices](https://walker-data.com/spopt-r/articles/travel-time-matrices.md) -
  Use real-world travel times with r5r and other routing engines
