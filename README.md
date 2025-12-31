# spopt <img src="man/figures/logo.png" align="right" height="139" />
<!-- badges: start -->
<!-- badges: end -->
Spatial Optimization for R. An R-native implementation of spatial optimization algorithms, inspired by [Python spopt](https://pysal.org/spopt/).

## Features

### Regionalization (Spatial Clustering)
- `skater()` - Spatial K'luster Analysis by Tree Edge Removal
- `ward_spatial()` - Spatially constrained Ward hierarchical clustering
- `max_p_regions()` - Maximize number of regions with threshold constraint
- `azp()` - Automatic Zoning Procedure (basic, tabu, simulated annealing)
- `spenc()` - Spatially-encouraged spectral clustering

### Facility Location
- `lscp()` - Location Set Covering Problem
- `mclp()` - Maximum Coverage Location Problem
- `p_median()` - P-Median (minimize total weighted distance)
- `p_center()` - P-Center (minimize maximum distance)
- `p_dispersion()` - P-Dispersion (maximize minimum inter-facility distance)
- `cflp()` - Capacitated Facility Location (facilities with capacity limits)
- `frlm()` - Flow Refueling Location Model

### Market Analysis
- `huff()` - Huff Model for market share and trade area analysis

### Utilities
- `sp_weights()` - Compute spatial contiguity weights (queen/rook)
- `distance_matrix()` - Compute distance matrices between sf objects

## Installation

### System Requirements

This package requires a Rust toolchain to compile the backend.

**macOS/Linux:**
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

**Windows:**
Windows installation uses pre-built binaries and does not require Rust. Just install directly:
```r
pak::pak("walkerke/spopt-r")
```

If pre-built binaries are not yet available for your version, you can use [WSL](https://learn.microsoft.com/en-us/windows/wsl/) with the Linux instructions above.

**macOS/Linux:** After Rust installation, restart your terminal/R session so `cargo` is on your PATH.

### Install Package

```r
# install.packages("pak")
pak::pak("walkerke/spopt-r")
```

The first installation will compile the Rust code, which takes 1-2 minutes.

## Usage

### Regionalization

```r
library(spopt)
library(sf)

# Load example data
nc <- st_read(system.file("shape/nc.shp", package = "sf"))

# Cluster into 8 regions based on SIDS rates
result <- skater(nc, attrs = c("SID74", "SID79"), n_regions = 8)

# Result is sf with .region column
plot(result[".region"])

# Max-P: maximize regions where each has >= 100,000 births
result <- max_p_regions(
 nc,
 attrs = c("SID74", "SID79"),
 threshold_var = "BIR74",
 threshold = 100000
)
attr(result, "spopt")$n_regions
```

### Facility Location

```r
library(spopt)
library(sf)

# Create demand and facility points
set.seed(42)
demand <- st_as_sf(
  data.frame(x = runif(100), y = runif(100), population = rpois(100, 500)),
  coords = c("x", "y")
)
facilities <- st_as_sf(
  data.frame(x = runif(20), y = runif(20)),
  coords = c("x", "y")
)

# P-Median: minimize total weighted distance with 5 facilities
result <- p_median(demand, facilities, n_facilities = 5, weight_col = "population")

# MCLP: maximize population covered within 0.3 units with 3 facilities
result <- mclp(demand, facilities, service_radius = 0.3,
               n_facilities = 3, weight_col = "population")
attr(result, "spopt")$coverage_pct

# View selected facilities
result$facilities[result$facilities$.selected, ]
```

### Network Distances

All facility location functions accept a custom `cost_matrix` for network-based analysis:

```r
library(dodgr)  # for road network distances

# Get street network and compute travel times
net <- dodgr_streetnet_sc(pts = st_coordinates(demand))
graph <- weight_streetnet(net, wt_profile = "motorcar")

# Compute travel time matrix (minutes)
cost_matrix <- dodgr_times(graph,
                           from = st_coordinates(demand),
                           to = st_coordinates(facilities)) / 60

# Use network times instead of Euclidean distance
result <- p_median(demand, facilities, n_facilities = 5,
                   weight_col = "population", cost_matrix = cost_matrix)
```

## Design Principles

- **sf first**: All functions accept sf objects and return sf objects
- **Snake case API**: `max_p_regions()`, `p_median()`, `sp_weights()`
- **R-native**: No Python/reticulate dependency
- **Fast**: Rust backend via extendr for performance-critical algorithms
- **HiGHS solver**: Open-source MIP solver for facility location problems

## License

MIT

## Acknowledgements

- [PySAL spopt](https://pysal.org/spopt/) for algorithm implementations and inspiration
- [extendr](https://extendr.github.io/) for R/Rust bindings
- [HiGHS](https://highs.dev/) for the optimization solver
