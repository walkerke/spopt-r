# spopt <img src="man/figures/logo.png" align="right" height="139" />
<!-- badges: start -->
<!-- badges: end -->
Spatial Optimization for R. An R-native implementation of spatial optimization algorithms, inspired by [Python spopt](https://pysal.org/spopt/).

## Features

### Regionalization (Spatial Clustering)
- `skater()` - Spatial K'luster Analysis by Tree Edge Removal
- `ward_spatial()` - Spatially constrained Ward hierarchical clustering
- `max_p_regions()` - Maximize number of regions with threshold constraint (coming soon)
- `azp()` - Automatic Zoning Procedure (coming soon)
- `spenc()` - Spectral clustering with spatial constraints (coming soon)

### Facility Location
- `lscp()` - Location Set Covering Problem
- `mclp()` - Maximum Coverage Location Problem
- `p_median()` - P-Median (minimize total weighted distance)
- `p_center()` - P-Center (minimize maximum distance)
- `p_dispersion()` - P-Dispersion (maximize minimum inter-facility distance)
- `frlm()` - Flow Refueling Location Model (coming soon)

## Installation

**Requirements:** Rust toolchain (install from https://rustup.rs)

```r
# Install from GitHub
# install.packages("pak")
pak::pak("kylewalker/spopt-r")
```

## Usage

### Regionalization

```r
library(spopt)
library(sf)

# Load example data
nc <- st_read(system.file("shape/nc.shp", package = "sf"))

# Cluster into 8 regions based on SIDS rates
result <- skater(nc, attrs = ~ SID74 + SID79, n_regions = 8)

# Result is sf with .region column
plot(result[".region"])
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

# View selected facilities
result$facilities[result$facilities$.selected, ]

# View metadata
attr(result, "spopt")
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
