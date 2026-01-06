# Script to generate tarrant_travel_times example dataset
# Creates a list with origin tracts, candidate facilities, and travel-time matrix
#
# Prerequisites:
#   - osmium-tool installed (brew install osmium-tool on macOS)
#   - Java 21 for r5r (installed via rJavaEnv)
#   - Sufficient disk space for OSM download (~3.7GB)

library(tidycensus)
library(tidyverse)
library(sf)

# ============================================================================
# Step 1: Download and clip OSM data
# ============================================================================

# Create directory for OSM data
osm_dir <- file.path(tempdir(), "osm_data")
dir.create(osm_dir, showWarnings = FALSE)

us_south_pbf <- file.path(osm_dir, "us-south-latest.osm.pbf")
tarrant_pbf <- file.path(osm_dir, "tarrant.osm.pbf")

# Download US South extract from GeoFabrik (~700 MB)
if (!file.exists(us_south_pbf)) {
  message("Downloading US South OSM data from GeoFabrik...")
  options(timeout = 600)
  download.file(
    url = "https://download.geofabrik.de/north-america/us-south-latest.osm.pbf",
    destfile = us_south_pbf,
    mode = "wb"
  )
}

# Clip to Tarrant County bounding box using osmium
# bbox: xmin, ymin, xmax, ymax (approx Tarrant County with buffer)
if (!file.exists(tarrant_pbf)) {
  message("Clipping OSM data to Tarrant County...")
  osmium_cmd <- sprintf(
    'osmium extract -b -97.55,32.55,-97.0,33.05 "%s" -o "%s" --overwrite',
    us_south_pbf,
    tarrant_pbf
  )
  system(osmium_cmd)
}

# Remove the large US South file so r5r only uses the clipped file
if (file.exists(us_south_pbf)) {
  message("Removing large US South OSM file...")
  file.remove(us_south_pbf)
}

# ============================================================================
# Step 2: Set up r5r routing network
# ============================================================================

# Set Java memory before loading r5r
options(java.parameters = "-Xmx8G")

library(r5r)

# Install Java 21 if needed (only needs to be done once)
# rJavaEnv::java_quick_install(version = 21)

# Build routing network from clipped OSM data
message("Building r5r routing network...")
r5r_core <- build_network(data_path = osm_dir, verbose = FALSE)

# ============================================================================
# Step 3: Prepare Census tract data and candidate facilities
# ============================================================================

# Get tract-level population data for Tarrant County, TX
tarrant_tracts <- get_acs(
  geography = "tract",
  variables = "B01003_001",
  state = "TX",
  county = "Tarrant",
  geometry = TRUE,
  year = 2023
) |>
  filter(estimate > 0) |>
  rename(population = estimate) |>
  select(GEOID, NAME, population, geometry)

# Demand points: tract centroids in WGS84
demand_pts <- tarrant_tracts |>
  st_centroid() |>
  st_transform(4326) |>
  mutate(id = row_number())

# Candidate facilities: sample 30 locations across the county
set.seed(1983)
n_candidates <- 30

county_boundary <- tarrant_tracts |> st_union()
candidate_pts <- st_sample(county_boundary, n_candidates) |>
  st_as_sf() |>
  st_transform(4326) |>
  mutate(id = row_number())

# ============================================================================
# Step 4: Generate travel-time matrix with r5r
# ============================================================================

# Prepare points for r5r (requires id, lon, lat columns)
demand_r5r <- demand_pts |>
  st_coordinates() |>
  as_tibble() |>
  rename(lon = X, lat = Y) |>
  mutate(id = as.character(demand_pts$id))

candidates_r5r <- candidate_pts |>
  st_coordinates() |>
  as_tibble() |>
  rename(lon = X, lat = Y) |>
  mutate(id = as.character(candidate_pts$id))

# Calculate travel times (driving)
message("Calculating travel-time matrix (this may take a few minutes)...")
ttm <- travel_time_matrix(
  r5r_core,
  origins = demand_r5r,
  destinations = candidates_r5r,
  mode = "CAR",
  departure_datetime = as.POSIXct("2025-03-15 08:00:00"),
  max_trip_duration = 120,
  progress = TRUE
)

# Clean up r5r
stop_r5(r5r_core)

# ============================================================================
# Step 5: Convert to matrix format and bundle
# ============================================================================

# Pivot to wide format matrix (rows = demand, cols = facilities)
ttm_matrix <- ttm |>
  select(from_id, to_id, travel_time_p50) |>
  pivot_wider(
    names_from = to_id,
    values_from = travel_time_p50
  ) |>
  arrange(as.numeric(from_id)) |>
  select(-from_id) |>
  select(as.character(1:n_candidates)) |>
  as.matrix()

# Handle any missing values (unreachable pairs)
ttm_matrix[is.na(ttm_matrix)] <- Inf

# Verify dimensions
stopifnot(nrow(ttm_matrix) == nrow(demand_pts))
stopifnot(ncol(ttm_matrix) == nrow(candidate_pts))

message(sprintf(
  "Matrix dimensions: %d x %d",
  nrow(ttm_matrix),
  ncol(ttm_matrix)
))

# Bundle as a list
tarrant_travel_times <- list(
  tracts = tarrant_tracts,
  demand = demand_pts,
  candidates = candidate_pts,
  matrix = ttm_matrix
)

# ============================================================================
# Step 6: Save to data/
# ============================================================================

usethis::use_data(tarrant_travel_times, overwrite = TRUE)

message("Done! tarrant_travel_times saved to data/")
