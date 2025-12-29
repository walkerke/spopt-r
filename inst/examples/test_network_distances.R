# =============================================================================
# Testing Network Distances with spopt
# =============================================================================
# This script demonstrates using network distances instead of Euclidean
# for facility location problems.
#
# Two approaches:
#   1. dodgr - Road network distances (OSM-based, local computation)
#   2. r5r - Multimodal routing (car, transit, walk, bike)
#
# Install if needed:
#   install.packages("dodgr")
#   install.packages("r5r")  # Requires Java 21+

library(spopt)
library(sf)

# =============================================================================
# PART 1: dodgr - Road Network Distances
# =============================================================================
# dodgr downloads OSM data and computes distances locally
# Very fast: 1M pairwise distances in ~1.5 seconds

if (requireNamespace("dodgr", quietly = TRUE)) {
  library(dodgr)

  cat("\n========== DODGR WORKFLOW ==========\n\n")

  # --- Step 1: Download street network for Portland, OR ---
  cat("Downloading OSM street network (Portland, OR)...\n")

  # Define bounding box for area of interest
  bbox <- c(
    xmin = -122.70,
    ymin = 45.50,
    xmax = -122.64,
    ymax = 45.54
  )

  # Download network (silicate format for speed)
  net <- dodgr_streetnet_sc(bbox)
  cat("  Network downloaded:", nrow(net$edge), "edges\n")

  # --- Step 2: Create weighted graph ---
  cat("Creating weighted graph...\n")

  # Weight by travel time (motorcar profile)
  graph <- weight_streetnet(net, wt_profile = "motorcar")
  cat("  Graph edges:", nrow(graph), "\n")

  # Contract graph for faster routing (optional but recommended)
  graph_contracted <- dodgr_contract_graph(graph)
  cat("  Contracted to:", nrow(graph_contracted), "edges\n")

  # --- Step 3: Sample points FROM the network ---
  # This ensures all points are on roads (not in lakes, buildings, etc.)
  cat("\nSampling points from network vertices...\n")

  # Use vertices from the CONTRACTED graph - these are guaranteed connected
  # (contraction removes disconnected components)
  vertices <- dodgr_vertices(graph_contracted)
  cat("  Contracted network has", nrow(vertices), "vertices\n")

  set.seed(42)
  n_demand <- 20
  n_facilities <- 8

  # Sample vertices for demand and facilities (no overlap)
  sampled_idx <- sample(nrow(vertices), n_demand + n_facilities)
  demand_idx <- sampled_idx[1:n_demand]
  facility_idx <- sampled_idx[(n_demand + 1):(n_demand + n_facilities)]

  # Create sf objects from sampled vertices
  demand <- st_as_sf(
    data.frame(
      id = 1:n_demand,
      population = rpois(n_demand, 1000),
      x = vertices$x[demand_idx],
      y = vertices$y[demand_idx]
    ),
    coords = c("x", "y"),
    crs = 4326
  )

  facilities <- st_as_sf(
    data.frame(
      id = 1:n_facilities,
      x = vertices$x[facility_idx],
      y = vertices$y[facility_idx]
    ),
    coords = c("x", "y"),
    crs = 4326
  )

  cat("  Demand points:", nrow(demand), "(sampled from network)\n")
  cat("  Candidate facilities:", nrow(facilities), "(sampled from network)\n")

  # --- Step 4: Compute cost matrix ---
  cat("\nComputing travel time matrix...\n")

  from_coords <- st_coordinates(demand)
  to_coords <- st_coordinates(facilities)

  start_time <- Sys.time()

  # Travel times in seconds
  cost_matrix_time <- dodgr_times(
    graph_contracted,
    from = from_coords,
    to = to_coords
  )

  elapsed <- as.numeric(Sys.time() - start_time)
  cat("  Matrix computed in", round(elapsed, 3), "seconds\n")
  cat("  Dimensions:", nrow(cost_matrix_time), "x", ncol(cost_matrix_time), "\n")

  # Convert to minutes for interpretability
  cost_matrix_minutes <- cost_matrix_time / 60

  # Check for NA values (indicates unreachable points)
  n_na <- sum(is.na(cost_matrix_minutes))
  if (n_na > 0) {
    cat("  WARNING:", n_na, "NA values in cost matrix (unreachable point pairs)\n")
    cat("  This shouldn't happen when sampling from contracted graph vertices.\n")
    cat("  Checking connectivity...\n")
  } else {
    cat("  No NA values - all points are mutually reachable\n")
  }

  cat("  Travel time range:",
      round(min(cost_matrix_minutes, na.rm = TRUE), 1), "to",
      round(max(cost_matrix_minutes, na.rm = TRUE), 1), "minutes\n\n")

  # --- Step 4: Compare Euclidean vs Network ---
  cat("Comparing Euclidean vs Network distances...\n\n")

  # Euclidean distance matrix (what spopt uses by default)
  cost_matrix_euclidean <- as.matrix(st_distance(demand, facilities))
  # Convert from meters to km
  cost_matrix_euclidean <- matrix(
    as.numeric(cost_matrix_euclidean) / 1000,
    nrow = nrow(cost_matrix_euclidean)
  )

  # P-median with Euclidean distances
  cat("P-median with Euclidean distances:\n")
  result_euclidean <- p_median(
    demand, facilities,
    n_facilities = 3,
    weight_col = "population"
  )
  selected_euc <- which(result_euclidean$facilities$.selected)
  cat("  Selected facilities:", paste(selected_euc, collapse = ", "), "\n")
  # cat("  Objective:", round(attr(result_euclidean, "spopt")$total_cost, 2), "\n\n")

  # P-median with network travel times
  cat("P-median with network travel times:\n")
  result_network <- p_median(
    demand, facilities,
    n_facilities = 3,
    weight_col = "population",
    cost_matrix = cost_matrix_minutes
  )
  selected_net <- which(result_network$facilities$.selected)
  cat("  Selected facilities:", paste(selected_net, collapse = ", "), "\n")
  cat("  Objective:", round(attr(result_network, "spopt")$total_cost, 2), "pop-minutes\n\n")

  # Did the selection change?
  if (identical(selected_euc, selected_net)) {
    cat("Same facilities selected (Euclidean approximation was good here)\n")
  } else {
    cat("DIFFERENT facilities selected!\n")
    cat("  Euclidean:", paste(selected_euc, collapse = ", "), "\n")
    cat("  Network:", paste(selected_net, collapse = ", "), "\n")
  }

  # --- MCLP example ---
  cat("\n--- MCLP with Network Distances ---\n\n")

  # What's covered within 5 minutes drive?
  result_mclp <- mclp(
    demand, facilities,
    service_radius = 360,  # 5 minutes
    n_facilities = 1,
    weight_col = "population",
    cost_matrix = cost_matrix_minutes
  )

  meta <- attr(result_mclp, "spopt")
  cat("MCLP (5-minute service radius, 2 facilities):\n")
  cat("  Coverage:", round(meta$coverage_pct, 1), "%\n")
  cat("  Population covered:", meta$covered_weight, "/", meta$total_weight, "\n")

} else {
  cat("Install dodgr to run this section: install.packages('dodgr')\n")
}


# =============================================================================
# PART 2: r5r - Multimodal Routing
# =============================================================================
# r5r uses the R5 routing engine (Conveyal)
# Supports: car, bike, walk, transit (with GTFS)
# Requires: Java 21+, OSM PBF file, optionally GTFS

if (requireNamespace("r5r", quietly = TRUE)) {
  library(r5r)

  cat("\n\n========== R5R WORKFLOW ==========\n\n")

  cat("r5r requires:\n")
  cat("  1. Java 21+ installed\
")
  cat("  2. OSM PBF file for your region\n")
  cat("  3. GTFS file(s) for transit (optional)\n\n")

  cat("Download OSM PBF from: https://download.geofabrik.de/\n")
  cat("Find GTFS feeds at: https://transitfeeds.com/\n\n")

  # Example workflow (requires data files)
  cat("Example r5r workflow:\n\n")

  cat('
# --- Setup r5r ---
# Point to folder with OSM PBF and GTFS files
data_path <- "path/to/data"

# Initialize r5r (builds routing graph, takes a minute first time)
r5r_core <- setup_r5(data_path, verbose = FALSE)

# --- Compute travel time matrix ---
ttm <- travel_time_matrix(
  r5r_core,
  origins = demand,
  destinations = facilities,
  mode = c("WALK", "TRANSIT"),
  departure_datetime = as.POSIXct("2024-01-15 08:00:00"),
  max_trip_duration = 60  # minutes
)

# Reshape to matrix format
cost_matrix <- ttm |>
  tidyr::pivot_wider(
    id_cols = from_id,
    names_from = to_id,
    values_from = travel_time_p50
  ) |>
  dplyr::select(-from_id) |>
  as.matrix()

# --- Use with spopt ---
result <- p_median(
  demand, facilities,
  n_facilities = 5,
  weight_col = "population",
  cost_matrix = cost_matrix
)

# Clean up
stop_r5(r5r_core)
')

} else {
  cat("\n\nInstall r5r to run multimodal routing: install.packages('r5r')\n")
  cat("Note: r5r requires Java 21+\n")
}


# =============================================================================
# PART 3: Bring Your Own Matrix
# =============================================================================
# You can use ANY source for the cost matrix:
#   - Google Distance Matrix API
#   - HERE Routing API
#   - OSRM (self-hosted or demo server)
#   - Valhalla
#   - Pre-computed matrices from any source

cat("\n\n========== CUSTOM COST MATRIX ==========\n\n")

cat("spopt accepts any numeric matrix as cost_matrix.\n")
cat("The matrix should be:\n")
cat("  - Dimensions: n_demand x n_facilities\n")
cat("  - Units: consistent (all minutes, all meters, etc.)\n")
cat("  - No row/column names required\n\n")

cat("Example with manual matrix:\n\n")

# Small example
demand_small <- st_as_sf(data.frame(
  id = 1:5,
  pop = c(100, 200, 150, 300, 250),
  x = c(0, 1, 2, 3, 4),
  y = c(0, 0, 0, 0, 0)
), coords = c("x", "y"))

facilities_small <- st_as_sf(data.frame(
  id = 1:3,
  x = c(0.5, 2, 3.5),
  y = c(0, 0, 0)
), coords = c("x", "y"))

# Custom cost matrix (e.g., from Google API)
# Travel times in minutes
custom_costs <- matrix(c(
  2, 8, 15,   # demand 1 to facilities 1,2,3
  3, 5, 12,   # demand 2
  7, 2, 7,    # demand 3
  12, 5, 3,   # demand 4
  15, 8, 2    # demand 5
), nrow = 5, byrow = TRUE)

cat("Custom cost matrix (travel times in minutes):\n")
print(custom_costs)

result <- p_median(
  demand_small, facilities_small,
  n_facilities = 2,
  weight_col = "pop",
  cost_matrix = custom_costs
)

cat("\nP-median result:\n")
cat("  Selected:", which(result$facilities$.selected), "\n")
cat("  Total weighted travel time:",
    round(attr(result, "spopt")$total_cost, 1), "person-minutes\n")


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n\n========== SUMMARY ==========\n\n")
cat("Network distance options for spopt:\n\n")
cat("1. dodgr (recommended for driving distances)\n")
cat("   - Downloads OSM automatically\n")
cat("   - Fast local computation\n")
cat("   - dodgr_times() / dodgr_dists()\n\n")
cat("2. r5r (recommended for transit/multimodal)\n")
cat("   - Requires OSM PBF + GTFS files\n")
cat("   - Time-dependent routing\n")
cat("   - travel_time_matrix()\n\n")
cat("3. API services (Google, HERE, OSRM, etc.)\n")
cat("   - Build matrix yourself, pass to spopt\n")
cat("   - Watch rate limits/costs\n\n")
cat("All spopt facility location functions accept cost_matrix parameter:\n")
cat("  p_median(), mclp(), lscp(), p_center()\n")
