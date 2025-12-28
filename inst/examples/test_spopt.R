# =============================================================================
# spopt R Package - Test Script
# =============================================================================
# Run this after installing the package with devtools::install()

library(spopt)
library(sf)

# -----------------------------------------------------------------------------
# WHAT WORKS (Fully Implemented)
# -----------------------------------------------------------------------------

cat("\n========== TESTING WORKING FUNCTIONS ==========\n\n")

# --- Utility Functions ---

cat("1. Testing sp_weights()...\n")
nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
w_queen <- sp_weights(nc, type = "queen")
w_rook <- sp_weights(nc, type = "rook")
cat("   Queen neighbors for county 1:", paste(w_queen[[1]], collapse = ", "), "\n")
cat("   Rook neighbors for county 1:", paste(w_rook[[1]], collapse = ", "), "\n")
cat("   SUCCESS\n\n")

cat("2. Testing distance_matrix()...\n")
set.seed(42)
pts1 <- st_as_sf(data.frame(x = runif(10), y = runif(10)), coords = c("x", "y"))
pts2 <- st_as_sf(data.frame(x = runif(5), y = runif(5)), coords = c("x", "y"))
d <- distance_matrix(pts1, pts2)
cat("   Distance matrix dimensions:", nrow(d), "x", ncol(d), "\n")
cat("   Min distance:", round(min(d), 4), "Max distance:", round(max(d), 4), "\n")
cat("   SUCCESS\n\n")

# --- Regionalization ---

cat("3. Testing skater() - SKATER regionalization...\n")
# Using NC data with SIDS rates
result_skater <- skater(
  nc,
  attrs = ~ SID74 + SID79,
  n_regions = 8,
  seed = 42
)
cat("   Result class:", paste(class(result_skater), collapse = ", "), "\n")
cat("   Regions created:", length(unique(result_skater$.region)), "\n")
cat("   Metadata available:", !is.null(attr(result_skater, "spopt")), "\n")
cat("   SUCCESS\n\n")

cat("4. Testing ward_spatial() - Ward clustering...\n")
result_ward <- ward_spatial(
  nc,
  attrs = ~ SID74 + SID79,
  n_regions = 8
)
cat("   Result class:", paste(class(result_ward), collapse = ", "), "\n")
cat("   Regions created:", length(unique(result_ward$.region)), "\n")
cat("   SUCCESS\n\n")

cat("5. Testing max_p_regions() - Maximize regions with threshold constraint...\n")
result_maxp <- max_p_regions(
  nc,
  attrs = ~ SID74 + SID79,
  threshold_var = "BIR74",
  threshold = 50000,
  n_iterations = 100,
  n_sa_iterations = 100,
  seed = 42
)
meta_maxp <- attr(result_maxp, "spopt")
cat("   Regions created:", meta_maxp$n_regions, "\n")
cat("   Objective (within-region SSD):", round(meta_maxp$objective, 2), "\n")
cat("   Solve time:", round(meta_maxp$solve_time, 4), "seconds\n")

# Verify all regions meet threshold
all_meet <- all(sapply(unique(result_maxp$.region), function(r) {
  sum(nc$BIR74[result_maxp$.region == r]) >= 50000
}))
cat("   All regions meet threshold:", all_meet, "\n")
cat("   SUCCESS\n\n")

# --- Facility Location ---

cat("6. Testing lscp() - Location Set Covering...\n")
set.seed(42)
demand <- st_as_sf(
  data.frame(x = runif(30), y = runif(30)),
  coords = c("x", "y")
)
facilities <- st_as_sf(
  data.frame(x = runif(10), y = runif(10)),
  coords = c("x", "y")
)

result_lscp <- lscp(demand, facilities, service_radius = 0.4)
cat("   Facilities selected:", sum(result_lscp$facilities$.selected), "\n")
cat("   Coverage:", round(attr(result_lscp, "spopt")$coverage_pct, 1), "%\n")
cat("   SUCCESS\n\n")

cat("7. Testing mclp() - Maximum Coverage Location...\n")
demand$population <- rpois(30, 500)
result_mclp <- mclp(
  demand, facilities,
  service_radius = 0.3,
  n_facilities = 3,
  weight_col = "population"
)
cat("   Facilities selected:", sum(result_mclp$facilities$.selected), "\n")
cat("   Coverage:", round(attr(result_mclp, "spopt")$coverage_pct, 1), "%\n")
cat("   SUCCESS\n\n")

cat("8. Testing p_median() - Minimize total weighted distance...\n")
result_pmedian <- p_median(
  demand, facilities,
  n_facilities = 4,
  weight_col = "population"
)
cat("   Facilities selected:", sum(result_pmedian$facilities$.selected), "\n")
cat("   Mean distance:", round(attr(result_pmedian, "spopt")$mean_distance, 4), "\n")
cat("   All demand assigned:", all(!is.na(result_pmedian$demand$.facility)), "\n")
cat("   SUCCESS\n\n")

cat("9. Testing p_center() - Minimize maximum distance...\n")
result_pcenter <- p_center(
  demand, facilities,
  n_facilities = 4
)
cat("   Facilities selected:", sum(result_pcenter$facilities$.selected), "\n")
cat("   Max distance:", round(attr(result_pcenter, "spopt")$max_distance, 4), "\n")
cat("   SUCCESS\n\n")

cat("10. Testing p_dispersion() - Maximize facility spread...\n")
result_pdisp <- p_dispersion(
  facilities,
  n_facilities = 4
)
cat("   Facilities selected:", sum(result_pdisp$.selected), "\n")
cat("   Min inter-facility distance:", round(attr(result_pdisp, "spopt")$min_distance, 4), "\n")
cat("   SUCCESS\n\n")

# -----------------------------------------------------------------------------
# MAX-P REGIONS BENCHMARK
# -----------------------------------------------------------------------------

cat("========== MAX-P REGIONS BENCHMARK ==========\n\n")

cat("Testing max_p_regions with varying thresholds...\n\n")
thresholds <- c(30000, 50000, 75000, 100000)
for (thresh in thresholds) {
  start <- Sys.time()
  result <- max_p_regions(
    nc,
    attrs = ~ SID74 + SID79,
    threshold_var = "BIR74",
    threshold = thresh,
    n_iterations = 100,
    seed = 42
  )
  elapsed <- as.numeric(Sys.time() - start, units = "secs")
  n_regions <- attr(result, "spopt")$n_regions

  # Verify threshold constraint
  all_meet <- all(sapply(unique(result$.region), function(r) {
    sum(nc$BIR74[result$.region == r]) >= thresh
  }))

  cat(sprintf("   threshold=%6.0f: %2d regions, valid=%s, time=%.3fs\n",
      thresh, n_regions, all_meet, elapsed))
}

# -----------------------------------------------------------------------------
# TEXAS TRACTS BENCHMARK (Large-scale test)
# -----------------------------------------------------------------------------

cat("\n========== TEXAS TRACTS BENCHMARK ==========\n\n")

if (requireNamespace("tidycensus", quietly = TRUE)) {
  library(tidycensus)

  cat("Fetching TX tracts from Census API...\n")
  options(tigris_use_cache = TRUE)

  tx <- get_acs(
    geography = "tract",
    variables = c(pop = "B01003_001", income = "B19013_001"),
    state = "TX",
    geometry = TRUE,
    year = 2022
  ) |>
    tidyr::pivot_wider(
      id_cols = c(GEOID, NAME, geometry),
      names_from = variable,
      values_from = c(estimate, moe)
    ) |>
    sf::st_as_sf() |>
    dplyr::filter(!sf::st_is_empty(geometry)) |>
    dplyr::filter(!is.na(estimate_pop), !is.na(estimate_income)) |>
    dplyr::filter(estimate_pop > 0)

  cat(sprintf("Dataset: %d tracts, total pop = %s\n\n",
      nrow(tx), format(sum(tx$estimate_pop), big.mark = ",")))

  # Benchmark with different thresholds
  thresholds <- c(50000, 100000, 200000)

  for (thresh in thresholds) {
    cat(sprintf("Threshold: %s population\n", format(thresh, big.mark = ",")))

    start <- Sys.time()
    result_tx <- max_p_regions(
      tx,
      attrs = ~ estimate_income,
      threshold_var = "estimate_pop",
      threshold = thresh,
      n_iterations = 100,
      n_sa_iterations = 100,
      seed = 42
    )
    elapsed <- as.numeric(Sys.time() - start, units = "secs")

    meta <- attr(result_tx, "spopt")
    all_meet <- all(meta$region_stats$meets_threshold)

    cat(sprintf("  Regions: %d\n", meta$n_regions))
    cat(sprintf("  All meet threshold: %s\n", all_meet))
    cat(sprintf("  Time: %.2f seconds\n\n", elapsed))
  }

  # Detailed stats for last run
  cat("Region size distribution (last run):\n")
  sizes <- table(result_tx$.region)
  cat(sprintf("  Tracts per region: min=%d, median=%.0f, max=%d\n",
      min(sizes), median(sizes), max(sizes)))

  pops <- sapply(unique(result_tx$.region), function(r) {
    sum(tx$estimate_pop[result_tx$.region == r])
  })
  cat(sprintf("  Population per region: min=%s, median=%s, max=%s\n",
      format(min(pops), big.mark = ","),
      format(median(pops), big.mark = ","),
      format(max(pops), big.mark = ",")))

} else {
  cat("Skipping TX tracts benchmark (tidycensus not installed)\n")
  cat("Install with: install.packages('tidycensus')\n")
}

# -----------------------------------------------------------------------------
# VISUALIZATION EXAMPLES
# -----------------------------------------------------------------------------

cat("\n========== VISUALIZATION ==========\n\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  cat("Plotting SKATER results...\n")
  p1 <- ggplot(result_skater) +
    geom_sf(aes(fill = factor(.region)), color = "white", linewidth = 0.2) +
    scale_fill_viridis_d(name = "Region") +
    labs(title = "SKATER Regionalization of NC Counties",
         subtitle = "8 regions based on SID74 + SID79") +
    theme_minimal()
  print(p1)
  cat("   SKATER plot displayed\n\n")

  cat("Plotting Max-P results...\n")
  p2 <- ggplot(result_maxp) +
    geom_sf(aes(fill = factor(.region)), color = "white", linewidth = 0.2) +
    scale_fill_viridis_d(name = "Region") +
    labs(title = "Max-P Regionalization of NC Counties",
         subtitle = sprintf("%d regions, each with BIR74 >= 50,000",
                            attr(result_maxp, "spopt")$n_regions)) +
    theme_minimal()
  print(p2)
  cat("   Max-P plot displayed\n\n")

} else {
  cat("Plotting with base R (install ggplot2 for better plots)...\n")
  par(mfrow = c(1, 2))
  plot(result_skater[".region"], main = "SKATER Regionalization")
  plot(result_maxp[".region"], main = "Max-P Regionalization")
  par(mfrow = c(1, 1))
  cat("   Plots displayed\n\n")
}

# -----------------------------------------------------------------------------
# WHAT'S STILL TODO / NOT YET IMPLEMENTED
# -----------------------------------------------------------------------------

cat("========== TODO / NOT IMPLEMENTED ==========\n\n")

cat("1. azp() - Automatic Zoning Procedure\n")
cat("   Status: Rust stub exists, returns empty. R wrapper not created.\n")
cat("   What's needed: Implement local search with tabu/SA variants in Rust\n\n")

cat("2. spenc() - Spectral clustering with spatial constraints\n")
cat("   Status: Not started\n")
cat("   What's needed: Eigendecomposition + k-means in Rust\n\n")

cat("3. frlm() - Flow Refueling Location Model\n")
cat("   Status: Not started\n")
cat("   What's needed: Network flow optimization, most complex algorithm\n\n")

cat("4. ward_spatial() spatial constraint\n")
cat("   Status: Works but does NOT enforce spatial contiguity\n")
cat("   What's needed: Implement constrained agglomerative clustering\n\n")

cat("========== SUMMARY ==========\n\n")
cat("WORKING (10 functions):\n")
cat("  - sp_weights(), distance_matrix()\n")
cat("  - skater(), ward_spatial(), max_p_regions()\n")
cat("  - lscp(), mclp(), p_median(), p_center(), p_dispersion()\n\n")

cat("TODO (4 functions):\n")
cat("  - azp(), spenc()\n")
cat("  - frlm()\n")
cat("  - ward_spatial() spatial constraint enforcement\n\n")

cat("All facility location algorithms use HiGHS MIP solver via Rust.\n")
cat("All regionalization algorithms use parallel Rust implementations.\n")
cat("All working functions return sf objects as first-class citizens.\n")
