#' Ward Spatial Clustering
#'
#' Performs spatially-constrained hierarchical clustering using Ward's
#' minimum variance method. Only spatially contiguous areas can be merged,
#' ensuring all resulting regions are spatially connected.
#'
#' @param data An sf object with polygon or point geometries.
#' @param attrs Character vector of column names to use for clustering
#'   (e.g., `c("var1", "var2")`). If NULL, uses all numeric columns.
#' @param n_regions Integer. Number of regions (clusters) to create.
#' @param weights Spatial weights specification. One of "queen" (default),
#'   "rook", or an nb object from spdep.
#' @param scale Logical. If TRUE (default), standardize attributes before clustering.
#' @param verbose Logical. Print progress messages.
#'
#' @return An sf object with a `.region` column containing cluster assignments.
#'   Metadata is stored in the "spopt" attribute.
#'
#' @details
#' This function implements spatially-constrained agglomerative hierarchical
#' clustering using Ward's minimum variance criterion. Unlike standard Ward
#' clustering, this version enforces spatial contiguity by only allowing
#' clusters that share a border to be merged.
#'
#' The algorithm:
#' 1. Starts with each observation as its own cluster
#' 2. At each step, finds the pair of \strong{adjacent} clusters with minimum
#'    Ward distance (increase in total within-cluster variance)
#' 3. Merges them into a single cluster
#' 4. Repeats until the desired number of regions is reached
#'
#' The result guarantees that all regions are spatially contiguous.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#'
#' # Cluster into 8 spatially-contiguous regions
#' result <- ward_spatial(nc, attrs = c("SID74", "SID79"), n_regions = 8)
#' plot(result[".region"])
#' }
#'
#' @export
ward_spatial <- function(data,
                         attrs = NULL,
                         n_regions,
                         weights = "queen",
                         scale = TRUE,
                         verbose = FALSE) {
  # Input validation
  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object", call. = FALSE)
  }

  if (!is.numeric(n_regions) || n_regions < 2) {
    stop("`n_regions` must be an integer >= 2", call. = FALSE)
  }

  # Determine which columns to check for NAs
  check_cols <- if (!is.null(attrs)) attrs else character(0)

  # Validate data: remove empty geometries, check for NAs
  validated <- validate_regionalization_data(data, check_cols, call_name = "ward_spatial")
  data <- validated$data

  n <- nrow(data)
  if (n_regions >= n) {
    stop("`n_regions` must be less than number of observations", call. = FALSE)
  }

  # Extract attributes
  attr_matrix <- extract_attrs(data, attrs)

  if (scale) {
    attr_matrix <- scale(attr_matrix)
  }

  # Prepare spatial weights
  nb <- prepare_weights(data, weights)

  # Convert nb to adjacency indices
  adj <- nb_to_adj_indices(nb)

  if (verbose) {
    message(sprintf(
      "Ward Spatial: n=%d, n_regions=%d, attrs=%d",
      n, n_regions, ncol(attr_matrix)
    ))
  }

  # Call Rust implementation
  start_time <- Sys.time()

  result_list <- rust_ward_constrained(
    attrs = attr_matrix,
    n_regions = as.integer(n_regions),
    adj_i = adj$i,
    adj_j = adj$j
  )

  end_time <- Sys.time()

  # Extract results
  labels <- result_list$labels
  objective <- result_list$objective
  actual_n_regions <- result_list$n_regions

  # Attach results to sf object
  result <- data
  result$.region <- as.integer(labels)

  if (verbose) {
    message(sprintf(
      "  Result: %d regions, objective=%.4f, time=%.3fs",
      actual_n_regions, objective,
      as.numeric(difftime(end_time, start_time, units = "secs"))
    ))
  }

  # Attach metadata
  metadata <- list(
    algorithm = "ward_spatial",
    n_regions = actual_n_regions,
    objective = objective,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    scaled = scale,
    contiguity_enforced = TRUE
  )

  attach_spopt_metadata(result, metadata)
}
