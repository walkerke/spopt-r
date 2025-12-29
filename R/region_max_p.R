#' Max-P Regions
#'
#' Perform Max-P regionalization to maximize the number of spatially contiguous
#' regions such that each region satisfies a minimum threshold constraint on a
#' specified attribute. This is useful for creating regions that meet minimum
#' population or sample size requirements.
#'
#' @param data An sf object with polygon or point geometries.
#' @param attrs Character vector of column names to use for computing
#'   within-region dissimilarity (e.g., `c("var1", "var2")`). If NULL,
#'   uses all numeric columns.
#' @param threshold_var Character. Name of the column containing the threshold
#'   variable (e.g., population, income). Each region must have a sum of this
#'   variable >= `threshold`.
#' @param threshold Numeric. Minimum sum of `threshold_var` required per region.
#' @param weights Spatial weights specification. One of "queen" (default),
#'   "rook", or an nb object from spdep.
#' @param n_iterations Integer. Number of construction phase iterations (default 100).
#'   Higher values explore more random starting solutions.
#' @param n_sa_iterations Integer. Number of simulated annealing iterations (default 100).
#'   Set to 0 to skip the SA refinement phase.
#' @param cooling_rate Numeric. SA cooling rate between 0 and 1 (default 0.99).
#'   Smaller values cool faster, larger values allow more exploration.
#' @param tabu_length Integer. Length of tabu list for SA phase (default 10).
#' @param scale Logical. If TRUE (default), standardize attributes before
#'   computing dissimilarity.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical. Print progress messages.
#'
#' @return An sf object with a `.region` column containing region assignments.
#'   Metadata is stored in the "spopt" attribute, including:
#'   \itemize{
#'     \item algorithm: "max_p"
#'     \item n_regions: Number of regions created (the "p" in max-p)
#'     \item objective: Total within-region sum of squared deviations
#'     \item threshold_var: Name of threshold variable
#'     \item threshold: Threshold value used
#'     \item solve_time: Time to solve in seconds
#'   }
#'
#' @details
#' The Max-P algorithm (Duque, Anselin & Rey, 2012; Wei, Rey & Knaap, 2020)
#' solves the problem of aggregating n geographic areas into the maximum number
#' of homogeneous regions while ensuring:
#'
#' 1. Each region is spatially contiguous (connected)
#' 2. Each region satisfies a minimum threshold on a specified attribute
#'
#' The algorithm has two phases:
#' \enumerate{
#'   \item Construction phase: Builds feasible solutions via randomized greedy
#'     region growing. Multiple random starts are explored in parallel.
#'   \item Simulated annealing phase: Refines solutions by moving border areas
#'     between regions to minimize within-region dissimilarity while respecting
#'     constraints.
#' }
#'
#' This implementation is optimized for speed using:
#' \itemize{
#'   \item Parallel construction with early termination
#'   \item Efficient articulation point detection for move eligibility
#'   \item Incremental threshold tracking
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#'
#' # Create regions where each has at least 100,000 in BIR74
#' result <- max_p_regions(
#'   nc,
#'   attrs = c("SID74", "SID79"),
#'   threshold_var = "BIR74",
#'   threshold = 100000
#' )
#'
#' # Check number of regions created
#' attr(result, "spopt")$n_regions
#'
#' # Plot results
#' plot(result[".region"])
#' }
#'
#' @references
#' Duque, J. C., Anselin, L., & Rey, S. J. (2012). The max-p-regions problem.
#' Journal of Regional Science, 52(3), 397-419.
#'
#' Wei, R., Rey, S., & Knaap, E. (2020). Efficient regionalization for spatially
#' explicit neighborhood delineation. International Journal of Geographical
#' Information Science, 35(1), 135-151.
#'
#' @export
max_p_regions <- function(data,
                          attrs = NULL,
                          threshold_var,
                          threshold,
                          weights = "queen",
                          n_iterations = 100L,
                          n_sa_iterations = 100L,
                          cooling_rate = 0.99,
                          tabu_length = 10L,
                          scale = TRUE,
                          seed = NULL,
                          verbose = FALSE) {
  # Input validation
  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object", call. = FALSE)
  }

  n <- nrow(data)
  if (n < 2) {
    stop("Data must have at least 2 observations", call. = FALSE)
  }

  # Validate threshold variable
  if (missing(threshold_var) || !is.character(threshold_var) || length(threshold_var) != 1) {
    stop("`threshold_var` must be a single column name", call. = FALSE)
  }
  if (!threshold_var %in% names(data)) {
    stop(paste0("Threshold variable '", threshold_var, "' not found in data"), call. = FALSE)
  }
  threshold_values <- as.numeric(data[[threshold_var]])
  if (any(is.na(threshold_values))) {
    stop("Threshold variable contains NA values", call. = FALSE)
  }

  # Validate threshold
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0) {
    stop("`threshold` must be a positive number", call. = FALSE)
  }

  total_threshold <- sum(threshold_values)
  if (total_threshold < threshold) {
    stop(paste0(
      "Total of threshold variable (", round(total_threshold, 2),
      ") is less than threshold (", threshold,
      "). No valid regionalization possible."
    ), call. = FALSE)
  }

  # Validate other parameters
  if (!is.numeric(n_iterations) || n_iterations < 1) {
    stop("`n_iterations` must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(n_sa_iterations) || n_sa_iterations < 0) {
    stop("`n_sa_iterations` must be a non-negative integer", call. = FALSE)
  }
  if (!is.numeric(cooling_rate) || cooling_rate <= 0 || cooling_rate >= 1) {
    stop("`cooling_rate` must be between 0 and 1 (exclusive)", call. = FALSE)
  }
  if (!is.numeric(tabu_length) || tabu_length < 0) {
    stop("`tabu_length` must be a non-negative integer", call. = FALSE)
  }

  # Extract attributes for dissimilarity
  attr_matrix <- extract_attrs(data, attrs)

  if (scale) {
    attr_matrix <- scale(attr_matrix)
  }

  # Prepare spatial weights
  nb <- prepare_weights(data, weights)

  # Convert nb to sparse matrix indices
  adj <- nb_to_adj_indices(nb)

  if (verbose) {
    message(sprintf(
      "Max-P: n=%d, threshold=%.2f (var=%s), attrs=%d, edges=%d",
      n, threshold, threshold_var, ncol(attr_matrix), length(adj$i)
    ))
    message(sprintf(
      "  Construction iterations: %d, SA iterations: %d",
      n_iterations, n_sa_iterations
    ))
  }

  # Call Rust implementation
  start_time <- Sys.time()

  result_list <- rust_max_p(
    attrs = attr_matrix,
    threshold_var = threshold_values,
    threshold = as.numeric(threshold),
    adj_i = adj$i,
    adj_j = adj$j,
    n_iterations = as.integer(n_iterations),
    n_sa_iterations = as.integer(n_sa_iterations),
    cooling_rate = as.numeric(cooling_rate),
    tabu_length = as.integer(tabu_length),
    seed = if (!is.null(seed)) as.integer(seed) else NULL
  )

  end_time <- Sys.time()

  # Extract results
  labels <- result_list$labels
  n_regions <- result_list$n_regions
  objective <- result_list$objective

  # Attach results to sf object
  result <- data
  result$.region <- as.integer(labels)

  if (verbose) {
    message(sprintf(
      "  Result: %d regions, objective=%.4f, time=%.3fs",
      n_regions, objective,
      as.numeric(difftime(end_time, start_time, units = "secs"))
    ))
  }

  # Compute region statistics
  region_stats <- compute_region_stats(result, threshold_var, threshold)

  # Attach metadata
  metadata <- list(
    algorithm = "max_p",
    n_regions = n_regions,
    objective = objective,
    threshold_var = threshold_var,
    threshold = threshold,
    region_stats = region_stats,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    scaled = scale,
    n_iterations = n_iterations,
    n_sa_iterations = n_sa_iterations
  )

  attach_spopt_metadata(result, metadata)
}

#' Compute region statistics
#' @keywords internal
compute_region_stats <- function(result, threshold_var, threshold) {
  # Compute stats for each region
  regions <- unique(result$.region)
  stats <- lapply(regions, function(r) {
    idx <- result$.region == r
    list(
      region = r,
      n_areas = sum(idx),
      threshold_sum = sum(result[[threshold_var]][idx]),
      meets_threshold = sum(result[[threshold_var]][idx]) >= threshold
    )
  })

  do.call(rbind, lapply(stats, as.data.frame))
}
