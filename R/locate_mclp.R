#' Maximum Coverage Location Problem (MCLP)
#'
#' Solves the Maximum Coverage Location Problem: maximize total weighted
#' demand covered by locating exactly p facilities.
#'
#' @param demand An sf object representing demand points.
#' @param facilities An sf object representing candidate facility locations.
#' @param service_radius Numeric. Maximum distance for coverage.
#' @param n_facilities Integer. Number of facilities to locate (p).
#' @param weight_col Character. Column name in `demand` containing demand weights.
#' @param cost_matrix Optional. Pre-computed distance matrix.
#' @param distance_metric Distance metric: "euclidean" (default) or "manhattan".
#' @param verbose Logical. Print solver progress.
#'
#' @return A list with two sf objects:
#'   \itemize{
#'     \item `$demand`: Original demand sf with `.covered` and `.facility` columns
#'     \item `$facilities`: Original facilities sf with `.selected` column
#'   }
#'   Metadata is stored in the "spopt" attribute.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Create demand with weights
#' demand <- st_as_sf(data.frame(
#'   x = runif(50), y = runif(50), population = rpois(50, 100)
#' ), coords = c("x", "y"))
#' facilities <- st_as_sf(data.frame(x = runif(10), y = runif(10)), coords = c("x", "y"))
#'
#' # Maximize population coverage with 3 facilities
#' result <- mclp(demand, facilities, service_radius = 0.3,
#'                n_facilities = 3, weight_col = "population")
#'
#' attr(result, "spopt")$coverage_pct
#' }
#'
#' @export
mclp <- function(demand,
                 facilities,
                 service_radius,
                 n_facilities,
                 weight_col,
                 cost_matrix = NULL,
                 distance_metric = "euclidean",
                 verbose = FALSE) {
  # Input validation
  if (!inherits(demand, "sf")) {
    stop("`demand` must be an sf object", call. = FALSE)
  }
  if (!inherits(facilities, "sf")) {
    stop("`facilities` must be an sf object", call. = FALSE)
  }
  if (!weight_col %in% names(demand)) {
    stop(paste0("Weight column '", weight_col, "' not found in demand"), call. = FALSE)
  }

  weights <- as.numeric(demand[[weight_col]])
  if (any(is.na(weights))) {
    stop("Weight column contains NA values", call. = FALSE)
  }

  # Compute distance matrix if not provided
  if (is.null(cost_matrix)) {
    cost_matrix <- distance_matrix(demand, facilities, type = distance_metric)
  }

  n_demand <- nrow(demand)
  n_fac <- nrow(facilities)

  if (n_facilities > n_fac) {
    stop("`n_facilities` cannot exceed number of candidate facilities", call. = FALSE)
  }

  start_time <- Sys.time()

  # Call Rust solver
  result <- rust_mclp(cost_matrix, weights, service_radius, as.integer(n_facilities))

  end_time <- Sys.time()

  # Build result sf objects
  demand_result <- demand
  facilities_result <- facilities

  selected_indices <- result$selected

  # Determine coverage
  covered <- rep(FALSE, n_demand)
  for (i in seq_len(n_demand)) {
    for (j in selected_indices) {
      if (cost_matrix[i, j] <= service_radius) {
        covered[i] <- TRUE
        break
      }
    }
  }

  demand_result$.covered <- covered
  demand_result$.facility <- NA_integer_

  for (i in seq_len(n_demand)) {
    if (covered[i]) {
      dists <- cost_matrix[i, selected_indices]
      demand_result$.facility[i] <- selected_indices[which.min(dists)]
    }
  }

  facilities_result$.selected <- seq_len(n_fac) %in% selected_indices
  facilities_result$.n_assigned <- 0L

  for (j in selected_indices) {
    facilities_result$.n_assigned[j] <- sum(demand_result$.facility == j, na.rm = TRUE)
  }

  output <- list(
    demand = demand_result,
    facilities = facilities_result
  )

  metadata <- list(
    algorithm = "mclp",
    n_selected = result$n_selected,
    objective = result$objective,
    service_radius = service_radius,
    n_facilities = n_facilities,
    covered_weight = result$covered_weight,
    total_weight = result$total_weight,
    coverage_pct = result$coverage_pct,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs"))
  )

  attr(output, "spopt") <- metadata
  class(output) <- c("spopt_mclp", "spopt_locate", "list")

  output
}
