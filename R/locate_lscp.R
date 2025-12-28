#' Location Set Covering Problem (LSCP)
#'
#' Solves the Location Set Covering Problem: find the minimum number of
#' facilities needed to cover all demand points within a given service radius.
#'
#' @param demand An sf object representing demand points (or polygons, using centroids).
#' @param facilities An sf object representing candidate facility locations.
#' @param service_radius Numeric. Maximum distance for a facility to cover a demand point.
#' @param cost_matrix Optional. Pre-computed distance matrix (demand x facilities).
#'   If NULL, computed from geometries.
#' @param distance_metric Distance metric: "euclidean" (default) or "manhattan".
#' @param verbose Logical. Print solver progress.
#'
#' @return A list with two sf objects:
#'   \itemize{
#'     \item `$demand`: Original demand sf with `.covered` column (logical)
#'     \item `$facilities`: Original facilities sf with `.selected` column (logical)
#'   }
#'   Metadata is stored in the "spopt" attribute.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Create demand and facility points
#' demand <- st_as_sf(data.frame(x = runif(50), y = runif(50)), coords = c("x", "y"))
#' facilities <- st_as_sf(data.frame(x = runif(10), y = runif(10)), coords = c("x", "y"))
#'
#' # Find minimum facilities to cover all demand within 0.3 units
#' result <- lscp(demand, facilities, service_radius = 0.3)
#'
#' # View selected facilities
#' result$facilities[result$facilities$.selected, ]
#' }
#'
#' @export
lscp <- function(demand,
                 facilities,
                 service_radius,
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
  if (!is.numeric(service_radius) || service_radius <= 0) {
    stop("`service_radius` must be a positive number", call. = FALSE)
  }

  # Compute distance matrix if not provided
  if (is.null(cost_matrix)) {
    cost_matrix <- distance_matrix(demand, facilities, type = distance_metric)
  }

  n_demand <- nrow(demand)
  n_facilities <- nrow(facilities)

  start_time <- Sys.time()

  # Call Rust solver
  result <- rust_lscp(cost_matrix, service_radius)

  end_time <- Sys.time()

  # Build result sf objects
  demand_result <- demand
  facilities_result <- facilities

  # Determine which demand points are covered
  selected_indices <- result$selected  # 1-based
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

  # Assign each demand to nearest selected facility
  for (i in seq_len(n_demand)) {
    if (covered[i]) {
      dists <- cost_matrix[i, selected_indices]
      demand_result$.facility[i] <- selected_indices[which.min(dists)]
    }
  }

  # Mark selected facilities
  facilities_result$.selected <- seq_len(n_facilities) %in% selected_indices
  facilities_result$.n_assigned <- 0L

  for (j in selected_indices) {
    facilities_result$.n_assigned[j] <- sum(demand_result$.facility == j, na.rm = TRUE)
  }

  # Build output list
  output <- list(
    demand = demand_result,
    facilities = facilities_result
  )

  # Attach metadata
  metadata <- list(
    algorithm = "lscp",
    n_selected = result$n_selected,
    objective = result$objective,
    service_radius = service_radius,
    covered_demand = result$covered_demand,
    total_demand = result$total_demand,
    coverage_pct = result$coverage_pct,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    solver_status = result$status,
    uncoverable_demand = result$uncoverable_demand
  )

  attr(output, "spopt") <- metadata
  class(output) <- c("spopt_lscp", "spopt_locate", "list")

  output
}
