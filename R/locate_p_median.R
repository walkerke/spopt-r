#' P-Median Problem
#'
#' Solves the P-Median problem: minimize total weighted distance from demand
#' points to their assigned facilities by locating exactly p facilities.
#'
#' @param demand An sf object representing demand points.
#' @param facilities An sf object representing candidate facility locations.
#' @param n_facilities Integer. Number of facilities to locate (p).
#' @param weight_col Character. Column name in `demand` containing demand weights.
#' @param cost_matrix Optional. Pre-computed distance matrix.
#' @param distance_metric Distance metric: "euclidean" (default) or "manhattan".
#' @param verbose Logical. Print solver progress.
#'
#' @return A list with two sf objects:
#'   \itemize{
#'     \item `$demand`: Original demand sf with `.facility` column (assigned facility)
#'     \item `$facilities`: Original facilities sf with `.selected` and `.n_assigned` columns
#'   }
#'   Metadata is stored in the "spopt" attribute.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' demand <- st_as_sf(data.frame(
#'   x = runif(100), y = runif(100), population = rpois(100, 500)
#' ), coords = c("x", "y"))
#' facilities <- st_as_sf(data.frame(x = runif(20), y = runif(20)), coords = c("x", "y"))
#'
#' # Locate 5 facilities minimizing total weighted distance
#' result <- p_median(demand, facilities, n_facilities = 5, weight_col = "population")
#'
#' # Mean distance to assigned facility
#' attr(result, "spopt")$mean_distance
#' }
#'
#' @export
p_median <- function(demand,
                     facilities,
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

  if (is.null(cost_matrix)) {
    cost_matrix <- distance_matrix(demand, facilities, type = distance_metric)
  }

  n_demand <- nrow(demand)
  n_fac <- nrow(facilities)

  if (n_facilities > n_fac) {
    stop("`n_facilities` cannot exceed number of candidate facilities", call. = FALSE)
  }

  start_time <- Sys.time()

  result <- rust_p_median(cost_matrix, weights, as.integer(n_facilities))

  end_time <- Sys.time()

  # Build output
  demand_result <- demand
  facilities_result <- facilities

  demand_result$.facility <- result$assignments  # 1-based facility index

  selected_indices <- result$selected
  facilities_result$.selected <- seq_len(n_fac) %in% selected_indices
  facilities_result$.n_assigned <- 0L

  for (j in selected_indices) {
    facilities_result$.n_assigned[j] <- sum(result$assignments == j)
  }

  output <- list(
    demand = demand_result,
    facilities = facilities_result
  )

  metadata <- list(
    algorithm = "p_median",
    n_selected = result$n_selected,
    n_facilities = n_facilities,
    objective = result$objective,
    mean_distance = result$mean_distance,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs"))
  )

  attr(output, "spopt") <- metadata
  class(output) <- c("spopt_pmedian", "spopt_locate", "list")

  output
}
