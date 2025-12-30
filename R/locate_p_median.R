#' P-Median Problem
#'
#' Solves the P-Median problem: minimize total weighted distance from demand
#' points to their assigned facilities by locating exactly p facilities.
#' This is an efficiency-focused objective that minimizes overall travel burden.
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
#'
#'   \itemize{
#'     \item `$demand`: Original demand sf with `.facility` column (assigned facility)
#'     \item `$facilities`: Original facilities sf with `.selected` and `.n_assigned` columns
#'   }
#'   Metadata is stored in the "spopt" attribute.
#'
#' @details
#' The p-median problem minimizes the total weighted distance (or travel cost)
#' between demand points and their nearest assigned facility. It is the most
#' widely used location model for efficiency-oriented facility siting.
#'
#' The integer programming formulation is:
#' \deqn{\min \sum_i \sum_j w_i d_{ij} x_{ij}}
#' Subject to:
#' \deqn{\sum_j y_j = p}
#' \deqn{\sum_j x_{ij} = 1 \quad \forall i}
#' \deqn{x_{ij} \leq y_j \quad \forall i,j}
#' \deqn{x_{ij}, y_j \in \{0,1\}}
#'
#' Where \eqn{w_i} is the demand weight at location i, \eqn{d_{ij}} is the
#' distance from demand i to facility j, \eqn{x_{ij} = 1} if demand i is
#' assigned to facility j, and \eqn{y_j = 1} if facility j is selected.
#'
#' @section Use Cases:
#' P-median is appropriate when minimizing total travel cost or distance:
#' \itemize{
#'   \item **Public facilities**: Schools, libraries, or community centers where
#'     the goal is to minimize total student/patron travel
#'   \item **Warehouses and distribution**: Locating distribution centers to
#'     minimize total shipping costs to customers
#'   \item **Healthcare**: Positioning clinics to minimize aggregate patient
#'     travel time across a population
#'   \item **Service depots**: Locating maintenance facilities to minimize
#'     total technician travel to service calls
#' }
#'
#' For equity-focused objectives where no demand point should be too far,
#' consider [p_center()] instead.
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
#' @references
#' Hakimi, S. L. (1964). Optimum Locations of Switching Centers and the
#' Absolute Centers and Medians of a Graph. Operations Research, 12(3), 450-459.
#' \doi{10.1287/opre.12.3.450}
#'
#' @seealso [p_center()] for minimizing maximum distance (equity objective)
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

  # Validate cost matrix
  if (any(is.na(cost_matrix))) {
    n_na <- sum(is.na(cost_matrix))
    warning(sprintf(
      "cost_matrix contains %d NA values (unreachable points). Replacing with large value.",
      n_na
    ))
    max_cost <- max(cost_matrix, na.rm = TRUE)
    cost_matrix[is.na(cost_matrix)] <- max_cost * 100
  }
  if (any(is.infinite(cost_matrix))) {
    n_inf <- sum(is.infinite(cost_matrix))
    warning(sprintf(
      "cost_matrix contains %d Inf values. Replacing with large value.",
      n_inf
    ))
    finite_max <- max(cost_matrix[is.finite(cost_matrix)])
    cost_matrix[is.infinite(cost_matrix)] <- finite_max * 100
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
