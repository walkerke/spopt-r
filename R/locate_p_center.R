#' P-Center Problem
#'
#' Solves the P-Center problem: minimize the maximum distance from any demand
#' point to its nearest facility by locating exactly p facilities.
#' This is an equity-focused (minimax) objective that ensures no demand point
#' is too far from service.
#'
#' @param demand An sf object representing demand points.
#' @param facilities An sf object representing candidate facility locations.
#' @param n_facilities Integer. Number of facilities to locate (p).
#' @param cost_matrix Optional. Pre-computed distance matrix.
#' @param distance_metric Distance metric: "euclidean" (default) or "manhattan".
#' @param verbose Logical. Print solver progress.
#'
#' @return A list with two sf objects:
#'   \itemize{
#'     \item `$demand`: Original demand sf with `.facility` column
#'     \item `$facilities`: Original facilities sf with `.selected` column
#'   }
#'   Metadata includes `max_distance` (the objective value).
#'
#' @details
#' The p-center problem minimizes the maximum distance between any demand point
#' and its assigned facility. This "minimax" objective ensures equitable access
#' by focusing on the worst-served location rather than average performance.
#'
#' The integer programming formulation is:
#' \deqn{\min W}
#' Subject to:
#' \deqn{\sum_j y_j = p}
#' \deqn{\sum_j x_{ij} = 1 \quad \forall i}
#' \deqn{x_{ij} \leq y_j \quad \forall i,j}
#' \deqn{\sum_j d_{ij} x_{ij} \leq W \quad \forall i}
#' \deqn{x_{ij}, y_j \in \{0,1\}}
#'
#' Where W is the maximum distance to minimize, \eqn{d_{ij}} is the distance
#' from demand i to facility j, \eqn{x_{ij} = 1} if demand i is assigned to
#' facility j, and \eqn{y_j = 1} if facility j is selected.
#'
#' @section Use Cases:
#' P-center is appropriate when equity and worst-case performance matter:
#' \itemize{
#'   \item **Emergency services**: Fire stations or ambulance depots where
#'     response time standards must be met for all residents
#'   \item **Equity-focused planning**: Ensuring no community is underserved,
#'     even if it increases average travel distance
#'   \item **Critical infrastructure**: Backup facilities or emergency shelters
#'     where everyone must be within reach
#'   \item **Service level guarantees**: When contracts or regulations specify
#'     maximum acceptable distance or response time
#' }
#'
#' For efficiency-focused objectives that minimize total travel, consider
#' [p_median()] instead.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' demand <- st_as_sf(data.frame(x = runif(50), y = runif(50)), coords = c("x", "y"))
#' facilities <- st_as_sf(data.frame(x = runif(15), y = runif(15)), coords = c("x", "y"))
#'
#' # Minimize maximum distance with 4 facilities
#' result <- p_center(demand, facilities, n_facilities = 4)
#'
#' # Maximum distance any demand point must travel
#' attr(result, "spopt")$max_distance
#' }
#'
#' @references
#' Hakimi, S. L. (1965). Optimum Distribution of Switching Centers in a
#' Communication Network and Some Related Graph Theoretic Problems.
#' Operations Research, 13(3), 462-475. \doi{10.1287/opre.13.3.462}
#'
#' @seealso [p_median()] for minimizing total weighted distance (efficiency objective)
#'
#' @export
p_center <- function(demand,
                     facilities,
                     n_facilities,
                     cost_matrix = NULL,
                     distance_metric = "euclidean",
                     verbose = FALSE) {
  if (!inherits(demand, "sf")) {
    stop("`demand` must be an sf object", call. = FALSE)
  }
  if (!inherits(facilities, "sf")) {
    stop("`facilities` must be an sf object", call. = FALSE)
  }

  if (is.null(cost_matrix)) {
    cost_matrix <- distance_matrix(demand, facilities, type = distance_metric)
  }

  # Validate cost matrix for NA/Inf values
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

  start_time <- Sys.time()

  # Call Rust MIP solver
  result <- rust_p_center(cost_matrix, as.integer(n_facilities))

  end_time <- Sys.time()

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
    algorithm = "p_center",
    n_selected = length(result$selected),
    n_facilities = n_facilities,
    objective = result$max_distance,
    max_distance = result$max_distance,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs"))
  )

  attr(output, "spopt") <- metadata
  class(output) <- c("spopt_pcenter", "spopt_locate", "list")

  output
}

# Greedy heuristic for p-center (internal)
greedy_p_center <- function(cost_matrix, p) {
  n_demand <- nrow(cost_matrix)
  n_fac <- ncol(cost_matrix)

  selected <- integer(0)
  remaining <- seq_len(n_fac)

  # Start with facility that minimizes max distance if only one chosen
  max_dists <- apply(cost_matrix, 2, max)
  first <- which.min(max_dists)
  selected <- c(selected, first)
  remaining <- setdiff(remaining, first)

  # Greedily add facilities
  while (length(selected) < p && length(remaining) > 0) {
    # Current max distance for each demand
    current_min_dists <- apply(cost_matrix[, selected, drop = FALSE], 1, min)

    # For each remaining facility, compute max distance if added
    best_fac <- remaining[1]
    best_max <- Inf

    for (fac in remaining) {
      new_min_dists <- pmin(current_min_dists, cost_matrix[, fac])
      new_max <- max(new_min_dists)
      if (new_max < best_max) {
        best_max <- new_max
        best_fac <- fac
      }
    }

    selected <- c(selected, best_fac)
    remaining <- setdiff(remaining, best_fac)
  }

  # Compute assignments
  min_dists <- apply(cost_matrix[, selected, drop = FALSE], 1, min)
  assignments <- sapply(seq_len(n_demand), function(i) {
    selected[which.min(cost_matrix[i, selected])]
  })

  list(
    selected = selected,
    assignments = assignments,
    max_distance = max(min_dists)
  )
}
