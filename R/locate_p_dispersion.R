#' P-Dispersion Problem
#'
#' Solves the P-Dispersion problem: maximize the minimum distance between
#' any two selected facilities. This "maximin" objective ensures facilities
#' are spread out as much as possible.
#'
#' @param facilities An sf object representing candidate facility locations.
#'   Note: This problem does not use demand points.
#' @param n_facilities Integer. Number of facilities to locate (p).
#' @param cost_matrix Optional. Pre-computed inter-facility distance matrix.
#' @param distance_metric Distance metric: "euclidean" (default) or "manhattan".
#' @param verbose Logical. Print solver progress.
#'
#' @return An sf object (the facilities input) with a `.selected` column.
#'   Metadata includes `min_distance` (the objective value).
#'
#' @details
#' The p-dispersion problem selects p facilities from a set of candidates such
#' that the minimum pairwise distance between any two selected facilities is
#' maximized. Unlike p-median or p-center, this problem does not consider
#' demand points---it focuses solely on spreading facilities apart.
#'
#' The mixed integer programming formulation uses a Big-M approach:
#' \deqn{\max D}
#' Subject to:
#' \deqn{\sum_j y_j = p}
#' \deqn{D \leq d_{ij} + M(2 - y_i - y_j) \quad \forall i < j}
#' \deqn{y_j \in \{0,1\}, \quad D \geq 0}
#'
#' Where D is the minimum separation distance to maximize, \eqn{d_{ij}} is the
#' distance between facilities i and j, \eqn{y_j = 1} if facility j is selected,
#' and M is a large constant. When both facilities i and j are selected
#' (\eqn{y_i = y_j = 1}), the constraint reduces to \eqn{D \leq d_{ij}},
#' ensuring D is at most the distance between any pair of selected facilities.
#'
#' @section Use Cases:
#' P-dispersion is appropriate when facilities should be spread apart:
#' \itemize{
#'   \item **Obnoxious facilities**: Hazardous waste sites, prisons, or other
#'     undesirable facilities that should be separated from each other
#'   \item **Franchise territories**: Retail locations where stores should not
#'     cannibalize each other's market
#'   \item **Redundant systems**: Backup servers or emergency caches that should
#'     be geographically distributed for resilience
#'   \item **Monitoring networks**: Air quality sensors or seismic monitors that
#'     should cover distinct areas without overlap
#'   \item **Spatial sampling**: Selecting representative sample locations that
#'     are well-distributed across a study area
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' facilities <- st_as_sf(data.frame(x = runif(20), y = runif(20)), coords = c("x", "y"))
#'
#' # Select 5 facilities maximally dispersed
#' result <- p_dispersion(facilities, n_facilities = 5)
#'
#' # Minimum distance between any two selected facilities
#' attr(result, "spopt")$min_distance
#' }
#'
#' @references
#' Kuby, M. J. (1987). Programming Models for Facility Dispersion: The
#' p-Dispersion and Maxisum Dispersion Problems. Geographical Analysis,
#' 19(4), 315-329. \doi{10.1111/j.1538-4632.1987.tb00133.x}
#'
#' @export
p_dispersion <- function(facilities,
                         n_facilities,
                         cost_matrix = NULL,
                         distance_metric = "euclidean",
                         verbose = FALSE) {
  if (!inherits(facilities, "sf")) {
    stop("`facilities` must be an sf object", call. = FALSE)
  }

  n_fac <- nrow(facilities)

  if (n_facilities > n_fac) {
    stop("`n_facilities` cannot exceed number of facilities", call. = FALSE)
  }
  if (n_facilities < 2) {
    stop("`n_facilities` must be at least 2", call. = FALSE)
  }

  # Compute inter-facility distance matrix
  if (is.null(cost_matrix)) {
    cost_matrix <- distance_matrix(facilities, facilities, type = distance_metric)
  }

  start_time <- Sys.time()

  # Call Rust MIP solver
  result <- rust_p_dispersion(cost_matrix, as.integer(n_facilities))

  end_time <- Sys.time()

  facilities_result <- facilities
  selected_indices <- result$selected
  facilities_result$.selected <- seq_len(n_fac) %in% selected_indices

  metadata <- list(
    algorithm = "p_dispersion",
    n_selected = length(result$selected),
    n_facilities = n_facilities,
    objective = result$min_distance,
    min_distance = result$min_distance,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs"))
  )

  attr(facilities_result, "spopt") <- metadata

  facilities_result
}

# Greedy heuristic for p-dispersion (internal)
greedy_p_dispersion <- function(cost_matrix, p) {
  n <- nrow(cost_matrix)

  # Start with the two facilities furthest apart
  max_dist <- 0
  best_pair <- c(1, 2)

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (cost_matrix[i, j] > max_dist) {
        max_dist <- cost_matrix[i, j]
        best_pair <- c(i, j)
      }
    }
  }

  selected <- best_pair
  remaining <- setdiff(seq_len(n), selected)

  # Greedily add facilities that maximize min distance to selected set
  while (length(selected) < p && length(remaining) > 0) {
    best_fac <- remaining[1]
    best_min <- 0

    for (fac in remaining) {
      min_dist_to_selected <- min(cost_matrix[fac, selected])
      if (min_dist_to_selected > best_min) {
        best_min <- min_dist_to_selected
        best_fac <- fac
      }
    }

    selected <- c(selected, best_fac)
    remaining <- setdiff(remaining, best_fac)
  }

  # Compute minimum distance between any pair of selected
  min_distance <- Inf
  for (i in seq_along(selected)) {
    for (j in seq_along(selected)) {
      if (i < j) {
        d <- cost_matrix[selected[i], selected[j]]
        if (d < min_distance) {
          min_distance <- d
        }
      }
    }
  }

  list(
    selected = selected,
    min_distance = min_distance
  )
}
