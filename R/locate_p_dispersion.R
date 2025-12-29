#' P-Dispersion Problem
#'
#' Solves the P-Dispersion problem: maximize the minimum distance between
#' any two selected facilities. This ensures facilities are spread out
#' as much as possible.
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
