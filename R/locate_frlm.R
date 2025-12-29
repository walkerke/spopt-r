#' Flow Refueling Location Model (FRLM)
#'
#' Solves the Flow Refueling Location Model to optimally place refueling
#' facilities along network paths. This model maximizes the volume of
#' origin-destination flows that can be served given vehicle range constraints.
#'
#' @param flows A data frame or sf object containing flow information with columns:
#'   \itemize{
#'     \item \code{origin}: Origin identifier
#'     \item \code{destination}: Destination identifier
#'     \item \code{volume}: Flow volume (e.g., number of trips)
#'   }
#' @param candidates An sf object with candidate facility locations (points).
#' @param network A network representation. Can be:
#'   \itemize{
#'     \item An sfnetwork object
#'     \item An igraph object with edge weights
#'     \item A distance matrix between candidates
#'   }
#' @param vehicle_range Numeric. Maximum vehicle range (same units as network distances).
#' @param n_facilities Integer. Number of facilities to place.
#' @param method Character. Optimization method: "greedy" (default and currently only option).
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with class "spopt_frlm" containing:
#'   \itemize{
#'     \item \code{facilities}: The candidates sf object with a \code{.selected} column
#'     \item \code{selected_indices}: 1-based indices of selected facilities
#'     \item \code{coverage}: Coverage statistics
#'   }
#'   Metadata is stored in the "spopt" attribute.
#'
#' @details
#' The Flow Refueling Location Model (Kuby & Lim, 2005) addresses the problem
#' of locating refueling stations for range-limited vehicles (e.g., electric
#' vehicles, hydrogen fuel cell vehicles) along travel paths.
#'
#' A flow (origin-destination path) is "covered" if a vehicle can complete
#' the trip with refueling stops at the selected facilities, never exceeding
#' its maximum range between stops.
#'
#' This implementation uses a greedy heuristic that iteratively selects the
#' facility providing the greatest marginal increase in covered flow volume.
#'
#' @section Input Format:
#' For simple cases, you can provide:
#' \itemize{
#'   \item \code{flows}: Data frame with origin, destination, volume
#'   \item \code{candidates}: sf points for potential facility locations
#'   \item \code{network}: Distance matrix or network object
#' }
#'
#' @examples
#' \dontrun{
#' # Simple example with distance matrix
#' library(sf)
#'
#' # Create candidate locations
#' candidates <- st_as_sf(data.frame(
#'   id = 1:10,
#'   x = runif(10, 0, 100),
#'   y = runif(10, 0, 100)
#' ), coords = c("x", "y"))
#'
#' # Create flows (using candidate indices as origins/destinations)
#' flows <- data.frame(
#'   origin = c(1, 1, 3, 5),
#'   destination = c(8, 10, 7, 9),
#'   volume = c(100, 200, 150, 300)
#' )
#'
#' # Solve with vehicle range of 50 units
#' result <- frlm(flows, candidates, vehicle_range = 50, n_facilities = 3)
#'
#' # View selected facilities
#' result$facilities[result$facilities$.selected, ]
#' }
#'
#' @references
#' Kuby, M., & Lim, S. (2005). The flow-refueling location problem for
#' alternative-fuel vehicles. Socio-Economic Planning Sciences, 39(2), 125-145.
#'
#' @export
frlm <- function(flows,
                 candidates,
                 network = NULL,
                 vehicle_range,
                 n_facilities,
                 method = c("greedy"),
                 verbose = FALSE) {

  method <- match.arg(method)

  # Input validation
  if (!is.numeric(vehicle_range) || vehicle_range <= 0) {
    stop("`vehicle_range` must be a positive number", call. = FALSE)
  }

  if (!is.numeric(n_facilities) || n_facilities < 1) {
    stop("`n_facilities` must be a positive integer", call. = FALSE)
  }

  # Extract candidate coordinates
  if (!inherits(candidates, "sf")) {
    stop("`candidates` must be an sf object with point geometries", call. = FALSE)
  }

  n_candidates <- nrow(candidates)
  if (n_facilities > n_candidates) {
    stop("`n_facilities` cannot exceed number of candidates", call. = FALSE)
  }

  # Compute distance matrix between candidates if not provided
  if (is.null(network)) {
    dist_matrix <- as.matrix(sf::st_distance(candidates))
    # Convert to numeric (removes units)
    dist_matrix <- matrix(as.numeric(dist_matrix), nrow = nrow(dist_matrix))
  } else if (is.matrix(network)) {
    dist_matrix <- network
  } else {
    stop("Currently only distance matrix networks are supported", call. = FALSE)
  }

  # Validate flows
  if (!all(c("origin", "destination", "volume") %in% names(flows))) {
    stop("`flows` must have columns: origin, destination, volume", call. = FALSE)
  }

  n_flows <- nrow(flows)

  if (verbose) {
    message(sprintf(
      "FRLM: %d candidates, %d flows, range=%.1f, n_facilities=%d",
      n_candidates, n_flows, vehicle_range, n_facilities
    ))
  }

  # Build paths for each flow
  # For simplicity, assume origin and destination are indices into candidates
  # and the path goes through all candidates between them ordered by distance

  path_candidates <- integer(0)
  path_distances <- numeric(0)
  path_offsets <- integer(n_flows)

  for (i in seq_len(n_flows)) {
    path_offsets[i] <- length(path_candidates)

    origin <- flows$origin[i]
    destination <- flows$destination[i]

    # Get distances from origin to all candidates
    if (origin <= n_candidates && destination <= n_candidates) {
      origin_dists <- dist_matrix[origin, ]
      dest_dists <- dist_matrix[destination, ]

      # Simple path: candidates ordered by distance from origin
      # that are between origin and destination
      total_dist <- dist_matrix[origin, destination]

      # Find candidates roughly along the path
      # (within a corridor of total_dist from both endpoints)
      on_path <- which(origin_dists + dest_dists <= total_dist * 1.5)

      # Sort by distance from origin
      on_path <- on_path[order(origin_dists[on_path])]

      # Add to path arrays (0-based for Rust)
      path_candidates <- c(path_candidates, on_path - 1L)
      path_distances <- c(path_distances, origin_dists[on_path])
    }
  }

  # Ensure we have at least one candidate per path
  if (length(path_candidates) == 0) {
    stop("No valid paths could be constructed from flows and candidates", call. = FALSE)
  }

  start_time <- Sys.time()

  # Call Rust implementation
  result_list <- rust_frlm_greedy(
    n_candidates = as.integer(n_candidates),
    path_candidates = as.integer(path_candidates),
    path_offsets = as.integer(path_offsets),
    path_distances = as.numeric(path_distances),
    flow_volumes = as.numeric(flows$volume),
    vehicle_range = as.numeric(vehicle_range),
    n_facilities = as.integer(n_facilities)
  )

  end_time <- Sys.time()

  # Build result
  selected_indices <- result_list$selected
  candidates$.selected <- seq_len(nrow(candidates)) %in% selected_indices

  coverage <- list(
    covered_volume = result_list$covered_volume,
    total_volume = result_list$total_volume,
    coverage_pct = result_list$coverage_pct
  )

  result <- list(
    facilities = candidates,
    selected_indices = selected_indices,
    coverage = coverage
  )

  class(result) <- c("spopt_frlm", "spopt_locate", "list")

  # Attach metadata
  metadata <- list(
    algorithm = "frlm",
    method = method,
    n_candidates = n_candidates,
    n_facilities = result_list$n_selected,
    n_flows = n_flows,
    vehicle_range = vehicle_range,
    covered_volume = coverage$covered_volume,
    total_volume = coverage$total_volume,
    coverage_pct = coverage$coverage_pct,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs"))
  )

  attr(result, "spopt") <- metadata

  if (verbose) {
    message(sprintf(
      "  Selected %d facilities, coverage=%.1f%%, time=%.3fs",
      result_list$n_selected, coverage$coverage_pct,
      metadata$solve_time
    ))
  }

  result
}
