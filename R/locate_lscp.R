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
#' @details
#' The LSCP minimizes the number of facilities required to ensure that every
#' demand point is within the service radius of at least one facility. This is
#' a mandatory coverage model where full coverage is required.
#'
#' The integer programming formulation is:
#' \deqn{\min \sum_j y_j}
#' Subject to:
#' \deqn{\sum_j a_{ij} y_j \geq 1 \quad \forall i}
#' \deqn{y_j \in \{0,1\}}
#'
#' Where \eqn{y_j = 1} if facility j is selected, and \eqn{a_{ij} = 1} if
#' facility j can cover demand point i (distance \eqn{\leq} service radius).
#'
#' @section Use Cases:
#' LSCP is appropriate when complete coverage is mandatory:
#' \itemize{
#'   \item **Emergency services**: Fire stations, ambulance depots, or hospitals
#'     where every resident must be reachable within a response time standard
#'   \item **Public services**: Schools, polling places, or post offices where
#'     universal access is required by law or policy
#'   \item **Infrastructure**: Cell towers or utility substations where gaps
#'     in coverage are unacceptable
#'   \item **Retail/logistics**: Warehouse locations to ensure all customers
#'     can receive same-day or next-day delivery
#' }
#'
#' For situations where complete coverage is not required or not feasible
#' within budget constraints, consider [mclp()] instead.
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
#' @references
#' Toregas, C., Swain, R., ReVelle, C., & Bergman, L. (1971). The Location of
#' Emergency Service Facilities. Operations Research, 19(6), 1363-1373.
#' \doi{10.1287/opre.19.6.1363}
#'
#' @seealso [mclp()] for maximizing coverage with a fixed number of facilities
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
