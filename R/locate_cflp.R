#' Capacitated Facility Location Problem (CFLP)
#'
#' Solves the Capacitated Facility Location Problem: minimize total weighted
#' distance from demand points to facilities, subject to capacity constraints
#' at each facility. Unlike standard p-median, facilities have limited capacity
#' and demand may need to be split across multiple facilities.
#'
#' @param demand An sf object representing demand points.
#' @param facilities An sf object representing candidate facility locations.
#' @param n_facilities Integer. Number of facilities to locate. Set to 0 if
#'   using `facility_cost_col` to determine optimal number.
#' @param weight_col Character. Column name in `demand` containing demand weights
#'   (e.g., population, customers, volume).
#' @param capacity_col Character. Column name in `facilities` containing capacity
#'   of each facility.
#' @param facility_cost_col Optional character. Column name in `facilities`
#'   containing fixed cost to open each facility. If provided and `n_facilities = 0`,
#'   the solver determines the optimal number of facilities to minimize total cost.
#' @param cost_matrix Optional. Pre-computed distance/cost matrix (demand x facilities).
#' @param distance_metric Distance metric: "euclidean" (default) or "manhattan".
#' @param verbose Logical. Print solver progress.
#'
#' @return A list with two sf objects:
#'   \itemize{
#'     \item `$demand`: Original demand sf with `.facility` column (primary assignment)
#'       and `.split` column (TRUE if demand is split across facilities)
#'     \item `$facilities`: Original facilities sf with `.selected`, `.n_assigned`,
#'       and `.utilization` columns
#'   }
#'   Metadata is stored in the "spopt" attribute, including:
#'   \itemize{
#'     \item `objective`: Total cost (transportation + facility costs if applicable)
#'     \item `mean_distance`: Mean weighted distance
#'     \item `n_split_demand`: Number of demand points split across facilities
#'     \item `allocation_matrix`: Full allocation matrix (n_demand x n_facilities)
#'   }
#'
#' @details
#' The CFLP extends the p-median problem by adding capacity constraints. Each
#' facility \eqn{j} has a maximum capacity \eqn{Q_j}, and the total demand
#' assigned to it cannot exceed this capacity.
#'
#' When demand exceeds available capacity at the nearest facility, the solver
#' may split demand across multiple facilities. The `.split` column indicates
#' which demand points have been split, and the `allocation_matrix` in metadata
#' shows the exact fractions.
#'
#' Two modes of operation:
#' \enumerate{
#'   \item **Fixed number**: Set `n_facilities` to select exactly that many facilities
#'   \item **Cost-based**: Set `n_facilities = 0` and provide `facility_cost_col` to
#'     let the solver determine the optimal number based on fixed + variable costs
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Demand points with population
#' demand <- st_as_sf(data.frame(
#'   x = runif(100), y = runif(100), population = rpois(100, 500)
#' ), coords = c("x", "y"))
#'
#' # Facilities with varying capacities
#' facilities <- st_as_sf(data.frame(
#'   x = runif(15), y = runif(15),
#'   capacity = c(rep(5000, 5), rep(10000, 5), rep(20000, 5)),
#'   fixed_cost = c(rep(100, 5), rep(200, 5), rep(400, 5))
#' ), coords = c("x", "y"))
#'
#' # Fixed number of facilities
#' result <- cflp(demand, facilities, n_facilities = 5,
#'                weight_col = "population", capacity_col = "capacity")
#'
#' # Check utilization
#' result$facilities[result$facilities$.selected, c("capacity", ".utilization")]
#'
#' # Cost-based (optimal number of facilities)
#' result <- cflp(demand, facilities, n_facilities = 0,
#'                weight_col = "population", capacity_col = "capacity",
#'                facility_cost_col = "fixed_cost")
#' attr(result, "spopt")$n_selected
#' }
#'
#' @references
#' Cornuejols, G., Fisher, M. L., & Nemhauser, G. L. (1977). Location of bank
#' accounts to optimize float: An analytic study of exact and approximate
#' algorithms. Management Science, 23(8), 789-810.
#'
#' @export
cflp <- function(demand,
                 facilities,
                 n_facilities,
                 weight_col,
                 capacity_col,
                 facility_cost_col = NULL,
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
  if (!capacity_col %in% names(facilities)) {
    stop(paste0("Capacity column '", capacity_col, "' not found in facilities"), call. = FALSE)
  }

  weights <- as.numeric(demand[[weight_col]])
  capacities <- as.numeric(facilities[[capacity_col]])

  if (any(is.na(weights))) {
    stop("Weight column contains NA values", call. = FALSE)
  }
  if (any(is.na(capacities))) {
    stop("Capacity column contains NA values", call. = FALSE)
  }
  if (any(capacities <= 0)) {
    stop("All capacities must be positive", call. = FALSE)
  }

  # Check if problem is feasible

  total_demand <- sum(weights)
  total_capacity <- sum(capacities)
  if (total_capacity < total_demand) {
    stop(sprintf(
      "Total capacity (%.2f) is less than total demand (%.2f). Problem is infeasible.",
      total_capacity, total_demand
    ), call. = FALSE)
  }

  # Handle facility costs
  facility_costs <- NULL
  if (!is.null(facility_cost_col)) {
    if (!facility_cost_col %in% names(facilities)) {
      stop(paste0("Facility cost column '", facility_cost_col, "' not found in facilities"),
           call. = FALSE)
    }
    facility_costs <- as.numeric(facilities[[facility_cost_col]])
    if (any(is.na(facility_costs))) {
      stop("Facility cost column contains NA values", call. = FALSE)
    }
  }

  # Validate n_facilities
  n_fac <- nrow(facilities)
  if (n_facilities > n_fac) {
    stop("`n_facilities` cannot exceed number of candidate facilities", call. = FALSE)
  }
  if (n_facilities == 0 && is.null(facility_costs)) {
    stop("When n_facilities = 0, must provide facility_cost_col", call. = FALSE)
  }

  # Compute cost matrix if needed
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

  start_time <- Sys.time()

  # Call Rust solver
  result <- rust_cflp(
    cost_matrix,
    weights,
    capacities,
    as.integer(n_facilities),
    facility_costs
  )

  end_time <- Sys.time()

  # Check for errors
  if (!is.null(result$error)) {
    stop(result$error, call. = FALSE)
  }

  # Build output
  demand_result <- demand
  facilities_result <- facilities

  demand_result$.facility <- result$assignments  # Primary assignment (1-based)

  # Determine which demands are split
  allocation_matrix <- matrix(
    result$allocation_matrix,
    nrow = n_demand,
    ncol = n_fac,
    byrow = TRUE
  )
  primary_allocation <- sapply(seq_len(n_demand), function(i) {
    allocation_matrix[i, result$assignments[i]]
  })
  demand_result$.split <- primary_allocation < 0.999

  # Mark selected facilities
  selected_indices <- result$selected
  facilities_result$.selected <- seq_len(n_fac) %in% selected_indices
  facilities_result$.n_assigned <- 0L
  facilities_result$.utilization <- 0.0

  for (j in seq_len(n_fac)) {
    facilities_result$.n_assigned[j] <- sum(result$assignments == j)
    facilities_result$.utilization[j] <- result$utilizations[j]
  }

  output <- list(
    demand = demand_result,
    facilities = facilities_result
  )

  metadata <- list(
    algorithm = "cflp",
    n_selected = result$n_selected,
    n_facilities = n_facilities,
    objective = result$objective,
    mean_distance = result$mean_distance,
    n_split_demand = result$n_split_demand,
    allocation_matrix = allocation_matrix,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    solver_status = result$status
  )

  attr(output, "spopt") <- metadata
  class(output) <- c("spopt_cflp", "spopt_locate", "list")

  output
}
