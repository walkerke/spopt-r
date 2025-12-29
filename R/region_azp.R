#' Automatic Zoning Procedure (AZP)
#'
#' Performs regionalization using the Automatic Zoning Procedure algorithm.
#' AZP uses local search to minimize within-region heterogeneity while
#' maintaining spatial contiguity. Three variants are available: basic
#' (greedy), tabu search, and simulated annealing.
#'
#' @param data An sf object with polygon or point geometries.
#' @param attrs Character vector of column names to use for clustering
#'   (e.g., `c("var1", "var2")`). If NULL, uses all numeric columns.
#' @param n_regions Integer. Number of regions (clusters) to create.
#' @param weights Spatial weights specification. One of "queen" (default),
#'   "rook", or an nb object from spdep.
#' @param method Character. Optimization method: "basic" (greedy local search),
#'   "tabu" (tabu search), or "sa" (simulated annealing). Default is "tabu".
#' @param max_iterations Integer. Maximum number of iterations (default 100).
#' @param tabu_length Integer. Length of tabu list for tabu method (default 10).
#' @param cooling_rate Numeric. Cooling rate for SA method, between 0 and 1
#'   (default 0.99).
#' @param initial_temperature Numeric. Initial temperature for SA method.
#'   If 0 (default), automatically set based on initial objective.
#' @param scale Logical. If TRUE (default), standardize attributes before clustering.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical. Print progress messages.
#'
#' @return An sf object with a `.region` column containing cluster assignments.
#'   Metadata is stored in the "spopt" attribute, including:
#'   \itemize{
#'     \item algorithm: "azp"
#'     \item method: The optimization method used
#'     \item n_regions: Number of regions created
#'     \item objective: Total within-region sum of squared deviations
#'     \item solve_time: Time to solve in seconds
#'   }
#'
#' @details
#' The Automatic Zoning Procedure (AZP) was introduced by Openshaw (1977) and
#' refined by Openshaw & Rao (1995). It is a local search algorithm that:
#'
#' 1. Starts with an initial random partition into n_regions
#' 2. Iteratively moves border areas between regions to reduce heterogeneity
#' 3. Maintains spatial contiguity throughout
#' 4. Terminates when no improving moves are found
#'
#' Three variants are available:
#' \itemize{
#'   \item \strong{basic}: Greedy local search that only accepts improving moves
#'   \item \strong{tabu}: Tabu search that can accept non-improving moves to escape
#'     local optima, with a tabu list preventing cycling
#'   \item \strong{sa}: Simulated annealing that accepts worse moves with decreasing
#'     probability as temperature cools
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#'
#' # Basic AZP with 8 regions
#' result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 8)
#'
#' # Tabu search (often finds better solutions)
#' result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 8,
#'               method = "tabu", tabu_length = 15)
#'
#' # Simulated annealing
#' result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 8,
#'               method = "sa", cooling_rate = 0.95)
#'
#' # View results
#' plot(result[".region"])
#' }
#'
#' @references
#' Openshaw, S. (1977). A geographical solution to scale and aggregation
#' problems in region-building, partitioning and spatial modelling.
#' Transactions of the Institute of British Geographers, 2(4), 459-472.
#'
#' Openshaw, S., & Rao, L. (1995). Algorithms for reengineering 1991 Census
#' geography. Environment and Planning A, 27(3), 425-446.
#'
#' @export
azp <- function(data,
                attrs = NULL,
                n_regions,
                weights = "queen",
                method = c("tabu", "basic", "sa"),
                max_iterations = 100L,
                tabu_length = 10L,
                cooling_rate = 0.99,
                initial_temperature = 0,
                scale = TRUE,
                seed = NULL,
                verbose = FALSE) {
  # Input validation
  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object", call. = FALSE)
  }

  method <- match.arg(method)

  if (!is.numeric(n_regions) || n_regions < 2) {
    stop("`n_regions` must be an integer >= 2", call. = FALSE)
  }

  # Determine which columns to check for NAs
  check_cols <- if (!is.null(attrs)) attrs else character(0)

  # Validate data: remove empty geometries, check for NAs
  validated <- validate_regionalization_data(data, check_cols, call_name = "azp")
  data <- validated$data

  n <- nrow(data)
  if (n_regions >= n) {
    stop("`n_regions` must be less than number of observations", call. = FALSE)
  }

  if (!is.numeric(max_iterations) || max_iterations < 1) {
    stop("`max_iterations` must be a positive integer", call. = FALSE)
  }

  if (method == "tabu" && (!is.numeric(tabu_length) || tabu_length < 1)) {
    stop("`tabu_length` must be a positive integer", call. = FALSE)
  }

  if (method == "sa" && (!is.numeric(cooling_rate) || cooling_rate <= 0 || cooling_rate >= 1)) {
    stop("`cooling_rate` must be between 0 and 1 (exclusive)", call. = FALSE)
  }

  # Extract attributes
  attr_matrix <- extract_attrs(data, attrs)

  if (scale) {
    attr_matrix <- scale(attr_matrix)
  }

  # Prepare spatial weights
  nb <- prepare_weights(data, weights)

  # Convert nb to adjacency indices
  adj <- nb_to_adj_indices(nb)

  if (verbose) {
    message(sprintf(
      "AZP: n=%d, n_regions=%d, method=%s, attrs=%d",
      n, n_regions, method, ncol(attr_matrix)
    ))
  }

  # Call Rust implementation
  start_time <- Sys.time()

  result_list <- rust_azp(
    attrs = attr_matrix,
    n_regions = as.integer(n_regions),
    adj_i = adj$i,
    adj_j = adj$j,
    method = method,
    max_iterations = as.integer(max_iterations),
    tabu_length = as.integer(tabu_length),
    cooling_rate = as.numeric(cooling_rate),
    initial_temperature = as.numeric(initial_temperature),
    seed = if (!is.null(seed)) as.integer(seed) else NULL
  )

  end_time <- Sys.time()

  # Extract results
  labels <- result_list$labels
  objective <- result_list$objective

  # Attach results to sf object
  result <- data
  result$.region <- as.integer(labels)

  if (verbose) {
    message(sprintf(
      "  Result: %d regions, objective=%.4f, time=%.3fs",
      length(unique(labels)), objective,
      as.numeric(difftime(end_time, start_time, units = "secs"))
    ))
  }

  # Attach metadata
  metadata <- list(
    algorithm = "azp",
    method = method,
    n_regions = length(unique(labels)),
    objective = objective,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    scaled = scale,
    max_iterations = max_iterations
  )

  attach_spopt_metadata(result, metadata)
}
